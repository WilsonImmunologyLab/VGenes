import os
import shutil
import subprocess
import sys
import time
from collections import Counter

import VGenesSQL


DEFAULT_DEFINE_CLONES_DISTANCE = "0.15"
DEFAULT_DEFINE_CLONES_MODE = "gene"
DEFAULT_DEFINE_CLONES_LINK = "single"
DEFAULT_DEFINE_CLONES_MODEL = "ham"


def looks_like_nucleotide(sequence):
    seq = str(sequence or "").upper()
    if not seq:
        return False
    return all(base in {"A", "T", "C", "G", "N"} for base in seq)


def species_paths(species):
    if species == 'Human':
        return {
            'db_v': 'IG/Human/HumanVGenes.nt',
            'db_j': 'IG/Human/HumanJGenes.nt',
            'db_d': 'IG/Human/HumanDGenes.nt',
            'repo_v': '../IgBlast/IG/Human/HumanVGenes.fasta',
            'repo_d': '../IgBlast/IG/Human/HumanDGenes.fasta',
            'repo_j': '../IgBlast/IG/Human/HumanJGenes.fasta',
            'organism': 'human',
            'aux': 'optional_file/human_gl.aux',
        }
    if species == 'Mouse':
        return {
            'db_v': 'IG/Mouse/MouseVGenes.nt',
            'db_j': 'IG/Mouse/MouseJGenes.nt',
            'db_d': 'IG/Mouse/MouseDGenes.nt',
            'repo_v': '../IgBlast/IG/Mouse/MouseVGenes.fasta',
            'repo_d': '../IgBlast/IG/Mouse/MouseDGenes.fasta',
            'repo_j': '../IgBlast/IG/Mouse/MouseJGenes.fasta',
            'organism': 'mouse',
            'aux': 'optional_file/mouse_gl.aux',
        }
    raise ValueError('Unsupported species for Change-O cloning: ' + str(species))


def write_fasta(data_rows, temp_folder):
    time_stamp = str(int(time.time() * 100))
    seq_pathname = os.path.join(temp_folder, time_stamp + '.fasta')
    seq_index = 2 if looks_like_nucleotide(data_rows[0][2]) else 1
    with open(seq_pathname, 'w') as current_file:
        for row in data_rows:
            current_file.write('>' + str(row[0]) + '\n')
            current_file.write(str(row[seq_index]) + '\n')
    return time_stamp, seq_pathname


def resolve_tool_paths(working_prefix):
    make_db_script = os.path.join(working_prefix, 'MakeDb.py')
    if not os.path.isfile(make_db_script):
        raise FileNotFoundError('MakeDb.py was not found in the application directory.')

    define_clones_script = shutil.which('DefineClones.py')
    if not define_clones_script:
        raise FileNotFoundError('DefineClones.py was not found on PATH. Please install the Change-O tools in this environment.')

    return make_db_script, define_clones_script


def run_igblast_fmt7(seq_pathname, output_path, species, working_prefix, igblast_path):
    species_info = species_paths(species)
    workingdir = os.path.join(working_prefix, 'IgBlast')
    cmd = [
        igblast_path,
        '-germline_db_V', species_info['db_v'],
        '-germline_db_J', species_info['db_j'],
        '-germline_db_D', species_info['db_d'],
        '-organism', species_info['organism'],
        '-domain_system', 'imgt',
        '-ig_seqtype', 'Ig',
        '-query', seq_pathname,
        '-auxiliary_data', species_info['aux'],
        '-outfmt', '7 std qseq sseq btop',
        '-out', output_path,
    ]
    subprocess.run(cmd, cwd=workingdir, check=True)


def run_make_db(working_prefix, temp_folder, fmt7_path, fasta_path, species):
    make_db_script, _ = resolve_tool_paths(working_prefix)
    species_info = species_paths(species)
    output_path = os.path.join(temp_folder, os.path.splitext(os.path.basename(fmt7_path))[0] + '_db-pass.tsv')
    cmd = [
        sys.executable,
        make_db_script,
        'igblast',
        '-i', fmt7_path,
        '-s', fasta_path,
        '-r',
        species_info['repo_v'],
        species_info['repo_d'],
        species_info['repo_j'],
        '--extended',
        '--format', 'airr',
        '-o', output_path,
    ]
    subprocess.run(cmd, cwd=temp_folder, check=True)
    return output_path


def run_define_clones(
    working_prefix,
    temp_folder,
    db_path,
    distance=DEFAULT_DEFINE_CLONES_DISTANCE,
    mode=DEFAULT_DEFINE_CLONES_MODE,
    link=DEFAULT_DEFINE_CLONES_LINK,
    model=DEFAULT_DEFINE_CLONES_MODEL,
):
    _, define_clones_script = resolve_tool_paths(working_prefix)
    output_path = os.path.join(temp_folder, os.path.splitext(os.path.basename(db_path))[0] + '_clone-pass.tsv')
    cmd = [
        sys.executable,
        define_clones_script,
        '-d', db_path,
        '--format', 'airr',
        '--act', 'set',
        '--mode', str(mode),
        '--model', str(model),
        '--norm', 'len',
        '--link', str(link),
        '--dist', str(distance),
        '-o', output_path,
    ]
    subprocess.run(cmd, cwd=temp_folder, check=True)
    return output_path


def parse_changeo_output(pathname):
    name_index = 0
    clone_index = 48
    line_index = 0
    err = 0
    err_msg = 'Please double check your input file, make sure it is the output of Change-O DefineClones.py!\n'
    result = []
    clone_ids = []

    with open(pathname, 'r') as readfile:
        for line in readfile:
            line = line.rstrip('\r\n')
            if line_index == 0:
                tmp_list = line.split('\t')
                try:
                    name_index = tmp_list.index('sequence_id')
                except ValueError:
                    err_msg += "Can not find sequence_id information from your file!\n"
                    err = 1
                try:
                    clone_index = tmp_list.index('clone_id')
                except ValueError:
                    err_msg += "Can not find clone_id information from your file!\n"
                    err = 1
                if err == 1:
                    break
            else:
                tmp_list = line.split('\t')
                result.append([tmp_list[name_index], tmp_list[clone_index]])
                clone_ids.append(tmp_list[clone_index])
            line_index += 1

    if err == 1:
        return err, err_msg

    clone_dict = Counter(clone_ids)
    for record in result:
        if clone_dict[record[1]] < 2:
            record[1] = '0'
    return err, result


def apply_changeo_results(db_filename, seq_names, parsed_rows):
    VGenesSQL.UpdateFieldWhereIn(db_filename, 'ClonalPool', '0', 'SeqName', seq_names)
    clone_map = {}
    for seq_name, clone_id in parsed_rows:
        clone_map.setdefault(str(clone_id), []).append(seq_name)
    for clone_id, names in clone_map.items():
        VGenesSQL.UpdateFieldWhereIn(db_filename, 'ClonalPool', str(clone_id), 'SeqName', names)


def build_clone_items(db_filename):
    data_in = VGenesSQL.RunSQL(db_filename, 'SELECT GeneType,ClonalPool FROM vgenesDB WHERE ClonalPool <> "0"')
    clone_dict = {}
    for gene_type, pool_id in data_in:
        clone_name = str(gene_type) + '|' + 'Clone' + str(pool_id)
        clone_dict[clone_name] = clone_dict.get(clone_name, 0) + 1

    clone_items = []
    for key, value in sorted(clone_dict.items(), key=lambda x: x[1], reverse=True):
        clone_items.append(key + '|Num of seq: ' + str(value))
    clone_items.sort(key=lambda x: x[0] if x else '')
    return clone_items


def run_integrated_changeo(
    db_filename,
    data_rows,
    working_prefix,
    temp_folder,
    igblast_path,
    current_record='',
    define_clones_distance=DEFAULT_DEFINE_CLONES_DISTANCE,
    define_clones_mode=DEFAULT_DEFINE_CLONES_MODE,
    define_clones_link=DEFAULT_DEFINE_CLONES_LINK,
    define_clones_model=DEFAULT_DEFINE_CLONES_MODEL,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None
    if not data_rows:
        raise ValueError('No records selected for Change-O clone analysis.')

    species = data_rows[0][3]
    seq_names = [row[0] for row in data_rows]

    emit_progress(5, 'Fetching data ...')
    _, fasta_path = write_fasta(data_rows, temp_folder)

    emit_progress(20, 'Running IgBlast ...')
    fmt7_path = os.path.join(temp_folder, os.path.splitext(os.path.basename(fasta_path))[0] + '.fmt7')
    run_igblast_fmt7(fasta_path, fmt7_path, species, working_prefix, igblast_path)

    emit_progress(40, 'Parsing IgBlast results ...')
    airr_db_path = run_make_db(working_prefix, temp_folder, fmt7_path, fasta_path, species)

    emit_progress(65, 'Running DefineClones ...')
    clone_output_path = run_define_clones(
        working_prefix,
        temp_folder,
        airr_db_path,
        distance=define_clones_distance,
        mode=define_clones_mode,
        link=define_clones_link,
        model=define_clones_model,
    )

    emit_progress(82, 'Importing Change-O results ...')
    err, parsed_rows = parse_changeo_output(clone_output_path)
    if err == 1:
        raise ValueError(parsed_rows)
    apply_changeo_results(db_filename, seq_names, parsed_rows)

    emit_progress(92, 'Collecting clone summary ...')
    clone_items = build_clone_items(db_filename)

    emit_progress(100, 'Change-O clonal calling finished.')
    return {
        'integrated': True,
        'clone_items': clone_items,
        'current_record': current_record,
        'output_path': clone_output_path,
        'summary_message': 'Successfully identified clones with Change-O!',
    }
