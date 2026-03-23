import math
import os
import shutil
import traceback

import pandas as pd
from PyQt5.QtCore import QObject, QRunnable, Qt, pyqtSignal, pyqtSlot

import VGenesCloneCaller
import VGenesSQL
from changeo_runner import run_igblast_fmt7, run_integrated_changeo, write_fasta


class TaskSignals(QObject):
    started = pyqtSignal(str)
    progress = pyqtSignal(int, str)
    result = pyqtSignal(object)
    error = pyqtSignal(str, str)
    finished = pyqtSignal()


class FunctionTask(QRunnable):
    def __init__(self, fn, *args, task_name="Task", **kwargs):
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.task_name = task_name
        self.signals = TaskSignals()

    @pyqtSlot()
    def run(self):
        try:
            self.signals.started.emit(self.task_name)
            result = self.fn(*self.args, emit_progress=self.signals.progress.emit, **self.kwargs)
        except Exception as exc:
            self.signals.error.emit(str(exc), traceback.format_exc())
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()


def build_seq_table_order_clause(sort_field, sort_order, default_sort_fields):
    order_fields = []
    if sort_field:
        direction = "ASC" if sort_order == Qt.AscendingOrder else "DESC"
        order_fields.append(f"{sort_field} {direction}")

    for field in default_sort_fields:
        if field in ("", None, "None"):
            continue
        if any(entry.split()[0] == field for entry in order_fields):
            continue
        order_fields.append(field)

    if not order_fields:
        order_fields.append("SeqName")
    return " ORDER BY " + ", ".join(order_fields)


def fetch_seq_table_page(
    db_filename,
    page_size,
    current_page_index,
    sort_field,
    sort_order,
    default_sort_fields,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None

    if not db_filename or db_filename == "none":
        return {"status": "no_db"}

    emit_progress(5, "Counting records ...")
    total_records = VGenesSQL.RunSQL(db_filename, "SELECT COUNT(*) FROM vgenesDB")[0][0]
    total_pages = max(1, math.ceil(total_records / page_size))

    if total_records == 0:
        return {
            "status": "empty",
            "total_records": 0,
            "total_pages": total_pages,
            "page_size": page_size,
            "current_page": 1,
        }

    emit_progress(15, "Loading display fields ...")
    header_rows = VGenesSQL.RunSQL(
        db_filename,
        'SELECT Field,FieldNickName FROM fieldsname WHERE display = "yes" ORDER BY display_priority,ID LIMIT 0,200',
    )
    current_field_list = [row[0] for row in header_rows]
    current_nickname_list = [row[1] for row in header_rows]

    fields = ",".join(current_field_list)
    if fields == "":
        fields = "*"

    emit_progress(40, "Loading page records ...")
    record_limit_statement = f" LIMIT {current_page_index * page_size},{page_size}"
    sql_statement = "select " + fields + " from vgenesDB"
    sql_statement += build_seq_table_order_clause(sort_field, sort_order, default_sort_fields)
    sql_statement += record_limit_statement
    data_rows = VGenesSQL.RunSQL(db_filename, sql_statement)

    emit_progress(90, "Preparing table view ...")
    return {
        "status": "ok",
        "total_records": total_records,
        "total_pages": total_pages,
        "page_size": page_size,
        "current_page": current_page_index + 1,
        "field_list": current_field_list,
        "nickname_list": current_nickname_list,
        "rows": data_rows,
        "row_names": [str(row[0]) for row in data_rows],
    }


class _ProgressEmitter:
    def __init__(self, emit_progress, start_pct, end_pct):
        self.emit_progress = emit_progress
        self.start_pct = start_pct
        self.end_pct = end_pct

    def emit(self, pct, label):
        span = max(1, self.end_pct - self.start_pct)
        mapped_pct = self.start_pct + int((pct / 100.0) * span)
        self.emit_progress(mapped_pct, label)


def _read_seq_fields(db_filename, seq_name, fields):
    sql = "SELECT " + ",".join(fields) + " FROM vgenesDB WHERE SeqName = ?"
    rows = VGenesSQL.run_sql_query(db_filename, sql, [seq_name])
    if len(rows) == 0:
        return tuple("" for _ in fields)
    return rows[0]


def run_conventional_clone_calling(
    db_filename,
    clonal_pools,
    duplicates,
    remove,
    total_sequences,
    error_log_text,
    error_count,
    pool_names,
    current_record,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None

    emit_progress(2, "Preparing clonal calling ...")
    if len(pool_names) == 0:
        pool_names = ["Current Clone pool"]

    cp_seqs = 0
    cp_count = 0

    existing_clone_rows = VGenesSQL.RunSQL(db_filename, 'SELECT DISTINCT(ClonalPool) FROM vgenesDB')
    existing_clone_list = [row[0] for row in existing_clone_rows]
    clone_id = 1
    while str(clone_id) in existing_clone_list:
        clone_id += 1

    total_pools = max(1, len(clonal_pools))
    for pool_index, pool in enumerate(clonal_pools):
        pool_start = 5 + int((pool_index / total_pools) * 55)
        pool_end = 5 + int(((pool_index + 1) / total_pools) * 55)
        progress_adapter = _ProgressEmitter(emit_progress, pool_start, pool_end)

        cp_list = VGenesCloneCaller.CloneCaller(list(pool), duplicates, progress_adapter, pool_names[pool_index])

        for record in cp_list:
            cp_count += 1
            duplicate_count = 1
            duplicate_label = "Sequences identical: "
            first_seq_name = None

            for seq_name in record:
                if not duplicates:
                    VGenesSQL.UpdateFieldbySeqName(seq_name, str(clone_id), 'ClonalPool', db_filename)
                    existing_clone_list.append(str(clone_id))
                else:
                    if duplicate_count == 1:
                        first_seq_name = seq_name
                        duplicate_seq_label = 'Duplicate of:  ' + seq_name
                    else:
                        if remove:
                            VGenesSQL.UpdateFieldbySeqName(seq_name, 'Duplicate', 'Quality', db_filename)
                            VGenesSQL.UpdateFieldbySeqName(seq_name, 'Delete', 'Project', db_filename)
                        else:
                            VGenesSQL.UpdateFieldbySeqName(seq_name, duplicate_seq_label, 'Quality', db_filename)
                        duplicate_label += seq_name + ', '
                    duplicate_count += 1
                cp_seqs += 1

            depth = 'Depth = ' + str(duplicate_count - 1)
            if duplicates and first_seq_name:
                existing_comments, existing_quality = _read_seq_fields(db_filename, first_seq_name, ['Comments', 'Quality'])
                if duplicate_label.endswith(', '):
                    duplicate_label = duplicate_label[:-2]
                if existing_comments not in ('', ' ', 'Comments', None):
                    duplicate_label = duplicate_label + ', ' + str(existing_comments)
                if existing_quality not in ('', ' ', 'Quality', None):
                    depth = depth + '  ' + str(existing_quality)
                VGenesSQL.UpdateFieldbySeqName(first_seq_name, duplicate_label, 'Comments', db_filename)
                VGenesSQL.UpdateFieldbySeqName(first_seq_name, depth, 'Quality', db_filename)

            while str(clone_id) in existing_clone_list:
                clone_id += 1

    emit_progress(70, 'Integrating cell barcode info ...')
    sql = 'SELECT GeneType,ClonalPool,Blank10 FROM vgenesDB WHERE ClonalPool != 0'
    clone_df = pd.DataFrame(VGenesSQL.RunSQL(db_filename, sql), columns=['GeneType', 'ClonalPool', 'CellBarcode'])
    if not clone_df.empty:
        has_barcode = ~clone_df['CellBarcode'].isin(['Blank10', ''])
        data_with_barcode = clone_df[has_barcode].copy()
        grouped = data_with_barcode.groupby('CellBarcode')
        for barcode, group in grouped:
            hc_pools = []
            kappa_pools = []
            lambda_pools = []

            for _, row in group.iterrows():
                gene_type = row['GeneType']
                pool = row['ClonalPool']
                if gene_type == 'Heavy':
                    hc_pools.append(f'H{pool}')
                elif gene_type == 'Kappa':
                    kappa_pools.append(f'K{pool}')
                elif gene_type == 'Lambda':
                    lambda_pools.append(f'L{pool}')

            lc_pools = kappa_pools + lambda_pools
            pairs = ",".join([f"{hc}_{lc}" for hc in hc_pools for lc in lc_pools]) if hc_pools and lc_pools else "0"
            VGenesSQL.UpdateFieldWhere(db_filename, 'ClonalRank', str(pairs), 'Blank10', barcode)

    emit_progress(85, 'Collecting clone summary ...')
    data_in = VGenesSQL.RunSQL(db_filename, 'SELECT GeneType,ClonalPool FROM vgenesDB WHERE ClonalPool <> "0"')
    clone_dict = {}
    for gene_type, pool_id in data_in:
        clone_name = str(gene_type) + '|' + 'Clone' + str(pool_id)
        clone_dict[clone_name] = clone_dict.get(clone_name, 0) + 1

    clone_items = []
    for key, value in sorted(clone_dict.items(), key=lambda x: x[1], reverse=True):
        clone_items.append(key + '|Num of seq: ' + str(value))
    clone_items.sort(key=lambda x: x[0] if x else '')

    error_file = ""
    summary = str(cp_count) + ' clonal pools containing ' + str(cp_seqs) + ' sequences were identified from ' + str(total_sequences) + ' total sequences analyzed.\n'
    if len(summary) > 0:
        combined_log = summary + 'The following ' + str(error_count) + ' sequences could not be anaylzed for\nclonality because no CDR3s are indicated:\n' + error_log_text
        error_file = os.path.join(os.path.dirname(db_filename), 'Temp', 'ErLog.txt')
        try:
            os.makedirs(os.path.dirname(error_file), exist_ok=True)
            with open(error_file, 'w') as current_file:
                current_file.write(combined_log)
        except Exception:
            error_file = ""

    emit_progress(100, 'Clonal calling finished.')
    return {
        'clone_items': clone_items,
        'summary_message': 'Integtated cell barcode info\nSuccessfully identified clones!',
        'error_file': error_file,
        'current_record': current_record,
        'remove_duplicates': remove,
    }


def run_changeo_igblast_export(
    data_rows,
    output_path,
    working_prefix,
    temp_folder,
    igblast_path,
    preset_id=None,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None
    if not data_rows:
        raise ValueError('No records selected for Change-O IgBlast export.')

    emit_progress(15, 'Fetching data ...')
    _, seq_pathname = write_fasta(data_rows, temp_folder)

    emit_progress(45, 'Running IgBlast ...')
    from igblast_presets import resolve_preset

    species = data_rows[0][3]
    preset = resolve_preset(preset_id, species, "IG")
    run_igblast_fmt7(seq_pathname, output_path, preset, working_prefix, igblast_path)

    emit_progress(100, 'IgBlast export finished.')
    return {
        'summary_message': 'Successfully finished IgBlast, the results have been saved!',
        'output_path': output_path,
    }


def run_changeo_clone_pipeline(
    db_filename,
    data_rows,
    working_prefix,
    temp_folder,
    igblast_path,
    current_record='',
    preset_id=None,
    define_clones_distance='0.15',
    define_clones_mode='gene',
    define_clones_link='single',
    define_clones_model='ham',
    emit_progress=None,
):
    return run_integrated_changeo(
        db_filename,
        data_rows,
        working_prefix,
        temp_folder,
        igblast_path,
        current_record=current_record,
        preset_id=preset_id,
        define_clones_distance=define_clones_distance,
        define_clones_mode=define_clones_mode,
        define_clones_link=define_clones_link,
        define_clones_model=define_clones_model,
        emit_progress=emit_progress,
    )


def run_hclc_pairing(
    db_filename,
    checked_records,
    temp_folder,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None

    if len(checked_records) == 0:
        raise ValueError('Please check some sequences to start!')

    progress = 0
    new_checks = []
    err_msg_type1 = ''
    err_msg_type2 = ''
    err_msg_type3 = ''
    total = len(checked_records)

    for item in checked_records:
        safe_item = str(item).replace('"', '""')
        rows = VGenesSQL.RunSQL(
            db_filename,
            'SELECT SeqName,Blank10 FROM vgenesDB WHERE SeqName = "' + safe_item + '"',
        )
        if len(rows) == 0:
            progress += 1
            emit_progress(int(progress / total * 100), 'Matching HC/LC pairs ...')
            continue

        barcode = rows[0][1]
        if barcode == 'Blank10' or barcode == '':
            err_msg_type1 += 'Sequence ' + item + ' does not have barcode information!\n'
        else:
            safe_barcode = str(barcode).replace('"', '""')
            pair_rows = VGenesSQL.RunSQL(
                db_filename,
                'SELECT SeqName,Blank10 FROM vgenesDB WHERE Blank10 = "' + safe_barcode + '"',
            )
            pair_count = len(pair_rows) - 1
            if pair_count == 0:
                err_msg_type2 += 'For ' + item + ', did not find any Heavy/Light chain using same barcode!\n'
            else:
                err_msg_type3 += 'For ' + item + ', find ' + str(pair_count) + ' Heavy/Light chain using same barcode!\n'
                for record in pair_rows:
                    seq_name = record[0]
                    if seq_name not in checked_records and seq_name not in new_checks:
                        new_checks.append(seq_name)

        progress += 1
        emit_progress(int(progress / total * 100), 'Matching HC/LC pairs ...')

    error_file = os.path.join(temp_folder, 'ErLog.txt')
    with open(error_file, 'w') as current_file:
        current_file.write('Running finished!\n')
        current_file.write('\nThe following records have paired HC/LC:\n')
        current_file.write(err_msg_type3)
        current_file.write('\nThe following records do not have barcode information:\n')
        current_file.write(err_msg_type1)
        current_file.write('\nThe following records do not have any paired HC/LC:\n')
        current_file.write(err_msg_type2)

    return {
        'checked_records': list(checked_records),
        'new_checks': new_checks,
        'error_file': error_file,
    }


def run_hclc_export_db(
    db_filename,
    output_path,
    checked_records,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None

    emit_progress(2, 'Creating paired HC/LC database ...')
    try:
        shutil.copy(db_filename, output_path)
    except Exception:
        return {
            'sign': 1,
            'message': 'Can not save file in this path! You do not have write permission!',
            'pathname': output_path,
        }

    if len(checked_records) == 0:
        sql_statement = 'SELECT DISTINCT(Blank10) FROM vgenesDB WHERE GeneType = "Heavy"'
        data_in = VGenesSQL.RunSQL(db_filename, sql_statement)
    else:
        list_str = '("' + '","'.join(checked_records) + '")'
        sql_statement = 'SELECT DISTINCT(Blank10) FROM vgenesDB WHERE SeqName IN ' + list_str
        data_in = VGenesSQL.RunSQL(db_filename, sql_statement)

    if len(data_in) < 2:
        return {
            'sign': 1,
            'message': 'Your VGene DB do not have any barcode information!',
            'pathname': output_path,
        }

    good_list = []
    total = len(data_in)
    progress = 0
    for record in data_in:
        barcode = record[0]
        if barcode not in ('Blank10', ''):
            sql_statement1 = (
                'SELECT SeqName FROM vgenesDB WHERE Blank10 = "' + str(barcode).replace('"', '""') +
                '" AND GeneType IN ("Heavy","Beta","Delta")'
            )
            data_in1 = VGenesSQL.RunSQL(db_filename, sql_statement1)
            sql_statement2 = (
                'SELECT SeqName FROM vgenesDB WHERE Blank10 = "' + str(barcode).replace('"', '""') +
                '" AND GeneType NOT IN ("Heavy","Beta","Delta")'
            )
            data_in2 = VGenesSQL.RunSQL(db_filename, sql_statement2)
            if len(data_in1) == 1 and len(data_in2) == 1:
                good_list.append(barcode)

        progress += 1
        emit_progress(int(progress / total * 100), 'Filtering complete HC/LC pairs ...')

    seq_num = len(good_list)
    if seq_num > 0:
        list_str = '("' + '","'.join(good_list) + '")'
        sql_statement = 'DELETE FROM vgenesDB WHERE Blank10 NOT IN ' + list_str
        VGenesSQL.RunUpdateSQL(output_path, sql_statement)
        return {
            'sign': 0,
            'message': 'Total ' + str(seq_num) + ' HC/LC pairs were found!',
            'pathname': output_path,
        }

    return {
        'sign': 1,
        'message': 'Did not find any HC/LC pair in your current DB!',
        'pathname': output_path,
    }


def run_copy_records_db(
    db_filename,
    output_path,
    checked_records,
    emit_progress=None,
):
    if emit_progress is None:
        emit_progress = lambda pct, label: None

    emit_progress(5, 'Creating database copy ...')
    try:
        shutil.copy(db_filename, output_path)
    except Exception:
        return {
            'sign': 1,
            'message': 'Can not save file in this path! You do not have write permission!',
            'pathname': output_path,
        }

    list_str = '("' + '","'.join(checked_records) + '")'
    sql_statement = 'DELETE FROM vgenesDB WHERE SeqName NOT IN ' + list_str
    emit_progress(60, 'Removing unchecked records ...')
    VGenesSQL.RunUpdateSQL(output_path, sql_statement)

    emit_progress(100, 'Selected-record database finished.')
    return {
        'sign': 0,
        'message': 'New DB with selected records have been created!',
        'pathname': output_path,
    }
