import json
import os
import re
import shutil
import subprocess
from copy import deepcopy


CONFIG_FILENAME = "igblast_reference_presets.json"


def _working_prefix():
    return os.path.dirname(os.path.realpath(__file__))


def config_path():
    return os.path.join(_working_prefix(), "Conf", CONFIG_FILENAME)


def custom_reference_root():
    return os.path.join(_working_prefix(), "IgBlast", "Custom")


def makeblastdb_path():
    bundled = os.path.join(_working_prefix(), "Tools", "makeblastdb")
    if os.path.isfile(bundled):
        return bundled
    if os.path.isfile(bundled + ".exe"):
        return bundled + ".exe"
    found = shutil.which("makeblastdb")
    if found:
        return found
    raise FileNotFoundError("makeblastdb was not found in Tools/ or on PATH.")


def _preset(
    preset_id,
    name,
    organism_label,
    sequence_type,
    species_label,
    domain_system,
    ig_seqtype,
    root_folder,
    v_nt,
    d_nt,
    j_nt,
    aux_file,
    v_fasta,
    d_fasta,
    j_fasta,
    enabled=True,
    default=False,
    builtin=True,
):
    return {
        "id": preset_id,
        "name": name,
        "organism_label": organism_label,
        "sequence_type": sequence_type,
        "species_label": species_label,
        "domain_system": domain_system,
        "ig_seqtype": ig_seqtype,
        "root_folder": root_folder,
        "db_v": v_nt,
        "db_d": d_nt,
        "db_j": j_nt,
        "repo_v": v_fasta,
        "repo_d": d_fasta,
        "repo_j": j_fasta,
        "auxiliary_data": aux_file,
        "enabled": enabled,
        "default": default,
        "builtin": builtin,
    }


def built_in_presets():
    root = _working_prefix()
    ig_root = os.path.join(root, "IgBlast")
    return [
        _preset(
            "human_ig",
            "Human IG",
            "human",
            "IG",
            "Human",
            "kabat",
            "Ig",
            os.path.join(ig_root, "IG", "Human"),
            os.path.join(ig_root, "IG", "Human", "HumanVGenes.nt"),
            os.path.join(ig_root, "IG", "Human", "HumanDGenes.nt"),
            os.path.join(ig_root, "IG", "Human", "HumanJGenes.nt"),
            os.path.join(ig_root, "optional_file", "human_gl.aux"),
            os.path.join(ig_root, "IG", "Human", "HumanVGenes.fasta"),
            os.path.join(ig_root, "IG", "Human", "HumanDGenes.fasta"),
            os.path.join(ig_root, "IG", "Human", "HumanJGenes.fasta"),
            default=True,
        ),
        _preset(
            "mouse_ig",
            "Mouse IG",
            "mouse",
            "IG",
            "Mouse",
            "kabat",
            "Ig",
            os.path.join(ig_root, "IG", "Mouse"),
            os.path.join(ig_root, "IG", "Mouse", "MouseVGenes.nt"),
            os.path.join(ig_root, "IG", "Mouse", "MouseDGenes.nt"),
            os.path.join(ig_root, "IG", "Mouse", "MouseJGenes.nt"),
            os.path.join(ig_root, "optional_file", "mouse_gl.aux"),
            os.path.join(ig_root, "IG", "Mouse", "MouseVGenes.fasta"),
            os.path.join(ig_root, "IG", "Mouse", "MouseDGenes.fasta"),
            os.path.join(ig_root, "IG", "Mouse", "MouseJGenes.fasta"),
        ),
        _preset(
            "human_tr",
            "Human TR",
            "human",
            "TR",
            "Human",
            "imgt",
            "TCR",
            os.path.join(ig_root, "TR", "Human"),
            os.path.join(ig_root, "TR", "Human", "HumanVGenes.nt"),
            os.path.join(ig_root, "TR", "Human", "HumanDGenes.nt"),
            os.path.join(ig_root, "TR", "Human", "HumanJGenes.nt"),
            os.path.join(ig_root, "optional_file", "human_gl.aux"),
            os.path.join(ig_root, "TR", "Human", "HumanVGenes.fasta"),
            os.path.join(ig_root, "TR", "Human", "HumanDGenes.fasta"),
            os.path.join(ig_root, "TR", "Human", "HumanJGenes.fasta"),
        ),
        _preset(
            "mouse_tr",
            "Mouse TR",
            "mouse",
            "TR",
            "Mouse",
            "imgt",
            "TCR",
            os.path.join(ig_root, "TR", "Mouse"),
            os.path.join(ig_root, "TR", "Mouse", "MouseVGenes.nt"),
            os.path.join(ig_root, "TR", "Mouse", "MouseDGenes.nt"),
            os.path.join(ig_root, "TR", "Mouse", "MouseJGenes.nt"),
            os.path.join(ig_root, "optional_file", "mouse_gl.aux"),
            os.path.join(ig_root, "TR", "Mouse", "MouseVGenes.fasta"),
            os.path.join(ig_root, "TR", "Mouse", "MouseDGenes.fasta"),
            os.path.join(ig_root, "TR", "Mouse", "MouseJGenes.fasta"),
        ),
    ]


def _sanitize_preset_paths(preset):
    preset = deepcopy(preset)
    for key in ("db_v", "db_d", "db_j"):
        value = str(preset.get(key, "")).strip()
        if not value:
            continue
        if value.lower().endswith(".fasta"):
            candidate = value[:-6] + ".nt"
            if os.path.isfile(candidate) or os.path.isfile(candidate + ".nhr"):
                preset[key] = candidate
    return preset


def _seed_config():
    os.makedirs(os.path.dirname(config_path()), exist_ok=True)
    data = {
        "version": 1,
        "default_preset_id": "human_ig",
        "presets": built_in_presets(),
    }
    save_config(data)
    return data


def load_config():
    path = config_path()
    if not os.path.isfile(path):
        return _seed_config()
    try:
        with open(path, "r") as handle:
            data = json.load(handle)
    except Exception:
        return _seed_config()

    builtin_map = {preset["id"]: preset for preset in built_in_presets()}
    presets = {preset.get("id"): preset for preset in data.get("presets", []) if preset.get("id")}
    for preset_id, preset in builtin_map.items():
        merged = deepcopy(preset)
        if preset_id in presets:
            merged.update(presets[preset_id])
            merged["builtin"] = True
        presets[preset_id] = merged

    ordered = []
    seen = set()
    for preset in built_in_presets():
        ordered.append(_sanitize_preset_paths(presets[preset["id"]]))
        seen.add(preset["id"])
    for preset in data.get("presets", []):
        preset_id = preset.get("id")
        if preset_id and preset_id not in seen:
            ordered.append(_sanitize_preset_paths(preset))
            seen.add(preset_id)

    data["presets"] = ordered
    if data.get("default_preset_id") not in seen:
        data["default_preset_id"] = "human_ig"
    return data


def save_config(data):
    os.makedirs(os.path.dirname(config_path()), exist_ok=True)
    with open(config_path(), "w") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def list_presets(sequence_type=None, enabled_only=False):
    presets = load_config().get("presets", [])
    results = []
    for preset in presets:
        if sequence_type and preset.get("sequence_type") != sequence_type:
            continue
        if enabled_only and not preset.get("enabled", True):
            continue
        results.append(preset)
    return results


def get_preset(preset_id):
    for preset in load_config().get("presets", []):
        if preset.get("id") == preset_id:
            return preset
    return None


def get_default_preset_id(sequence_type=None):
    data = load_config()
    default_preset = get_preset(data.get("default_preset_id"))
    if default_preset and (sequence_type is None or default_preset.get("sequence_type") == sequence_type):
        return default_preset["id"]
    presets = list_presets(sequence_type=sequence_type, enabled_only=True)
    if presets:
        return presets[0]["id"]
    return None


def normalize_preset_id(name):
    raw = "".join(ch.lower() if ch.isalnum() else "_" for ch in str(name or "").strip())
    raw = "_".join(part for part in raw.split("_") if part)
    return raw or "custom_preset"


def upsert_preset(preset):
    data = load_config()
    presets = data.get("presets", [])
    preset = deepcopy(preset)
    preset_id = preset.get("id") or normalize_preset_id(preset.get("name"))
    preset["id"] = preset_id
    replaced = False
    for index, current in enumerate(presets):
        if current.get("id") == preset_id:
            if current.get("builtin") and not preset.get("builtin", current.get("builtin")):
                preset["builtin"] = True
            presets[index] = preset
            replaced = True
            break
    if not replaced:
        presets.append(preset)
    data["presets"] = presets
    if preset.get("default"):
        data["default_preset_id"] = preset_id
        for current in data["presets"]:
            current["default"] = current.get("id") == preset_id
    save_config(data)
    return preset_id


def delete_preset(preset_id):
    data = load_config()
    new_presets = []
    deleted = False
    for preset in data.get("presets", []):
        if preset.get("id") != preset_id:
            new_presets.append(preset)
            continue
        if preset.get("builtin"):
            raise ValueError("Built-in presets cannot be deleted.")
        deleted = True
    data["presets"] = new_presets
    if deleted and data.get("default_preset_id") == preset_id:
        data["default_preset_id"] = get_default_preset_id()
    save_config(data)
    return deleted


def set_default_preset(preset_id):
    data = load_config()
    found = False
    for preset in data.get("presets", []):
        is_default = preset.get("id") == preset_id
        preset["default"] = is_default
        if is_default:
            found = True
    if not found:
        raise ValueError("Preset was not found: " + str(preset_id))
    data["default_preset_id"] = preset_id
    save_config(data)


def detect_reference_files(root_folder):
    detected = {}
    if not root_folder:
        return detected
    base = os.path.abspath(root_folder)
    if not os.path.isdir(base):
        return detected
    names = os.listdir(base)
    for kind in ("V", "D", "J"):
        nt_match = None
        fasta_match = None
        for entry in names:
            lower = entry.lower()
            if kind.lower() not in lower:
                continue
            if lower.endswith(".nt") and nt_match is None:
                nt_match = os.path.join(base, entry)
            if lower.endswith(".fasta") and fasta_match is None:
                fasta_match = os.path.join(base, entry)
        if nt_match:
            detected[f"db_{kind.lower()}"] = nt_match
        if fasta_match:
            detected[f"repo_{kind.lower()}"] = fasta_match
    return detected


def validate_preset(preset, igblast_executable=None):
    errors = []
    if igblast_executable:
        if not os.path.isfile(igblast_executable) and shutil.which(igblast_executable) is None:
            errors.append("IgBlast executable was not found: " + str(igblast_executable))
    required = [
        ("name", "Preset name"),
        ("sequence_type", "Sequence type"),
        ("organism_label", "Organism label"),
        ("domain_system", "Domain system"),
        ("db_v", "V reference"),
        ("db_j", "J reference"),
        ("auxiliary_data", "Auxiliary data"),
    ]
    for key, label in required:
        if not str(preset.get(key, "")).strip():
            errors.append(label + " is required.")
    for key in ("db_v", "db_j", "db_d", "repo_v", "repo_d", "repo_j", "auxiliary_data"):
        value = str(preset.get(key, "")).strip()
        if not value:
            continue
        if key.startswith("db_"):
            if value.lower().endswith(".fasta"):
                errors.append("IgBlast DB path must point to the generated BLAST database base name, not a FASTA file: " + value)
                continue
            index_candidates = [value + suffix for suffix in (".nhr", ".nin", ".nsq")]
            if os.path.isfile(value):
                if not all(os.path.isfile(path) for path in index_candidates):
                    errors.append("IgBlast DB index files are missing for " + value)
            elif not all(os.path.isfile(path) for path in index_candidates):
                errors.append("IgBlast DB index files are missing for " + value)
            continue
        if not os.path.isfile(value):
            errors.append("Missing file for " + key + ": " + value)
    return errors


def legacy_preset_id(species_label, sequence_type):
    species = str(species_label or "").strip().lower()
    seq_type = "TR" if str(sequence_type or "").upper() == "TR" else "IG"
    if species == "mouse":
        return "mouse_tr" if seq_type == "TR" else "mouse_ig"
    return "human_tr" if seq_type == "TR" else "human_ig"


def resolve_preset(preset_id=None, species_label=None, sequence_type="IG"):
    if preset_id:
        preset = get_preset(preset_id)
        if preset:
            return preset
    fallback_id = legacy_preset_id(species_label, sequence_type)
    preset = get_preset(fallback_id)
    if preset:
        return preset
    default_id = get_default_preset_id(sequence_type=sequence_type)
    preset = get_preset(default_id)
    if preset:
        return preset
    raise ValueError("No IgBlast preset is available for sequence type " + str(sequence_type))


def build_igblast_command(
    igblast_executable,
    preset,
    query_path,
    outfmt,
    output_path=None,
    show_translation=True,
    extra_args=None,
):
    def _db_base(path_value):
        path_value = os.path.abspath(path_value)
        if path_value.lower().endswith(".fasta"):
            candidate = path_value[:-6] + ".nt"
            if os.path.isfile(candidate) or os.path.isfile(candidate + ".nhr"):
                return candidate
        return path_value

    cmd = [
        igblast_executable,
        "-germline_db_V", _db_base(preset["db_v"]),
        "-germline_db_J", _db_base(preset["db_j"]),
        "-organism", str(preset["organism_label"]).lower(),
        "-domain_system", preset["domain_system"],
        "-query", query_path,
        "-auxiliary_data", os.path.abspath(preset["auxiliary_data"]),
        "-outfmt", str(outfmt),
    ]
    db_d = str(preset.get("db_d", "")).strip()
    if db_d:
        cmd.extend(["-germline_db_D", _db_base(db_d)])
    if show_translation:
        cmd.append("-show_translation")
    ig_seqtype = str(preset.get("ig_seqtype", "")).strip()
    if ig_seqtype:
        cmd.extend(["-ig_seqtype", ig_seqtype])
    if output_path:
        cmd.extend(["-out", output_path])
    if extra_args:
        cmd.extend(list(extra_args))
    return cmd


def build_make_db_repo_args(preset):
    repos = []
    for key in ("repo_v", "repo_d", "repo_j"):
        value = str(preset.get(key, "")).strip()
        if value:
            repos.append(os.path.abspath(value))
    return repos


def _copy_file(source_path, dest_path):
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    shutil.copy2(source_path, dest_path)
    return dest_path


def _normalize_fasta(source_path, dest_path):
    allowed = set("ACGTURYSWKMBDHVN-")

    def _safe_header(raw_header, index):
        raw_header = raw_header[1:].strip()
        parts = [part.strip() for part in raw_header.split("|")]
        if len(parts) > 1 and parts[1]:
            candidate = parts[1]
        elif parts and parts[0]:
            candidate = parts[0]
        else:
            candidate = "seq_" + str(index)
        candidate = re.sub(r"[^A-Za-z0-9_.*\-]+", "_", candidate).strip("_")
        return ">" + (candidate or ("seq_" + str(index)))

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(source_path, "r") as infile, open(dest_path, "w") as outfile:
        header = None
        seq_parts = []
        record_index = 0
        for raw_line in infile:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).upper()
                    seq = "".join(ch for ch in seq if ch in allowed)
                    if seq:
                        outfile.write(header + "\n")
                        outfile.write(seq + "\n")
                record_index += 1
                header = _safe_header(line, record_index)
                seq_parts = []
            else:
                cleaned = re.sub(r"\s+", "", line).upper()
                seq_parts.append(cleaned)
        if header is not None:
            seq = "".join(seq_parts).upper()
            seq = "".join(ch for ch in seq if ch in allowed)
            if seq:
                outfile.write(header + "\n")
                outfile.write(seq + "\n")
    if os.path.getsize(dest_path) == 0:
        raise ValueError("No valid FASTA sequences were found after normalization: " + source_path)


def _build_db_from_fasta(fasta_path, output_base):
    cmd = [
        makeblastdb_path(),
        "-parse_seqids",
        "-dbtype", "nucl",
        "-in", fasta_path,
        "-out", output_base,
    ]
    result = subprocess.run(
        cmd,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or ("makeblastdb failed with exit code " + str(result.returncode))
        raise RuntimeError(detail)


def build_custom_preset_assets(preset):
    preset = deepcopy(preset)
    preset_id = preset.get("id") or normalize_preset_id(preset.get("name"))
    preset["id"] = preset_id

    if not str(preset.get("repo_v", "")).strip():
        raise ValueError("V FASTA is required for a custom preset.")
    if not str(preset.get("repo_j", "")).strip():
        raise ValueError("J FASTA is required for a custom preset.")
    if not str(preset.get("auxiliary_data", "")).strip():
        raise ValueError("Auxiliary data file is required for a custom preset.")

    root_folder = os.path.join(custom_reference_root(), preset.get("sequence_type", "IG"), preset_id)
    os.makedirs(root_folder, exist_ok=True)

    v_fasta = os.path.join(root_folder, "VGenes.fasta")
    _normalize_fasta(preset["repo_v"], v_fasta)
    j_fasta = os.path.join(root_folder, "JGenes.fasta")
    _normalize_fasta(preset["repo_j"], j_fasta)
    d_fasta = ""
    if str(preset.get("repo_d", "")).strip():
        d_fasta = os.path.join(root_folder, "DGenes.fasta")
        _normalize_fasta(preset["repo_d"], d_fasta)
    aux_path = _copy_file(preset["auxiliary_data"], os.path.join(root_folder, preset_id + ".aux"))

    db_v = os.path.join(root_folder, "VGenes.nt")
    db_j = os.path.join(root_folder, "JGenes.nt")
    db_d = os.path.join(root_folder, "DGenes.nt") if d_fasta else ""

    _build_db_from_fasta(v_fasta, db_v)
    _build_db_from_fasta(j_fasta, db_j)
    if d_fasta:
        _build_db_from_fasta(d_fasta, db_d)

    preset["root_folder"] = root_folder
    preset["repo_v"] = v_fasta
    preset["repo_j"] = j_fasta
    preset["repo_d"] = d_fasta
    preset["db_v"] = db_v
    preset["db_j"] = db_j
    preset["db_d"] = db_d
    preset["auxiliary_data"] = aux_path
    return preset
