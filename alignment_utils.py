from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices


def load_substitution_matrix(name):
    matrix_name = str(name or "").strip().upper()
    if matrix_name == "PAM60":
        return substitution_matrices.load("PAM60")
    if matrix_name == "BENNER22":
        return substitution_matrices.load("BENNER22")
    return substitution_matrices.load("BLOSUM62")


def global_alignment_strings(seq_a, seq_b, matrix_name, gap_open, gap_extend):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = load_substitution_matrix(matrix_name)
    aligner.open_gap_score = float(gap_open)
    aligner.extend_gap_score = float(gap_extend)

    alignment = aligner.align(seq_a, seq_b)[0]
    return str(alignment[0]), str(alignment[1])
