# Set of Global Variables for editing things:

# 2 variables for user experience
global TASK_TITLE, OPEN_XLSX_AUTO

# 8 variables for changing results
global MAT, MIS, GAP_OPEN, GAP_EXTEND, PAM_RANGE_MAX, ERR_MAX, ERR_PADDING, PHRED_MEANINGFUL_MIN

TASK_TITLE = ""
OPEN_XLSX_AUTO = False
MAT = 2
MIS = -1
GAP_OPEN = -50
GAP_EXTEND = -4
PAM_RANGE_MAX = 5
ERR_MAX = 0.03
ERR_PADDING = 1
PHRED_MEANINGFUL_MIN = 30


def get_align_matrix_for_subsequence_positioning():
    # 'X' = for both end of main sequence, meaning the subsequence must be between the sequence.
    matrix = {
        ('A', 'A'): MAT, ('A', 'T'): MIS, ('A', 'G'): MIS, ('A', 'C'): MIS, ('A', '-'): MIS, ('A', 'X'): -1000,
        ('T', 'A'): MIS, ('T', 'T'): MAT, ('T', 'G'): MIS, ('T', 'C'): MIS, ('T', '-'): MIS, ('T', 'X'): -1000,
        ('G', 'A'): MIS, ('G', 'T'): MIS, ('G', 'G'): MAT, ('G', 'C'): MIS, ('G', '-'): MIS, ('G', 'X'): -1000,
        ('C', 'A'): MIS, ('C', 'T'): MIS, ('C', 'G'): MIS, ('C', 'C'): MAT, ('C', '-'): MIS, ('C', 'X'): -1000,
        ('-', 'A'): MIS, ('-', 'T'): MIS, ('-', 'G'): MIS, ('-', 'C'): MIS, ('-', '-'): MAT, ('-', 'X'): -1000,
        ('X', 'A'): -1000, ('X', 'T'): -1000, ('X', 'G'): -1000, ('X', 'C'): -1000, ('X', '-'): -1000, ('X', 'X'): MAT,
    }
    return matrix
