import os
import datetime
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO
from InDel_Counter import PAM_MAX, ERR_MAX, MATCH_LETTER, get_result_text_set
from InDel_Counter import RESULT_LOG_ADDRESS, ALIGN_LOG_ADDRESS, DATA_ADDRESS, GUIDE_RNA_ADDRESS, \
    GUIDE_RNA_SET_ADDRESS, REF_SET_ADDRESS, TASK_TITLE
from InDel_Counter import read_g_rna, get_indel_count_result_list, find_final_position, get_indel_shape_text
from InDel_Counter import write_result_sub_log, write_result_csv

from log_writer import write_gene_seq_log_for_file, write_error_seq_log_for_file

MATCH = MAT = 2
MISMATCH = MIS = -1
GAP_OPEN = -50
GAP_EXTEND = -4

ALIGN_MIN = 50
PHRED_MEANINGFUL_MIN = 30
MULTIPROCESSING_MAX = 1
FILIAL_NO = 1

ERR_PADDING = 1

# 'X' = for both end of main sequence, meaning the subsequence must be between this
ALIGN_MATRIX_FOR_SUBSEQUENCE_POSITION = {
    ('A', 'A'): MAT,    ('A', 'T'): MIS,    ('A', 'G'): MIS,    ('A', 'C'): MIS,    ('A', 'X'): -1000,
    ('T', 'A'): MIS,    ('T', 'T'): MAT,    ('T', 'G'): MIS,    ('T', 'C'): MIS,    ('T', 'X'): -1000,
    ('G', 'A'): MIS,    ('G', 'T'): MIS,    ('G', 'G'): MAT,    ('G', 'C'): MIS,    ('G', 'X'): -1000,
    ('C', 'A'): MIS,    ('C', 'T'): MIS,    ('C', 'G'): MIS,    ('C', 'C'): MAT,    ('C', 'X'): -1000,
    ('X', 'A'): -1000,  ('X', 'T'): -1000,  ('X', 'G'): -1000,  ('X', 'C'): -1000,  ('X', 'X'): MAT,
}


def get_result_text_set(sorted_result: list, f_no: int = 1):
    res_text = "hetero(+/-)"
    err_list = []
    n = n1 = n2 = ne = 0
    s1 = s2 = 0
    p1 = p2 = pe = -1
    if len(sorted_result) < 1:
        err_list.append("WARNING: no result")
        return res_text, err_list, s1, s2

    for i, res in enumerate(sorted_result):
        n += res[1]
        if n1 == 0:
            s1, n1 = res[0], res[1]
            p1 = i
        elif n2 == 0:
            s2, n2 = res[0], res[1]
            p2 = i

    if ne > n * 0.1:
        err_list.append("error rate is higher than 10%")
    if f_no == 0:
        return res_text, err_list, s1, s2

    if n1 + n2 < ((n - ne) * 0.9):
        err_list.append("two largest genotypes are less than 90%(without error)")

    if n2 > n1 * 0.5:
        if "WT" in (s1, s2):
            res_text = "hetero(+/-)"
        else:
            res_text = "hetero(1/2)"
    else:
        if "WT" in (s1,):
            res_text = "homo(-/-)"
        else:
            res_text = "homo(+/+)"
    return res_text, err_list, s1, s2


def get_align_line_set(ref_str: str, seq_str: str):
    ref_str = "X" + ref_str + "X"

    aligner = PairwiseAligner()

    matrix = substitution_matrices.Array(data=ALIGN_MATRIX_FOR_SUBSEQUENCE_POSITION)

    aligner.substitution_matrix = matrix

    aligner.open_gap_score = GAP_OPEN
    aligner.extend_gap_score = GAP_EXTEND

    alignments = aligner.align(ref_str, seq_str)

    ref_line_untrimmed = alignments[0][0]
    seq_line_untrimmed = alignments[0][1]

    ref_line = match_line = seq_line = ""
    count_base = 0
    for i, c in enumerate(seq_line_untrimmed):
        if count_base >= len(seq_str):
            break
        if len(seq_line) == 0 and c == "-":
            continue
        seq_line += c
        ref_line += ref_line_untrimmed[i]
        count_base += 1
        if seq_line_untrimmed[i] == '-':
            match_line += 'x'
            count_base -= 1
        elif ref_line_untrimmed[i] == '-':
            match_line += '+'
        elif seq_line_untrimmed[i] == ref_line_untrimmed[i]:
            match_line += '|'
        else:
            match_line += '.'

    line_set = {'ref_line': ref_line, 'match_line': match_line, 'seq_line': seq_line}
    return line_set


def get_align_score(indel_length: int, ref_line: str, match_line: str, seq_line: str, **kwargs):
    if len(match_line) < 5:
        return 0
    a = 0
    for c in match_line[1:-2]:
        if c == '|':
            a += 1
    score = a / (len(match_line[ERR_PADDING:-ERR_PADDING])-indel_length)
    return score


def get_phred_line(seq_raw):
    # only from fastq

    quality = seq_raw.letter_annotations['phred_quality']

    if len(quality) < 1:
        return "# phred score not available"

    phred = ""
    for a in quality:
        phred += str(chr(a+33))
    return phred


def get_aligned_phred_line(phred_line: str, seq_line: str, **kwargs):
    aligned_phred_line = ""

    for c in seq_line:
        if c == '-':
            aligned_phred_line += " "
        else:
            aligned_phred_line += phred_line[0]
            phred_line = phred_line[1:]

    return aligned_phred_line


def get_best_align_line_set_list_dict(ref_raw_list: list, seq_raw_list: list, g_rna_raw_list: list):

    result_line_set_list_dict = {"err": []}

    for ref_raw in ref_raw_list:
        result_line_set_list_dict[str(ref_raw.name)] = []

    for i, seq_raw in enumerate(seq_raw_list):
        seq_str = str(seq_raw.seq).upper()

        phred_line = get_phred_line(seq_raw)

        best_line_set = {'ref_line': "", 'match_line': "", 'seq_line': "", 'align_score': -1000,
                         'name': seq_raw.name, 'phred_line': "", "indel_type": "", "pos_line": "", "reason": "."}
        test_line_set = {'ref_line': "", 'match_line': "", 'seq_line': "", 'name': seq_raw.name, 'phred_line': ""}
        best_ref = "err"
        best_indel_length = 0

        for j, ref_raw in enumerate(ref_raw_list):
            ref_str = str(ref_raw.seq).upper()
            test_line_set = get_align_line_set(ref_str, seq_str)

            test_line_set['phred_line'] = get_aligned_phred_line(phred_line, **test_line_set)
            # test_line_set['name'] = seq_raw.name

            g_rna = str(g_rna_raw_list[j].seq).upper()

            indel_dict = get_main_indel(test_line_set, g_rna)
            test_indel_length = indel_dict["indel_length"]

            if get_align_score(test_indel_length, **test_line_set) > get_align_score(best_indel_length, **best_line_set):
                best_line_set['ref_line'] = test_line_set['ref_line']
                best_line_set['match_line'] = test_line_set['match_line']
                best_line_set['seq_line'] = test_line_set['seq_line']
                best_line_set['phred_line'] = test_line_set['phred_line']

                best_line_set['align_score'] = int((get_align_score(best_indel_length, **best_line_set)-1)*1000)

                if len(test_line_set['match_line']) > ALIGN_MIN:
                    best_ref = ref_raw.name
                    best_indel_length = test_indel_length

                    best_line_set["pos_line"] = indel_dict["pos_line"]
                    best_line_set["indel_type"] = indel_dict["indel_type"]
                    best_line_set["reason"] = indel_dict["reason"]

            if get_align_score(best_indel_length, **best_line_set) < (1-ERR_MAX):
                best_line_set["reason"] = 'too many mismatch/error at the best matching sequence of ' + best_ref
                best_ref = 'err'
                best_line_set["indel_type"] = 'err'

        print(f"\r{round(i/len(seq_raw_list), 3)}", end="")
        best_line_set['name'] += f" at {best_ref}"
        result_line_set_list_dict[best_ref].append(best_line_set)

        # print()
        # print(best_line_set["indel_type"])
        # print(best_line_set["pos_line"])
        # print(best_line_set["ref_line"])
        # print(best_line_set["match_line"])
        # print(best_line_set["seq_line"])
        # print(best_line_set["phred_line"])

    return result_line_set_list_dict


def get_indel_pos_line(line_set: dict, pri: int, pre: int):

    if pri < 0 or pre < 0:
        return ""
    ref_line_buffer = line_set["ref_line"] + "X"*PAM_MAX

    pos_line = " "*pri
    for i in range(pri, pre):
        if ref_line_buffer[i] == '-':
            pos_line += "-"
        else:
            pos_line += ">"

    for i in range(PAM_MAX):
        if ref_line_buffer[pre + i] == ref_line_buffer[pre + i + 1] == "G":
            pos_line += '<<<'
            break
        else:
            pos_line += ' '

    if len(pos_line) < len(line_set["ref_line"]):
        pos_line += " "*(len(line_set["ref_line"]) - len(pos_line))
    else:
        pos_line = pos_line[:len(line_set["ref_line"])]

    return pos_line


def get_main_indel(line_set: dict, g_rna: str = ""):
    ref_line = line_set['ref_line']
    match_line = line_set['match_line']
    seq_line = line_set['seq_line']

    p1 = p2 = -1
    indel_i, indel_d, indel_length = 0, 0, 0

    indel_dict = {
        "pos_line": "", "indel_type": "WT", "indel_length": 0, "reason": "."
    }

    pri, pre = find_final_position(ref_line, g_rna)

    indel_dict["pos_line"] = get_indel_pos_line(line_set, pri, pre)

    for i, m in enumerate(match_line + str("X"*(PAM_MAX*2))):
        if i < 2:
            continue

        if m in MATCH_LETTER[1:]:
            p2 = i + 1
            indel_length += 1
            if p1 < 0:
                p1 = i
            if m == 'x':
                indel_d += 1
            if m == '+':
                indel_i += 1
            if m == '.':
                indel_d += 1
                indel_i += 1

        elif p2 >= 0:
            if indel_d == indel_i == indel_length == 1 and (ord(line_set["phred_line"][p1])-33) < PHRED_MEANINGFUL_MIN:
                    p1 = p2 = -1
                    indel_i = indel_d = indel_length = 0
                    continue

            if ((pre - PAM_MAX) < p2) and (p1 < (pre + PAM_MAX)):
                indel_dict["indel_type"] = get_indel_shape_text(indel_i, indel_d, p2 - pre)
                indel_dict["reason"] = "Position confirmed by guide RNA sequence"
                indel_dict["indel_length"] = indel_length
                break
            p1 = p2 = -1
            indel_i = indel_d = indel_length = 0

        if i > (len(ref_line)):
            p1 = p2 = -1
            break

    return indel_dict


def get_indel_count_for_file_dict(line_set_list_dict: dict):

    indel_count_for_file = {
        'err': {'err':0}
    }

    for key in line_set_list_dict.keys():

        indel_count_dict = {}

        line_set_list = line_set_list_dict[key]

        for line_set in line_set_list:
            if line_set['indel_type'] in indel_count_dict.keys():
                indel_count_dict[line_set['indel_type']] += 1
            else:
                indel_count_dict[line_set['indel_type']] = 1

        indel_count_list = list(sorted(indel_count_dict.items(), key=lambda k: k[1], reverse=True))

        indel_count_for_file[key] = indel_count_list

    return indel_count_for_file


if __name__ == '__main__':

    # <TEST>

    # </TEST>

    # Get the addresses of testing files, reference sequences, and the guide RNA sequence.
    address_list = [file_name for file_name in os.listdir(DATA_ADDRESS)
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_name)) and file_name[-6:] == '.fastq']
    print(address_list)
    ref_raw_iter = SeqIO.parse(REF_SET_ADDRESS, "fasta")

    ref_raw_list = []
    for item in ref_raw_iter:
        ref_raw_list.append(item)

    g_rna_seq_iter = SeqIO.parse(GUIDE_RNA_SET_ADDRESS, "fasta")

    g_rna_seq_list = []
    for item in g_rna_seq_iter:
        g_rna_seq_list.append(item)

    g_rna_seq = read_g_rna(GUIDE_RNA_ADDRESS)

    # For each file, this will read, align and build the line set.
    # Also, the line set will be used to count the indels, and making the result.
    
    indel_count_total = {}
    
    for file_name in address_list:
        seq_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")

        seq_raw_list = []
        for item in seq_raw_iter:
            seq_raw_list.append(item)

        print(f"for {file_name}:")
        line_set_list_dict = get_best_align_line_set_list_dict(ref_raw_list, seq_raw_list, g_rna_seq_list)

        # line_set = {name, (pos, ref, match, seq, phred)_line, indel_type, align_score, reason}

        indel_count_for_file_dict = get_indel_count_for_file_dict(line_set_list_dict)

        indel_count_total[file_name] = indel_count_for_file_dict

        write_gene_seq_log_for_file(file_name, line_set_list_dict, indel_count_for_file_dict)
        write_error_seq_log_for_file(file_name, line_set_list_dict, indel_count_for_file_dict)

    # write_result_main_log(indel_count_total)


        # TODO
