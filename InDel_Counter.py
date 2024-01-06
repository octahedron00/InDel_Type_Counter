import csv
import datetime
import json
import os

import numpy as np
import pandas as pd


DIV_LETTER = ('@',)
SEQ_LETTER = ('A', 'T', 'G', 'C', '-')
MATCH_LETTER = ('|', 'x', '+', '\\', '/', '.')

PAM_MAX = 5
ERR_MAX = 0.05
DATA_ADDRESS = "./data/"
ALIGN_LOG_ADDRESS = "./log/"
RESULT_LOG_ADDRESS = "./log/"
GUIDE_RNA_ADDRESS = "./ref/guide_RNA.txt"
GUIDE_RNA_SET_ADDRESS = "./ref/guide_RNA_set.fasta"
REF_SET_ADDRESS = "./ref/reference_seq_set.fasta"
TASK_TITLE = ""


# variables / functions naming rule
# text_ = english message
# line_ = sequence or string to perform with
# is_ = boolean

# etc. = int / char / etc.

# The sequence must have at least 5 NT span from the guide RNA sequence


def get_reformat_line_set_list(match_result_list: list):
    line_set_list = []

    # Making multiple-line seq to one long line
    # p_1: data / p_2: match / p_3: ref_seq position,
    p11, p12, p13, p21, p22, p23 = -1, -1, -1, -1, -1, -1
    for i, line in enumerate(match_result_list):
        if line[0] in DIV_LETTER:
            p11, p12, p13, p21, p22, p23 = -1, -1, -1, -1, -1, -1
        if line[0] in SEQ_LETTER:
            if p12 < 0:
                p11 = i
            elif p12 >= 0 and p13 < 0:
                p13 = i
            elif p13 >= 0 and p21 < 0:
                p21 = i
            elif p22 >= 0 and p23 < 0:
                p23 = i
                # print(f"{p21} to {p11}, {p22} to {p12}, {p23} to {p13}")
                match_result_list[p11] = match_result_list[p11][:-1] + match_result_list[p21]
                match_result_list[p12] = match_result_list[p12][:-1] + match_result_list[p22]
                match_result_list[p13] = match_result_list[p13][:-1] + match_result_list[p23]
                match_result_list[p21] = match_result_list[p22] = match_result_list[p23] = '.'
                p21 = p22 = p23 = -1

        if line[0] in MATCH_LETTER:
            if p12 < 0:
                p12 = i
            else:
                p22 = i

    # Making a list
    p11 = p12 = p13 = -1
    for i, line in enumerate(match_result_list):
        if line[0] in DIV_LETTER:
            line_set_list.append([match_result_list[p11], match_result_list[p12], match_result_list[p13]])
            p11 = p12 = p13 = -1
        if line[0] in SEQ_LETTER:
            if p11 < 0:
                p11 = i
            elif p12 >= 0 and p13 < 0:
                p13 = i
        if line[0] in MATCH_LETTER:
            p12 = i

    return line_set_list


def get_indel_shape_text(indel_i: int, indel_d: int, pos: int):
    if indel_i < 0:
        return 'err'
    if indel_i == 0:
        if indel_d == 0:
            return 'WT'
        else:
            return f"{indel_d}D{pos}"
    else:
        if indel_d == 0:
            return f"{indel_i}I{pos}"
        else:
            return f"{indel_i}I{indel_d}D{pos}"


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
        if n1 == 0 and res[0] != "err":
            s1, n1 = res[0], res[1]
            p1 = i
        elif n2 == 0 and res[0] != "err":
            s2, n2 = res[0], res[1]
            p2 = i
        elif res[0] == 'err':
            ne = res[1]
            pe = i

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


def find_final_position(ref_line: str, g_rna: str):
    pri = -1
    pre = -1
    if len(g_rna) > 0:
        for i in range(len(ref_line) - len(g_rna)):
            k, is_good = 0, False
            pri = i
            for j, n in enumerate(g_rna):
                while ref_line[i + j + k] == '-':
                    k += 1
                if ref_line[i + j + k] != n:
                    break
                if j == (len(g_rna) - 1):
                    is_good = True
                    pre = (i + j + k + 1)
            if is_good:
                break
    return pri, pre


def get_indel_count_result_list(no: int, seq_line: str, match_line: str, ref_line: str, g_rna: str = ""):
    p1 = p2 = -1
    pr = 0
    is_by_cas9 = False
    indel_i, indel_d, count_e = 0, 0, 0
    result = [no, "WT", [seq_line, match_line, ref_line], "", "."]

    pr = find_final_position(ref_line, g_rna)

    for i, m in enumerate(match_line[1:-3]):
        if m in MATCH_LETTER[1:]:
            count_e += 1

    for i, m in enumerate(match_line):
        if i < 2:
            continue

        if m in ('x', '+'):
            p2 = i + 1
            if p1 < 0: p1 = i
            if m == 'x':    indel_d += 1
            if m == '+':    indel_i += 1

        elif p2 >= 0:
            if (pr > 0) and (abs(p2 - 1 - pr) < PAM_MAX):
                result[1] = get_indel_shape_text(indel_i, indel_d, p2 - 1 - pr)
                result[4] = "Position confirmed by guide RNA sequence"
                break
            for j in range(PAM_MAX):
                if ref_line[i + j] == ref_line[i + j + 1] == 'G':
                    is_by_cas9 = True
            if is_by_cas9:
                result[1] = get_indel_shape_text(indel_i, indel_d, p2 - 1 - pr)
                result[4] = "Position confirmed by NGG PAM sequence"
            p1 = p2 = -1
            indel_i = indel_d = 0

        if i > (len(ref_line) - (PAM_MAX + 2)):
            p1 = p2 = -1
            break

    if ((count_e - p2 + p1) / (len(ref_line) - p2 + p1)) > ERR_MAX:
        result[1] = "err"
        result[4] = "Too many mismatch"

    result[3] = f"{int(-(count_e - p2 + p1) / (len(ref_line) - p2 + p1) * 1000)}"

    return result


def InDel_count_result(file_address: str, g_rna_seq: str = "", f_no: int = 1):
    return_result = {
        'data': '',
        'err': '',
        'result': ''
    }

    result_dict = {
        'err': 0
    }

    match_result_file = open(os.path.join(DATA_ADDRESS, file_address), "r")
    match_result_list = match_result_file.readlines()
    match_result_file.close()

    # line_set_builder
    line_set_list = get_reformat_line_set_list(match_result_list)
    result_list = []

    # Compare the length of each sequence
    for i, line_set in enumerate(line_set_list):
        if len(line_set[0]) is len(line_set[1]) and len(line_set[1]) is len(line_set[2]):
            pass
        else:
            err_m = f"SOMETHING WRONG: Length of seq/match is not the same at around line no. {i} of {file_address}"
            print(err_m)
            print(f"PLEASE CHECK the missing letter or accidentally added \\n around line no. {i} of {file_address}")
            return_result['err'] = err_m
            return return_result

    for k, line_set in enumerate(line_set_list):
        result_list.append(get_indel_count_result_list(k, line_set[0], line_set[1], line_set[2], g_rna_seq))

    for result_row in result_list:
        if result_row[1] in result_dict:
            result_dict[result_row[1]] += 1
        else:
            result_dict[result_row[1]] = 1

    return_result['data'] = result_list
    return_result['result'] = result_dict
    return return_result


def read_g_rna(g_rna_address: str = ""):
    if not os.path.isfile(g_rna_address):
        return ""
    file = open(g_rna_address, 'r')
    file_lines = file.readlines()
    file.close()

    g_rna_name, g_rna_seq = "", ""
    for file_line in file_lines:
        if len(file_line.strip()) == 0 or file_line[0] == "#":
            continue
        g_rna_map = json.JSONDecoder().decode(file_line)
        g_rna_name = g_rna_map['name']
        g_rna_seq = g_rna_map['seq']
        print("gene name of guide RNA : ", g_rna_name)
        print("guide RNA sequence : ", g_rna_seq)
    return g_rna_seq


def write_result_csv(list_file_result: iter, g_rna_seq: str = "", f_no: int = 1):
    file_csv = open("Count_result.csv", 'w', newline="")
    file_csv_writer = csv.writer(file_csv)

    file_csv_writer.writerow(["<InDel_Counter Main Log>"])
    file_csv_writer.writerow(
        [f"Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})"])
    file_csv_writer.writerow(["PAM_MAX", PAM_MAX, "ERR_MAX", ERR_MAX])
    file_csv_writer.writerow(["Guide RNA sequence", g_rna_seq])
    file_csv_writer.writerow([])
    file_csv_writer.writerow(["file", "T filial No", "error", "genotype", "_of", "# total", "# of genotype"])

    for file_result in list_file_result:
        row = [file_result[0], f_no]

        res_text, err_list, s1, s2 = get_result_text_set(file_result[2])

        row.append("")
        for err_text in err_list:
            row[-1] += ("# " + err_text)

        row.append(res_text)
        if res_text[:2] == "ho":
            row.append(f"{s1} / {s1}")
        else:
            row.append(f"{s1} / {s2}")

        row.append(len(file_result[1]['data']))

        for item in file_result[2]:
            if item[0] == 'err':
                continue
            row.append(f"{item[0]} {item[1]}")
        row.append(f"err {file_result[1]['result']['err']}")

        file_csv_writer.writerow(row)

    file_csv.close()


def write_result_main_log(list_file_result: iter, g_rna_seq: str = "", f_no: int = 1):
    file_log = open("Count_result.txt", 'w')

    file_log.write("# <InDel_Counter Main Log>\n")
    file_log.write(
        f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n")
    file_log.write(f"# PAM_MAX = {PAM_MAX}, ERR_MAX = {ERR_MAX}\n\n")
    file_log.write(f"# Guide RNA sequence : {g_rna_seq}\n\n\n")

    for file_result in list_file_result:
        file_log.write(f"# {file_result[0]} as a data for Transgenic filial {f_no}\n")
        res_text, err_list, s1, s2 = get_result_text_set(file_result[2])

        for err_text in err_list:
            file_log.write(f"# WARNING : {err_text}\n")
        file_log.write("\n[Result]\n")
        for item in file_result[2]:
            if item[0] == 'err': continue
            file_log.write(f"{item[0]}\t{item[1]}\n")
        file_log.write(f"err\t{file_result[1]['result']['err']}\n")
        if res_text[:2] == "ho":
            file_log.write(f"\n{res_text} of "
                           f"{s1} ({file_result[1]['result'][s1]} of {len(file_result[1]['data'])})\n")
        else:
            file_log.write(f"\n{res_text} of "
                           f"{s1} ({file_result[1]['result'][s1]} of {len(file_result[1]['data'])}) and "
                           f"{s2} ({file_result[1]['result'][s2]} of {len(file_result[1]['data'])})\n")
        file_log.write("\n\n@\n")
    file_log.close()


def write_result_sub_log(file_result: tuple, g_rna_seq: str = "", f_no: int = 1):
    file_sub = open(f"log/log_{file_result[0]}" , "w")
    file_sub.write(f"# <InDel_Counter Side Log for {file_result[0]}>\n")
    file_sub.write(
        f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n")
    file_sub.write(f"# PAM_MAX = {PAM_MAX}, ERR_MAX = {ERR_MAX}\n\n")
    file_sub.write(f"# Guide RNA sequence : {g_rna_seq}\n\n\n")

    file_sub.write(f"# {file_result[0]} as a data for Transgenic filial {f_no}\n")
    res_text, err_list, s1, s2 = get_result_text_set(file_result[2])

    for err_text in err_list:
        file_sub.write(f"# WARNING : {err_text}\n")
    file_sub.write("\n[Result]\n")
    for item in file_result[2]:
        if item[0] == 'err': continue
        file_sub.write(f"{item[0]}\t{item[1]}\n")
    file_sub.write(f"err\t{file_result[1]['result']['err']}\n")
    if res_text[:2] == "ho":
        file_sub.write(f"\n{res_text} of "
                       f"{s1} ({file_result[1]['result'][s1]} of {len(file_result[1]['data'])})\n")
    else:
        file_sub.write(f"\n{res_text} of "
                       f"{s1} ({file_result[1]['result'][s1]} of {len(file_result[1]['data'])}) and "
                       f"{s2} ({file_result[1]['result'][s2]} of {len(file_result[1]['data'])})\n")
    file_sub.write("\n\n@\n")

    for result_set in file_result[1]['data']:
        file_sub.write(f"# {result_set[0]}\n")
        file_sub.write(f"# {result_set[1]}  \t{result_set[3]}\t\t{result_set[4]}\n")
        file_sub.write(f"{result_set[2][0]}")
        file_sub.write(f"{result_set[2][1]}")
        file_sub.write(f"{result_set[2][2]}")
        file_sub.write("\n\n@\n")
    file_sub.close()


if __name__ == '__main__':
    g_rna_seq = read_g_rna(GUIDE_RNA_ADDRESS)
    list_file_result = []

    file_address_list = [file_address for file_address in os.listdir(DATA_ADDRESS)
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_address))]

    for file_address in file_address_list:
        count_result = InDel_count_result(file_address, g_rna_seq, f_no)
        sorted_result = list(sorted(count_result['result'].items(), key=lambda k: k[1], reverse=True))
        list_file_result.append((file_address, count_result, sorted_result))

    write_result_main_log(list_file_result, g_rna_seq, f_no)
    write_result_csv(list_file_result, g_rna_seq, f_no)

    for file_result in list_file_result:
        write_result_sub_log(file_result, g_rna_seq, f_no)

