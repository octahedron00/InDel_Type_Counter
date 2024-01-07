import datetime
import os
from pathlib import Path

from InDel_Counter import PAM_MAX, ERR_MAX, MATCH_LETTER, get_result_text_set
from InDel_Counter import RESULT_LOG_ADDRESS, ALIGN_LOG_ADDRESS, DATA_ADDRESS, GUIDE_RNA_ADDRESS, \
    GUIDE_RNA_SET_ADDRESS, REF_SET_ADDRESS, TASK_TITLE

g_rna_seq = "g_rna_seq_,,,issue"
FILIAL_NO = 1


def get_file_txt_name(file_name: str):
    return Path(file_name).stem + ".txt"


def write_gene_seq_log_for_file(file_name: str, line_set_list_dict: dict, indel_count_for_file_dict: dict):
    file_log = open(os.path.join(ALIGN_LOG_ADDRESS, "gene_"+get_file_txt_name(file_name)), 'w')

    # print(line_set_list_dict['err'])
    # print(sorted(line_set_list_dict['err'], key=lambda x: len(x['match_line'])))

    file_log.write(f""
                   f"# <InDel_Counter Side Log for {file_name}>\n"
                   f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n"
                   f"# Guide RNA sequence : {g_rna_seq}\n"
                   f"# \n"
                   f"# \n"
                   f"# {file_name} as a data for Transgenic filial {FILIAL_NO}\n")

    for key in indel_count_for_file_dict.keys():
        indel_count_list = indel_count_for_file_dict[key]

        res_text, err_list, s1, s2 = get_result_text_set(indel_count_list, FILIAL_NO)

        file_log.write("\n\n")
        file_log.write(f"For {key}:\n")

        for err_text in err_list:
            file_log.write(f"# WARNING : {err_text}\n")
        file_log.write("\n[Result]\n")
        for item in indel_count_list:
            file_log.write(f"{item[0]}\t{item[1]}\n")
        print(indel_count_list)
        if len(indel_count_list) == 0:
            continue
        if res_text[:2] == "ho":
            file_log.write(f"\n{res_text} of "
                           f"{indel_count_list[0][0]} ({indel_count_list[0][1]} of {len(line_set_list_dict[key])})\n")
        else:
            file_log.write(f"\n{res_text} of "
                           f"{indel_count_list[0][0]} ({indel_count_list[0][1]} of {len(line_set_list_dict[key])}) and "
                           f"{indel_count_list[1][0]} ({indel_count_list[1][1]} of {len(line_set_list_dict[key])})\n")
        file_log.write("\n\n@\n")

    for key in line_set_list_dict.keys():
        if key == 'err':
            continue
        # sort by the length: to put the same sequences in the similar spot
        line_set_list = sorted(line_set_list_dict[key], key=lambda x: len(x['match_line']))
        for line_set in line_set_list:
            file_log.write(f"@\t{key}\t{line_set['name']}\n")
            file_log.write(f"@\t{line_set['indel_type']}\t{line_set['align_score']}\t{line_set['reason']}\n")
            file_log.write(f"{line_set['pos_line']}\n")
            file_log.write(f"{line_set['ref_line']}\n")
            file_log.write(f"{line_set['match_line']}\n")
            file_log.write(f"{line_set['seq_line']}\n")
            file_log.write(f"{line_set['phred_line']}\n\n\n")
    file_log.close()


def write_error_seq_log_for_file(file_name: str, line_set_list_dict: dict, indel_count_for_file_dict: dict):
    file_log = open(os.path.join(ALIGN_LOG_ADDRESS, "error_" + get_file_txt_name(file_name)), 'w')

    # print(line_set_list_dict['err'])
    # print(sorted(line_set_list_dict['err'], key=lambda x: len(x['match_line'])))

    file_log.write(f"# <InDel_Counter Side Error Log for {file_name}>\n")
    file_log.write(
        f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n")
    file_log.write(f"# PAM_MAX = {PAM_MAX}, ERR_MAX = {ERR_MAX}\n\n")
    file_log.write(f"# Guide RNA sequence : {g_rna_seq}\n\n\n")

    file_log.write(f"# {file_name} as a data for Transgenic filial {FILIAL_NO}\n")

    #
    for key in indel_count_for_file_dict.keys():

        file_log.write("\n\n")
        file_log.write(f"For {key}:\n")

        indel_count_list = indel_count_for_file_dict[key]

        res_text, err_list, s1, s2 = get_result_text_set(indel_count_list, FILIAL_NO)

        for err_text in err_list:
            file_log.write(f"# WARNING : {err_text}\n")
        file_log.write("\n[Result]\n")
        for item in indel_count_list:
            file_log.write(f"{item[0]}\t{item[1]}\n")
        if len(indel_count_list) == 0:
            continue
        if res_text[:2] == "ho":
            file_log.write(f"\n{res_text} of "
                                 f"{indel_count_list[0][0]} ({indel_count_list[0][1]} of {len(line_set_list_dict[key])})\n")
        else:
            file_log.write(f"\n{res_text} of "
                                 f"{indel_count_list[0][0]} ({indel_count_list[0][1]} of {len(line_set_list_dict[key])}) and "
                                 f"{indel_count_list[1][0]} ({indel_count_list[1][1]} of {len(line_set_list_dict[key])})\n")
        file_log.write("\n\n@\n")

    # sort by the length: to put the same sequences in the similar spot

    key = 'err'

    print("lsld", key)
    print(line_set_list_dict.keys())
    print(len(line_set_list_dict[key]))

    line_set_list = sorted(line_set_list_dict[key], key=lambda x: len(x['match_line']))
    for line_set in line_set_list:
        file_log.write(f"@\t{key}\t{line_set['name']}\n")
        file_log.write(f"@\t{line_set['indel_type']}\t{line_set['align_score']}\t{line_set['reason']}\n")
        file_log.write(f"{line_set['pos_line']}\n")
        file_log.write(f"{line_set['ref_line']}\n")
        file_log.write(f"{line_set['match_line']}\n")
        file_log.write(f"{line_set['seq_line']}\n")
        file_log.write(f"{line_set['phred_line']}\n\n\n")
    file_log.close()


def write_result_main_log():
    pass
    # file_log = open("Count_result.txt", 'w')
#
#     # print(line_set_list_dict['err'])
#     # print(sorted(line_set_list_dict['err'], key=lambda x: len(x['match_line'])))
#
#     file_log.write(f"# <InDel_Counter Main Log>\n")
#     file_log.write(
#         f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n")
#     file_log.write(f"# PAM_MAX = {PAM_MAX}, ERR_MAX = {ERR_MAX}\n\n")
#     file_log.write(f"# Guide RNA sequence : {g_rna_seq}\n\n\n")
#
#     file_log.write(f"# {file_name} as a data for Transgenic filial {FILIAL_NO}\n")
#
#     #
#     for key in indel_count_for_file_dict.keys():
#         indel_count_list = indel_count_for_file_dict[key]
#
#         res_text, err_list, s1, s2 = get_result_text_set(indel_count_list, FILIAL_NO)
#
#         file_log.write("\n\n")
#         file_log.write(f"For {key}:\n")
#
#         for err_text in err_list:
#             file_log.write(f"# WARNING : {err_text}\n")
#         file_log.write("\n[Result]\n")
#         for item in indel_count_list:
#             file_log.write(f"{item[0]}\t{item[1]}\n")
#         if res_text[:2] == "ho":
#             file_log.write(f"\n{res_text} of "
#                            f"{indel_count_list[s1][0]} ({indel_count_list[s1][1]} of {len(line_set_list_dict[key])})\n")
#         else:
#             file_log.write(f"\n{res_text} of "
#                            f"{indel_count_list[s1][0]} ({indel_count_list[s1][1]} of {len(line_set_list_dict[key])}) and "
#                            f"{indel_count_list[s2][0]} ({indel_count_list[s2][1]} of {len(line_set_list_dict[key])})\n")
#         file_log.write("\n\n@\n")
#
#     file_log.close()
#
