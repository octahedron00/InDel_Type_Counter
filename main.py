import datetime
import os
import re

from Bio import SeqIO

from line_set import Line_Set, Reference_Set, InDel_Counter_for_Ref
from log_writer import write_main_log, write_sub_log, write_main_csv_log

DATA_ADDRESS = "./data/"
GUIDE_RNA_ADDRESS = "./ref/guide_RNA.txt"
GUIDE_RNA_SET_ADDRESS = "./ref/guide_RNA_set.fasta"
REF_SET_ADDRESS = "./ref/reference_seq_set.fasta"
TASK_TITLE = ""
FILIAL_NO = 1

# with spans, 3 second for 1000 lines: 20000 for a minute, 600000: 30 minutes


def get_best_line_set(read: SeqIO.SeqRecord, ref_set_list: list):
    if len(ref_set_list) == 0:
        return None
    best_line_set = Line_Set(read_raw=read, ref_set=ref_set_list[0])

    for ref_set in ref_set_list:
        test_line_set = Line_Set(read_raw=read, ref_set=ref_set)

        if test_line_set.score > best_line_set.score:
            best_line_set = test_line_set

    return best_line_set


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    address_list = [file_name for file_name in os.listdir(DATA_ADDRESS)
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_name)) and file_name[-6:] == '.fastq']
    # print(address_list)

    address_list.sort(key=lambda f: int(''.join(filter(str.isdigit, f)) + '0'))
    print(address_list)

    ref_raw_iter = SeqIO.parse(REF_SET_ADDRESS, "fasta")

    ref_raw_list = []
    for item in ref_raw_iter:
        ref_raw_list.append(item)

    g_rna_seq_iter = SeqIO.parse(GUIDE_RNA_SET_ADDRESS, "fasta")

    g_rna_raw_list = []
    for item in g_rna_seq_iter:
        g_rna_raw_list.append(item)

    ref_set_list = []
    for i, ref_raw in enumerate(ref_raw_list):
        if i >= len(g_rna_raw_list):
            i = len(g_rna_raw_list) - 1
        ref_set = Reference_Set(ref_raw=ref_raw, guide_rna_raw=g_rna_raw_list[i])
        ref_set_list.append(ref_set)

    log = open("log_test.txt", 'w')

    all_indel_counter_list_list = []

    print()
    total_length = 0
    for i, file_name in enumerate(address_list):
        read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
        print(f"\r({i+1}/{len(address_list)}) reading {file_name}", end="")
        for item in read_raw_iter:
            total_length += 1
    print(f"\rTotal reads :{total_length} for {len(address_list)} files")
    finished_length = 0

    print()
    start_time = datetime.datetime.now()

    for file_no, file_name in enumerate(address_list):

        start_time_for_file = datetime.datetime.now()
        indel_counter_list = []
        for ref_set in ref_set_list:
            indel_counter = InDel_Counter_for_Ref(ref_set=ref_set)
            indel_counter.set_file_name(file_name=file_name)
            indel_counter_list.append(indel_counter)

        read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")

        read_raw_list = []
        for item in read_raw_iter:
            read_raw_list.append(item)

        line_set_list = []
        for i, read_raw in enumerate(read_raw_list):
            finished_length += 1
            now_time = datetime.datetime.now()
            delta_time = now_time - start_time

            print(f"\r({file_no + 1}/{len(address_list)}) "
                  f"for {file_name}: {((i+1)/len(read_raw_list)):.3f} / "
                  f"remaining: {(delta_time/finished_length)*(total_length-finished_length)} "
                  f"(for this file: {(delta_time/finished_length)*(len(read_raw_list)-(i+1))})", end="")

            best_line_set = get_best_line_set(read_raw, ref_set_list)
            best_line_set.set_file_name(file_name=file_name)

            line_set_list.append(best_line_set)

        print(f"\r({file_no + 1}/{len(address_list)}) for {file_name}: Complete / Writing log files", end="")
        for line_set in line_set_list:
            for indel_counter in indel_counter_list:
                indel_counter.count(line_set)

        for indel_counter in indel_counter_list:
            write_sub_log(line_set_list=line_set_list, indel_counter=indel_counter, file_name=file_name)

        all_indel_counter_list_list.append(indel_counter_list)

        end_time_for_file = datetime.datetime.now()

        print(f"\r({file_no + 1}/{len(address_list)}) for {file_name}: Complete / Log written / "
              f"{end_time_for_file - start_time_for_file} is used (length: {len(line_set_list)})")

    write_main_log(indel_counter_list_list=all_indel_counter_list_list)
    write_main_csv_log(indel_counter_list_list=all_indel_counter_list_list, ref_set_list=ref_set_list)
    log.close()
    print(f"Work Completed! (total time: {datetime.datetime.now() - start_time})")

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
