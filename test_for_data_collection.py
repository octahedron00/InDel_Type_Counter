#!/bin/sh

import datetime
import os
import click
from Bio import SeqIO

from src.reference import Reference
from src.line_set import Line_Set
import src.globals as glv

DATA_ADDRESS = "./data_raw/"
GUIDE_RNA_ADDRESS = "./ref/guide_RNA.txt"
GUIDE_RNA_SET_ADDRESS = "./ref/guide_RNA_set.fasta"
REF_SET_ADDRESS = "./ref/reference_seq_set.fasta"


# with spans, 3 second for 1000 lines: 20000 for a minute, 600000: 30 minutes
# > total 800 nt of ref, 150 nt for a line: 40,000,000 for a second.


def get_best_line_set(read: SeqIO.SeqRecord, reference_list: list):
    if len(reference_list) == 0:
        return None
    best_line_set = Line_Set(read_raw=read, ref_set=reference_list[0])

    for ref_set in reference_list:
        test_line_set = Line_Set(read_raw=read, ref_set=ref_set)

        if test_line_set.score > best_line_set.score:
            best_line_set = test_line_set

    return best_line_set


def get_reference_list_from_file():
    # Get reference and guide RNA sequences,
    # and match them as a 'Reference' class

    reference_list = []

    ref_raw_iter = SeqIO.parse(REF_SET_ADDRESS, "fasta")
    ref_raw_list = []
    for item in ref_raw_iter:
        ref_raw_list.append(item)

    g_rna_seq_iter = SeqIO.parse(GUIDE_RNA_SET_ADDRESS, "fasta")
    g_rna_raw_list = []
    for item in g_rna_seq_iter:
        g_rna_raw_list.append(item)

    for i, ref_raw in enumerate(ref_raw_list):
        if i >= len(g_rna_raw_list):
            i = len(g_rna_raw_list) - 1
        ref_set = Reference(ref_raw=ref_raw, guide_rna_raw=g_rna_raw_list[i])
        reference_list.append(ref_set)

    return reference_list


def get_file_address_list():
    address_list = [file_name for file_name in os.listdir(DATA_ADDRESS)
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_name)) and file_name[-4:] == '.raw']
    address_list.sort(key=lambda f: int(''.join(filter(str.isdigit, f)) + '0'))
    return address_list


def get_total_number_of_reads(address_list: list[str]):
    total_reads_count = 0
    for i, file_name in enumerate(address_list):
        print(f"\r({i + 1}/{len(address_list)}) reading {file_name}", end="")

        read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
        for _ in read_raw_iter:
            total_reads_count += 1
    print(f"\rTotal reads :{total_reads_count} for {len(address_list)} files                     ")
    print()
    return total_reads_count


def key_for_sorting_err(line_set: Line_Set):
    if line_set.indel_type == 'err':
        return len(line_set)
    return 0


@click.command()
@click.option('-e', '--err_ratio_max', default=0.03,
              help=glv.EXPLANATION_MAP['err_ratio_max'])
@click.option('-p', '--err_padding_for_seq', default=1,
              help=glv.EXPLANATION_MAP['err_padding_for_seq'])
@click.option('-x', '--cut_pos_from_pam', default=-3,
              help=glv.EXPLANATION_MAP['cut_pos_from_pam'])
@click.option('-r', '--cut_pos_radius', default=5,
              help=glv.EXPLANATION_MAP['cut_pos_radius'])
#
@click.option('-s', '--phred_meaningful_score_min', default=30,
              help=glv.EXPLANATION_MAP['phred_meaningful_score_min'])
@click.option('-d', '--pam_distance_max', default=5,
              help=glv.EXPLANATION_MAP['pam_distance_max'])
@click.option('--score_match', default=2,
              help=glv.EXPLANATION_MAP['score_match'])
@click.option('--score_mismatch', default=-1,
              help=glv.EXPLANATION_MAP['score_mismatch'])
@click.option('--score_gap_open', default=-50,
              help=glv.EXPLANATION_MAP['score_gap_open'])
@click.option('--score_gap_extend', default=-4,
              help=glv.EXPLANATION_MAP['score_gap_extend'])
#
@click.option('-t', '--task_title', default="Task at " + str(datetime.datetime.now()),
              help=glv.EXPLANATION_MAP['task_title'])
@click.option('-o', '--open_xlsx_auto', default=False, is_flag=True,
              help=glv.EXPLANATION_MAP['open_xlsx_auto'])
def main(err_ratio_max, err_padding_for_seq, cut_pos_from_pam, cut_pos_radius,
         phred_meaningful_score_min, pam_distance_max,
         score_match, score_mismatch, score_gap_open, score_gap_extend, task_title, open_xlsx_auto):
    # set global variables
    glv.ERR_RATIO_MAX = err_ratio_max
    glv.ERR_PADDING_FOR_SEQ = err_padding_for_seq
    glv.CUT_POS_FROM_PAM = cut_pos_from_pam
    glv.CUT_POS_RADIUS = cut_pos_radius

    glv.PHRED_MEANINGFUL_MIN = phred_meaningful_score_min
    glv.PAM_DISTANCE_MAX = pam_distance_max

    glv.MAT = score_match
    glv.MIS = score_mismatch
    glv.GAP_OPEN = score_gap_open
    glv.GAP_EXTEND = score_gap_extend

    glv.TASK_TITLE = task_title
    glv.OPEN_XLSX_AUTO = open_xlsx_auto

    # Get sorted file address list from a folder, list[str]
    address_list = get_file_address_list()
    print("File list:", address_list)
    print()

    # get a list[Reference]
    # reference_list = get_reference_list_from_file()

    # for total result, making list[list[InDel_Counter]]
    # all_indel_counter_list_list = []

    # # # for counting expected time left,
    # # # count total number of reads,
    # # # count total number of finished number of reads,
    # # # and check the time of initiation
    # '''This function will make a text print: opening large file takes some time'''
    # total_reads_count = get_total_number_of_reads(address_list=address_list)
    # finish_reads_count = 0
    # start_time = datetime.datetime.now()

    for file_no, file_name in enumerate(address_list):
        with open(os.path.join(DATA_ADDRESS, file_name), 'r') as raw_log:

            lines = raw_log.readlines()

            indel_list = []

            for line in lines:
                line = line.strip()
                info_list = line.split('\t')

                if len(info_list) < 3:
                    continue

                type = info_list[0]
                indel_i = int(info_list[1])
                indel_d = int(info_list[2])
                indel_length = int(info_list[3])
                indel_pos = int(info_list[4])
                count = int(info_list[5])

                indel_list.append(IndelInfo(type, indel_i, indel_d, indel_length, indel_pos, count))

    pos_average = [0]*40
    pos_end = [0]*40

    phred = (',', ':', 'F')

    one_mismatch_total = {
        ',': 0, ':': 0, 'F': 0
    }

    for item in indel_list:

        for i in range(-20, 20):
            if (item.indel_pos-item.indel_length) < i <= item.indel_pos:
                pos_average[i] += int(item.count / item.indel_length)
            if item.indel_pos == i and item.indel_length > 0:
                pos_end[item.indel_pos] += item.count

        # print(item.type, pos_average)

        if len(item.type) > 4 and item.type[:4] == '1I1D' and item.type[-1] in phred:

            one_mismatch_total[item.type[-1]] += item.count

    print(pos_average)
    print(pos_end)
    print(one_mismatch_total)


class IndelInfo:
    type = ""
    indel_i = 0
    indel_d = 0
    indel_length = 0
    indel_pos = 0
    count = 0

    def __init__(self, type: str, indel_i: int, indel_d: int, indel_length: int, indel_pos: int, count: int):

        self.type = type
        self.indel_i = indel_i
        self.indel_d = indel_d
        self.indel_length = indel_length
        self.indel_pos = indel_pos
        self.count = count
        # print(type, indel_i, indel_d, indel_length, indel_pos, count)


if __name__ == '__main__':
    print("InDel Type Counter test for data collection. " + glv.VERSION)
    print()
    main()

    # this is how 'click' works...
    print("Things never works below here... << error message")




