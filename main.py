#!/bin/sh

import datetime
import os
import click
from Bio import SeqIO

from src.reference import Reference
from src.line_set import Line_Set
from src.indel_counter_for_genotype import InDel_Counter_for_Genotype
from src.log_writer import write_main_log, write_sub_log, write_main_csv_log, XLSX_LOG_NAME
import src.globals as glv

DATA_ADDRESS = "./data/"
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
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_name)) and file_name[-6:] == '.fastq']
    address_list.sort(key=lambda f: int(''.join(filter(str.isdigit, f)) + '0'))
    return address_list


def get_total_number_of_reads(address_list: list[str]):
    total_reads_count = 0
    for i, file_name in enumerate(address_list):
        print(f"\r({i + 1}/{len(address_list)}) reading {file_name}", end="")

        read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
        for _ in read_raw_iter:
            total_reads_count += 1
    print(f"\rTotal reads :{total_reads_count} for {len(address_list)} files")
    print()
    return total_reads_count


@click.command()
@click.option('-e', '--err_ratio_max', default=0.03,
              help=glv.EXPLANATION_MAP['err_ratio_max'])
@click.option('-p', '--err_padding_for_seq', default=1,
              help=glv.EXPLANATION_MAP['err_padding_for_seq'])
@click.option('-x', '--cut_pos_from_pam', default=-3,
              help=glv.EXPLANATION_MAP['cut_pos_from_pam'])
@click.option('-r', '--cut_pos_radius', default=3,
              help=glv.EXPLANATION_MAP['cut_pos_radius'])
#
@click.option('-d', '--pam_distance_max', default=5,
              help=glv.EXPLANATION_MAP['pam_distance_max'])
@click.option('-s', '--phred_meaningful_score_min', default=30,
              help=glv.EXPLANATION_MAP['phred_meaningful_score_min'])
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
         pam_distance_max, phred_meaningful_score_min,
         score_match, score_mismatch, score_gap_open, score_gap_extend, task_title, open_xlsx_auto):
    # set global variables
    glv.ERR_RATIO_MAX = err_ratio_max
    glv.PAM_DISTANCE_MAX = pam_distance_max
    glv.ERR_PADDING_FOR_SEQ = err_padding_for_seq
    glv.PHRED_MEANINGFUL_MIN = phred_meaningful_score_min
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
    reference_list = get_reference_list_from_file()

    # for total result, making list[list[InDel_Counter]]
    all_indel_counter_list_list = []

    # # # for counting expected time left,
    # # # count total number of reads,
    # # # count total number of finished number of reads,
    # # # and check the time of initiation
    '''This function will make a text print: opening large file takes some time'''
    total_reads_count = get_total_number_of_reads(address_list=address_list)
    finish_reads_count = 0
    start_time = datetime.datetime.now()

    for file_no, file_name in enumerate(address_list):
        '''
        For each file, This happens:
            get a list of 'reads' from file (NGS-fastq only),
            
            for each read:
                for each reference:
                    (make line_set from 1 read and 1 reference)
                    make possible aligns with each reference sequences,
                    align the phred score line,
                    get the position of guide_RNA and PAM sequence,
                    get the indel type by checking around the PAM starting point,
                    get the score(mismatch ratio without main indel) and check if it is error seq or not.
                
                pick the best aligned line_set by score, 
                and add it to the list.
                
            for each line_set in list:
                for each indel_counter in list:
                    if line_set.ref == indel_counter.ref:
                        indel_counter.count(line_set) < count the each indel type
                
            
        
        '''
        # # # for counting expected time left
        start_time_for_file = datetime.datetime.now()

        # build list[InDel_Counter_For_Ref]
        indel_counter_list = []
        for reference in reference_list:
            indel_counter = InDel_Counter_for_Genotype(ref_set=reference)
            indel_counter.set_file_name(file_name=file_name)
            indel_counter_list.append(indel_counter)

        # build list[Bio.SeqRecord]
        read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
        read_raw_list = [read_raw for read_raw in read_raw_iter]

        # build list[Line_Set]
        line_set_list = []
        for i, read_raw in enumerate(read_raw_list):

            # # # for showing expected time left
            finish_reads_count += 1
            if (i % 100) == 0:
                now_time = datetime.datetime.now()
                delta_time = now_time - start_time
                print(f"\r({file_no + 1}/{len(address_list)}) "
                      f"for {file_name}: {((i+1)/len(read_raw_list)):.3f} / "
                      f"remaining: {(delta_time/finish_reads_count)*(total_reads_count-finish_reads_count)} "
                      f"(for this file: {(delta_time/finish_reads_count)*(len(read_raw_list)-(i+1))}) "
                      f"(length: {len(read_raw_list)})", end="")

            best_line_set = get_best_line_set(read_raw, reference_list)
            best_line_set.set_file_name(file_name=file_name)

            line_set_list.append(best_line_set)

        # # # for showing expected time left / while log writing
        print(f"\r({file_no + 1}/{len(address_list)}) for {file_name}: Complete / "
              f"Writing log files (length: {len(line_set_list)})", end="")

        # count the number of each indel type,
        # also setting the best line_set for each indel type happens here
        for line_set in line_set_list:
            for indel_counter in indel_counter_list:
                if indel_counter.ref_name == line_set.ref_name:
                    indel_counter.count(line_set)

        # # save the indel type count to each line set
        # # just for the log file
        for line_set in line_set_list:
            for indel_counter in indel_counter_list:
                if indel_counter.ref_name == line_set.ref_name:
                    line_set.set_indel_same_type_count(indel_counter.count_map)

        # add the file result to the total result
        all_indel_counter_list_list.append(indel_counter_list)

        # # Sorting Line Set List!
        # # err for the last / biggest indel type first / higher score / higher phred score / longer one first
        # # just for the log file
        indel_counter_map = {}
        for indel_counter in indel_counter_list:
            indel_counter_map[indel_counter.ref_name] = indel_counter
        line_set_list.sort(key=lambda l: len(l))
        line_set_list.sort(key=lambda l: l.phred_score, reverse=True)
        line_set_list.sort(key=lambda l: l.score, reverse=True)
        line_set_list.sort(key=lambda l: indel_counter_map[l.ref_name].count_map[l.indel_type], reverse=True)
        line_set_list.sort(key=lambda l: l.indel_type == 'err')

        # Writing sub log
        for indel_counter in indel_counter_list:
            write_sub_log(line_set_list=[l for l in line_set_list if l.ref_name == indel_counter.ref_name],
                          indel_counter=indel_counter, file_name=file_name)

        # # for showing time used
        end_time_for_file = datetime.datetime.now()
        print(f"\r({file_no + 1}/{len(address_list)}) for {file_name}: Complete / Log written / "
              f"{end_time_for_file - start_time} ({end_time_for_file - start_time_for_file} for this file) is passed "
              f"(length: {len(line_set_list)})")
        # end.

    # Writing total log
    write_main_log(indel_counter_list_list=all_indel_counter_list_list)
    write_main_csv_log(indel_counter_list_list=all_indel_counter_list_list, ref_set_list=reference_list)

    # # for showing time used
    print(f"Work Completed! (total time: {datetime.datetime.now() - start_time})")

    #
    if glv.OPEN_XLSX_AUTO:
        os.system(f"start EXCEL.EXE {XLSX_LOG_NAME}")


if __name__ == '__main__':
    print("InDel Type Counter ver. "+glv.VERSION)
    print()
    main()

    # this is how 'click' works...
    print("Things never works below here... << error message")




