import datetime
import os
import click
from Bio import SeqIO

from line_set import Line_Set, Reference, InDel_Counter_for_Ref
from log_writer import write_main_log, write_sub_log, write_main_csv_log, XLSX_LOG_NAME
import globals

DATA_ADDRESS = "./data/"
GUIDE_RNA_ADDRESS = "./ref/guide_RNA.txt"
GUIDE_RNA_SET_ADDRESS = "./ref/guide_RNA_set.fasta"
REF_SET_ADDRESS = "./ref/reference_seq_set.fasta"

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


@click.command()
@click.option('-e', '--err_max', default=0.03,
              help="The threshold of mismatch ratio in aligned line set, without the main indel.")
@click.option('-r', '--pam_range_max', default=5,
              help="The max range of main indel position from PAM sequence")
@click.option('-p', '--err_padding', default=1,
              help="The mismatch in this padding length from both end will not be counted")
@click.option('-s', '--phred_meaningful_score_min', default=30,
              help="One mismatch will be recognised as an indel, only if the phred score of the nucleotide is higher; "
                   "For ignoring all 'one mismatch', make it higher than 100")
@click.option('--score_match', default=2,
              help="Score for align: for Match")
@click.option('--score_mismatch', default=-1,
              help="Score for align: for Mismatch")
@click.option('--score_gap_open', default=-50,
              help="Score for align: for Gap Open; low penalty will make fake indels more often")
@click.option('--score_gap_extend', default=-4,
              help="Score for align: for Gap Extension; low penalty will make fake indels more often")
@click.option('-t', '--task_title', default="Task at " + str(datetime.datetime.now()), help="Title for this task")
@click.option('-o', '--open_xlsx_auto', default=False, help="Open the excel log file automatically if finished")
def main(err_max, pam_range_max, err_padding, phred_meaningful_score_min,
         score_match, score_mismatch, score_gap_open, score_gap_extend, task_title, open_xlsx_auto):
    #
    # set global variables
    globals.ERR_MAX = err_max
    globals.PAM_RANGE_MAX = pam_range_max
    globals.ERR_PADDING = err_padding
    globals.PHRED_MEANINGFUL_MIN = phred_meaningful_score_min
    globals.MAT = score_match
    globals.MIS = score_mismatch
    globals.GAP_OPEN = score_gap_open
    globals.GAP_EXTEND = score_gap_extend
    globals.TASK_TITLE = task_title
    globals.OPEN_XLSX_AUTO = open_xlsx_auto
    
    #
    # Get addresses of files from fixed file location, and sort
    address_list = [file_name for file_name in os.listdir(DATA_ADDRESS)
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_name)) and file_name[-6:] == '.fastq']
    address_list.sort(key=lambda f: int(''.join(filter(str.isdigit, f)) + '0'))
    print("File list:", address_list)

    #
    # Get reference and guide RNA sequences, 
    # and match them as a 'Reference' class
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
        ref_set = Reference(ref_raw=ref_raw, guide_rna_raw=g_rna_raw_list[i])
        ref_set_list.append(ref_set)

    #
    # for total result, making list[list[InDel_Counter]]
    all_indel_counter_list_list = []

    #
    # for counting expected time left, count total number of reads
    total_reads_count = 0
    for i, file_name in enumerate(address_list):
        print(f"\r({i+1}/{len(address_list)}) reading {file_name}", end="")

        read_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")
        for item in read_raw_iter:
            total_reads_count += 1
    print(f"\rTotal reads :{total_reads_count} for {len(address_list)} files")
    
    #
    # for counting expected time left, count total number of finished number of reads
    finish_reads_count = 0

    # for counting expected time left, check the time of start
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
            finish_reads_count += 1
            now_time = datetime.datetime.now()
            delta_time = now_time - start_time

            # to prevent print function's time delay problem:
            # This will show the update at about every 0.3 seconds
            if (i % 100) == 0:
                print(f"\r({file_no + 1}/{len(address_list)}) "
                      f"for {file_name}: {((i+1)/len(read_raw_list)):.3f} / "
                      f"remaining: {(delta_time/finish_reads_count)*(total_reads_count-finish_reads_count)} "
                      f"(for this file: {(delta_time/finish_reads_count)*(len(read_raw_list)-(i+1))}) "
                      f"(length: {len(read_raw_list)})", end="")

            best_line_set = get_best_line_set(read_raw, ref_set_list)
            best_line_set.set_file_name(file_name=file_name)

            line_set_list.append(best_line_set)

        print(f"\r({file_no + 1}/{len(address_list)}) for {file_name}: Complete / "
              f"Writing log files (length: {len(line_set_list)})", end="")

        for line_set in line_set_list:
            for indel_counter in indel_counter_list:
                if indel_counter.ref_name == line_set.ref_name:
                    indel_counter.count(line_set)

        # Sorting Line Set List.
        indel_counter_map = {}
        for indel_counter in indel_counter_list:
            indel_counter_map[indel_counter.ref_name] = indel_counter
        line_set_list.sort(key=lambda l: len(l))
        line_set_list.sort(key=lambda l: indel_counter_map[l.ref_name].count_map[l.indel_type])
        line_set_list.sort(key=lambda l: l.indel_type == 'err')

        # writing sub log
        for indel_counter in indel_counter_list:
            write_sub_log(line_set_list=[l for l in line_set_list if l.ref_name == indel_counter.ref_name],
                          indel_counter=indel_counter, file_name=file_name)

        # finish: add it to the list, end time calc, print.
        all_indel_counter_list_list.append(indel_counter_list)
        end_time_for_file = datetime.datetime.now()
        print(f"\r({file_no + 1}/{len(address_list)}) for {file_name}: Complete / Log written / "
              f"{end_time_for_file - start_time} ({end_time_for_file - start_time_for_file} for this file) is passed "
              f"(length: {len(line_set_list)})")

    write_main_log(indel_counter_list_list=all_indel_counter_list_list)
    write_main_csv_log(indel_counter_list_list=all_indel_counter_list_list, ref_set_list=ref_set_list)

    print(f"Work Completed! (total time: {datetime.datetime.now() - start_time})")
    if globals.OPEN_XLSX_AUTO:
        os.system(f"start EXCEL.EXE {XLSX_LOG_NAME}")


if __name__ == '__main__':
    print("InDel Type Counter ver. "+globals.VERSION)
    print("Initiating...")
    main()
    print("Things never works below here... << error message")




