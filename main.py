# This is a sample Python script.
import os

from Bio import SeqIO

from InDel_Counter import DATA_ADDRESS, REF_SET_ADDRESS, GUIDE_RNA_SET_ADDRESS
from line_set import Line_Set, Reference_Set, InDel_Counter_for_Ref

# Press the green button in the gutter to run the script.
if __name__ == '__main__':


    address_list = [file_name for file_name in os.listdir(DATA_ADDRESS)
                    if os.path.isfile(os.path.join(DATA_ADDRESS, file_name)) and file_name[-6:] == '.fastq']
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
            i = len(g_rna_raw_list) -1
        ref_set = Reference_Set(ref_raw=ref_raw, guide_rna_raw=g_rna_raw_list[i])
        ref_set_list.append(ref_set)

    log = open("log_test.txt", 'w')

    for file_name in address_list:

        indel_counter_list = []
        for ref_set in ref_set_list:
            indel_counter = InDel_Counter_for_Ref(ref_set=ref_set)
            indel_counter.set_file_name(file_name=file_name)
            indel_counter_list.append(indel_counter)

        seq_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")

        seq_raw_list = []
        for item in seq_raw_iter:
            seq_raw_list.append(item)

        line_set_list = []
        print()
        for i, seq_raw in enumerate(seq_raw_list):
            print(f"\r{round((i+1)/len(seq_raw_list), 3)}", end="")

            best_line_set = Line_Set(read_raw=seq_raw, ref_set=ref_set_list[0])
            best_line_set.set_file_name(file_name=file_name)

            for ref_set in ref_set_list:
                test_line_set = Line_Set(read_raw=seq_raw, ref_set=ref_set)
                test_line_set.set_file_name(file_name=file_name)

                if best_line_set.score < test_line_set.score:
                    best_line_set = test_line_set

            line_set_list.append(best_line_set)

        for line_set in line_set_list:
            for indel_counter in indel_counter_list:
                indel_counter.append(line_set)

        for indel_counter in indel_counter_list:
            log.write(str(indel_counter))
            log.write("\n\n")

        for line_set in line_set_list:
            log.write(str(line_set))
            log.write("\n\n")




    log.close()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
