# This is a sample Python script.
import os

from Bio import SeqIO

from InDel_Counter import DATA_ADDRESS, REF_SET_ADDRESS, GUIDE_RNA_SET_ADDRESS
from line_set import Line_Set, Reference_Set

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

    ref_set = Reference_Set(ref_raw=ref_raw_list[0], guide_rna_raw=g_rna_raw_list[0])

    log = open("log_test.txt", 'w')

    for file_name in address_list[:-1]:
        seq_raw_iter = SeqIO.parse(os.path.join(DATA_ADDRESS, file_name), "fastq")

        seq_raw_list = []
        for item in seq_raw_iter:
            line_set = Line_Set(read_raw=item, ref_set=ref_set)

            log.write(str(line_set))
            log.write("\n"
                      "\n")
            seq_raw_list.append(item)

    log.close()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
