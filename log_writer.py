import datetime
import glob
import os
from pathlib import Path
import csv
import pandas as pd
from xlsxwriter.workbook import Workbook


from line_set import InDel_Counter_for_Ref, Line_Set, Reference_Set, Genotype
from line_set import MAT, MIS, GAP_OPEN, GAP_EXTEND, PAM_MAX, ERR_MAX, ERR_PADDING

SUB_LOG_ADDRESS = "./log/"
RESULT_LOG_ADDRESS = "./log/"

MAIN_LOG_NAME = "./Count_result.txt"
CSV_LOG_NAME = "./Count_result.csv"
XLSX_LOG_NAME = "./Count_result.xlsx"

FILIAL_NO = 1


def write_sub_log(line_set_list: list[Line_Set], indel_counter: InDel_Counter_for_Ref, file_name: str):
    file_log = open(SUB_LOG_ADDRESS + file_name[:-6] + "---" + indel_counter.ref_name + ".txt", "w")

    file_log.write(f""
                   f"# <InDel_Counter Side Log for {file_name}>\n"
                   f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n"
                   f"# \n"
                   f"# {file_name} as a data for Transgenic filial {FILIAL_NO}\n"
                   f"\n"
                   f"\n")

    file_log.write(str(indel_counter))
    file_log.write("\n\n\n")

    for line_set in line_set_list:
        if line_set.ref_name == indel_counter.ref_name:
            file_log.write(str(line_set))
            file_log.write("\n\n")

    file_log.close()


def write_main_log(indel_counter_list_list: list[list[InDel_Counter_for_Ref]]):

    file_log = open("Count_result.txt", "w")

    file_log.write(f""
                   f"# <InDel_Counter Main Log>\n"
                   f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n"
                   f"\n"
                   f"\n"
                   f"\n")

    for indel_counter_list in indel_counter_list_list:
        file_log.write(f"<{indel_counter_list[0].file_name}>\n"
                       f"\n")
        for indel_counter in indel_counter_list:
            file_log.write(str(indel_counter))
            file_log.write("\n"
                           "\n"
                           "\n")

    file_log.close()


def write_main_csv_log(indel_counter_list_list: list[list[InDel_Counter_for_Ref]], ref_set_list: list[Reference_Set]):

    file_csv = open(CSV_LOG_NAME, 'w', newline="")
    file_csv_writer = csv.writer(file_csv)

    file_csv_writer.writerow(["<InDel_Counter Main Log>"])
    file_csv_writer.writerow(
        [f"Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})"])
    file_csv_writer.writerow(["PAM_MAX", PAM_MAX, "ERR_MAX", ERR_MAX])
    file_csv_writer.writerow([])
    file_csv_writer.writerow(["References:"])
    file_csv_writer.writerow(["Name", "seq", "Guide RNA name", "Guide RNA seq"])
    for ref_set in ref_set_list:
        file_csv_writer.writerow([ref_set.ref_name, ref_set.ref_seq, ref_set.guide_rna_name, ref_set.guide_rna_seq])

    file_csv_writer.writerow([])
    file_csv_writer.writerow(["file", "reference", "warning", "genotype", "_of",
                              "# total", "# error", "# not err", "# of genotype"])

    for indel_counter_list in indel_counter_list_list:
        total_lines_for_file = 0
        for indel_counter in indel_counter_list:
            total_lines_for_file += len(indel_counter)
        for i, indel_counter in enumerate(indel_counter_list):
            row = [""]
            genotype = indel_counter.get_genotype()
            sorted_count_map_list = indel_counter.get_sorted_count_map_list()
            len_without_err = indel_counter.get_len(with_err=False) + 0.000001
            len_total = len(indel_counter) + 0.000001
            if i == 0:
                row[0] = indel_counter.file_name
            if i == 1:
                row[0] = f"{total_lines_for_file} lines"
            row += [indel_counter.ref_name,
                    genotype.warning.strip().replace("\n", " + "),
                    genotype.allele_set_text,
                    str(genotype).splitlines()[0].strip(),
                    len(indel_counter),
                    f"err {indel_counter.count_map['err']}({round(indel_counter.count_map['err']/len_total, 3)})",
                    indel_counter.get_len(with_err=False)]

            for key, value in sorted_count_map_list:
                if key == 'err':
                    continue
                row.append(f"{key} {value}({round(value/len_without_err, 3)})")
            file_csv_writer.writerow(row)
        file_csv_writer.writerow([])

    file_csv.close()

    # Write Excel file: CSV is not working well in Excel software
    workbook = Workbook(XLSX_LOG_NAME)
    worksheet = workbook.add_worksheet()
    with open(CSV_LOG_NAME, 'rt', encoding='utf8') as f:
        reader = csv.reader(f)
        for r, row in enumerate(reader):
            for c, col in enumerate(row):
                worksheet.write(r, c, col)
        worksheet.set_column_pixels(0, 0, 255/1.4275)
        worksheet.set_column_pixels(2, 2, 255/1.4275)
        worksheet.set_column_pixels(4, 4, 660/1.4275)
        worksheet.set_column_pixels(6, 6, 127/1.4275)
        worksheet.set_column_pixels(8, 9, 150/1.4275)
        worksheet.set_column_pixels(10, 80, 136/1.4275)
    workbook.close()


