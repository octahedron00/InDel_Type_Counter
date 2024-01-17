import datetime
import csv
from xlsxwriter.workbook import Workbook

from src.line_set import Line_Set
from src.reference import Reference
from src.indel_counter_for_genotype import InDel_Counter_for_Genotype
import src.globals as glv

SUB_LOG_ADDRESS = "./log/"
RESULT_LOG_ADDRESS = "./log/"

MAIN_LOG_NAME = "./count_result.txt"
CSV_LOG_NAME = "./count_result.csv"
XLSX_LOG_NAME = "./count_result.xlsx"

# Variables for uh... zero division
Z = 0.000000001


def write_sub_log(line_set_list: list[Line_Set], indel_counter: InDel_Counter_for_Genotype, file_name: str):
    file_log = open(SUB_LOG_ADDRESS + file_name[:-6] + "---" + indel_counter.ref_name + ".txt", "w")

    file_log.write(f""
                   f"# <InDel_Type_Counter {glv.VERSION} Side Log for {file_name}>\n"
                   f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n"
                   f"# \n"
                   f"# Task Title: {glv.TASK_TITLE}\n"
                   f"# {file_name} as a data / {indel_counter.ref_name} as a reference sequence\n"
                   f"\n"
                   f"{glv.get_text_of_global_variables()}"
                   f"\n"
                   f"{_showing_selected_area_to_text(guide_rna_seq=indel_counter.guide_rna_seq)}"
                   f"\n"
                   f"\n")

    file_log.write(indel_counter.get_simple_example_text())
    file_log.write("\n"
                   "\n"
                   "----------------------\n")

    # file_log.write(indel_counter.get_examples_text())
    # file_log.write("\n"
    #                "\n"
    #                "----------------------\n")

    file_log.write("[Raw Data]\n"
                   "\n"
                   "\n")

    for line_set in line_set_list:
        if line_set.ref_name == indel_counter.ref_name:
            file_log.write(str(line_set))
            file_log.write("\n"
                           "\n")

    file_log.close()


def write_main_log(indel_counter_list_list: list[list[InDel_Counter_for_Genotype]]):

    file_log = open(MAIN_LOG_NAME, "w")

    file_log.write(f"# <InDel_Type_Counter {glv.VERSION} Main Log>\n"
                   f"# Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})\n"
                   f"# \n"
                   f"# Task Title: {glv.TASK_TITLE}\n"
                   f"\n"
                   f"{glv.get_text_of_global_variables()}"
                   f"\n"
                   f"\n")

    for indel_counter_list in indel_counter_list_list:
        file_log.write(f"<{indel_counter_list[0].file_name}>\n"
                       f"\n")
        for indel_counter in indel_counter_list:
            file_log.write(indel_counter.get_simple_example_text())
            file_log.write("\n"
                           "\n"
                           "\n")

    file_log.close()


def write_main_csv_log(indel_counter_list_list: list[list[InDel_Counter_for_Genotype]], ref_set_list: list[Reference]):

    file_csv = open(CSV_LOG_NAME, 'w', newline="")
    file_csv_writer = csv.writer(file_csv)

    file_csv_writer.writerow([f"<InDel_Type_Counter {glv.VERSION} Main Log>"])
    file_csv_writer.writerow(
        [f"Log at {datetime.datetime.now()} (UTC {datetime.datetime.now() - datetime.datetime.utcnow()})"])
    file_csv_writer.writerow(["Task Title", glv.TASK_TITLE])
    file_csv_writer.writerows(glv.get_row_of_global_variables())
    file_csv_writer.writerow([])

    file_csv_writer.writerow(["References:"])
    file_csv_writer.writerow(["Name", "seq", "Guide RNA name", "Guide RNA seq"])
    for ref_set in ref_set_list:
        file_csv_writer.writerow([ref_set.ref_name, ref_set.ref_seq, ref_set.guide_rna_name, ref_set.guide_rna_seq])

    file_csv_writer.writerow([])
    file_csv_writer.writerow(["Results:"])
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
            len_without_err = indel_counter.get_len(with_err=False) + Z
            len_total = len(indel_counter) + Z
            if i == 0:
                row[0] = indel_counter.file_name
            if i == 1:
                row[0] = f"{total_lines_for_file} lines"
            row += [indel_counter.ref_name,
                    genotype.warning.strip().replace("\n", ", "),
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

    #
    # Write Excel file: CSV is not working well in Excel software
    workbook = Workbook(XLSX_LOG_NAME)
    worksheet = workbook.add_worksheet()

    italic_format = workbook.add_format({'italic': True, 'font_color': 'silver'})
    with open(CSV_LOG_NAME, 'rt', encoding='utf8') as f:
        reader = csv.reader(f)
        for r, row in enumerate(reader):
            cell_format = workbook.add_format({})
            if len(row) > 2 and len(str(row[2])) > 8:
                if str(row[2])[:5] in (str("Reads not enough")[:5], str("Error only")[:5]):
                    cell_format = italic_format

            for c, col in enumerate(row):
                worksheet.write(r, c, col)
                if c > 0:
                    worksheet.write(r, c, col, cell_format)
        worksheet.set_column_pixels(0, 0, 255/1.4275)
        worksheet.set_column_pixels(2, 2, 255/1.4275)
        worksheet.set_column_pixels(4, 4, 660/1.4275)
        worksheet.set_column_pixels(6, 6, 127/1.4275)
        worksheet.set_column_pixels(8, 9, 150/1.4275)
        worksheet.set_column_pixels(10, 80, 136/1.4275)
    workbook.close()


def _showing_selected_area_to_text(guide_rna_seq: str):

    pos_line = ""
    ref_line = guide_rna_seq
    selected_area_line = ""

    # pre = position rna ending
    pre = len(ref_line)
    std_pos = len(ref_line)
    cut_pos = len(ref_line) + glv.CUT_POS_FROM_PAM

    ref_line += "NGG_______"

    for i, a in enumerate(ref_line):
        if i < pre:
            pos_line += '>'
        elif std_pos <= i < std_pos+3:
            pos_line += '<'
        else:
            pos_line += ' '

        if (cut_pos - glv.CUT_POS_RADIUS) <= i < cut_pos:
            selected_area_line += '('
        elif cut_pos <= i < (cut_pos + glv.CUT_POS_RADIUS):
            selected_area_line += ')'
        else:
            selected_area_line += ' '

    pos_line += '    '
    ref_line += ' or '
    selected_area_line += '    '

    new_starting_point = len(ref_line)


    ref_line += guide_rna_seq

    pre = len(ref_line)
    std_pos = len(ref_line) + 1
    cut_pos = len(ref_line) + 1 + glv.CUT_POS_FROM_PAM

    ref_line += "_NGG_______"

    for i, a in enumerate(ref_line):
        if i < new_starting_point:
            continue

        if i < pre:
            pos_line += '>'
        elif std_pos <= i < std_pos + 3:
            pos_line += '<'
        else:
            pos_line += ' '

        if (cut_pos - glv.CUT_POS_RADIUS) <= i < cut_pos:
            selected_area_line += '('
        elif cut_pos <= i < (cut_pos + glv.CUT_POS_RADIUS):
            selected_area_line += ')'
        else:
            selected_area_line += ' '

    return f"Selected area for estimated indel position is: \n" \
           f"{pos_line}\n" \
           f"{ref_line}\n" \
           f"{selected_area_line}\n"
