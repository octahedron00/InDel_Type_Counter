from src.line_set import Line_Set
from src.reference import Reference
import src.globals as glv

# Variables for validation
HOMO_RATIO_MIN = 0.8
HETERO_RATIO_MIN = 0.35
THIRD_RATIO_MAX = 0.02
ERR_RATIO_MAX = 0.1
ERR_ONLY_WARNING_RATIO_MIN = 0.7
READ_MIN = 30

# Variables for showing awesome simple result!
MARGIN_FOR_SAMPLE = 5

# Variables for uh... zero division
Z = 0.000000001
DICT_ERR_NT = {
    'A': '!',
    'T': '@',
    'G': '#',
    'C': '$',
    '-': '%',
}
DICT_ERR_NT_REV = {
    '!': 'A',
    '@': 'T',
    '#': 'G',
    '$': 'C',
    '%': '-',
}


class InDel_Counter_for_Genotype:
    file_name = "(file_name not defined)"
    ref_name = ""
    guide_rna_name = ""
    guide_rna_seq = ""
    count_map = {}
    best_example_map = {}

    def __init__(self, reference: Reference):

        # Only init the initial values of variables:
        # The real value will be updated with each 'count' function

        self.ref_name = reference.ref_name
        self.guide_rna_name = reference.guide_rna_name
        self.guide_rna_seq = reference.guide_rna_seq
        self.count_map = dict({'err': 0})
        self.best_example_map = {}

    def __str__(self):
        genotype = self.get_genotype()
        str = f"for {self.ref_name} in {self.file_name}: \n" \
              f"guide_rna: {self.guide_rna_name} ({self.guide_rna_seq})\n" \
              f"\n" \
              f"[Result] \n" \
              f"{genotype}\n" \
              f"\n" \
              f"total {len(self)} (without err: {self.get_len(with_err=False)}) \n"
        sorted_count_tuple_list = self.get_sorted_count_map_list()
        for count_tuple in sorted_count_tuple_list:
            key, value = count_tuple[0], count_tuple[1]
            if key == 'err':
                continue
            str += f"{key}: \t{value} ({round(int(value) / self.get_len(with_err=False), 3)} without err)\n"
        if len(self) > 0:
            str += f"err: \t{self.count_map['err']} ({round(self.count_map['err'] / self.get_len(with_err=True), 3)})\n"
        return str

    def __len__(self, with_err: bool = True):
        length = 0
        for key in self.count_map.keys():
            if with_err or key != 'err':
                length += self.count_map[key]
        return length

    def get_len(self, with_err: bool = True):
        length = 0
        for key in self.count_map.keys():
            if with_err or key != 'err':
                length += self.count_map[key]
        return length

    def set_file_name(self, file_name: str):
        self.file_name = file_name

    def count(self, line_set: Line_Set):
        if line_set.ref_name != self.ref_name:
            return
        if line_set.indel_type in self.count_map.keys():
            self.count_map[line_set.indel_type] += 1
        else:
            self.count_map[line_set.indel_type] = 1

        if line_set.indel_type in self.best_example_map.keys():
            # check the highest score, the highest phred_score, and the largest length.

            if self.best_example_map[line_set.indel_type].score < line_set.score:
                self.best_example_map[line_set.indel_type] = line_set
            elif self.best_example_map[line_set.indel_type].score == line_set.score:

                if self.best_example_map[line_set.indel_type].phred_score < line_set.phred_score:
                    self.best_example_map[line_set.indel_type] = line_set
                elif self.best_example_map[line_set.indel_type].phred_score == line_set.phred_score:

                    if len(self.best_example_map[line_set.indel_type]) < len(line_set):
                        self.best_example_map[line_set.indel_type] = line_set
        else:
            # initial set for the indel type
            self.best_example_map[line_set.indel_type] = line_set

    def get_sorted_count_map_list(self):
        sorted_count_tuple_list = list(sorted(self.count_map.items(), key=lambda k: k[1], reverse=True))
        return sorted_count_tuple_list

    def get_genotype(self):
        sorted_count_tuple_list = self.get_sorted_count_map_list()
        return Genotype(indel_counter=self, sorted_count_tuple_list=sorted_count_tuple_list)

    def get_abstract_text(self, is_html=False):

        genotype = self.get_genotype()
        text = f"for {self.ref_name} in {self.file_name}: \n" \
               f"guide_rna: {self.guide_rna_name} ({self.guide_rna_seq})\n" \
               f"\n" \
               f"[Result] \n" \
               f"{genotype}\n"

        if is_html:
            html_genotype = str(genotype).replace("\n", '<br>')

            if html_genotype[:3] == str('err(err)')[:3]:
                text = f"<h2>for <div class=ref_name>{self.ref_name}</div> in <div class=file_name>{self.file_name}</div>:</h2>" \
                       f"<h2><b class=not_important>{html_genotype}</b></h2>" \
                       f"<div class=not_important>guide_rna: {self.guide_rna_name} ({self.guide_rna_seq})\n"
            else:
                text = f"<h2>for <div class=ref_name>{self.ref_name}</div> in <div class=file_name>{self.file_name}</div>:</h2>" \
                       f"<h2><b class=important>{html_genotype}</b></h2>" \
                       f"<div>guide_rna: {self.guide_rna_name} ({self.guide_rna_seq})\n"

        sorted_best_example_tuple = sorted(self.best_example_map.items(),
                                           key=lambda f: self.count_map[f[0]], reverse=True)

        if len(sorted_best_example_tuple) < 2:
            text += f"err: {self.count_map['err']} ({(self.count_map['err'] / (self.get_len(with_err=True) + Z)):.3f})\n"
            if is_html:
                text += '</div> <here code ends />\n'
            return text

        text += get_simple_example_lines(sorted_best_example_tuple, self, is_html=is_html)
        if is_html:
            text += '</div> <here code ends />\n'
        return text


def get_simple_example_lines(sorted_best_example_tuple: list, self: InDel_Counter_for_Genotype, is_html=False):
    wt_seq = ""
    wt_pos = ""
    text = ""

    sample_set = dict()

    key, line_set = sorted_best_example_tuple[0]
    if key == 'err':
        key, line_set = sorted_best_example_tuple[1]

    wt_set = get_simple_example_line_set(line_set=line_set)
    wt_basic = wt_set["ref_line"]
    pos_basic = wt_set["pos_line"]

    for i, a in enumerate(wt_basic):
        if a == '-':
            continue
        wt_seq += a
        wt_pos += pos_basic[i]

    for key, line_set in sorted_best_example_tuple:
        if key == 'err':
            continue
        sample_set[key] = get_simple_example_line_set(line_set)

    i = 0
    while i < len(wt_seq):
        is_added = False
        for key in sample_set.keys():
            if sample_set[key]["match_line"][i] == '+':
                is_added = True
        if is_added:
            wt_seq = wt_seq[:i] + ' ' + wt_seq[i:]
            wt_pos = wt_pos[:i] + ' ' + wt_pos[i:]
            for key in sample_set.keys():
                if sample_set[key]["match_line"][i] == '+':
                    if is_html:
                        sample_set[key]["read_line"] = sample_set[key]["read_line"][:i] + \
                                                       DICT_ERR_NT[sample_set[key]["read_line"][i]] + \
                                                       sample_set[key]["read_line"][i + 1:]
                    continue
                sample_set[key]["ref_line"] = sample_set[key]["ref_line"][:i] + ' ' + sample_set[key]["ref_line"][
                                                                                      i:]
                sample_set[key]["read_line"] = sample_set[key]["read_line"][:i] + ' ' + sample_set[key][
                                                                                            "read_line"][i:]
                sample_set[key]["match_line"] = sample_set[key]["match_line"][:i] + ' ' + sample_set[key][
                                                                                              "match_line"][i:]
                sample_set[key]["pos_line"] = sample_set[key]["pos_line"][:i] + ' ' + sample_set[key]["pos_line"][
                                                                                      i:]
        else:
            for key in sample_set.keys():
                if sample_set[key]["match_line"][i] == '|':
                    continue
                sample_set[key]["read_line"] = sample_set[key]["read_line"][:i] + \
                                               DICT_ERR_NT[sample_set[key]["read_line"][i]] + \
                                               sample_set[key]["read_line"][i + 1:]
        i += 1

    for key in sample_set.keys():
        for nt_key in DICT_ERR_NT_REV:
            if is_html:
                sample_set[key]['read_line'] = \
                    sample_set[key]['read_line'].replace(nt_key, f"<div class=point>{DICT_ERR_NT_REV[nt_key]}</div>")
            else:
                sample_set[key]['read_line'] = \
                    sample_set[key]['read_line'].replace(nt_key, DICT_ERR_NT_REV[nt_key].lower())

    text += f"\n" \
            f"{wt_pos}\n" \
            f"{wt_seq} : Reference      {self.get_len(with_err=False)}: Total without err\n" \
            f"\n"
    for key in sample_set.keys():
        text += f"{sample_set[key]['read_line']} : {key:<8}" \
                f"({self.count_map[key]:>5}/{self.get_len(with_err=False):}, " \
                f"{self.count_map[key] / (self.get_len(with_err=False) + Z):.3f})\n"

    text += f"\n"
    if is_html:
        text += f"<div class=point>{'-' * len(wt_seq)}</div> : err     " \
                f"({self.count_map['err']:>5}/{self.get_len(with_err=True)}, " \
                f"{self.count_map['err'] / (self.get_len(with_err=True) + Z):.3f})\n"
    else:
        text += f"{'-' * len(wt_seq)} : err     ({self.count_map['err']:>5}/{self.get_len(with_err=True)}, " \
                f"{self.count_map['err'] / (self.get_len(with_err=True) + Z):.3f})\n"
    return text


def get_simple_example_line_set(line_set: Line_Set):
    std_pos = line_set.std_pos
    rna_pos = line_set.rna_pos

    pam_len = 3
    rna_len = len(line_set.guide_rna_seq)

    ins_up = 0
    for i in range(std_pos, rna_pos - rna_len - MARGIN_FOR_SAMPLE, -1):
        while line_set.ref_line[i - ins_up] == '-':
            # no worries of 'out of range'
            ins_up += 1
    ins_down = 0
    for i in range(std_pos + 1, std_pos + pam_len + MARGIN_FOR_SAMPLE, 1):
        while line_set.ref_line[i + ins_down] == '-':
            ins_down += 1
    # print(line_set.read_name)
    # print(ins_up, ins_down)
    # print(line_set.ref_line[(std_pos - rna_len - ins_up - MARGIN_FOR_SAMPLE):
    #                         (std_pos + pam_len + ins_down + MARGIN_FOR_SAMPLE)])
    # print(line_set.read_line[(std_pos - rna_len - ins_up - MARGIN_FOR_SAMPLE):
    #                          (std_pos + pam_len + ins_down + MARGIN_FOR_SAMPLE)])

    simple_ref_line = line_set.ref_line[(std_pos - rna_len - ins_up - MARGIN_FOR_SAMPLE):
                                        (std_pos + pam_len + ins_down + MARGIN_FOR_SAMPLE)]

    simple_read_line = line_set.read_line[(std_pos - rna_len - ins_up - MARGIN_FOR_SAMPLE):
                                          (std_pos + pam_len + ins_down + MARGIN_FOR_SAMPLE)]

    simple_match_line = line_set.match_line[(std_pos - rna_len - ins_up - MARGIN_FOR_SAMPLE):
                                            (std_pos + pam_len + ins_down + MARGIN_FOR_SAMPLE)]

    simple_pos_line = line_set.pos_line[(std_pos - rna_len - ins_up - MARGIN_FOR_SAMPLE):
                                        (std_pos + pam_len + ins_down + MARGIN_FOR_SAMPLE)]
    #
    # if glv.DEBUG:
    #     print(line_set.read_name)
    #     print(std_pos, rna_pos, rna_len, ins_up, ins_down)
    #     print(simple_pos_line)
    #     print(simple_ref_line)
    #     print(simple_match_line)
    #     print(simple_read_line)
    #     print()

    return {"ref_line": simple_ref_line,
            "read_line": simple_read_line,
            "match_line": simple_match_line,
            "pos_line": simple_pos_line}


class Genotype:
    name = ""
    warning = ""
    allele1_name = ""
    allele2_name = ""
    allele1_ratio = 0
    allele2_ratio = 0
    allele_set_text = ""
    allele_set_shape = ""

    def __init__(self, indel_counter: InDel_Counter_for_Genotype, sorted_count_tuple_list: list):

        if indel_counter.get_len(with_err=False) < READ_MIN:
            self.append_warning("Reads not enough")

        if (indel_counter.count_map['err'] / (len(indel_counter) + Z)) > ERR_ONLY_WARNING_RATIO_MIN:
            self.append_warning("Error ratio too high")

        # now, set the main genotype
        if len(sorted_count_tuple_list) < 2:
            self.name = "err"
            self.allele1_name = self.allele2_name = 'err'
            self.append_warning("Error only")

        elif len(sorted_count_tuple_list) < 3:
            self.name = "homo"
            for key, value in sorted_count_tuple_list:
                if key != 'err':
                    self.allele1_name = self.allele2_name = key
                    self.allele1_ratio = 1.000

        else:
            for key, value in sorted_count_tuple_list:
                if key != 'err':
                    if self.allele1_ratio == 0:
                        self.allele1_name = key
                        self.allele1_ratio = round(value / indel_counter.get_len(with_err=False), 3)
                    elif self.allele2_ratio == 0:
                        self.allele2_name = key
                        self.allele2_ratio = round(value / indel_counter.get_len(with_err=False), 3)
                    else:
                        break
            if self.allele1_ratio > HOMO_RATIO_MIN:
                self.name = "homo"
            elif self.allele2_ratio > HETERO_RATIO_MIN:
                self.name = "hetero"
            else:
                self.name = "ambiguous"
                self.append_warning("Genotype is ambiguous")

        # setting the allele set
        if self.name == 'err':
            self.allele_set_text = 'err/err'
            self.allele_set_shape = 'err'

        if self.name == 'homo':
            if self.allele1_name == 'WT':
                self.allele_set_shape = '+/+'
            else:
                self.allele_set_shape = '-/-'
            self.allele_set_text = self.allele1_name + "/" + self.allele1_name

        if self.name == 'ambiguous':
            self.allele_set_text = self.allele1_name + "/" + self.allele2_name
            self.allele_set_shape = 'err'

        if self.name == 'hetero':
            if 'WT' in (self.allele1_name, self.allele2_name):
                self.allele_set_shape = '-/+'
            else:
                self.allele_set_shape = '1/2'
            self.allele_set_text = self.allele1_name + "/" + self.allele2_name

        # appending warning messages
        # err ratio
        if (indel_counter.count_map['err'] / (len(indel_counter) + Z)) > ERR_RATIO_MAX:
            self.append_warning("Too many err")

        # other allele ratio test
        if self.name == 'homo' and len(sorted_count_tuple_list) > 2:
            pos = 1
            for i in range(2):
                key, value = sorted_count_tuple_list[i]
                if key == 'err':
                    pos = 2
            key, value = sorted_count_tuple_list[pos]
            if value / indel_counter.get_len(with_err=False) > THIRD_RATIO_MAX:
                self.append_warning("Other allele ratio is too large")
        if self.name == 'hetero' and len(sorted_count_tuple_list) > 3:
            pos = 2
            for i in range(3):
                key, value = sorted_count_tuple_list[i]
                if key == 'err':
                    pos = 3
            key, value = sorted_count_tuple_list[pos]
            if value / indel_counter.get_len(with_err=False) > THIRD_RATIO_MAX:
                self.append_warning("Other allele ratio is too large")

    def __str__(self):
        string = ""
        if self.name in ('hetero', 'ambiguous'):
            string = f"{self.name}({self.allele_set_shape}) of " \
                     f"{self.allele1_name}({self.allele1_ratio:.3f} without err) and " \
                     f"{self.allele2_name}({self.allele2_ratio:.3f} without err) " \
                     f"(sum: {(self.allele1_ratio + self.allele2_ratio):.3f})"
        else:
            string = f"{self.name}({self.allele_set_shape}) of " \
                     f"{self.allele1_name}({self.allele1_ratio:.3f} without err)"
        if len(self.warning) > 0:
            string += "\n"
            string += self.warning
        return string

    def append_warning(self, new_warning: str):
        if len(self.warning) > 0:
            self.warning += "\n"
        self.warning += new_warning
