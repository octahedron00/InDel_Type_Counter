from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO

MATCH = MAT = 2
MISMATCH = MIS = -1
GAP_OPEN = -50
GAP_EXTEND = -4

ALIGN_MIN = 50

PHRED_MEANINGFUL_MIN = 30

MATCH_LETTER = ('|', 'x', '+', '\\', '/', '.')
MATCH_ERR_LETTER = MATCH_LETTER[1:]

PAM_MAX = 5
ERR_MAX = 0.05

FILIAL_NO = 1

ERR_PADDING = 1
HOMO_RATIO_MIN = 0.9
HETERO_RATIO_SUM_MIN = 0.8
ERR_RATIO_MAX = 0.1

# 'X' = for both end of main sequence, meaning the subsequence must be between this
ALIGN_MATRIX_FOR_SUBSEQUENCE_POSITION = {
    ('A', 'A'): MAT, ('A', 'T'): MIS, ('A', 'G'): MIS, ('A', 'C'): MIS, ('A', '-'): MIS, ('A', 'X'): -1000,
    ('T', 'A'): MIS, ('T', 'T'): MAT, ('T', 'G'): MIS, ('T', 'C'): MIS, ('T', '-'): MIS, ('T', 'X'): -1000,
    ('G', 'A'): MIS, ('G', 'T'): MIS, ('G', 'G'): MAT, ('G', 'C'): MIS, ('G', '-'): MIS, ('G', 'X'): -1000,
    ('C', 'A'): MIS, ('C', 'T'): MIS, ('C', 'G'): MIS, ('C', 'C'): MAT, ('C', '-'): MIS, ('C', 'X'): -1000,
    ('-', 'A'): MIS, ('-', 'T'): MIS, ('-', 'G'): MIS, ('-', 'C'): MIS, ('-', '-'): MAT, ('-', 'X'): -1000,
    ('X', 'A'): -1000, ('X', 'T'): -1000, ('X', 'G'): -1000, ('X', 'C'): -1000, ('X', '-'): -1000, ('X', 'X'): MAT,
}


class Reference_Set:
    ref_raw = None
    ref_seq = ""
    ref_name = ""
    guide_rna_raw = None
    guide_rna_seq = ""
    guide_rna_name = ""

    def __init__(self, ref_raw: SeqIO.SeqRecord, guide_rna_raw: SeqIO.SeqRecord):
        self.ref_raw = ref_raw
        self.ref_seq = str(ref_raw.seq).upper()
        self.ref_name = ref_raw.name
        self.guide_rna_raw = guide_rna_raw
        self.guide_rna_seq = str(guide_rna_raw.seq).upper()
        self.guide_rna_name = guide_rna_raw.name

    def __len__(self):
        return len(self.ref_seq)

    def __str__(self):
        if len(self.ref_seq) < 30:
            return f"<Class Reference_set>\n" \
                   f"ref_raw: {self.ref_raw}\n" \
                   f"ref_seq: {self.ref_seq}\n" \
                   f"ref_name: {self.ref_name}\n" \
                   f"guide_rna_raw: {self.guide_rna_raw}\n" \
                   f"guide_rna_seq: {self.guide_rna_seq}\n" \
                   f"guide_rna_name: {self.guide_rna_name}\n"

        return f"<Class Reference_set>\n" \
               f"ref_raw: {self.ref_raw}\n" \
               f"ref_seq: {self.ref_seq[:20]}...{self.ref_seq[-5:]}\n" \
               f"guide_rna_raw: {self.guide_rna_raw}\n" \
               f"guide_rna_seq: {self.guide_rna_seq}\n" \
               f"guide_rna_name: {self.guide_rna_name}\n"


class Line_Set:
    ref_set = None
    pos_line = ""
    ref_line = ""
    match_line = ""
    read_line = ""
    phred_line = ""

    guide_rna_seq = ""

    read_name = ""
    ref_name = ""
    guide_rna_name = ""
    file_name = "(file_name not defined)"

    indel_type = ""
    indel_length = 0
    indel_reason = ""
    cut_pos = 0

    score = 0
    int_score = 0

    def __init__(self, read_raw: SeqIO.SeqRecord, ref_set: Reference_Set):

        self.ref_set = ref_set
        self.ref_name = ref_set.ref_name
        self.guide_rna_seq = ref_set.guide_rna_seq
        self.guide_rna_name = ref_set.guide_rna_name

        self.read_name = read_raw.name
        self._set_align_line_set(ref_seq=ref_set.ref_seq, read_seq=str(read_raw.seq).upper())
        self._set_aligned_phred_line(read_raw=read_raw)
        self._set_indel_pos_line(guide_rna_seq=ref_set.guide_rna_seq)
        self._set_main_indel(cut_pos=self.cut_pos)

        self._set_score_and_check_err()

    def __len__(self):
        return len(self.ref_line)

    def __str__(self):
        return f"@\tref: {self.ref_name}\t/read: {self.read_name}\t/file: {self.file_name}\n" \
               f"@\tindel: {self.indel_type}\t/score: {self.int_score}\t/reason: {self.indel_reason}\n" \
               f"position: {self.pos_line}\n" \
               f"ref_line: {self.ref_line}\n" \
               f"match   : {self.match_line}\n" \
               f"readline: {self.read_line}\n" \
               f"phred   : {self.phred_line}\n"

    def set_file_name(self, file_name: str):
        self.file_name = file_name

    def _set_align_line_set(self, ref_seq, read_seq):
        ref_seq = "X" + ref_seq + "X"

        aligner = PairwiseAligner()
        matrix = substitution_matrices.Array(data=ALIGN_MATRIX_FOR_SUBSEQUENCE_POSITION)
        aligner.substitution_matrix = matrix
        aligner.open_gap_score = GAP_OPEN
        aligner.extend_gap_score = GAP_EXTEND

        alignments = aligner.align(ref_seq, read_seq)

        ref_line_untrimmed = alignments[0][0]
        read_line_untrimmed = alignments[0][1]

        ref_line = match_line = read_line = ""
        count_base = 0
        for i, c in enumerate(read_line_untrimmed):
            if count_base >= len(read_seq):
                break
            if len(read_line) == 0 and c == "-":
                continue
            read_line += c
            ref_line += ref_line_untrimmed[i]
            count_base += 1
            if read_line_untrimmed[i] == '-':
                match_line += 'x'
                count_base -= 1
            elif ref_line_untrimmed[i] == '-':
                match_line += '+'
            elif read_line_untrimmed[i] == ref_line_untrimmed[i]:
                match_line += '|'
            else:
                match_line += '.'

        self.ref_line = ref_line
        self.match_line = match_line
        self.read_line = read_line

    def _set_aligned_phred_line(self, read_raw: SeqIO.SeqRecord):

        quality = read_raw.letter_annotations['phred_quality']

        if len(quality) < 1:
            return "# phred score not available"

        phred_line = ""
        for a in quality:
            phred_line += str(chr(a + 33))

        aligned_phred_line = ""

        for c in self.read_line:
            if c == '-':
                aligned_phred_line += " "
            else:
                aligned_phred_line += phred_line[0]
                phred_line = phred_line[1:]

        self.phred_line = aligned_phred_line

    def _set_indel_pos_line(self, guide_rna_seq: str):
        pri = -1
        pre = -1
        if len(self.guide_rna_seq) > 0:
            for i in range(len(self.ref_line) - len(self.guide_rna_seq)):
                k, is_good = 0, False
                pri = i
                for j, n in enumerate(self.guide_rna_seq):
                    while self.ref_line[i + j + k] == '-':
                        k += 1
                    if self.ref_line[i + j + k] != n:
                        break
                    if j == (len(self.guide_rna_seq) - 1):
                        is_good = True
                        pre = (i + j + k + 1)
                if is_good:
                    break

        if pri < 0 or pre < 0:
            self.cut_pos = -1000
            return

        ref_line_buffer = self.ref_line + "X" * (PAM_MAX * 2)

        pos_line = " " * pri
        for i in range(pri, pre):
            if ref_line_buffer[i] == '-':
                pos_line += "-"
            else:
                pos_line += ">"

        for i in range(PAM_MAX):
            if ref_line_buffer[pre + i] == ref_line_buffer[pre + i + 1] == "G":
                pos_line += '<<<'
                break
            else:
                pos_line += ' '

        if len(pos_line) < len(self.ref_line):
            pos_line += " " * (len(self.ref_line) - len(pos_line))
        else:
            pos_line = pos_line[:len(self.ref_line)]

        self.pos_line = pos_line
        self.cut_pos = pre

    def _set_main_indel(self, cut_pos: int):
        ref_line = self.ref_line
        match_line = self.match_line
        phred_line = self.phred_line
        pre = self.cut_pos

        p1 = p2 = -1
        indel_i, indel_d, indel_length = 0, 0, 0

        indel_type = "WT"
        indel_length = 0
        indel_reason = "(WT)"

        for i, m in enumerate(match_line + str("X" * (PAM_MAX * 2))):
            if i < 2:
                continue
            if i >= (len(phred_line)):
                break

            if m in MATCH_ERR_LETTER:
                p2 = i + 1
                indel_length += 1
                if p1 < 0:
                    p1 = i
                if m == 'x':
                    indel_d += 1
                if m == '+':
                    indel_i += 1
                if m == '.':
                    indel_d += 1
                    indel_i += 1

            elif p2 >= 0:
                if indel_d == indel_i == indel_length == 1 and ((pre - PAM_MAX) < p2) and (p1 < (pre + PAM_MAX)):
                    if (ord(phred_line[p1]) - 33) > PHRED_MEANINGFUL_MIN:
                        indel_type = self._get_indel_shape_text(indel_i, indel_d, p2 - pre)
                        indel_reason = "Indel position confirmed by guide RNA and signal score"
                        indel_length = indel_length
                elif ((pre - PAM_MAX) < p2) and (p1 < (pre + PAM_MAX)):
                    indel_type = self._get_indel_shape_text(indel_i, indel_d, p2 - pre)
                    indel_reason = "Indel position confirmed by guide RNA sequence"
                    indel_length = indel_length
                    break
                p1 = p2 = -1
                indel_i = indel_d = indel_length = 0

        self.indel_type = indel_type
        self.indel_length = indel_length
        self.indel_reason = indel_reason

    def _get_indel_shape_text(self, indel_i: int, indel_d: int, pos: int):
        if indel_i < 0:
            return 'err'
        if indel_i == 0:
            if indel_d == 0:
                return 'WT'
            else:
                return f"{indel_d}D{pos}"
        else:
            if indel_d == 0:
                return f"{indel_i}I{pos}"
            else:
                return f"{indel_i}I{indel_d}D{pos}"

    def _set_score_and_check_err(self):
        if len(self) < 5:
            return 0
        a = 0
        for c in self.match_line[ERR_PADDING:-ERR_PADDING]:
            if c == '|':
                a += 1
        score = a / (len(self) - 2 * ERR_PADDING - self.indel_length)

        self.score = score

        if score < (1 - ERR_MAX):
            self.indel_type = 'err'
            self.indel_reason = 'Too many mismatch and error'
        self.int_score = int((score - 1) * 1000)


class InDel_Counter_for_Ref:
    file_name = "(file_name not defined)"
    ref_name = ""
    guide_rna_name = ""
    guide_rna_seq = ""

    def __init__(self, ref_set: Reference_Set):
        self.ref_name = ref_set.ref_name
        self.guide_rna_name = ref_set.guide_rna_name
        self.guide_rna_seq = ref_set.guide_rna_seq
        self.count_map = dict({'err': 0})

    def __str__(self):
        genotype = self.get_genotype()
        str = f"for {self.ref_name} in {self.file_name}: \n" \
              f"guide_rna: {self.guide_rna_name} ({self.guide_rna_seq})\n" \
              f"\n" \
              f"[Result] \n" \
              f"{genotype}\n" \
              f"\n" \
              f"total {len(self)} (without err: {self.get_len(with_err=False)}) \n"
        sorted_count_tuple_list = self._get_sorted_count_map_list()
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
        if line_set.ref_name == self.ref_name:
            if line_set.indel_type in self.count_map.keys():
                self.count_map[line_set.indel_type] += 1
            else:
                self.count_map[line_set.indel_type] = 1

    def _get_sorted_count_map_list(self):
        sorted_count_tuple_list = list(sorted(self.count_map.items(), key=lambda k: k[1], reverse=True))
        return sorted_count_tuple_list

    def get_genotype(self):
        sorted_count_tuple_list = self._get_sorted_count_map_list()
        return Genotype(indel_counter=self, sorted_count_tuple_list=sorted_count_tuple_list)


class Genotype:
    name = ""
    warning = ""
    allele1_name = ""
    allele2_name = ""
    allele1_ratio = 0
    allele2_ratio = 0
    allele_set_text = ""
    allele_set_shape = ""

    def __init__(self, indel_counter: InDel_Counter_for_Ref, sorted_count_tuple_list: list):
        if len(sorted_count_tuple_list) < 2:
            self.name = "err"
            self.warning = "error only"
            self.allele1_name = self.allele2_name = 'err'
            self.allele_set_text = "err/err"
            self.allele_set_shape = "err"

        elif len(sorted_count_tuple_list) < 3:
            self.name = "homo"
            for key, value in sorted_count_tuple_list:
                if key != 'err':
                    self.allele1_name = key
                    self.allele1_ratio = 1
                    self.allele_set_text = key + "/" + key
        else:
            for key, value in sorted_count_tuple_list:
                if key != 'err':
                    if 0 < self.allele1_ratio:
                        self.allele2_name = key
                        self.allele2_ratio = round(value / indel_counter.get_len(with_err=False), 3)
                        break
                    else:
                        self.allele1_name = key
                        self.allele1_ratio = round(value / indel_counter.get_len(with_err=False), 3)
            if self.allele1_ratio > HOMO_RATIO_MIN:
                self.name = "homo"
                if self.allele1_name == "WT":
                    self.allele_set_shape = "+/+"
                    self.allele_set_text = "WT/WT"
                else:
                    self.allele_set_shape = "-/-"
                    self.allele_set_text = self.allele1_name + "/" + self.allele1_name

            elif self.allele1_ratio + self.allele2_ratio > HETERO_RATIO_SUM_MIN:
                self.name = "hetero"
                self.allele_set_text = self.allele1_name + "/" + self.allele2_name
                if "WT" in (self.allele1_name, self.allele2_name):
                    self.allele_set_shape = "-/+"
                else:
                    self.allele_set_shape = "1/2"
            else:
                self.name = "ambiguous"
                self.warning = "No genotype set is dominant enough"
                self.allele_set_text = self.allele1_name + "/" + self.allele2_name
                self.allele_set_shape = "err"
            if (indel_counter.count_map['err'] / len(indel_counter)) > ERR_RATIO_MAX:
                if len(self.warning) > 0:
                    self.warning += '\n'
                    self.warning += "error ratio is high"
                else:
                    self.warning = "error ratio is high"

    def __str__(self):
        if self.name in ('hetero', 'ambiguous'):
            string = f"{self.name}({self.allele_set_shape}) of " \
                     f"{self.allele1_name}({self.allele1_ratio} without err) and " \
                     f"{self.allele2_name}({self.allele2_ratio} without err) " \
                     f"(sum: {round(self.allele1_ratio+self.allele2_ratio, 3)})"
            if len(self.warning) > 0:
                string += "\n"
                string += self.warning
            return string
        else:
            string = f"{self.name}({self.allele_set_shape}) of " \
                     f"{self.allele1_name}({self.allele1_ratio} without err)"
            if len(self.warning) > 0:
                string += "\n"
                string += self.warning
            return string
