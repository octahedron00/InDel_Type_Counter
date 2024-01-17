from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO

# Variables for uh... making the results:
# Can be changed by command line interface; All in globals.py
import src.globals as glv

# Variables for letter recognition
MATCH_LETTER = ('|', 'x', '+', '\\', '/', '.')
MATCH_ERR_LETTER = MATCH_LETTER[1:]

# Variables for validation
HOMO_RATIO_MIN = 0.8
HETERO_RATIO_MIN = 0.35
THIRD_RATIO_MAX = 0.02
ERR_RATIO_MAX = 0.1
READ_MIN = 30


class Reference:
    '''
    class Reference
    a set of reference for alignment and indel validation.

    The class has:
        a reference sequence to align,
        a guide_RNA (spacer) sequence to validate the main indels,
        a name for both sequence to show.
    '''

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
        ref_seq_for_print = self.ref_seq
        if len(ref_seq_for_print) < 30:
            ref_seq_for_print = ref_seq_for_print[:20] + "..." + ref_seq_for_print[-5:]

        return f"<Class Reference>\n" \
               f"ref_raw: {self.ref_raw}\n" \
               f"ref_seq: {ref_seq_for_print}\n" \
               f"ref_name: {self.ref_name}\n" \
               f"guide_rna_raw: {self.guide_rna_raw}\n" \
               f"guide_rna_seq: {self.guide_rna_seq}\n" \
               f"guide_rna_name: {self.guide_rna_name}\n"


def get_guide_rna_seq_position(ref_line: str, guide_rna_seq: str):
    pri = -1
    pre = -1
    if len(guide_rna_seq) > 0:
        for i in range(len(ref_line) - len(guide_rna_seq)):
            k, is_good = 0, False
            pri = i
            for j, n in enumerate(guide_rna_seq):
                while ref_line[i + j + k] == '-':
                    k += 1
                if ref_line[i + j + k] != n:
                    break
                if j == (len(guide_rna_seq) - 1):
                    is_good = True
                    pre = (i + j + k + 1)
            if is_good:
                break
    return pri, pre


class Line_Set:
    '''
    class Line_Set
    a set of line(sequence and string) to show the alignment and indel type,
        with a validated indel type.

    The class has:
        name of reference and read(id with number from sequencer), guide_RNA
            to calculate / show / debug,

        ref(erence)_line
        read_line
        match_line
            which are aligned & trimmed to show the mismatch and indel.
        pos(ition)_line
            to show the position of guide_RNA with '>',
            and the possible PAM sequence with '<'.
        phred_line
            aligned to show the signal strength(score) of each NT of read,

        cut_pos(ition)
            a starting point of PAM sequence: starting position to look forward.
        indel_type
        indel_reason
            which are set by calculation,
        indel_length
            for the calculation of the score, like...
            2I2D means two mismatch, but the indel_length is 2

        score
            calculated by the perfect match rate, without the main indel mismatch.
        int_score
            (score-1) * 1000 : <more big 'minus' value means more error>
        phred_score
            average phred score value, 'the quality' of the read.
            all score is used to find the best line_set,
            to check the best reference for a Read to match, and
            to get the best line set for an InDel_Type


    '''

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
    indel_type_pos = 0
    cut_pos = 0

    # # only for the beauty of the log
    indel_same_type_count = 0

    score = 0
    int_score = 0
    phred_score = 0

    def __init__(self, read_raw: SeqIO.SeqRecord, ref_set: Reference):

        self.ref_set = ref_set
        self.ref_name = ref_set.ref_name
        self.guide_rna_seq = ref_set.guide_rna_seq
        self.guide_rna_name = ref_set.guide_rna_name

        self.read_name = read_raw.name

        # try aligning, and get ref, match, seq line here
        self._set_align_line_set(ref_seq=ref_set.ref_seq, read_seq=str(read_raw.seq).upper())

        # get aligned phred line here
        self._set_aligned_phred_line_and_score(read_raw=read_raw)

        # get pos line and cut position (endpoint of guide_RNA / startpoint of PAM)
        self._set_indel_pos_line(guide_rna_seq=ref_set.guide_rna_seq)

        # get indel type from the position... only if guide RNA sequence matches.
        self._set_main_indel(cut_pos=self.cut_pos)

        # check whether the type is error or not... and set align score with ratio of mismatch without main indel
        self._set_score_and_check_err()

    def __len__(self):
        return len(self.ref_line)

    def __str__(self):
        return f"@\tref: {self.ref_name}\t /read: {self.read_name}\t /file: {self.file_name}\n" \
               f"@\tscore: {self.int_score}\t /indel: {self.indel_type}({self.indel_same_type_count})\t /reason: \"{self.indel_reason}\"\n" \
               f"position: {self.pos_line}\n" \
               f"ref_line: {self.ref_line}\n" \
               f"match   : {self.match_line}\n" \
               f"readline: {self.read_line}\n" \
               f"phred   : {self.phred_line}\n"

    # # only for the beauty of the sub-log
    def get_str_simple(self):
        return f"| read: {self.read_name}\t /score: {self.int_score}\n" \
               f"| indel: {self.indel_type}\t /reason: \"{self.indel_reason}\"\n" \
               f"| position: {self.pos_line}\n" \
               f"| ref_line: {self.ref_line}\n" \
               f"| match   : {self.match_line}\n" \
               f"| readline: {self.read_line}\n" \
               f"L phred   : {self.phred_line}\n"

    def set_file_name(self, file_name: str):
        self.file_name = file_name

    def _set_align_line_set(self, ref_seq, read_seq):
        # 'X' means that the read must be between the reference
        # = Aligning the subsequence in the sequence
        ref_seq = "X" + ref_seq + "X"

        # align with the pre-set global variables for alignments.
        aligner = PairwiseAligner()
        matrix = substitution_matrices.Array(data=glv.get_align_matrix_for_subsequence_positioning())
        aligner.substitution_matrix = matrix
        aligner.open_gap_score = glv.GAP_OPEN
        aligner.extend_gap_score = glv.GAP_EXTEND
        alignments = aligner.align(ref_seq, read_seq)

        # get the untrimmed alignment result, based on the reference sequence
        ref_line_untrimmed = alignments[0][0]
        read_line_untrimmed = alignments[0][1]

        # trim the sequence, and build the match_line
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

        # set the lines
        self.ref_line = ref_line
        self.match_line = match_line
        self.read_line = read_line

    def _set_aligned_phred_line_and_score(self, read_raw: SeqIO.SeqRecord):
        # get quality score list[int] from SeqRecord: fastq files
        quality_list = read_raw.letter_annotations['phred_quality']
        if len(quality_list) < 1:
            return

        # Build unaligned phred line from the quality score list, use 33-encoding
        phred_line = ""
        phred_score_sum = 0
        for a in quality_list:
            phred_line += str(chr(a + glv.PHRED_ENCODING))
            phred_score_sum += a

        # Align the phred line
        aligned_phred_line = ""
        for c in self.read_line:
            if c == '-':
                aligned_phred_line += " "
            else:
                aligned_phred_line += phred_line[0]
                phred_line = phred_line[1:]

        self.phred_score = phred_score_sum / (len(phred_line) + 0.0000001)
        self.phred_line = aligned_phred_line

    def _set_indel_pos_line(self, guide_rna_seq: str):

        # get the guide_rna_seq alignment, from the function outside
        pri, pre = get_guide_rna_seq_position(self.ref_line, self.guide_rna_seq)

        # if guide RNA is not aligned: the sequence will be considered as an error,
        # since no indel is close the cut_pos(ition)
        if pri < 0 or pre < 0:
            self.cut_pos = -1000
            return

        # buffer: for finding the PAM sequence
        ref_line_buffer = self.ref_line + "X" * (glv.PAM_RANGE_MAX * 2)

        # build the pos(ition)_line
        # add ' ' to the starting point of guide RNA
        pos_line = " " * pri
        # add '>' at the position of guide RNA
        for i in range(pri, pre):
            if ref_line_buffer[i] == '-':
                pos_line += "-"
            else:
                pos_line += ">"
        # find the first possible PAM sequence by 'NGG' (range: PAM_RANGE_MAX)
        for i in range(glv.PAM_RANGE_MAX):
            if ref_line_buffer[pre + i] == ref_line_buffer[pre + i + 1] == "G":
                pos_line += '<<<'
                break
            else:
                pos_line += ' '
        # cut or add more ' ' to make the both length match
        if len(pos_line) < len(self.ref_line):
            pos_line += " " * (len(self.ref_line) - len(pos_line))
        else:
            pos_line = pos_line[:len(self.ref_line)]

        # set the lines
        self.pos_line = pos_line
        self.cut_pos = pre

    def _set_main_indel(self, cut_pos: int):
        '''
        set main indel, type and reason.

        How it works:
            by using the end position of 'RNA'

        :param cut_pos:
        '''

        ref_line = self.ref_line
        match_line = self.match_line
        phred_line = self.phred_line
        pre = self.cut_pos

        p1 = p2 = -1
        indel_i, indel_d, indel_length = 0, 0, 0

        indel_type = "WT"
        indel_length = 0
        indel_reason = "(WT)"
        indel_pos = -99999

        for i, m in enumerate(match_line + str("X" * (glv.PAM_RANGE_MAX * 2))):
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
                if indel_d == indel_i == indel_length == 1 and ((pre - glv.PAM_RANGE_MAX) < p2) and (p1 < (pre + glv.PAM_RANGE_MAX)):
                    if (ord(phred_line[p1]) - 33) > glv.PHRED_MEANINGFUL_MIN:
                        indel_type = self._get_indel_shape_text(indel_i, indel_d, p2 - pre)
                        indel_pos = p2-pre
                        indel_reason = "Indel position confirmed by guide RNA and signal score"
                        indel_length = indel_length
                elif ((pre - glv.PAM_RANGE_MAX) < p2) and (p1 < (pre + glv.PAM_RANGE_MAX)):
                    indel_type = self._get_indel_shape_text(indel_i, indel_d, p2 - pre)
                    indel_pos = p2-pre
                    indel_reason = "Indel position confirmed by guide RNA sequence"
                    indel_length = indel_length
                    break
                p1 = p2 = -1
                indel_i = indel_d = indel_length = 0

        self.indel_type_pos = indel_pos
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
        for c in self.match_line[glv.ERR_PADDING:-glv.ERR_PADDING]:
            if c == '|':
                a += 1
        score = a / (len(self) - 2 * glv.ERR_PADDING - self.indel_length)

        self.score = score

        if score < (1 - glv.ERR_MAX):
            self.indel_type = 'err'
            self.indel_reason = 'Too many mismatch and error'
        self.int_score = int((score - 1) * 1000)

    def set_indel_same_type_count(self, counter_dict: dict):
        self.indel_same_type_count = counter_dict[self.indel_type]


class InDel_Counter_for_Ref:
    file_name = "(file_name not defined)"
    ref_name = ""
    guide_rna_name = ""
    guide_rna_seq = ""
    count_map = {}
    best_example_map = {}

    def __init__(self, ref_set: Reference):
        self.ref_name = ref_set.ref_name
        self.guide_rna_name = ref_set.guide_rna_name
        self.guide_rna_seq = ref_set.guide_rna_seq
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
        if line_set.ref_name == self.ref_name:
            if line_set.indel_type in self.count_map.keys():
                self.count_map[line_set.indel_type] += 1
            else:
                self.count_map[line_set.indel_type] = 1

            if line_set.indel_type in self.best_example_map.keys():
                # check the highest score, highest phred_score, largest length.
                if self.best_example_map[line_set.indel_type].score < line_set.score:
                    self.best_example_map[line_set.indel_type] = line_set
                elif self.best_example_map[line_set.indel_type].score == line_set.score:
                    if self.best_example_map[line_set.indel_type].phred_score < line_set.phred_score:
                        self.best_example_map[line_set.indel_type] = line_set
                    if self.best_example_map[line_set.indel_type].phred_score == line_set.phred_score:
                        if len(self.best_example_map[line_set.indel_type]) < len(line_set):
                            self.best_example_map[line_set.indel_type] = line_set
            else:
                self.best_example_map[line_set.indel_type] = line_set

    def get_sorted_count_map_list(self):
        sorted_count_tuple_list = list(sorted(self.count_map.items(), key=lambda k: k[1], reverse=True))
        return sorted_count_tuple_list

    def get_genotype(self):
        sorted_count_tuple_list = self.get_sorted_count_map_list()
        return Genotype(indel_counter=self, sorted_count_tuple_list=sorted_count_tuple_list)

    def get_examples_text(self):
        #
        # best_example_map = dict{str : Line_Set}
        sorted_best_example_tuple = sorted(self.best_example_map.items(),
                                           key=lambda f: self.count_map[f[0]], reverse=True)

        example_text = "[Best Examples for each InDel type]\n" \
                       "\n" \
                       "\n"
        for key, line_set in sorted_best_example_tuple:
            if key == 'err':
                continue
            example_text += f"<< {key} ({self.count_map[key]}/{self.get_len(with_err=False)}, " \
                            f"{self.count_map[key]/self.get_len(with_err=False):.3f}, without err) >>\n" \
                            f"{line_set.get_str_simple()}\n" \
                            f"\n"
        return example_text


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

        if len(indel_counter) < READ_MIN:
            if len(self.warning) > 0:
                self.warning += '\n'
                self.warning += "Reads not enough"
            else:
                self.warning = "Reads not enough"
                # TODO: add warning_append function inside the class, rewrite warnings shorter
                # TODO: make the code easy to read
                # TODO: make another function for all these works... maybe?

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
                    if self.allele1_name == "WT":
                        self.allele_set_shape = "+/+"
                    else:
                        self.allele_set_shape = "-/-"
        else:
            is_third_ratio_good = True
            for key, value in sorted_count_tuple_list:
                if key != 'err':
                    if 0 < self.allele2_ratio:
                        if key == 'err':
                            continue
                        if value / indel_counter.get_len(with_err=False) > THIRD_RATIO_MAX:
                            is_third_ratio_good = False
                        break
                    elif 0 < self.allele1_ratio:
                        self.allele2_name = key
                        self.allele2_ratio = round(value / indel_counter.get_len(with_err=False), 3)
                    else:
                        self.allele1_name = key
                        self.allele1_ratio = round(value / indel_counter.get_len(with_err=False), 3)
            if self.allele1_ratio > HOMO_RATIO_MIN:
                self.name = "homo"
                if self.allele2_ratio > THIRD_RATIO_MAX:
                    is_third_ratio_good = False
                if self.allele1_name == "WT":
                    self.allele_set_shape = "+/+"
                    self.allele_set_text = "WT/WT"
                else:
                    self.allele_set_shape = "-/-"
                    self.allele_set_text = self.allele1_name + "/" + self.allele1_name

            elif self.allele2_ratio > HETERO_RATIO_MIN:
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

            if not is_third_ratio_good:
                if len(self.warning) > 0:
                    self.warning += '\n'
                    self.warning += "The ratio of next biggest allele is too large"
                else:
                    self.warning = "The ratio of next biggest allele is too large"

            if (indel_counter.count_map['err'] / len(indel_counter)) > ERR_RATIO_MAX:
                if len(self.warning) > 0:
                    self.warning += '\n'

                    self.warning += "Error ratio is high"
                else:
                    self.warning = "Error ratio is high"

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
