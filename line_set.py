from Bio.Align import PairwiseAligner, substitution_matrices
from Bio import SeqIO


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


    def __init__(self, ref_seq: str, ref_name: str, guide_rna_seq: str, guide_rna_name: str):
        self.ref_seq = ref_seq.upper()
        self.ref_name = ref_name
        self.guide_rna_seq = guide_rna_seq.upper()
        self.guide_rna_name = guide_rna_name


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
    seq_line = ""
    phred_line = ""

    guide_rna_seq = ""

    name = ""
    ref_name = ""
    guide_rna_name = ""

    indel_type = ""
    indel_length = 0
    indel_start = -1
    indel_end = -1




    def __init__(self, read_seq: SeqIO.SeqRecord, ref_set: Reference_Set):



