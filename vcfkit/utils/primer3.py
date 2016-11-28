from subprocess import Popen, PIPE
from vcfkit.utils import *
from vcfkit.utils.blastn import blast
import signal
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA as DNA_SET
from Bio.Restriction import AllEnzymes, CommOnly, RestrictionBatch
signal.signal(signal.SIGINT, lambda x,y: sys.exit(0))


class seqprimer:
    """
        An individual primer.
    """
    def __init__(self, primer_values, template, reference, left_primer = True):
        for k,v in primer_values.items():
            setattr(self, k, v)
        self.left = left_primer
        template = str(template)
        if left_primer:
            self.START = template.find(self.SEQUENCE)
        else:
            # Reverse - complement right primer to find its location
            pright_rc = str(Seq(self.SEQUENCE).reverse_complement())
            self.START = template.find(pright_rc)
        self.END = self.START + len(self.SEQUENCE)

        # Blast primer sequence
        b = blast(reference, num_alignments = 10, word_size = 14)
        self.unique_copies = b.check_primer(self.SEQUENCE)

    def unique(self):
        return self.unique_copies == 1

    def calc_position(template_start):
        pass

    def __str__(self):
        return self.SEQUENCE

    def __repr__(self):
        return "<primer:{self.SEQUENCE} - {self.START}-{self.END} ({self.unique_copies})>".format(**locals())


class primer_group:
    """
        A group of primers with a designated purpose.
        Can be a single primer or a pair.

        Used to calculate PCR product, start site, etc. 
    """
    def __init__(self, primer_values, template, reference):
        self.reference = reference
        pleft = {k.replace("PRIMER_LEFT_", ""):v for k,v in primer_values.items() if k.startswith("PRIMER_LEFT")}
        pright = {k.replace("PRIMER_RIGHT_", ""):v for k,v in primer_values.items() if k.startswith("PRIMER_RIGHT")}
        self.primer_left = seqprimer(pleft, template, reference)
        self.primer_right = seqprimer(pright, template, reference, False)
        self.amplicon = template[self.primer_left.START:self.primer_right.END]
        for k,v in primer_values.items():
            if not k.startswith("PRIMER_LEFT") and not k.startswith("PRIMER_RIGHT"):
                setattr(self, k, v)

    def unique_primer_group(self):
        """
            Return a set of filtered primers that are usable.
        """
        return all([self.primer_left.unique(), self.primer_right.unique()])

    def __repr__(self):
        return repr(self.primer_left) + "\t" + repr(self.primer_right)


class primer3:
    """
        A wrapper for primer3
    """
    def __init__(self, reference, method="pcr"):
        self.reference = reference
        self.PRIMER_OPT_SIZE = 20
        self.PRIMER_MIN_SIZE = 18  # Must be set
        self.PRIMER_MAX_SIZE = 20  # Must be set
        self.PRIMER_NUM_RETURN = 5
        self.PRIMER_PRODUCT_SIZE_RANGE = "600-800"
        
        if method == "pcr":
            seq_template_length = self.PRIMER_PRODUCT_SIZE_RANGE.split("-")[1]
            self.seq_template_length = int(seq_template_length)
            self.PRIMER_TASK = "pick_pcr_primers"
            self.generate_pcr_template = True

        # Global default
        thermo_path = "/usr/local/share/primer3_config/"
        self.PRIMER_THERMODYNAMIC_PARAMETERS_PATH = thermo_path

        # Set primer returns:
        if self.PRIMER_TASK == "pick_left_only":
            self.left_or_right = ["PRIMER_LEFT"]
        else:
            self.left_or_right = ["PRIMER_LEFT", "PRIMER_RIGHT"]

    def _generate_record(self):
        # Generates text record ready for input
        # into primer3
        attributes = [x for x in dir(self) if x.upper() == x and not x.startswith("_")]
        values = [str(getattr(self, x)) for x in attributes]
        att_val = zip(attributes, values)
        return '\n'.join(["=".join(x) for x in att_val]) + "\n=\n"

    def fetch_primers(self, sequence_template):
        # Runs primer3 with the generated record.
        self.SEQUENCE_TEMPLATE = sequence_template
        primer3_run = Popen(["primer3_core"], stdin=PIPE, stdout=PIPE)
        record = self._generate_record()
        resp, err = primer3_run.communicate(record)
        resp = resp.strip().split("\n")
        if err:
            exit(message(err))
        p3_results = dict([x.split("=") for x in resp
                               if x.split("=")[0] != ""])
        p3_results = {k: autoconvert(v) for k,v in p3_results.items()}

        if "PRIMER_LEFT_NUM_RETURNED" in p3_results:
            n_primers = p3_results["PRIMER_LEFT_NUM_RETURNED"]
            primer_return = []
            for primer_num in xrange(0, n_primers):
                pn = "_{n}_".format(n=primer_num)
                primer_dataset = {k.replace(pn,"_"):v for k,v in p3_results.items() if pn in k}
                primer_return.append(primer_group(primer_dataset,
                                                  self.SEQUENCE_TEMPLATE,
                                                  self.reference))
            return primer_return
