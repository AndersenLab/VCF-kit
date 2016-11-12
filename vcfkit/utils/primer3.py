from subprocess import Popen, PIPE
from vcfkit.utils import *



class primer3:
    """
        A wrapper for primer3
    """
    def __init__(self, method="pcr"):
        if method == "pcr":
            self.PRIMER_PRODUCT_SIZE_RANGE = "600-800"  # Must be set
            self.PRIMER_OPT_SIZE = 20
            self.PRIMER_MIN_SIZE = 18  # Must be set
            self.PRIMER_MAX_SIZE = 20  # Must be set
            self.PRIMER_NUM_RETURN = 5
            seq_template_length = self.PRIMER_PRODUCT_SIZE_RANGE.split("-")[1]
            self.seq_template_length = int(seq_template_length)
            self.PRIMER_TASK = "pick_pcr_primers"
            self.generate_pcr_template = True
        # Global default
        thermo_path = "/usr/local/Cellar/primer3/2.3.7/share/primer3/primer3_config/"
        self.PRIMER_THERMODYNAMIC_PARAMETERS_PATH = thermo_path

        # Set primer returns:
        if self.PRIMER_TASK == "pick_left_only":
            self.left_or_right = ["PRIMER_LEFT"]
        else:
            self.left_or_right = ["PRIMER_LEFT", "PRIMER_RIGHT"]

    def generate_record(self):
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
        record = self.generate_record()
        resp, err = primer3_run.communicate(record)
        resp = resp.strip().split("\n")
        if err:
            raise Exception(err)
        primer3_results = dict([x.split("=") for x in resp
                               if x.split("=")[0] != ""])
        self.results = {}
        for k, v in primer3_results.items():
            self.results[k] = autoconvert(v)
        return self.primers()

    def primers(self):
        # Iterates through primers
        if "PRIMER_LEFT_NUM_RETURNED" in self.results:
            n_primers = self.results["PRIMER_LEFT_NUM_RETURNED"]
            for primer_num in xrange(0, n_primers):
                primer_return = {}
                for side in self.left_or_right:
                    primer_prefix = side + "_" + str(primer_num)
                    if primer_prefix not in self.results:
                        break
                    primer_num = str(primer_num)
                    primer_return[side] = {}
                    for key, val in [(k.replace(primer_prefix + "_",""), v) 
                        for k, v in self.results.items() 
                            if k.startswith(primer_prefix)]:
                                if key == primer_prefix:
                                    start, length = val.split(",")
                                    primer_return[side]["PRIMER_START"] = int(start)
                                    primer_return[side]["PRIMER_LENGTH"] = int(length)
                                else:
                                    primer_return[side][key] = autoconvert(val)

                # Add PCR Template
                if self.generate_pcr_template:
                    primer_left_pos = self.results["PRIMER_LEFT_" + primer_num]
                    primer_right_pos = self.results["PRIMER_RIGHT_" + primer_num]
                    left_start, left_length = map(int, primer_left_pos.split(","))
                    right_start, right_length = map(int, primer_right_pos.split(","))
                    primer_return["PCR_PRODUCT_SEQ"] = self.results["SEQUENCE_TEMPLATE"][left_start:right_start+1]
                    primer_return["PCR_PRODUCT_LEN"] = len(primer_return["PCR_PRODUCT_SEQ"])
                yield primer_return
