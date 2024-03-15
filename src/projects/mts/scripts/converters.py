import re

class ProfileParser:
    def header(self, input):
        pass

    def parse(self, line):
        params = line.split("\t")
        return (params[0], params[1:])

    def run(self, var, file_in, file_out):
        input = open(file_in, "r")
        output = open(file_out, "w")
        self.header(input)
        for line in input:
            (contig, values) = self.parse(line.strip())
            print(contig, *values, sep="\t", file=output)

class ProfileFormatter:
    def header(self, first_line):
        pass

    def format(self, contig, profile):
        print(contig, *profile, sep="\t", file=self.out)

    def run(self, file_in, file_out):
        self.out = open(file_out, "w")
        input = open(file_in, "r")
        first_line = True
        for line in input:
            params = line.strip().split("\t")
            if first_line:
                first_line = False
                self.header(params[1:])
            self.format(params[0], params[1:])

extract_num = re.compile("\d+")

class BinningParser:
    def __init__(self, sep=",", filter=None):
        self.sep = sep
        self.filter = filter

    def parse(self, line):
        annotation_str = line.split(self.sep, 1)
        bin_id = annotation_str[1].strip()
        sample_contig = annotation_str[0].strip()
        return (sample_contig, bin_id)

    def run(self, file_in, file_out):
        input = open(file_in, "r")
        output = open(file_out, "w")
        for line in input:
            sample_contig, bin_id = self.parse(line)
            bin_num = extract_num.findall(bin_id)[0]
            if bin_num != self.filter:
                print(sample_contig, "BIN" + bin_num, sep="\t", file=output)
