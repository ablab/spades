import os.path
import subprocess
import sys
import gzip


class SequenceDistanceMeter(object):
    def __init__(self, path_to_bowtie, path_to_input, path_to_output, index, args="-a -c"):
        self.safe_initialize_path_to_bowtie(path_to_bowtie)
        self.safe_initialize_input(path_to_input)
        self.output = open(path_to_output, "w")
        self.index = index
        self.args = args
    
    def __del__(self):
        if hasattr(self, "input"):
                if self.input:
                    self.input.close()
        if hasattr(self, "output"):
                if self.output:
                    self.output.close()

    def unify_path(self, path_to_bowtie):
        if os.path.isdir(path_to_bowtie):
            if path_to_bowtie[-1] == "/":
                return path_to_bowtie
            else:
                return path_to_bowtie + "/"
        else:
            if os.path.isfile(path_to_bowtie):
                return self.join_with_delimiter(path_to_bowtie.split("/")[:-1])
            else:
                return None

    def join_with_delimiter(self, words, delimiter="/", start_with_delimiter=True, end_with_delimiter=True):
        if start_with_delimiter:
            string = delimiter
        else:
            string = ""
        for word in words:
            if word:
                if type(word) is int:
                    string = string + str(word) + delimiter
                else:
                    string = string + word + delimiter
        if end_with_delimiter:
            return string
        else:
            return string[:-len(delimiter)]

    def initialize_path_to_bowtie(self, path_to_bowtie):
        self.path_to_bowtie = self.unify_path(path_to_bowtie)
        if not self.path_to_bowtie:
            raise IOError

    def safe_initialize_path_to_bowtie(self, path_to_bowtie):
        try:
            self.initialize_path_to_bowtie(path_to_bowtie)
        except IOError:
            sys.exit("No bowtie found. Terminated.")

    def initialize_input(self, path_to_file, mode="r"):
        if ".gz" in path_to_file:
            openedFile = gzip.open(path_to_file, mode)
        else:
            openedFile = open(path_to_file, mode)
        return openedFile

    def safe_initialize_input(self, path_to_input):
        try:
            self.input = self.initialize_input(path_to_input)
        except IOError:
            sys.exit("No input file found. Terminated")

    def align_sequence(self, sequence):
        full_args = [self.path_to_bowtie + "./bowtie"]
        full_args += self.args.split()
        full_args += ["--suppress", "1,2,3,5,6,7,8", self.index, sequence]
        p = subprocess.Popen(full_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result_string, err = p.communicate()
        if p.returncode:
            raise IOError(err)
        result = []
        for line in result_string.split():
            result.append(int(line))
        return result

    def process_pair(self, pair):
        first_aligns = self.align_sequence(pair[0])
        second_aligns = self.align_sequence(pair[1])
        return self.compute_distances(first_aligns, second_aligns)

    def compute_distances(self, first_aligns, second_aligns):
        distances = []
        for align_from_first in first_aligns:
            for align_from_second in second_aligns:
                distances += [abs(align_from_second - align_from_first)]
        distances.sort()
        return distances

    def process_all_pairs(self):
        for line in self.input:
            pair = line.split()[:2]
            self.output.write(self.join_with_delimiter(self.process_pair(pair), ":", False, False) + "\n")






