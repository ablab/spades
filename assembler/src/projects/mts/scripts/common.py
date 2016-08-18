from __future__ import print_function
try:
    from future_builtins import zip
except:
    pass

import os
import os.path
try:
    import yaml
    def load_dict(input):
        return yaml.load(input)
    def dump_dict(dict, output):
        yaml.dump(dict, output)
except:
    def load_dict(input):
        def load_pairs():
            for line in input:
                params = line.split(":", 2)
                yield (params[0].strip(), params[1].strip())
        return dict(load_pairs())
    def dump_dict(dict, output):
        for k, v in dict.items():
            print(k, ": ", v, sep="", file=output)

FASTA_EXTS = {".fasta", ".fa", ".fna", ".fsa", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fna.gz"}
def gather_paths(path, basename=False):
    for filename in os.listdir(path):
        name = os.path.basename(filename)
        for ext in FASTA_EXTS:
            if not name.endswith(ext):
                continue
            filepath = os.path.join(path, filename)
            if basename:
                yield (name[0:-len(ext)], filepath)
            else:
                yield filepath

def detect_reads(dir):
    return sorted(list(gather_paths(dir)))[:2]

#Autodetect references
def gather_refs(data):
    if type(data) is list:
        for path in data:
            for ref in gather_refs(path):
                yield ref
    else:
        if data.startswith("@"):
            with open(data[1:]) as input:
                for ref in load_dict(input).items():
                    yield ref
        elif os.path.isdir(data):
            for ref in gather_paths(data, True):
                yield ref
        else:
            yield (os.path.splitext(os.path.basename(data))[0], data)

def get_id(internal_id, sample):
    res = internal_id.split("_", 2)[1]
    return sample + "-" + res

def load_annotation(file, normalize=True):
    res = dict()
    sample, _ = os.path.splitext(os.path.basename(file))
    with open(file) as input:
        for line in input:
            info = line.split(" : ")
            id = get_id(info[0], sample) if normalize else info[0]
            bins = info[1].split()
            if id in res:
                res[id].update(bins)
            else:
                res[id] = set(bins)
    return res

class Row:
    def __init__(self, data, colnames):
        self.data = data
        self.colnames = colnames

    def __getitem__(self, index):
        return self.data[self.colnames[index]]

class Table:
    def __init__(self):
        self.data = []
        self.colnames = None
        self.rownames = None
        self.rows = 0

    @staticmethod
    def read(filepath, sep="\t", headers=False):
        res = Table()
        with open(filepath) as input:
            for line in input:
                params = line.strip("\n").split(sep)
                if not res.colnames:
                    res.rownames = dict()
                    if headers:
                        res.colnames = dict(zip(params[1:], range(len(params))))
                        continue
                    else:
                        res.colnames = dict((i, i) for i in range(len(params)))
                if headers:
                    res.rownames[params[0]] = res.rows
                    res.data.append(params[1:])
                else:
                    res.rownames[res.rows] = res.rows
                    res.data.append(params)
                res.rows += 1
        return res

    def __getitem__(self, index):
        return Row(self.data[self.rownames[index]], self.colnames)

    def zip_with(self, other, method):
        for rowname, i in self.rownames.items():
            for colname, j in self.colnames.items():
                other_cell = other.data[other.rownames[rowname]][other.colnames[colname]]
                method(rowname, colname, self.data[i][j], other_cell)
