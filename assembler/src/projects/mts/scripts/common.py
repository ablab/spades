from __future__ import print_function
try:
    from future_builtins import zip
except:
    pass

import os
import os.path
import re

default_values = {
    "threads":     16,
    "assembly":    {"assembler": "spades", "k": 55, "groups": []},
    "profile":     {"profiler": "mts", "k": 21, "split": 10000, "min_mult": 2, "max_mult": 65535, "min_samples": 2},
    "binning":     {"binner": "canopy", "contig_length": 2000, "min_samples": 2, "bin_length": 500000, "max_clusters": 400},
    "propagation": {"enabled": True},
    "reassembly":  {"enabled": True}
}

# Taken from http://stackoverflow.com/questions/36831998/how-to-fill-default-parameters-in-yaml-file-using-python
def setdefault_recursively(tgt, default = default_values):
    for k in default:
        if isinstance(default[k], dict): # if the current item is a dict,
            # expand it recursively
            setdefault_recursively(tgt.setdefault(k, {}), default[k])
        else:
            # ... otherwise simply set a default value if it's not set before
            tgt.setdefault(k, default[k])

def fill_default_values(config):
    local_dir = config.get("LOCAL_DIR")
    if local_dir:
        default_values["bin"] = os.path.join(local_dir, "build/release/bin")
        default_values["scripts"] = os.path.join(local_dir, "src/projects/mts/scripts")
        default_values["assembly"]["dir"] = os.path.join(local_dir, "bin")
    setdefault_recursively(config)
    config["reassembly"].setdefault("dir", config["assembly"].get("dir"))

def sample_name(fullname):
    return os.path.splitext(os.path.basename(fullname))[0]

FASTA_EXTS = {".fasta", ".fasta.gz", ".fa", ".fna", ".fsa", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fna.gz"}
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
            yield (sample_name(data), data)

def get_id(internal_id, sample):
    res = internal_id.split("_", 2)[1]
    return sample + "-" + res

id_re = re.compile("\\d+")
split_format = re.compile("^([\w.-]+)_\(\d+_\d+\)$")

def extract_id(name):
    bin_id = None
    params = name.split("-", 1)
    if len(params) > 1:
        bin_id = int(id_re.findall(params[0])[0])
        name = params[1]
    contig_id = int(id_re.findall(name)[0])
    if bin_id is None:
        return contig_id
    else:
        return (bin_id, contig_id)

def load_annotation(file, normalize=True):
    res = dict()
    sample, _ = os.path.splitext(os.path.basename(file))
    with open(file) as input:
        for line in input:
            info = line.split("\t")
            id = get_id(info[0], sample) if normalize else info[0]
            bins = info[1].split()
            if id in res:
                res[id].update(bins)
            else:
                res[id] = set(bins)
    return res

def contig_length(name):
    # Length of contig split
    split = re.search("\((\d+)_(\d+)\)", name)
    if split:
        return int(split.group(2)) - int(split.group(1))
    #Default format
    else:
        return int(name.split("_")[3])

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
