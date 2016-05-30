from os import path

def get_id(internal_id, sample):
    res = internal_id.split("_", 2)[1]
    return sample + "-" + res

def load_annotation(file):
    res = dict()
    sample, _, _ = path.basename(file).partition(".")
    with open(file) as input:
        for line in input:
            info = line.split(" : ")
            id = get_id(info[0], sample)
            bins = info[1].split()
            if id in res:
                res[id].update(bins)
            else:
                res[id] = set(bins)
    return res
