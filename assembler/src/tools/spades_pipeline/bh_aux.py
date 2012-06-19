#!/usr/bin/env python

import os
import sys
import re

def verify(expr, message):
    if (not (expr)):
        print "Assertion failed. Message: " + message
        exit(0)


def determine_it_count(tmp_dir, prefix):
    import re

    files = os.listdir(tmp_dir)
    answer = 0;
    for f in files:
        m = re.match(r"^" + prefix + "\.(\d{2})\..*", f)
        if m:
            val = int(m.group(1))
            if (val > answer):
                answer = val
    return "%02d" % answer


def determine_read_files(folder, str_it_count, input_files, num_paired):
    answer = dict()
    answer["paired_reads"] = '"'
    answer["single_reads"] = '"'

    # paired files
    for id, input_file in enumerate(input_files):
        prefix = os.path.basename(input_file) + "." + str_it_count
        full_name = folder + prefix + ".cor.fastq"
        verify(os.path.isfile(full_name), "corrected file not found: " + full_name)
        if id < num_paired:
            answer["paired_reads"] += full_name + '  '
        else:
            answer["single_reads"] += full_name + '  '

    answer["paired_reads"] += '"'
    answer["single_reads"] += '"'

    if answer["paired_reads"] == '""':
        del answer["paired_reads"]
    if answer["single_reads"] == '""':
        del answer["single_reads"]

    return answer

# based on "hammer" function in ./src/tools/datasets.py
def hammer(given_props, output_dir, compress):
    def read_files():
        read_files = ["paired_reads", "single_reads"]
        return read_files

    cmd = []
    dataset_entry = []
    for prop in given_props.iterkeys():
        val = given_props[prop]
        if prop in read_files():
            new_val = '"'
            for oldfile in val[1:-1].strip().split("  "):
                newfile_relpath = os.path.basename(oldfile)
                newfile = output_dir + "/" + newfile_relpath
                cmd += ["cp " + oldfile + " " + output_dir + "/"]
                if compress:
                    cmd += ["gzip -9 -f " + newfile]
                    newfile_relpath += ".gz"
                    newfile += ".gz"
                new_val += newfile_relpath + '  '
            new_val += '"'
            val = new_val
        dataset_entry += [(prop, val)]
        #dataset_entry = map(lambda (a, b): a + "\t" + b, dataset_entry)
    #dataset_entry = reduce(lambda x, y: x + "\n" + y, dataset_entry)
    for c in cmd:
        os.system(c)
    return dataset_entry


def dataset_print(dataset):
    result = ""
    dataset_dict = dict(dataset)
    for key, value in dataset_dict.iteritems():
        result += key + "\t" + value + "\n"
    return result


def get_max_prefix(filename1, filename2):
    str1 = os.path.basename(filename1)
    str2 = os.path.basename(filename2)
    prefix = ""
    for i in range(min(len(str1), len(str2))):
        if str1[i] == str2[i]:
            prefix += str1[i]
        else:
            break
    return prefix


def generate_unpaired_basename(filename1, filename2):
    prefix = get_max_prefix(filename1, filename2)
    return prefix + "unpaired"


def generate_paired_basename(filename1, filename2, i):
    prefix = get_max_prefix(filename1, filename2)
    return prefix + "_paired_" + str(i + 1)


def generate_dataset(cfg):
    tmp_dir = cfg.working_dir
    if len(cfg.single_reads) == 0:
        cfg.single_reads = [generate_unpaired_basename(cfg.paired_reads[0], cfg.paired_reads[1])]
    input_files = cfg.paired_reads + cfg.single_reads

    str_it_count = determine_it_count(tmp_dir, os.path.basename(input_files[0]))

    dataset_cfg = determine_read_files(tmp_dir + r"/", str_it_count, input_files, len(cfg.paired_reads))

    import process_cfg

    dataset_cfg["single_cell"] = process_cfg.bool_to_str(cfg.single_cell)
    for key, value in cfg.__dict__.iteritems():
        if key.startswith("original_"):
            dataset_cfg[key] = value
        elif key == "reference_genome":
            dataset_cfg[key] = value

    return dataset_print(hammer(dataset_cfg, cfg.output_dir, cfg.gzip_output))


#### auxiliary function to manage input files 

def split_paired_file(input_filename, output_folder):
    ext = os.path.splitext(input_filename)[1]

    input_file = file
    out_basename = ""

    if ext == '.gz':
        import gzip

        input_file = gzip.open(input_filename, 'r')
        ungzipped = os.path.splitext(input_filename)[0]
        out_basename = os.path.splitext(os.path.basename(ungzipped))[0]
    else:
        input_file = open(input_filename, 'r')
        out_basename = os.path.splitext(os.path.basename(input_filename))[0]

    out_left_filename = os.path.join(output_folder, out_basename + "_1.fastq")
    out_right_filename = os.path.join(output_folder, out_basename + "_2.fastq")

    print("== Splitting " + input_filename + " into left and right reads")

    out_left_file = open(out_left_filename, 'w')
    out_right_file = open(out_right_filename, 'w')
    for id, line in enumerate(input_file):
        if id % 8 < 4:
            out_left_file.write(line)
        else:
            out_right_file.write(line)

    out_left_file.close()
    out_right_file.close()
    input_file.close()

    return [out_left_filename, out_right_filename]


def merge_paired_files(src_paired_reads, dst_paired_reads, output_folder):
    merged = []

    for i in [0, 1]:
        dst_basename = generate_paired_basename(dst_paired_reads[i], src_paired_reads[i], i)
        dst_filename = ""
        if dst_paired_reads[i].startswith(output_folder):
            dst_filename = os.path.join(os.path.dirname(dst_paired_reads[i]), dst_basename)
            os.rename(dst_paired_reads[i], dst_filename)
        else:
            import shutil

            dst_filename = os.path.join(output_folder, dst_basename)
            shutil.copy(dst_paired_reads[i], dst_filename)

        merged.append(dst_filename)

        print("== Merging " + src_paired_reads[i] + " and " + dst_paired_reads[
                                                              i] + " into one file with paired reads: " + dst_filename)

        src_file = open(src_paired_reads[i], 'r')
        dst_file = open(dst_filename, "a")
        dst_file.write(src_file.read())
        dst_file.close()
        src_file.close()

    return merged


def merge_single_files(src_single_read, dst_single_read, output_folder):
    dst_filename = ""
    if dst_single_read.startswith(output_folder):
        dst_filename = dst_single_read
    else:
        import shutil

        shutil.copy(dst_single_read, output_folder)
        dst_filename = os.path.join(output_folder, os.path.basename(dst_single_read))
    merged_filename = os.path.join(os.path.dirname(dst_filename),
        generate_unpaired_basename(src_single_read, dst_filename))

    print(
    "== Merging " + src_single_read + " and " + dst_single_read + " into one file with unpaired reads: " + merged_filename)

    src_file = open(src_single_read, 'r')
    dst_file = open(dst_filename, "a")
    dst_file.write(src_file.read())
    dst_file.close()
    src_file.close()
    os.rename(dst_filename, merged_filename)

    return merged_filename


def ungzip_if_needed(filename, output_folder):
    file_basename, file_extension = os.path.splitext(filename)
    if file_extension == ".gz":
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        ungzipped_filename = os.path.join(output_folder, os.path.basename(file_basename))
        ungzipped_file = open(ungzipped_filename, 'w')

        print("== Ungzipping " + filename + " into " + ungzipped_filename)

        import subprocess

        subprocess.call(['gunzip', filename, '-c'], stdout=ungzipped_file)
        ungzipped_file.close()
        filename = ungzipped_filename

    return filename

####
