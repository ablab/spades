#!/usr/bin/env python

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import os
import sys
import argparse
from traceback import print_exc
from datetime import datetime


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", "-o", type=str, help="file name or directory for recursive copyright change", required=True)
    parser.add_argument("--extensions", "-e", type=str, nargs="+", help="file extensions to be modified (all)")
    parser.add_argument('--quite', action='store_true', default=False, help="dry run")
    args = parser.parse_args(argv)
    return args

current_year = str(datetime.now().year)
spades_team_year = '2023'
copyright_str = 'Copyright (c)'
spades_team = 'SPAdes team'
spbu = 'Saint Petersburg State University'
spbau = 'Saint Petersburg Academic University'

spades_team_copyright = '%s %s-%s %s' % (copyright_str, spades_team_year, current_year, spades_team)
spbu_team_copyright = '%s 2015-2022 %s' % (copyright_str, spades_team)
all_rights_reseved = 'All Rights Reserved'
license_msg = 'See file LICENSE for details.'


py_comment = [
    '############################################################################',
    '# ' + spades_team_copyright,
    '# ' + all_rights_reseved,
    '# ' + license_msg,
    '############################################################################',
    '']

cpp_comment = [
    '//***************************************************************************',
    '//* ' + spades_team_copyright,
    '//* ' + all_rights_reseved,
    '//* ' + license_msg,
    '//***************************************************************************',
    '']


def update_year(line, final_year):
    dash_pos = line.find('-')
    if dash_pos == -1:
        year_pos = line.find('20')
        if year_pos == -1:
            return line
        elif line[year_pos:year_pos + 4] == final_year:
            return line
        else:
            return line[:year_pos + 4] + '-' + final_year + line[year_pos + 4:]
    else:
        return line[:dash_pos] + '-' + final_year + line[dash_pos + 5:]


def insert_in_script(filename, only_show, border_line, comment_sign, full_header):
    print(filename)
    if only_show:
        return

    lines = open(filename).readlines()
    if (spades_team_copyright + '\n') in lines:
        # up to date
        return

    has_copyright_header = any(license_msg in line for line in lines)

    modified = open(filename, 'w')

    if not has_copyright_header:
        # simply insert copyright header
        modified.write('\n')
        for com_line in full_header:
            modified.write(com_line + '\n')
        for line in lines:
            modified.write(line)
    else:
        # modify header accordingly
        header_done = False
        spades_team_inserted = False
        spbu_inserted = False
        bound_count = 0

        for line in lines:
            if header_done:
                # write after header as is
                modified.write(line)
            elif border_line in line:
                modified.write(line)
                bound_count += 1
                header_done = bound_count > 1
            elif spades_team in line:
                if bound_count == 0:
                    modified.write(comment_sign + (border_line * 9) + "\n")
                    bound_count += 1
                modified.write(update_year(line, current_year))
                spades_team_inserted = True
            elif copyright_str in line and spbu in line:
                if not spades_team_inserted:
                    if bound_count == 0:
                        modified.write(comment_sign + (border_line * 9) + "\n")
                        bound_count += 1
                    modified.write("%s %s\n" % (comment_sign, spades_team_copyright))
                    spades_team_inserted = True
                modified.write(update_year(line, '2022'))
                spbu_inserted = True
            elif copyright_str in line and spbau in line:
                if not spades_team_inserted:
                    if bound_count == 0:
                        modified.write(comment_sign + (border_line * 9) + "\n")
                        bound_count += 1
                    modified.write("%s %s\n" % (comment_sign, spades_team_copyright))
                    spades_team_inserted = True
                if not spbu_inserted:
                    modified.write("%s %s\n" % (comment_sign, spbu_team_copyright))
                    spbu_inserted = True
                modified.write(line)
            elif all_rights_reseved in line or license_msg in line:
                modified.write(line)
            elif bound_count == 0:
                # write before header as is
                modified.write(line)

    modified.close()


def visit(dirname, fnames, extensions, dryrun):
    for fname in fnames:
        path = os.path.join(dirname, fname)
        if not os.path.isfile(path):
            continue
        ext = os.path.splitext(fname)[1]
        if extensions and ext not in extensions:
            continue
        if (ext in ['.py', '.sh']) or fname.lower().startswith('cmake'):
            insert_in_script(path, dryrun, border_line="########", comment_sign="#", full_header=py_comment)
        elif ext in ['.hpp', '.cpp', '.h', '.c']:
            insert_in_script(path, dryrun, border_line="********", comment_sign="//*", full_header=cpp_comment)


def main():
    args = parse_args(sys.argv[1:])

    if not os.path.exists(args.input):
        print("Error! " + args.input + " does not exist!")
        return -1

    if os.path.isfile(args.input):
        visit('', [args.input], args.extensions, args.quite)
        return 0

    if not os.path.isdir(args.input):
        print("Please provide valid file or dir.")
        return -1

    for dir_tuple in os.walk(args.input):
        visit(dir_tuple[0], dir_tuple[2], args.extensions, args.quite)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
