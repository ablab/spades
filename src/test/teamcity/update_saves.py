#!/usr/bin/env python3

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

# script for testing SPAdes
# provide a path to .info file, path to spades.py and output folder
import sys
import os
import shutil
import datetime
import argparse
import logging
from traceback import print_exc
from io import StringIO

log = logging.getLogger('Saves')


def set_logger(args, logger_instance):
    output_level = logging.INFO
    logger_instance.setLevel(output_level)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)

    log.info("Running saves updater")


# Check etalon saves using detect_diffs.sh
def check_etalon_saves(etalon_saves, working_dir, output_dir):
    exit_code = 0
    log.info("Comparing etalon saves now")
    ecode = os.system(os.path.join(working_dir, "detect_diffs.sh") + " " + output_dir + " " + etalon_saves)
    return ecode == 0


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--force", default=False,
                        help="do not check saves, copy everything", action='store_true')

    parser.add_argument("--saves_dir", "-s", help="base folder with current etalon saves", type=str, required=True)
    parser.add_argument("--output_dir", "-o", help="folder with SPAdes output", type=str, required=True)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    set_logger(args, log)
    working_dir = os.path.dirname(os.path.join(os.getcwd(), __file__))
    latest_etalon = os.path.join(args.saves_dir, "etalon")

    if not args.force and check_etalon_saves(latest_etalon, working_dir, args.output_dir):
        log.info("Saves are identical, will not be copied")
        return 0

    new_etalon_base = os.path.join(args.saves_dir, "etalon_saves_" + datetime.date.today().strftime('%d.%m.%y'))
    new_etalon = new_etalon_base
    counter = 0
    while os.path.exists(new_etalon):
        new_etalon = new_etalon_base + "_%d" % counter
        counter += 1

    log.info("Copying files to %s" % new_etalon)
    shutil.copytree(args.output_dir, new_etalon)
    log.info("Copying done")
    if os.path.exists(latest_etalon):
        os.remove(latest_etalon)
    os.symlink(new_etalon, latest_etalon, target_is_directory=True)
    log.info("New symlinks created, exiting now")


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except KeyboardInterrupt:
        raise
    except:
        if log.handlers:
            strout = StringIO()
            print_exc(file=strout)
            s = strout.getvalue()
            if s:
                log.critical("SPAdes runner failed with the following error:\n" + s)
            else:
                print_exc()
        else:
            sys.stderr.write("SPAdes runner failed with the following error:")
            print_exc()
        sys.exit(-1)

