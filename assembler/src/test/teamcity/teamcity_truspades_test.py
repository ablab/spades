import os
import shutil
import sys


def run_test():
    shutil.rmtree('bin', True)
    shutil.rmtree('build', True)
    shutil.rmtree('build_spades', True)
    ecode = os.system('./prepare_cfg')
    if ecode != 0:
        print("Preparing configuration files finished abnormally with exit code " + str(ecode))
        sys.exit(2)
    ecode = os.system('./spades_compile.sh')
    if ecode != 0:
        print("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(3)
    cmd = "./truspades.py --disable-gzip-output --test"
    ecode = os.system(cmd)
    output_dir = "test_dataset_truspades"
    os.system("chmod -R 777 " + output_dir)
    if ecode != 0:
        print("SPAdes finished abnormally with exit code " + str(ecode))
        sys.exit(4)
    log = open(os.path.join(output_dir, "truspades.log")).readlines()
    if log[-1] != "TruSPAdes test passed corectly":
        sys.exit(5)

if __name__ == "__main__":
    run_test()