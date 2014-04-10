#!/usr/bin/env python
import json
import os
import subprocess
import sys

def fileExistsAndNonEmpty(filename):
    if not os.path.exists(filename):
        return False
    return os.stat(filename).st_size > 0

class AssemblerRunner(object):
    def __init__(self, sample_id, sample_seq, bam_file):
        with open("startplugin.json", "r") as fh:
            self.config = json.load(fh)
            self.params = self.config['pluginconfig']

        # launch.sh creates a symlink to the input BAM file in this directory
        self.output_dir = self.config['runinfo']['results_dir']
        self.sample_id = sample_id
        self.sample_seq = sample_seq
        self.sample_name = sample_id + "." + sample_seq
        self.sample_output_dir = os.path.join(self.output_dir, self.sample_name)
        self.bam_file = bam_file
        self.bam_rel_path = os.path.join(self.sample_name, self.bam_file)

        # relative path to the input bam file
        self.bam_to_assemble = os.path.join(self.output_dir, self.bam_rel_path)

        # how much to downsample (the step is skipped if it equals to 1)
        self.fraction_of_reads = float(self.params['fraction_of_reads'])

        # all executables are located in bin/ subdirectory
        self.assembler_path = os.path.join(os.environ['DIRNAME'], 'bin')

        # where to output HTML with results
        self.url_root = self.config['runinfo']['url_root']

        # skip assembly (and run only QUAST) if contigs exist
        self.quast_only = self.params.has_key('quastOnly')

        # information will be printed to "info.json"
        self.info = { 'params' : self.params, 'executedCommands' : [] }
        if sample_id != '' and sample_seq != '':
            self.info['sampleId'] = sample_id
            self.info['sampleSeq'] = sample_seq
            self.info['sampleName'] = self.sample_name

    # Prints 'pluginconfig' section of 'startplugin.json'
    def printAssemblyParameters(self):
        print("AssemblerPlus run parameters:")
        print(self.params)

    def writeInfo(self, json_filename):
        with open(json_filename, 'w+') as f:
            json.dump(self.info, f, indent=4)

    def runCommand(self, command, description=None):
        if description:
            print(description)
        else:
            print(command)
        sys.stdout.flush()
        os.system(command)
        self.info['executedCommands'].append(command)

    def runDownsampling(self):
        print("\nSubsampling using Picard")
        downsampler = os.path.join(self.assembler_path, 'DownsampleSam.jar')
        out = os.path.join(self.sample_output_dir, self.bam_file + "_scaled")

        cmd = ("java -Xmx2g -jar {downsampler} "
               "INPUT={self.bam_to_assemble} OUTPUT={out} "
               "PROBABILITY={self.fraction_of_reads}").format(**locals())
        self.runCommand(cmd)

        cmd = ("mv {out} {self.bam_to_assemble}").format(**locals())
        self.runCommand(cmd)

    def execute(self):
        self.printAssemblyParameters()
        read_count_cmd = "samtools view -c " + self.bam_rel_path
        read_count_process = subprocess.Popen(read_count_cmd, shell=True,
                                              stdout=subprocess.PIPE)
        num_reads = int(read_count_process.communicate()[0])

        def tooFewReads():
            if not self.params.has_key('min_reads'):
                return False
            self.min_reads = int(self.params['min_reads'])
            return num_reads <= self.min_reads

        print("%d reads in %s" % (num_reads, self.bam_file))
        if tooFewReads():
            print(("\tDoes not have more than %d reads. "
                   "Skipping this file") % (self.min_reads,))
            return

        if self.fraction_of_reads < 1:
            self.runDownsampling()

        if self.params.has_key('runMira'):
            self.runMira()

        if self.params.has_key('runSpades'):
            self.runSPAdes()

    def runMira(self):
        version = self.params['miraversion']
        assert(version >= "4.0")

        assembly_type = self.params['type']
        sff_extract_path = os.path.join(self.assembler_path, "sff_extract")
        rel_path = "mira-%s/bin/mira" % version
        mira_path = os.path.join(self.assembler_path, rel_path)
        mira_reference = self.params['agenome'] # FIXME: unused for some reason

        mira_info = { 'type' : assembly_type, 'reference' : mira_reference,
                      'version' : version }

        path_prefix = os.path.splitext(self.bam_to_assemble)[0]
        project_name = os.path.splitext(self.bam_file)[0]
        sff = path_prefix + ".sff"

        fastq = path_prefix + "_in.iontor.fastq"
        xml = path_prefix + "_traceinfo_in.iontor.xml"

        results_dir = os.path.join(self.sample_name, project_name + "_assembly",
                                   project_name + "_d_results")
        contigs_fn = os.path.join(results_dir,
                                  project_name + "_out.unpadded.fasta")
        log_fn = os.path.join(self.sample_output_dir, "mira.log")

        mira_info['contigs'] = contigs_fn
        mira_info['log'] = log_fn
        mira_info['wig'] = os.path.join(results_dir, project_name + "_out.wig")
        mira_info['ace'] = os.path.join(results_dir, project_name + "_out.ace")
        mira_info['info'] = os.path.join(self.sample_name,
                                         project_name + "_assembly",
                                         project_name + "_d_info",
                                         project_name + "_info_assembly.txt")

        skip_assembly = self.quast_only and fileExistsAndNonEmpty(contigs_fn)
        if not skip_assembly:

            cmd = "bam2sff -o {sff} {self.bam_to_assemble}".format(**locals())
            print("Running bam2sff")
            self.runCommand(cmd)

            cmd = "{sff_extract_path} -s {fastq} -x {xml} {sff}"\
                .format(**locals())
            print("Running sff_extract")
            self.runCommand(cmd)

            parameters = " ".join(["-DI:trt=/tmp", "-MI:IONTOR_SETTINGS",
                                   "-AS:mrpc=100"])
            manifest_content = ("""
#MIRA Manifest File

#Settings
project = {project_name}
job = denovo,genome,{assembly_type}
parameters = {parameters}
#Reads
readgroup = {self.sample_name}
data = {fastq} {xml}
technology = iontor

""").format(**locals())
            manifest_fn = os.path.join(self.sample_output_dir, "manifest.txt")
            with open(manifest_fn, "w+") as manifest:
                manifest.write(manifest_content)

            mira_info['manifestFilename'] = manifest_fn

            cmd = ("cd ./{self.sample_name}; "
                   "{mira_path} {manifest_fn} > {log_fn}").format(**locals())
            print("Running AssemblerPlus - Mira %s" % version)
            self.runCommand(cmd)

        output_dir = self.sample_name
        report_dir = self.createQuastReport(contigs_fn, output_dir)
        mira_info['quastReportDir'] = report_dir
        self.info['mira'] = mira_info

    def runSPAdes(self):
        version = self.params['spadesversion']
        assert(version >= "3.0.0")

        rel_path = os.path.join("SPAdes-%s-Linux" % version, "bin", "spades.py")
        spades_path = os.path.join(self.assembler_path, rel_path)

        output_dir = os.path.join(self.sample_name, "spades")
        contigs_fn = os.path.join(output_dir, "contigs.fasta")
        scaffolds_fn = os.path.join(output_dir, "scaffolds.fasta")
        log_fn = os.path.join(output_dir, "spades.log")
        skip_assembly = self.quast_only and fileExistsAndNonEmpty(contigs_fn)
        user_options = self.params['spadesOptions']

        spades_info = {'contigs' : contigs_fn,
                       'scaffolds' : scaffolds_fn,
                       'log' : log_fn,
                       'userOptions' : user_options,
                       'version' : version }

        pid = os.getpid()
        if not skip_assembly:
            cmd = ("{spades_path} --iontorrent --tmp-dir /tmp/{pid} "
                   "-s {self.bam_to_assemble} -o {output_dir} "
                   "{user_options} > /dev/null").format(**locals())
            print("Running AssemblerPlus - SPAdes %s" % version)
            self.runCommand(cmd)

        report_dir = self.createQuastReport(contigs_fn, output_dir)
        spades_info['quastReportDir'] = report_dir
        self.info['spades'] = spades_info

    def createQuastReport(self, contigs_fn, output_dir):
        version = "2.3"
        rel_path = os.path.join("quast-%s" % version, "quast.py")
        quast_path = os.path.join(self.assembler_path, rel_path)

        quast_reference = self.params['bgenome']
        quast_results_dir = os.path.join(output_dir, "quast_results")

        print("Running QUAST %s" % version)
        reference_param = ("-R " + quast_reference) if quast_reference else " "
        cmd = ("{quast_path} -o {quast_results_dir} "
               "{reference_param} {contigs_fn}").format(**locals())
        self.runCommand(cmd)

        try:
            if os.path.isfile(os.path.join(quast_results_dir, "report.html")):
                return os.path.abspath(quast_results_dir)
            else:
                return None
        except:
            return None

import sys
if __name__ == "__main__":
    if len(sys.argv) == 4:
        sample_id = sys.argv[1]
        sample_seq = sys.argv[2]
        bam_file = sys.argv[3]
        runner = AssemblerRunner(sample_id, sample_seq, bam_file)
        runner.execute()
        runner.writeInfo("info_%s.%s.json" % (sample_id, sample_seq))
    else:
        assert(len(sys.argv) == 2) # not a barcode run
        bam_file = sys.argv[1]

        # HACK: sample_name = '.' => essentially vanishes from all paths
        runner = AssemblerRunner('', '', bam_file)
        runner.execute()
        runner.writeInfo("info.json")
