#!/usr/bin/env python
import json
import os
import subprocess

def fileExistsAndNonEmpty(filename):
    if not os.path.exists(filename):
        return False
    return os.stat(filename).st_size > 0

class AssemblerRunner(object):
    def __init__(self, sample_name, bam_file):
        with open("startplugin.json", "r") as fh:
            self.config = json.load(fh)
            self.params = self.config['pluginconfig']

        self.sample_name = sample_name
        self.bam_file = bam_file

        # relative path to the input bam file
        self.bam_rel_path = os.path.join(self.sample_name, self.bam_file)

        # launch.sh creates a symlink to the input BAM file in this directory
        self.output_dir = self.config['runinfo']['results_dir']
        self.bam_to_assemble = os.path.join(self.output_dir, self.bam_rel_path)

        # how much to downsample (the step is skipped if it equals to 1)
        self.fraction_of_reads = float(self.params['fraction_of_reads'])

        # all executables are located in bin/ subdirectory
        self.assembler_path = os.path.join(os.environ['DIRNAME'], 'bin')

        # where to output HTML with results
        self.url_root = self.config['runinfo']['url_root']

        # skip assembly (and run only QUAST) if contigs exist
        self.quast_only = self.params.has_key('quastonly')

    # Prints 'pluginconfig' section of 'startplugin.json'
    def printAssemblyParameters(self):
        print("AssemblerPlus run parameters:")
        print(self.params)

    def runCommand(self, command, description=None):
        if description:
            print(description)
        else:
            print(command)
        os.system(command)
        
    def runDownsampling(self):
        print("\nSubsampling using Picard")
        downsampler_path = os.path.join(self.assembler_path, 'DownsampleSam.jar')
        output_file = os.path.join(self.output_dir, self.sample_name,
                                   self.bam_file + "_scaled")

        cmd = ("java -Xmx2g -jar {downsampler_path} "
               "INPUT={self.bam_to_assemble} OUTPUT={output_file} "
               "PROBABILITY={self.fraction_of_reads}").format(**locals())
        self.runCommand(cmd)

        cmd = ("mv {output_file} {self.bam_to_assemble}").format(**locals())
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

    #TODO
    def runMira(self):
        version = self.params['miraversion']
        assert(version >= "4.0")

        assembly_type = self.params['type']
        sff_extract_path = os.path.join(self.assembler_path, "sff_extract")
        rel_path = "mira-" + version + "/bin/mira"
        mira_path = os.path.join(self.assembler_path, rel_path)
        mira_reference = self.params['agenome']
        pass

    def runSPAdes(self):
        version = self.params['spadesversion']
        assert(version >= "3.0.0")

        rel_path = os.path.join("SPAdes-%s-Linux" % version, "bin", "spades.py")
        spades_path = os.path.join(self.assembler_path, rel_path)

        output_dir = os.path.join(self.sample_name, "spades")
        contigs_fn = os.path.join(output_dir, "contigs.fasta")
        skip_assembly = self.quast_only and fileExistsAndNonEmpty(contigs_fn)
        user_options = self.params['spadesOptions']
        pid = os.getpid()
        if not skip_assembly:
            cmd = ("{spades_path} --iontorrent --tmp-dir /tmp/{pid} "
                   "-s {self.bam_to_assemble} -o {output_dir} "
                   "{user_options} > /dev/null").format(**locals())
            print("Running AssemblerPlus - SPAdes %s" % version)
            self.runCommand(cmd)

        self.createQuastReport(contigs_fn, output_dir,
                               "spades_results.template", "SPAdes_QUAST.html")
        
    #TODO
    def createQuastReport(self, contigs_fn, output_dir,
                          template_filename, output_filename):
        version = "2.3"
        rel_path = os.path.join("quast-%s" % version, "quast.py")
        quast_path = os.path.join(self.assembler_path, rel_path)

        quast_reference = self.params['bgenome']
        if quast_reference == "None":
            return

        quast_results_dir = os.path.join(output_dir, "quast_results")
        cmd = ("{quast_path} -R {quast_reference} {contigs_fn};"
               "rm -rf {quast_results_dir};"
               "mv quast_results {output_dir}").format(**locals())
        print("Running QUAST %s" % version)
        self.runCommand(cmd)

        report = None

        try:
            # after the move latest/ symlink is broken
            link_path = os.path.join(quast_results_dir, "latest")
            latest = os.readlink(link_path)
            path_components = latest.rsplit("/", 2)
            fixed_path = os.path.join(path_components[0], output_dir,
                                      path_components[1], path_components[2])
            # we could leave the symlink broken but fixing it takes just 2 lines
            os.unlink(link_path)
            os.symlink(fixed_path, link_path)
            report = os.path.join(fixed_path, "report.html")
            if not os.path.exists(report):
                report = None
        except:
            pass
        
        #TODO use template_filename and output_filename
        html = ("<p align=left><a target=_blank "
                "href={self.sample_name}/spades/quast_results/latest/report.html"
                "> SPAdes-QUAST report</a></p>").format(**locals())
        with open('SPAdes_QUAST.html', 'w') as f:
            f.write(html)
            
import sys
if __name__ == "__main__":
    sample_name = sys.argv[1]
    bam_file = sys.argv[2]
    runner = AssemblerRunner(sample_name, bam_file)
    runner.execute()
