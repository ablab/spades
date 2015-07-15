############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import os
import sys
sys.path.append("src/spades_pipeline")
import SeqIO

__author__ = 'anton'

class BarcodeWriterFactory:
    def __init__(self, output_folder, barcode_suffix = ""):
        self.output_folder =  output_folder
        self.barcode_suffix = barcode_suffix

    def GenerateBarcodeWriter(self, code):
        return BarcodeWriter(self.output_folder, code, self.barcode_suffix)

class BarcodeWriter:
    def __init__(self, output_folder, barcode_name, barcode_suffix):
        if barcode_suffix != "":
            barcode_suffix = "_" + barcode_suffix
        self.handlers = (open(os.path.join(output_folder, "BC_" + barcode_name + "_R1" + barcode_suffix + ".fastq"), "w"),
                         open(os.path.join(output_folder, "BC_" + barcode_name + "_R2" + barcode_suffix + ".fastq"), "w"))
        self.current_handler = 0
        self.count = 0
        self.code = barcode_name

    def write(self, rec):
        SeqIO.write(rec, self.handlers[self.current_handler], "fastq")
        self.current_handler = 1 - self.current_handler
        self.count += 1


    def close(self):
        if self.count > 1000:
           sys.stdout.write(str(self.code) + " " + str(self.count) + "\n")
        assert(self.current_handler == 0)
        self.handlers[0].close()
        self.handlers[1].close()

class Demultiplexer:
    def __init__(self, factory, bardodes, extract_barcode):
        self.barcodes = barcodes
        self.factory = factory
        self.extract_barcode = extract_barcode
        self.barcode_mapping = dict()

    def write(self, rec):
        barcode = self.extract_barcode(rec)
        if self.barcodes == None or barcode in self.barcodes:
            if barcode not in self.barcode_mapping:
                self.barcode_mapping[barcode] = self.factory.GenerateBarcodeWriter(barcode)
            self.barcode_mapping[barcode].write(rec)

    def close(self):
        for barcode, writer in self.barcode_mapping.iteritems():
            writer.close()

def CollectBarcodes(f, extract_barcode):
    h = SeqIO.parse_fastq(SeqIO.Open(f, "r"))
    cnt = 0
    counter = dict()
    for record in h:
        barcode = extract_barcode(record)
        if barcode not in counter:
            counter[barcode] = 0
        counter[barcode] += 1
        cnt += 1
        if cnt == 100000:
            break
    return [barcode for barcode, count in counter.iteritems() if count >= 100]

def demultiplex(f, barcodes, extract_barcode):
    d = Demultiplexer(BarcodeWriterFactory(sys.argv[2]), barcodes, extract_barcode)
    f = SeqIO.parse_fastq(SeqIO.Open(input_file, "r"))
    cnt = 0
    for record in f:
        cnt += 1
        d.write(record)
    d.close()
    f.close()


if __name__ == '__main__':
    input_file = sys.argv[1]
    if os.path.exists(sys.argv[2]):
        os.rmdir(sys.argv[2])
    os.makedirs(sys.argv[2])
    barcodes = None
    if len(sys.argv) > 3:
        barcodes = [a.split(" ")[0] for a in open(sys.argv[3]).readlines()]
    else:
        barcodes = CollectBarcodes(input_file, lambda b: b.id[-6:])
    demultiplex(input_file, barcodes, lambda b: b.id[-6:])
