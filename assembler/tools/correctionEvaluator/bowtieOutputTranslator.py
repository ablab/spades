#!/usr/bin/python -OO

import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from common import safeOpenFile
import time
import argparse

class BowtieOutputTranslator:
    def __init__(self, bowtieOutputFile, alignedReadsFile):
        self.numberRegularExpression = re.compile("\d+")
        self.bowtieOutputFile = safeOpenFile(bowtieOutputFile, "r", "Error: no bowtie output file for translation found. Terminated.")
        self.alignedReadsFile = safeOpenFile(alignedReadsFile, "w", "Error: file can not be created. Terminated.")

    def __del__(self):
        if hasattr(self, "bowtieOutputFile"):
            if self.bowtieOutputFile:
                self.bowtieOutputFile.close()
        if hasattr(self, "alignedReadsFile"):
            if self.alignedReadsFile:
                self.alignedReadsFile.close()

    def applyRulesUsingRegExp(self, string, rules):
        """applies mismatch rules using regular expressions"""
        correctedString = ""
        position = -1
        for rule in rules:
            numberStringMatch = self.numberRegularExpression.match(rule)
            if numberStringMatch:
                numberString = numberStringMatch.group()
                previousPosition = position
                position = int(numberString)
                nucleotide = rule[len(numberString) + 1]
                correctedString += string[previousPosition + 1 : position] + nucleotide
        correctedString += string[position + 1 :]
        return correctedString

    def applyRulesStrictly(self, string, rules):
        """applies mismatch rules using string.split()"""
        correctedString = ""
        position = -1
        for rule in rules:
            numberString = rule.split(":")[0]
            previousPosition = position
            position = int(numberString)
            nucleotide = rule[len(numberString) + 1]
            correctedString += string[previousPosition + 1 : position] + nucleotide
        correctedString += string[position + 1 :]
        return correctedString

    def translate(self, shouldItUseRegexp, shouldItUseRam):
        """translating procedure"""
        if shouldItUseRam:
            output = ""
        for line in self.bowtieOutputFile:
            words = line.split()
            readID = words[0]
            isFromNegativeStrand = (words[1] == "-")
            sequence = words[4]
            phredQualities = words[5]
            if len(words) > 7:
                rules = words[7]
            else:
                rules = None
            if isFromNegativeStrand:
                sequence = sequence[::-1] #reverse string
            if rules:
                if shouldItUseRegexp:
                    sequence = self.applyRulesUsingRegExp(sequence, rules.split(","))
                else:
                    sequence = self.applyRulesStrictly(sequence, rules.split(","))
            if isFromNegativeStrand:
                sequence = str(Seq(sequence, IUPAC.unambiguous_dna).complement())
            if shouldItUseRam:
                output += ("@" + readID + "\n" + sequence + "\n+\n" + phredQualities + "\n")
            else:
                self.alignedReadsFile.write("@" + readID + "\n" + sequence + "\n+\n" + phredQualities + "\n")
        if shouldItUseRam:
            self.alignedReadsFile.write(output)
        return None

argParser = argparse.ArgumentParser(description='Translates bowtie output to .fastq format.')
argParser.add_argument('bowtieOutputFile', metavar='B', type=str, help='path to the bowtie output file')
argParser.add_argument('alignedReadsFile', metavar='A', type=str, help='path to the resulting aligned reads file')
argParser.add_argument('--useRam', action='store_true', help='use ram to make execution a little faster')
argParser.add_argument('--useRegExp', action='store_true', help='use RegExp for parsing')
args = argParser.parse_args()

startTime = time.clock()
print("Beginning translation...")
bowtieTranslator = BowtieOutputTranslator(args.bowtieOutputFile, args.alignedReadsFile)
bowtieTranslator.translate(args.useRegExp, args.useRam)
print("Translation finished. Took: %s." % (time.clock() - startTime))
print("Done.")