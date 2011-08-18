#!/usr/bin/python -O

import argparse
from common import safeOpenFile
from common import fileLinesCount
import numpy
import time
import re
import matplotlib.pyplot

class Comparator:

    def __init__(self, originalReadsFile, alignedReadsFile, correctedReadsFile, failedToAlignReadsFile, reportDir, reportTotal, noSmartGraphics):
        if reportDir[-1] == "/":
            self.reportDir = reportDir
        else:
            self.reportDir = reportDir + "/"

        if reportTotal:
            #overall reads counts (might be slow)
            self.alignedReadsCount = fileLinesCount(alignedReadsFile) / 4
            self.failedToAlignReadsCount = fileLinesCount(failedToAlignReadsFile) / 4
            self.originalReadsCount = fileLinesCount(originalReadsFile) / 4

        self.noSmartGraphics = noSmartGraphics    

        #file handlers
        self.originalReadsFile = safeOpenFile(originalReadsFile, "r", "Comparator error: no file with original reads found. Terminated.")
        self.alignedReadsFile = safeOpenFile(alignedReadsFile, "r", "Comparator error: no file with aligned reads found. Terminated.")
        self.correctedReadsFile = safeOpenFile(correctedReadsFile, "r", "Comparator error: no file with corrected reads found. Terminated.")
        self.failedToAlignReadsFile = safeOpenFile(failedToAlignReadsFile, "r", "Comparator error: no file with reads failed to align found. Terminated.")
        self.reportFile = safeOpenFile(self.reportDir + "report", "w", "Error: file can not be created. Terminated.")

        #constants
        self.readLength = self.defineReadLength()

        #temporary
        self.alignedReadsCache = {}
        self.failedToAlignReadsCache = {}

        #report
        self.notAlignedReadsCount = 0
        self.notFoundInAlignedReadsCount = 0
        self.notFoundInOriginalReadsCount = 0
        self.processedReadsCount = 0
        self.correctedAndTrimmedReadsCount = 0 #corrected and trimmed
        self.correctedOnlyReadsCount = 0
        self.trimmedOnlyReadsCount = 0

        self.wronglyCorrectedNucleotidesQualities = numpy.zeros((127,), dtype=numpy.int32)
        self.wronglyUncorrectedNucleotidesQualities = numpy.zeros((127,), dtype=numpy.int32)
        self.finelyCorrectedNucleotidesQualities = numpy.zeros((127,), dtype=numpy.int32)
        self.shouldNotHaveBeenCorrectedNucleotidesQualities = numpy.zeros((127,), dtype=numpy.int32)

        self.wronglyTrimmedNucleotidesQualities = numpy.zeros((127,), dtype=numpy.int32)
        self.finelyTrimmedNucleotidesQualities = numpy.zeros((127,), dtype=numpy.int32)

    def __del__(self):
        if hasattr(self, "originalReadsFile"):
            if self.originalReadsFile:
                self.originalReadsFile.close()
        if hasattr(self, "alignedReadsFile"):
            if self.alignedReadsFile:
                self.alignedReadsFile.close()
        if hasattr(self, "correctedReadsFile"):
            if self.correctedReadsFile:
                self.correctedReadsFile.close()
        if hasattr(self, "failedToAlignReadsFile"):
            if self.failedToAlignReadsFile:
                self.failedToAlignReadsFile.close()
        if hasattr(self, "reportFile"):
            if self.reportFile:
                self.reportFile.close()

    def compare(self):
        while True:
            correctedRead = self.getNextCorrectedRead()
            if not correctedRead:
                break

            self.processedReadsCount += 1

            idOfCorrectedRead = correctedRead[0]
            isReadActuallyCorrected = correctedRead[2]
            trimmed = correctedRead[3]

            alignedRead, isAligned = self.findAligned(idOfCorrectedRead)

            #if corrected is not actually corrected we can omit searching for original
            if alignedRead:
                if isAligned:
                    if not isReadActuallyCorrected and not trimmed: #i.e. the read is unchanged
                        self.compareCorrectedAndAligned(correctedRead, alignedRead)
                    else:
                        originalRead = self.findOriginal(idOfCorrectedRead)
                        if originalRead:
                            self.compareCorrectedAlignedAndOriginal(correctedRead, alignedRead, originalRead)
                    continue
                else:
                    self.notAlignedReadsCount += 1
                    continue

        self.writeReport()
        return None

    def writeReport(self):
        if hasattr(self, "alignedReadsCount"):
            self.reportFile.write("Total aligned reads count: %s\n" % self.alignedReadsCount)
            self.reportFile.write("Total failed to align reads count: %s\n" % self.failedToAlignReadsCount)
            self.reportFile.write("Total original reads count: %s\n\n" % self.originalReadsCount)

        width = 60

        self.reportFile.write(("Processed reads count: %s" % self.processedReadsCount).ljust(width))
        self.reportFile.write("[how many corrected reads were processed]\n")
        self.reportFile.write(("Corrected and trimmed reads count: %s" % self.correctedAndTrimmedReadsCount).ljust(width))
        self.reportFile.write("[how many reads (among processed) were corrected AND trimmed]\n")
        self.reportFile.write(("Corrected only reads count: %s" % self.correctedOnlyReadsCount).ljust(width))
        self.reportFile.write("[how many reads (among processed) were just corrected (and not trimmed)]\n")
        self.reportFile.write(("Trimmed only reads count: %s" % self.trimmedOnlyReadsCount).ljust(width))
        self.reportFile.write("[how many reads (among processed) were just trimmed (and not corrected)]\n")
        self.reportFile.write(("Not aligned reads count: %s" % self.notAlignedReadsCount).ljust(width))
        self.reportFile.write("[how many reads (among processed) were not found in aligned reads (but were found in failed to align)]\n")
        self.reportFile.write(("Not found in aligned reads count: %s" % self.notFoundInAlignedReadsCount).ljust(width))
        self.reportFile.write("[how many reads (among processed) were not found in aligned reads (neither in failed to align]\n")
        self.reportFile.write(("Not found in original reads count: %s" % self.notFoundInOriginalReadsCount).ljust(width))
        self.reportFile.write("[how many reads (among processed) were not found in original reads]\n\n")

        self.reportFile.write(("Total corrected nucleotides count: %s" % (sum(self.wronglyCorrectedNucleotidesQualities) + sum(self.finelyCorrectedNucleotidesQualities) + sum(self.shouldNotHaveBeenCorrectedNucleotidesQualities))).ljust(width))
        self.reportFile.write("[total amount of corrected nucleotides]\n")
        self.reportFile.write(("Wrongly corrected nucleotides count: %s" % sum(self.wronglyCorrectedNucleotidesQualities)).ljust(width))
        self.reportFile.write("[how many nucleotides were corrected but not the way they should've been]\n")
        self.reportFile.write(("Wrongly uncorrected nucleotides count: %s" % sum(self.wronglyUncorrectedNucleotidesQualities)).ljust(width))
        self.reportFile.write("[how many nucleotides weren't corrected but should've been]\n")
        self.reportFile.write(("Finely corrected nucleotides count: %s" % sum(self.finelyCorrectedNucleotidesQualities)).ljust(width))
        self.reportFile.write("[how many nucleotides were corrected the way they should've been]\n")
        self.reportFile.write(("Should not have been corrected nucleotides count: %s" % sum(self.shouldNotHaveBeenCorrectedNucleotidesQualities)).ljust(width))
        self.reportFile.write("[how many nucleotides were corrected but should've stayed unchanged]\n\n")

        self.reportFile.write(("Trimmed nucleotides count: %s" % (sum(self.finelyTrimmedNucleotidesQualities) + sum(self.wronglyTrimmedNucleotidesQualities))).ljust(width))
        self.reportFile.write("[total amount of trimmed nucleotides]\n")
        self.reportFile.write(("Finely trimmed nucleotides count: %s" % sum(self.finelyTrimmedNucleotidesQualities)).ljust(width))
        self.reportFile.write("[how many nucleotides were trimmed and appeared to be wrong in original read (so they had been trimmed for good)]\n")
        self.reportFile.write(("Wrongly trimmed nucleotides count: %s" % sum(self.wronglyTrimmedNucleotidesQualities)).ljust(width))
        self.reportFile.write("[how many nucleotides were trimmed and appeared to be fine in original read (so they should've stayed untrimmed)]\n")

        if self.noSmartGraphics:
            startNucleotideAsciiCode = 33
        else:
            startNucleotideAsciiCode = 36

        endNucleotideAsciiCode = 126

        self.savePlot(33, endNucleotideAsciiCode, self.wronglyCorrectedNucleotidesQualities, ',', "wronglyCorrectedNucleotidesQualities")
        self.savePlot(33, endNucleotideAsciiCode, self.wronglyUncorrectedNucleotidesQualities, ',', "wronglyUncorrectedNucleotidesQualities")
        self.savePlot(33, endNucleotideAsciiCode, self.finelyCorrectedNucleotidesQualities, ',', "finelyCorrectedNucleotidesQualities")
        self.savePlot(33, endNucleotideAsciiCode, self.shouldNotHaveBeenCorrectedNucleotidesQualities, ',', "shouldNotHaveBeenCorrectedNucleotidesQualities")
        self.savePlot(startNucleotideAsciiCode, endNucleotideAsciiCode, self.wronglyTrimmedNucleotidesQualities, ',', "wronglyTrimmedNucleotidesQualities")
        self.savePlot(startNucleotideAsciiCode, endNucleotideAsciiCode, self.finelyTrimmedNucleotidesQualities, ',', "finelyTrimmedNucleotidesQualities")

    def savePlot(self, startNucleotideAsciiCode, endNucleotideAsciiCode, data, line, name):
        matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(range(startNucleotideAsciiCode, endNucleotideAsciiCode), data[startNucleotideAsciiCode:endNucleotideAsciiCode], line)
        matplotlib.pyplot.xlabel("Ascii Code of Nucleotide Quality")
        matplotlib.pyplot.ylabel("Number of Nucleotides")
        matplotlib.pyplot.title((re.sub(r'(?<=.)([A-Z])', r' \1', name)).title())
        verticalRange = max(data[startNucleotideAsciiCode:endNucleotideAsciiCode])
        matplotlib.pyplot.axis([startNucleotideAsciiCode - 1, endNucleotideAsciiCode + 1, -0.01 * verticalRange, 1.01 * verticalRange])
        matplotlib.pyplot.savefig(self.reportDir + name + ".svg", format='svg')

    def compareCorrectedAlignedAndOriginal(self, correctedRead, alignedRead, originalRead):
        correctedSequence = correctedRead[1]
        isReadActuallyCorrected = correctedRead[2]
        trimmed = correctedRead[3]
        alignedSequence = alignedRead[1]
        originalSequence = originalRead[1]
        qualities = originalRead[2]

        if isReadActuallyCorrected and not trimmed:
            self.correctedOnlyReadsCount += 1
            for i in range(len(correctedSequence)):
                if correctedSequence[i] == alignedSequence[i]:
                    if correctedSequence[i] != alignedSequence[i]:
                        self.shouldNotHaveBeenCorrectedNucleotidesQualities[ord(qualities[i])] += 1
                else:
                    if correctedSequence[i] == alignedSequence[i]:
                        self.finelyCorrectedNucleotidesQualities[ord(qualities[i])] += 1
                    else:
                        self.wronglyCorrectedNucleotidesQualities[ord(qualities[i])] += 1
        elif not isReadActuallyCorrected and trimmed:
            self.trimmedOnlyReadsCount += 1
            for i in range(len(correctedSequence), len(correctedSequence) + trimmed):
                if originalSequence[i] != alignedSequence[i]:
                    self.finelyTrimmedNucleotidesQualities[ord(qualities[i])] += 1
                else:
                    self.wronglyTrimmedNucleotidesQualities[ord(qualities[i])] += 1
        else:
            self.correctedAndTrimmedReadsCount += 1
            for i in range(len(correctedSequence)):
                if originalSequence[i] == alignedSequence[i]:
                    if correctedSequence[i] != alignedSequence[i]:
                        self.shouldNotHaveBeenCorrectedNucleotidesQualities[ord(qualities[i])] += 1
                else:
                    if correctedSequence[i] == alignedSequence[i]:
                        self.finelyCorrectedNucleotidesQualities[ord(qualities[i])] += 1
                    else:
                        self.wronglyCorrectedNucleotidesQualities[ord(qualities[i])] += 1
            for i in range(len(correctedSequence), len(correctedSequence) + trimmed):
                if originalSequence[i] != alignedSequence[i]:
                    self.finelyTrimmedNucleotidesQualities[ord(qualities[i])] += 1
                else:
                    self.wronglyTrimmedNucleotidesQualities[ord(qualities[i])] += 1

    def compareCorrectedAndAligned(self, correctedRead, alignedRead):
        correctedSequence = correctedRead[1]
        alignedSequence = alignedRead[1]
        qualities = correctedRead[4]

        for i in range(len(correctedSequence)):
            if correctedSequence[i] != alignedSequence[i]:
                self.wronglyUncorrectedNucleotidesQualities[ord(qualities[i])] += 1

    def findOriginal(self, readID):
        while True:
            originalRead = self.getNextOriginalRead()
            if originalRead:
                idOfOriginalRead = originalRead[0]

                if idOfOriginalRead != readID:
                    if idOfOriginalRead in self.alignedReadsCache:
                        del self.alignedReadsCache[idOfOriginalRead]
                    elif idOfOriginalRead in self.failedToAlignReadsCache:
                        del self.failedToAlignReadsCache[idOfOriginalRead]
                    continue
                else:
                    return originalRead
            else:
                self.notFoundInOriginalReadsCount += 1
                break
        return None

    def findAligned(self, readID):
        if readID in self.alignedReadsCache:
            return [readID] + self.alignedReadsCache.pop(readID), True
        elif readID in self.failedToAlignReadsCache:
            return [readID] + self.failedToAlignReadsCache.pop(readID), False
        else:
            while True:
                alignedRead = self.getNextAlignedRead()
                if alignedRead:
                    idOfAlignedRead = alignedRead[0]

                    if readID != idOfAlignedRead:
                        self.alignedReadsCache[idOfAlignedRead] = list(alignedRead[1:])

                        failedRead = self.getNextFailedRead()

                        if failedRead:
                            idOfFailedRead = failedRead[0]

                            if readID != idOfFailedRead:
                                self.failedToAlignReadsCache[idOfFailedRead] = list(failedRead[1:])
                                continue
                            else:
                                return failedRead, False
                        else:
                            continue
                    else:
                        return alignedRead, True
                else:
                    self.notFoundInAlignedReadsCount += 1
                    break
        return None, None

    def defineReadLength(self):
        self.originalReadsFile.readline()
        readLength = len(self.originalReadsFile.readline()) - 1 #Because of the newline character
        self.originalReadsFile.seek(0)
        return readLength

    def getNextCorrectedRead(self):
        idString = self.correctedReadsFile.readline()[1:]
        if idString:
            idList = idString.split()

            read = self.correctedReadsFile.readline()[:-1]
            self.correctedReadsFile.seek(2, 1)
            qualities = self.correctedReadsFile.readline()[:-1]

            if len(idList) == 1:
                return idList[0], read, False, None, qualities
            elif len(idList) == 2:
                if idList[1][0] == "c": #begins with "c", i.e. "correct"
                    return idList[0], read, True, None, qualities
                else: #begins with "t", i.e. "trim"
                    trimmed = int(idList[1].split("=")[1])
                    return idList[0], read, False, trimmed, qualities
            elif len(idList) == 3:
                trimmed = int(idList[2].split("=")[1])
                return idList[0], read, True, trimmed, qualities
        else:
            return None

    def getNextOriginalRead(self):
        idString = self.originalReadsFile.readline()[1:]

        if idString:
            read = self.originalReadsFile.readline()[:-1]
            self.originalReadsFile.seek(2, 1)
            qualities = self.originalReadsFile.readline()[:-1]
            return idString[:-1], read, qualities
        else:
            return None

    def getNextAlignedRead(self):
        idString = self.alignedReadsFile.readline()[1:]
        if idString:
            read = self.alignedReadsFile.readline()[:-1]
            self.alignedReadsFile.seek(self.readLength + 3, 1)
            return idString[:-1], read
        else:
            return None

    def getNextFailedRead(self):
        idString = self.failedToAlignReadsFile.readline()[1:]
        if idString:
            read = self.failedToAlignReadsFile.readline()[:-1]
            self.failedToAlignReadsFile.seek(self.readLength + 3, 1)
            return idString[:-1], read
        else:
            return None

argParser = argparse.ArgumentParser(description='Compares reads to evaluate correction quality.')
argParser.add_argument('originalReadsFile', metavar='O', type=str, help='path to the original reads file')
argParser.add_argument('alignedReadsFile', metavar='A', type=str, help='path to the aligned reads file')
argParser.add_argument('correctedReadsFile', metavar='C', type=str, help='path to the corrected reads file')
argParser.add_argument('failedToAlignReadsFile', metavar='F', type=str, help='path to the file containing reads that failed to align')
argParser.add_argument('reportDir', metavar='R', type=str, help='path to the directory for report')
argParser.add_argument('--reportTotals', action='store_true', help='report total number of reads in all input files')
argParser.add_argument('--noSmartGraphics', action='store_true', help='turns off omitting unknown (\'N\') nucleotides from statistics')
args = argParser.parse_args()

t0 = time.clock()
print("Beginning comparison...")
comparator = Comparator(args.originalReadsFile, args.alignedReadsFile, args.correctedReadsFile, args.failedToAlignReadsFile, args.reportDir, args.reportTotals, args.noSmartGraphics)
comparator.compare()
print("Comparison finished. Took: %s." % (time.clock() - t0))