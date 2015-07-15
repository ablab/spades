#!/usr/bin/python -O

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import argparse
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches
import matplotlib.pyplot
import matplotlib.lines
import math


def similarThreshold(total):
	if total <= 2:
		return 2
	return total  / 2 + total % 2 



class Settings:
	def __init__(self, maxPos, maxCovPos, maxCov, assembliesNum):
		self.maxPos = maxPos
		self.maxCovPos = maxCovPos
		self.maxCov = maxCov

		#colors
		self.colorMisasembled = [ "#e41a1c" , "#b82525" ]
		self.colorMisasembledSimilar = [ "#ff7500" , "#e09110" ]
		self.colorCorrect = [ "#4daf4a" , "#40cf40" ]
		self.colorCorrectSimilar = [ "#377eb8" , "#576e88" ]

		#width		
		self.assemblyWidth = 800.0
		self.lastMargin = 20.0
		
		#scales
		self.scale = self.assemblyWidth / self.maxPos
		self.plotXScale = self.assemblyWidth / (self.maxCovPos + 1)

		#coverage plot 
		self.plotHeight = 130.0
		self.plotMargin = 40.0
		self.dotLength = self.plotXScale
		
		self.maxLogCov = math.ceil(math.log10(self.maxCov))
		self.plotYScale = self.plotHeight / self.maxLogCov

		ticNum = 5
		rawTicStep = self.maxPos / ticNum
		ticStepLog = math.pow(10, math.floor(math.log10(rawTicStep)))
		xStep = math.floor(rawTicStep / ticStepLog)
		if xStep >= 7:
			xStep = 10
		elif xStep >= 3:
			xStep = 5
		else:
			xStep = 1

		self.xTics = xStep * ticStepLog
		self.genomeAnnotation = "Genome, "
		
		ticStepLog = int(math.pow(10, math.floor(math.log10(self.xTics))))
		if xStep == 5:
			ticStepLog *= 10
			

		if ticStepLog == 1:
			self.genomeAnnotation += "bp"
		elif ticStepLog < 1000:
			self.genomeAnnotation += "x" + str(ticStepLog) + " bp"
		elif ticStepLog == 1000:
			self.genomeAnnotation += "Kbp"
		elif ticStepLog < 1000000:
			self.genomeAnnotation += "x" + str(ticStepLog / 1000) + " Kbp"
		elif ticStepLog == 1000000:
			self.genomeAnnotation += "Mbp"
		elif ticStepLog < 1000000000:
			self.genomeAnnotation += "x" + str(ticStepLog / 1000000) + " Mbp"
		elif ticStepLog == 1000000000:
			self.genomeAnnotation += "Gbp"

		if self.maxLogCov <= 3.0:
			self.yTics = 1.0
		else:
			self.yTics = 2.0
			self.maxLogCov -= int(self.maxLogCov) % 2
			self.plotYScale = self.plotHeight / self.maxLogCov

		self.genomeLength = self.maxPos
		self.genomeAnnotationScale = math.pow(10, math.ceil(math.log10(self.xTics)))

		self.zeroCovStep = -0.2
		self.dotWeight = 0.7

		self.zeroCoverageColor = "blue"
		self.coverageColor = "red"

		#dashed lines
		self.dashLines = True
		self.dashLineWeight = 0.2
		self.ticLength = 6
		self.axisWeight = 0.5

		#assembly display parameters
		self.assemblyStep = 70
		self.similarStep = 7
		self.goodStep = 7
		self.oddStep = [0, 4]
		self.contigHeight = 30
		self.simHeight = 6
		self.xOffset = 85.0

		#names parameters
		self.nameAnnotationXStep = 15
		self.nameAnnotationYStep = 15
		self.xticsStep = 6
		self.yticsStep = self.ticLength + 5
		self.xLabelStep = self.xticsStep + 25
		self.yLabelStep = self.yticsStep + 95

		self.totalHeight = self.plotHeight + self.plotMargin + assembliesNum * self.assemblyStep + self.lastMargin
		self.totalWidth = self.xOffset + self.assemblyWidth


		self.contigEdgeDelta = 3000
		self.minSimilarContig = 10000

		self.minConnectedBlock = 2000
		self.maxBlockGap = 10000
			
		self.drawArcs = False
		self.analyzeSimilar = False



class Alignment:
	def __init__(self, name, start, end, contigStart, rc):
		self.name = name
		self.start = start
		self.end = end
		self.contigPos = contigStart
		self.rc = rc
		
		self.order = 0
		self.similar = False
		self.misasembled = False
		self.color = "#000000"
		self.vPositionDelta = 0


	def Length(self):
		return self.end - self.start

	def Annotation(self):
		return self.name + "\n" +str(self.start) + "-" + str(self.end)


	def Center(self):
		return (self.end + self.start) / 2


	def CompareInexact(self, alignment, settings):
		return abs(alignment.start - self.start) <= settings.contigEdgeDelta and abs(alignment.end - self.end) <= settings.contigEdgeDelta



class Arc:
	def __init__(self, c1, c2):
		self.c1 = c1
		self.c2 = c2



class Contig:
	def __init__(self, name):
		self.name = name
		self.alignments = []
		self.arcs = []



class Assembly:
	def __init__(self, file_name, minVisualizedLength = 0):
		self.name, ext = os.path.splitext(file_name)

		self.blocks = []
		self.contigs = {}

		self.misassembled = set()

		mf = open(file_name + ".mis", 'r')
		for line in mf:
			self.misassembled.add(line.strip().split(' ', 1)[0])

		inf = open(file_name, 'r')
		i = 0

		for line in inf:
			l = line.strip().split(' ')
			
			if  int(l[1]) - int(l[0]) < minVisualizedLength:
				continue

			cid = l[2]
			rc = False
			if l[3] == "-":
				rc = True
			al = Alignment(cid, int(l[0]), int(l[1]), int(l[4]), rc)

			al.order = i % 2
			i += 1
			if cid in self.misassembled:
				al.misasembled = True
			
			self.blocks.append(al)

			if cid not in self.contigs:
				self.contigs[cid] = Contig(cid)

			self.contigs[cid].alignments.append(len(self.blocks) - 1)


	def Find(self, alignment, settings):
		if alignment.Length() < settings.minSimilarContig:
			return -1

		i = 0
		while i < len(self.blocks) and not alignment.CompareInexact(self.blocks[i], settings):
			i += 1

		if i == len(self.blocks):
			return -1
		else:
			return i



	def ApplyColor(self, settings):
		for al in self.blocks:
			al.vPositionDelta += settings.oddStep[al.order]
			if al.misasembled:
				if not al.similar:
					al.color = settings.colorMisasembled[al.order]
				else:
					al.color = settings.colorMisasembledSimilar[al.order]
			else:
				al.vPositionDelta + settings.goodStep
				if not al.similar:
					al.color = settings.colorCorrect[al.order]
				else:
					al.color = settings.colorCorrectSimilar[al.order]


	def DrawArcs(self, settings):
#		print (self.misassembled)
#		print (self.contigs.keys())
		for cid in self.misassembled:
			if cid not in self.contigs:
				continue

			contig = self.contigs[cid]
			sortedBlocks = sorted(contig.alignments, key=lambda x: self.blocks[x].contigPos)

			joinedAlignments = []
			currentStart = 0
			currentCStart = 0
			
			i = 0
			while i < len(sortedBlocks):
				block = sortedBlocks[i]
				
				currentStart = self.blocks[block].start
				currentCStart = self.blocks[block].contigPos

				while i < len(sortedBlocks) - 1 and abs(self.blocks[sortedBlocks[i]].end - self.blocks[sortedBlocks[i + 1]].start) < settings.maxBlockGap and self.blocks[sortedBlocks[i]].rc == self.blocks[sortedBlocks[i + 1]].rc:
					i += 1

				if (self.blocks[sortedBlocks[i]].end - currentStart < settings.minConnectedBlock ):
					i += 1
					continue

				joinedAlignments.append(Alignment("", currentStart, self.blocks[sortedBlocks[i]].end, self.blocks[sortedBlocks[i]].rc, currentCStart))
				i += 1

			i = 0
			while i < len(joinedAlignments) - 1:
				contig.arcs.append(Arc(joinedAlignments[i].Center(), joinedAlignments[i + 1].Center()))
				i += 1



class Assemblies:
	def __init__(self, file_list, max_pos, minVisualizedLength = 0):
		self.assemblies = []
		self.max_pos =  max_pos

		inf = open(file_list, 'r')
		
		for line in inf:
			if not inf == "":
				self.assemblies.append(Assembly(line.strip(), minVisualizedLength))

		inf.close()


	def FindSimilar(self, settings):
		for i in range(0, len(self.assemblies)):
#			print("processing assembly " + str(i))
			order = 0
			for block_num in range(0, len(self.assemblies[i].blocks)):
				al = self.assemblies[i].blocks[block_num]
				if al.similar:
					order = (al.order + 1) % 2
					continue

				total = 0
				sim_block_ids_within_asm = [-1 for jj in range(0, len(self.assemblies))]
				sim_block_ids_within_asm[i] = block_num

				for j in range(0, len(self.assemblies)):
					if i == j:
						continue

					block_id = self.assemblies[j].Find(al, settings)
					if block_id != -1 and al.misasembled == self.assemblies[j].blocks[block_id].misasembled:
						sim_block_ids_within_asm[j] = block_id
						total += 1

				if total < similarThreshold(len(self.assemblies)):
					continue

				for j in range(0, len(self.assemblies)):
					block_id = sim_block_ids_within_asm[j]
					if block_id == -1:
						continue
					self.assemblies[j].blocks[block_id].similar = True
					self.assemblies[j].blocks[block_id].order = order

				order = (order + 1) % 2


	def DrawArcs(self, settings):
		for a in self.assemblies:
			a.DrawArcs(settings)


	def ApplyColors(self, settings):
		for a in self.assemblies:
			a.ApplyColor(settings)


	def FindMaxPos(self):
		max_pos = 0
		for asm in self.assemblies:
			asm_max_pos = asm.blocks[len(asm.blocks) - 1].end
			if max_pos < asm_max_pos:
				max_pos = asm_max_pos
		self.max_pos = max_pos
		return max_pos			



class Visualizer:
	def __init__(self, assemblies, covHist, settings):
		self.assemblies = assemblies
		self.covHist = covHist
		self.settings = settings

		self.figure = matplotlib.pyplot.figure()
		self.subplot = self.figure.add_subplot(111)

	def __del__(self):
		pass

	def show(self):
		self.subplot.axis("equal")
		self.subplot.axis("off")
		matplotlib.pyplot.show()

	def save(self, fileName):
		self.subplot.axis("equal")
		self.subplot.axis("off")
		self.figure.savefig(fileName + ".svg", format='svg')



	def plot_genome_axis(self, offset):
		self.subplot.add_line(matplotlib.lines.Line2D((offset[0], offset[0] + self.settings.assemblyWidth), (offset[1], offset[1]), c="black", lw=self.settings.axisWeight))

		i = 0.0
		while i < self.settings.genomeLength - self.settings.xTics / 5.0:
			x = offset[0] + self.settings.assemblyWidth * float(i) / float(self.settings.genomeLength)

			if self.settings.dashLines:
				self.subplot.add_line(matplotlib.lines.Line2D((x, x), (offset[1] + self.settings.plotHeight , self.settings.lastMargin), c="grey", ls=':', lw=self.settings.dashLineWeight))

			self.subplot.add_line(matplotlib.lines.Line2D((x, x), (offset[1], offset[1] - self.settings.ticLength), c="black", lw=self.settings.axisWeight))

			self.subplot.annotate(str(round(float(i) / self.settings.genomeAnnotationScale ,1)), (x + self.settings.xticsStep, offset[1] - self.settings.xticsStep), fontsize=8, horizontalalignment='left', verticalalignment='top')

			i += self.settings.xTics

		
		x = offset[0] + self.settings.assemblyWidth	
		if self.settings.dashLines:
			self.subplot.add_line(matplotlib.lines.Line2D((x, x), (offset[1] + self.settings.plotHeight , self.settings.lastMargin), c="grey", ls=':', lw=self.settings.dashLineWeight))
		self.subplot.add_line(matplotlib.lines.Line2D((x, x), (offset[1], offset[1] - self.settings.ticLength), c="black", lw=self.settings.axisWeight))

		self.subplot.annotate(str(round(float(self.settings.genomeLength) / self.settings.genomeAnnotationScale ,2)), (x + self.settings.xticsStep, offset[1] - self.settings.xticsStep), fontsize=8, horizontalalignment='left', verticalalignment='top')

		self.subplot.annotate(self.settings.genomeAnnotation, (offset[0] + self.settings.assemblyWidth / 2.0, offset[1] - self.settings.xLabelStep), fontsize=11, horizontalalignment='center', verticalalignment='top')



	def plot_coverage(self, covHist, offset):
		self.subplot.add_line(matplotlib.lines.Line2D((offset[0], offset[0]), (offset[1], offset[1] + self.settings.plotHeight), c="black", lw=self.settings.axisWeight))

		cov = 0.0
		while (cov <= self.settings.maxLogCov):
			y = offset[1] + cov * self.settings.plotYScale
			self.subplot.add_line(matplotlib.lines.Line2D((offset[0] - self.settings.ticLength, offset[0]), (y, y), c="black", lw=self.settings.axisWeight))
			self.subplot.annotate(str(int(round(math.pow(10, cov)))), (offset[0] - self.settings.yticsStep, y), fontsize=8, horizontalalignment='right', verticalalignment='center')
			cov += self.settings.yTics

		self.subplot.annotate("Coverage", (offset[0] - self.settings.yLabelStep, offset[1] + self.settings.plotYScale * self.settings.maxLogCov / 2.0), fontsize=11, horizontalalignment='center', verticalalignment='center', rotation = "vertical")

		for pos in covHist:
			x = offset[0] + pos * self.settings.plotXScale
			if covHist[pos] != 0:
				y = offset[1] + math.log10(covHist[pos]) * self.settings.plotYScale
				color = self.settings.coverageColor
			else:
				y = offset[1] + self.settings.zeroCovStep * self.settings.plotYScale
				color = self.settings.zeroCoverageColor

		        self.subplot.add_line(matplotlib.lines.Line2D((x, x + self.settings.dotLength), (y, y), c=color, lw=self.settings.dotWeight))



	def plot_assembly(self, assembly, offset):
		for name in assembly.contigs:
			for arc in assembly.contigs[name].arcs:
				x = offset[0] + (arc.c1 + arc.c2) * self.settings.scale / 2
				y = offset[1]
				width = abs(arc.c1 - arc.c2) * self.settings.scale
				height = 0.1 * width
				if height < 20:
					height = 20
				if height > 90:
					height = 90
				
				self.subplot.add_patch(matplotlib.patches.Arc((x, y), width, height, angle=180.0, theta1=0.0, theta2=180.0,  ec="black", color="black", lw=0.2))
			

		for block in assembly.blocks:
			x = offset[0] + block.start * self.settings.scale 
			y = offset[1] + block.vPositionDelta
			height = self.settings.contigHeight
			width = block.Length() * self.settings.scale

			self.subplot.add_patch(matplotlib.patches.Rectangle((x,y), width, height,  ec="black", color=block.color, fill=True, lw=0.0))



	def visualize(self):
		self.subplot.add_patch(matplotlib.patches.Rectangle((-20,0), self.settings.totalWidth + 20 + self.settings.lastMargin, self.settings.totalHeight + 0, color="white", fill=True, lw=0))

		self.plot_genome_axis( (self.settings.xOffset, self.settings.totalHeight - self.settings.plotHeight) )

		if self.covHist is not None:
			self.plot_coverage(self.covHist, (self.settings.xOffset, self.settings.totalHeight - self.settings.plotHeight) )

		if self.assemblies is not None:
			offset = self.settings.plotHeight + self.settings.plotMargin + self.settings.assemblyStep

			for assembly in self.assemblies.assemblies:
				self.subplot.annotate(assembly.name, (self.settings.xOffset - self.settings.nameAnnotationXStep, self.settings.totalHeight - offset + self.settings.nameAnnotationYStep), fontsize=12, horizontalalignment='right', verticalalignment='bottom')

				self.plot_assembly(assembly, (self.settings.xOffset, self.settings.totalHeight - offset))
				offset += self.settings.assemblyStep




def readCoverage(fileName):
	inFile = open(fileName, 'r')

	hist = {}
	maxPos = 0
	maxCov = 0

	for line in inFile:
		pos = line.strip().split(' ')
		hist[int(pos[0])] = int(pos[1])
		if maxPos < int(pos[0]):
			maxPos = int(pos[0])
		if maxCov < int(pos[1]):
			maxCov = int(pos[1])

	inFile.close()
	return hist, maxPos, maxCov



if len(sys.argv) < 2:
	print("Usage: \n -a <file with list of assemblies files> \n --cov <file with coverage histogram> \n -g <genome length> \n -o <output file name> \n --arcs \n --similar")
	sys.exit()

inFileName = None
covFileName = None
genomeLength = 0
outputName = None
arcs = False
similar = False

i = 1
while i < len(sys.argv):
	if sys.argv[i] == "-a":
		inFileName = sys.argv[i + 1]
		i += 2	
	elif sys.argv[i] == "--cov":
		covFileName = sys.argv[i + 1]
		i += 2	
	elif sys.argv[i] == "-g":
		genomeLength =  int(sys.argv[i + 1])
		i += 2	
	elif sys.argv[i] == "-o":
		outputName = sys.argv[i + 1]
		i += 2	
	elif sys.argv[i] == "--arcs":
		arcs = True
		i += 1
	elif sys.argv[i] == "--similar":
		similar = True
		i += 1
	else:
		print("Unknown argument " + sys.argv[i])
		i += 1

minVisualizedLength = 200

hist, maxPos, maxCov = None, 10, 10
if covFileName is not None:
	hist, maxPos, maxCov = readCoverage(covFileName)

assemblies = None
if inFileName is not None:
	assemblies = Assemblies(inFileName, genomeLength, minVisualizedLength)
	asmNumber = len(assemblies.assemblies)

	if genomeLength == 0:
		genomeLength = assemblies.FindMaxPos()
else:
	asmNumber = 0

if genomeLength == 0:
	print("Set genome length")
	sys.exit()

settings = Settings(genomeLength, maxPos, maxCov, asmNumber)

if assemblies is not None:
	if arcs and assemblies is not None:
		settings.assemblyStep += 40
		assemblies.DrawArcs(settings)

	if similar and assemblies is not None:
		assemblies.FindSimilar(settings)

	assemblies.ApplyColors(settings)


v = Visualizer(assemblies, hist, settings)
v.visualize()

if outputName is None:
	outputName = inFileName
if outputName is None:
	outputName = "coverage"

v.save(outputName)

