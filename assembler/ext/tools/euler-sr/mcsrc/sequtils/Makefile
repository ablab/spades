include ../common.mak

#EXECS = joincontigs evolve extractseq excrep mskrep extractorthseq combmsk matchseq extractReads locateSeq extractQuals selectSeq extractBlastHits mapper sortSequences alignPairs spliceReads printOriginalReads printUnmaskedContigs matchShort splitFastaFile trimReads windowGC  chopReads countWords exciseNs validateSeq createRandomSequence
EXECS = joincontigs evolve extractseq excrep mskrep extractorthseq combmsk matchseq extractReads locateSeq extractQuals selectSeq extractBlastHits mapper sortSequences alignPairs spliceReads printOriginalReads printUnmaskedContigs matchShort splitFastaFile trimReads windowGC  chopReads countWords exciseNs validateSeq

.PHONY: arch_dir $(EXECS)

all: $(CURDIR)/make.dep $(EXECS) install

include ../make.rules
include $(CURDIR)/make.dep
chopReads: arch_dir $(MACHTYPE)/chopReads
combmsk: arch_dir $(MACHTYPE)/combmsk
extractBlastHits: arch_dir $(MACHTYPE)/extractBlastHits
joincontigs: arch_dir $(MACHTYPE)/joincontigs
evolve: arch_dir $(MACHTYPE)/evolve
extractseq: arch_dir $(MACHTYPE)/extractseq
extractorthseq: arch_dir $(MACHTYPE)/extractorthseq
selectSeq: arch_dir $(MACHTYPE)/selectSeq
excrep: arch_dir $(MACHTYPE)/excrep
exciseNs: arch_dir $(MACHTYPE)/exciseNs
mskrep: arch_dir $(MACHTYPE)/mskrep
extractReads: arch_dir $(MACHTYPE)/extractReads
extractQuals: arch_dir $(MACHTYPE)/extractQuals
matchseq: arch_dir $(MACHTYPE)/matchseq
locateSeq: arch_dir $(MACHTYPE)/locateSeq
mapper: arch_dir $(MACHTYPE)/mapper
sortSequences: arch_dir $(MACHTYPE)/sortSequences
alignPairs: arch_dir $(MACHTYPE)/alignPairs
spliceReads: arch_dir $(MACHTYPE)/spliceReads
printOriginalReads: arch_dir $(MACHTYPE)/printOriginalReads
printUnmaskedContigs: arch_dir $(MACHTYPE)/printUnmaskedContigs
matchShort: arch_dir $(MACHTYPE)/matchShort
splitFastaFile: arch_dir $(MACHTYPE)/splitFastaFile
trimReads: arch_dir $(MACHTYPE)/trimReads
windowGC: arch_dir $(MACHTYPE)/windowGC
countWords: arch_dir $(MACHTYPE)/countWords
validateSeq: arch_dir $(MACHTYPE)/validateSeq
createRandomSequence: arch_dir $(MACHTYPE)/createRandomSequence

LIBCOMMON = $(SRC_DIR)/lib/$(MACHTYPE)/libcommon.${LIB_EXT} 
LIBPB     = $(SRC_DIR)/lib/$(MACHTYPE)/libparseblast.${LIB_EXT} 
LIBALIGN  = $(SRC_DIR)/lib/$(MACHTYPE)/libalign.${LIB_EXT} 

$(MACHTYPE)/chopReads: ${MACHTYPE}/ChopReads.o
				$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@

$(MACHTYPE)/mapper: ${MACHTYPE}/PerfectSeqMapper.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/selectSeq: ${MACHTYPE}/SelectSequencesByLength.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/locateSeq: ${MACHTYPE}/LocateSequence.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/extractReads: ${MACHTYPE}/ExtractReads.o 
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/extractQuals: ${MACHTYPE}/ExtractQualities.o 
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/matchseq: ${MACHTYPE}/MatchSeq.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/combmsk: ${MACHTYPE}/CombineMasking.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/joincontigs: ${MACHTYPE}/joincontigs.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/evolve: ${MACHTYPE}/evolve.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@  $(LIBDIR) -lcommon $(LIBS)

$(MACHTYPE)/extractorthseq: $(MACHTYPE)/extractorthseq.o
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@   -lcommon $(LIBDIR) $(LIBS)

$(MACHTYPE)/extractseq: $(MACHTYPE)/ExtractSequence.o 
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@   -lcommon $(LIBDIR) $(LIBS)

$(MACHTYPE)/excrep: $(MACHTYPE)/ExciseRepeats.o 
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS)

$(MACHTYPE)/exciseNs: $(MACHTYPE)/ExciseNs.o 
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS)

$(MACHTYPE)/mskrep: $(MACHTYPE)/HardMaskRepeats.o  $(LIBCOMMON);
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS)

$(MACHTYPE)/extractBlastHits: $(MACHTYPE)/ExtractBlastHits.o  $(LIBCOMMON) $(LIBPB)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lparseblast -lcommon  $(LIBDIR) $(LIBS) 

$(MACHTYPE)/sortSequences: $(MACHTYPE)/SortSequencesByName.o  $(LIBCOMMON) $(LIBPB)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lparseblast -lcommon  $(LIBDIR) $(LIBS) 


$(MACHTYPE)/alignPairs: $(MACHTYPE)/AlignPairedSequences.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lparseblast -lcommon -lalign $(LIBDIR) $(LIBS) 

$(MACHTYPE)/spliceReads: $(MACHTYPE)/SpliceReads.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/printOriginalReads: $(MACHTYPE)/PrintOriginalReads.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/printUnmaskedContigs: $(MACHTYPE)/PrintUnmaskedContigs.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/matchShort: $(MACHTYPE)/MatchShortSequences.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/splitFastaFile: $(MACHTYPE)/SplitFastaFile.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/trimReads: $(MACHTYPE)/TrimReads.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/windowGC: $(MACHTYPE)/WindowGC.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/countWords: $(MACHTYPE)/CountWords.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 


$(MACHTYPE)/validateSeq: $(MACHTYPE)/ValidateSeq.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 

$(MACHTYPE)/createRandomSequence: $(MACHTYPE)/CreateRandomSequence.o  $(LIBCOMMON) $(LIBPB) $(LIBALIGN)
	$(CPP) $(CPPOPTS_lib) $(INCLUDEDIR) $< -o $@ -lcommon $(LIBDIR) $(LIBS) 


$(LIBPB):
	cd ../lib/blast && $(MAKE)
