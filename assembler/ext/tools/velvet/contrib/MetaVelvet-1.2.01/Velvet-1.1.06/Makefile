OBJS = readSet.o tightString.o kmer.o kmerOccurenceTable.o \
	roadMap.o splay.o splayTable.o \
	preGraph.o preGraphConstruction.o concatenatedPreGraph.o \
	graph.o graphReConstruction.o concatenatedGraph.o \
	graphStats.o correctedGraph.o \
	fib.o fibHeap.o dfib.o dfibHeap.o \
	passageMarker.o \
	scaffold.o shortReadPairs.o locallyCorrectedGraph.o readCoherentGraph.o \
	utility.o allocArray.o recycleBin.o
PROG =
INCLUDES = 
LIBS = 
SUBDIRS = . 

.SUFFIXES:.c .o

.c.o:
	$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(DEF) $(INCLUDES) $< -o $@

lib:libVelvet.a

libVelvet.a:$(OBJS)
	$(AR) -cru $@ $(OBJS)

cleanlocal:
	rm -f *.o $(PROG) *~ *.a

clean:cleanlocal
