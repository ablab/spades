OBJS = FastaSeq.o FastaWriter.o SplitStats.o SplitJudge.o InfiniteLoopChecker.o MetaGraph.o
PROG =
INCLUDES = 
LIBS = 
SUBDIRS = . 

.SUFFIXES:.c .o

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(DEF) $(INCLUDES) $< -o $@

lib:libCommon.a

libCommon.a:$(OBJS)
	$(AR) -cru $@ $(OBJS)

cleanlocal:
	rm -f *.o $(PROG) *~ *.a

clean:cleanlocal
