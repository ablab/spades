CC ?= gcc
CXX ?= g++
CFLAGS = -O3 -Wall
CXXFLAGS = -O2 -Wall
DFLAGS =
VelvetDir=Velvet-1.1.06
OBJS_V = $(VelvetDir)/splay.o $(VelvetDir)/splayTable.o \
	$(VelvetDir)/tightString.o $(VelvetDir)/graph.o \
	$(VelvetDir)/fibHeap.o $(VelvetDir)/fib.o \
	$(VelvetDir)/concatenatedGraph.o $(VelvetDir)/passageMarker.o \
	$(VelvetDir)/graphStats.o $(VelvetDir)/correctedGraph.o \
	$(VelvetDir)/dfib.o $(VelvetDir)/dfibHeap.o \
	$(VelvetDir)/recycleBin.o $(VelvetDir)/readSet.o \
	$(VelvetDir)/shortReadPairs.o $(VelvetDir)/scaffold.o \
	$(VelvetDir)/locallyCorrectedGraph.o $(VelvetDir)/graphReConstruction.o \
	$(VelvetDir)/roadMap.o $(VelvetDir)/preGraph.o \
	$(VelvetDir)/preGraphConstruction.o $(VelvetDir)/concatenatedPreGraph.o \
	$(VelvetDir)/readCoherentGraph.o $(VelvetDir)/utility.o \
	$(VelvetDir)/kmer.o $(VelvetDir)/kmerOccurenceTable.o $(VelvetDir)/allocArray.o
OBJS_H = Apps/meta-velveth.o $(OBJS_V)
OBJS_G = Apps/meta-velvetg.o $(OBJS_V)	
OBJS_IS = Apps/annoIS.o $(OBJS_V)
#PROG = meta-velvetg annoIS
PROG = meta-velvetg
INCLUDES = 
LIBS = -lz \
	-lVelvet -LVelvet-1.1.06 \
	-lVelvetAPI -LVelvetAPI \
	-lCommon -LCommon \
	-lPeak -LPeak \
	-lISGraph -L ISGraph
SUBDIRS = Velvet-1.1.06 VelvetAPI Common Peak ISGraph

MAXKMERLENGTH = 63
CATEGORIES = 2
DEF = -D MAXKMERLENGTH=$(MAXKMERLENGTH) -D CATEGORIES=$(CATEGORIES)

.SUFFIXES:.cc .o

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(DEF) $(INCLUDES) $< -o $@

all:lib-recur $(PROG)

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
	@target=`echo $@ | sed s/-recur//`; \
	wdir=`pwd`; \
	list='$(SUBDIRS)'; for subdir in $$list; do \
		cd $$subdir; \
		$(MAKE) CC="$(CC)" CXX="$(CXX)" \
			CFLAGS="$(CFLAGS)" CXXFLAGS="$(CXXFLAGS)" \
			DFLAGS="$(DFLAGS)" DEF="$(DEF)"\
			INCLUDES="$(INCLUDES)" $$target || exit 1; \
		cd $$wdir; \
	done;

lib:

meta-velveth:lib-recur $(OBJS_H)
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(OBJS_H) -o $@ $(LIBS)

meta-velvetg:lib-recur $(OBJS_G)
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(OBJS_G) -o $@ $(LIBS)

annoIS:lib-recur $(OBJS_IS)
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(OBJS_IS) -o $@ $(LIBS)

cleanlocal:
	rm -f *.o $(PROG) *~ *.a Apps/*.o

clean:cleanlocal clean-recur
