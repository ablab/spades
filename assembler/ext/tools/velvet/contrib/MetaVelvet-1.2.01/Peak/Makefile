OBJS = PeakDetectorParameters.o PeakDetectorFactory.o EMPeakDetector.o SimplePeakDetector.o
PROG =
INCLUDES = 
LIBS = 
SUBDIRS = . 

.SUFFIXES:.c .o

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(DEF) $(INCLUDES) $< -o $@

lib:libPeak.a

libPeak.a:$(OBJS)
	$(AR) -cru $@ $(OBJS)

cleanlocal:
	rm -f *.o $(PROG) *~ *.a

clean:cleanlocal
