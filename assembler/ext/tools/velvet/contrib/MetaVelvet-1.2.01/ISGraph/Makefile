OBJS = ISPath.o ISNode.o ISGraph.o
PROG =
INCLUDES = 
LIBS = 
SUBDIRS = . 

.SUFFIXES:.c .o

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(DEF) $(INCLUDES) $< -o $@

lib:libISGraph.a

libISGraph.a:$(OBJS)
	$(AR) -cru $@ $(OBJS)

cleanlocal:
	rm -f *.o $(PROG) *~ *.a

clean:cleanlocal
