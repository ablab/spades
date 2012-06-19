OBJS = VUtils.o VelvetUtils.o VelvetGraph.o
PROG =
INCLUDES = 
LIBS = 
SUBDIRS = . 

.SUFFIXES:.c .o

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(DEF) $(INCLUDES) $< -o $@

lib:libVelvetAPI.a

libVelvetAPI.a:$(OBJS)
	$(AR) -cru $@ $(OBJS)

cleanlocal:
	rm -f *.o $(PROG) *~ *.a

clean:cleanlocal
