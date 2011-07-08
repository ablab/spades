SRC=sequence/sequence.cpp sequence/sequence_tools.cpp
DIR=hammer
PROG=preproc
MYCFLAGS=-fPIC -std=c++0x
MYLDFLAGS=-fopenmp
builddir=../../build/$(DIR)

hammer: preproc cluster center
	@echo "== MADE HAMMER =="

center: sequence.o sequence_tools.o center.o
	$(CC) $(LDFLAGS) $(builddir)/sequence.o $(builddir)/sequence_tools.o $(builddir)/center.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made center =="

cluster: sequence.o sequence_tools.o cluster.o
	$(CC) $(LDFLAGS) $(builddir)/sequence.o $(builddir)/sequence_tools.o $(builddir)/cluster.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made cluster =="

preproc: sequence.o sequence_tools.o preproc.o
	$(CC) $(LDFLAGS) $(builddir)/sequence.o $(builddir)/sequence_tools.o $(builddir)/preproc.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made preproc =="

cluster.o: cluster.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c cluster.cpp -o $(builddir)/$@
center.o: center.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c center.cpp -o $(builddir)/$@
preproc.o: preproc.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c preproc.cpp -o $(builddir)/$@
sequence.o: ../common/sequence/sequence.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c ../common/sequence/sequence.cpp -o $(builddir)/$@
sequence_tools.o: ../common/sequence/sequence_tools.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c ../common/sequence/sequence_tools.cpp -o $(builddir)/$@


include ../makefile.inc

#preproc_exec: $(OBJS) 


#$(builddir)/preproc.o: preproc.cpp
#inner_preproc: preproc.o sequence.o sequence_tools.o
#	$(CC) preproc.o sequence.o sequence_tools.o $(LDFLAGS) -o preproc
#preproc.o : preproc.cpp
#	$(CC) $(CFLAGS) -c preproc.cpp -o preproc.o


