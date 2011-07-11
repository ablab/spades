SRC=sequence/sequence.hpp sequence/sequence_tools.hpp
DIR=hammer
PROG=preproc
MYCFLAGS=-fPIC -std=c++0x
MYLDFLAGS=-fopenmp
builddir=../../build/$(DIR)
basedir=../

hammer: preproc main reconstruct
	@echo "== MADE HAMMER =="

main:  hammer_tools.o main.o
	$(CC) $(LDFLAGS) $(builddir)/hammer_tools.o $(builddir)/main.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made main =="

center:  center.o
	$(CC) $(LDFLAGS) $(builddir)/center.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made center =="

cluster: cluster.o
	$(CC) $(LDFLAGS) $(builddir)/cluster.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made cluster =="

preproc: preproc.o
	$(CC) $(LDFLAGS) $(builddir)/preproc.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made preproc =="

preproc2: preproc2.o
	$(CC) $(LDFLAGS) $(builddir)/preproc2.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made preproc2 =="

reconstruct: reconstruct.o
	$(CC) $(LDFLAGS) $(builddir)/reconstruct.o $(LDFLAGS) -o $(builddir)/$@
	@echo "== made reconstruct =="


main.o: main.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c main.cpp -o $(builddir)/$@
hammer_tools.o: hammer_tools.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c hammer_tools.cpp -o $(builddir)/$@
cluster.o: cluster.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c cluster.cpp -o $(builddir)/$@
center.o: center.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c center.cpp -o $(builddir)/$@
reconstruct.o: reconstruct.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c reconstruct.cpp -o $(builddir)/$@
preproc.o: preproc.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c preproc.cpp -o $(builddir)/$@
preproc2.o: preproc2.cpp $(builddir)/.d
	$(CC) $(CFLAGS) -c preproc2.cpp -o $(builddir)/$@


include ../makefile.inc

#preproc_exec: $(OBJS) 


#$(builddir)/preproc.o: preproc.cpp
#inner_preproc: preproc.o sequence.o sequence_tools.o
#	$(CC) preproc.o sequence.o sequence_tools.o $(LDFLAGS) -o preproc
#preproc.o : preproc.cpp
#	$(CC) $(CFLAGS) -c preproc.cpp -o preproc.o


