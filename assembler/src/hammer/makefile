DIR=hammer
PROG=
MYCFLAGS=-fopenmp -fPIC -std=c++0x
MYLDFLAGS=-fopenmp
builddir=../../build/$(DIR)
basedir=../

hammer: $(builddir)/preproc $(builddir)/main
	@echo "== MADE HAMMER =="

$(builddir)/main:  $(builddir)/hammer_tools.o $(builddir)/main.o $(builddir)/kmer_cluster.o
	$(CC) $(LDFLAGS) $(builddir)/hammer_tools.o $(builddir)/main.o $(builddir)/kmer_cluster.o -o $(builddir)/$@
	@echo "== made main =="

$(builddir)/preproc: $(builddir)/preproc.o
	$(CC) $(LDFLAGS) $(builddir)/preproc.o -o $(builddir)/$@
	@echo "== made preproc =="

$(builddir)/kmer_cluster.o: kmer_cluster.cpp $(builddir)/.d
	$(CC) $(CFLAGS) kmer_cluster.cpp -o $(builddir)/$@
$(builddir)/main.o: main.cpp $(builddir)/.d
	$(CC) $(CFLAGS) main.cpp -o $(builddir)/$@
$(builddir)/hammer_tools.o: hammer_tools.cpp $(builddir)/.d
	$(CC) $(CFLAGS) hammer_tools.cpp -o $(builddir)/$@
$(builddir)/preproc.o: preproc.cpp $(builddir)/.d
	$(CC) $(CFLAGS) preproc.cpp -o $(builddir)/$@


include ../makefile.inc

#preproc_exec: $(OBJS) 


#$(builddir)/preproc.o: preproc.cpp
#inner_preproc: preproc.o sequence.o sequence_tools.o
#	$(CC) preproc.o sequence.o sequence_tools.o $(LDFLAGS) -o preproc
#preproc.o : preproc.cpp
#	$(CC) $(CFLAGS) -c preproc.cpp -o preproc.o


