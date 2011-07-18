DIR=hammer
PROG=
MYCFLAGS=-fopenmp -fPIC -std=c++0x
MYLDFLAGS=-fopenmp
builddir=../../build/$(DIR)
basedir=../

hammer: $(builddir)/preproc $(builddir)/main $(builddir)/reconstruct
	@echo "== MADE HAMMER =="

$(builddir)/main:  $(builddir)/hammer_tools.o $(builddir)/main.o $(builddir)/kmer_functions.o
	$(CC) $(LDFLAGS) $(builddir)/hammer_tools.o $(builddir)/main.o $(builddir)/kmer_functions.o -o $(builddir)/$@
	@echo "== made main =="

$(builddir)/center:  $(builddir)/center.o
	$(CC) $(LDFLAGS) $(builddir)/center.o -o $(builddir)/$@
	@echo "== made center =="

$(builddir)/cluster: $(builddir)/cluster.o
	$(CC) $(LDFLAGS) $(builddir)/cluster.o -o $(builddir)/$@
	@echo "== made cluster =="

$(builddir)/preproc: $(builddir)/preproc.o $(builddir)/kmer_functions.o $(builddir)/kmer_part_joiner.o
	$(CC) $(LDFLAGS) $(builddir)/kmer_functions.o $(builddir)/preproc.o $(builddir)/kmer_part_joiner.o -o $(builddir)/$@
	@echo "== made preproc =="

$(builddir)/preproc2: $(builddir)/preproc2.o
	$(CC) $(LDFLAGS) $(builddir)/preproc2.o -o $(builddir)/$@
	@echo "== made preproc2 =="

$(builddir)/reconstruct: $(builddir)/reconstruct.o $(builddir)/kmer_functions.o
	$(CC) $(LDFLAGS) $(builddir)/kmer_functions.o $(builddir)/reconstruct.o -o $(builddir)/$@
	@echo "== made reconstruct =="

$(builddir)/kmer_part_joiner.o: kmer_part_joiner.cpp $(builddir)/.d
	$(CC) $(CFLAGS) kmer_part_joiner.cpp -o $(builddir)/$@
$(builddir)/kmer_functions.o: kmer_functions.cpp  $(builddir)/.d
	$(CC) $(CFLAGS) kmer_functions.cpp -o $(builddir)/$@
$(builddir)/main.o: main.cpp $(builddir)/.d
	$(CC) $(CFLAGS) main.cpp -o $(builddir)/$@
$(builddir)/hammer_tools.o: hammer_tools.cpp $(builddir)/.d
	$(CC) $(CFLAGS) hammer_tools.cpp -o $(builddir)/$@
$(builddir)/cluster.o: cluster.cpp $(builddir)/.d
	$(CC) $(CFLAGS) cluster.cpp -o $(builddir)/$@
$(builddir)/center.o: center.cpp $(builddir)/.d
	$(CC) $(CFLAGS) center.cpp -o $(builddir)/$@
$(builddir)/reconstruct.o: reconstruct.cpp $(builddir)/.d
	$(CC) $(CFLAGS) reconstruct.cpp -o $(builddir)/$@
$(builddir)/preproc.o: preproc.cpp $(builddir)/.d
	$(CC) $(CFLAGS) preproc.cpp -o $(builddir)/$@
$(builddir)/preproc2.o: preproc2.cpp $(builddir)/.d
	$(CC) $(CFLAGS) preproc2.cpp -o $(builddir)/$@


include ../makefile.inc

#preproc_exec: $(OBJS) 


#$(builddir)/preproc.o: preproc.cpp
#inner_preproc: preproc.o sequence.o sequence_tools.o
#	$(CC) preproc.o sequence.o sequence_tools.o $(LDFLAGS) -o preproc
#preproc.o : preproc.cpp
#	$(CC) $(CFLAGS) -c preproc.cpp -o preproc.o


