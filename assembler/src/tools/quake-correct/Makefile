CC=g++
CFLAGS=-O3 -fopenmp -I/opt/local/var/macports/software/boost/1.46.1_0/opt/local/include -I.
LDFLAGS=-L. -lgzstream -lz
#INCLUDEDIR=
EXE_FILES = correct count-kmers count-qmers count_qmers reduce-kmers reduce-qmers trim build_bithash correct_stats
.PHONY: all clean

all: $(EXE_FILES)

clean:
	-rm $(EXE_FILES) *.o

correct: correct.cpp Read.o bithash.o edit.o libgzstream.a
	$(CC) $(CFLAGS) correct.cpp Read.o bithash.o edit.o -o correct $(LDFLAGS)

count-kmers: count-kmers.cpp count.o
	$(CC) $(CFLAGS) count-kmers.cpp count.o -o count-kmers

count-qmers: count-qmers.cpp count.o
	$(CC) $(CFLAGS) count-qmers.cpp count.o -o count-qmers

count_qmers: count_qmers.cpp count.o qmer_hash.o
	$(CC) $(CFLAGS) -o count_qmers count_qmers.cpp count.o qmer_hash.o

qmer_hash.o: qmer_hash.cpp qmer_hash.h
	$(CC) $(CFLAGS) -c qmer_hash.cpp

reduce-kmers: reduce-kmers.cpp
	$(CC) $(CFLAGS) reduce-kmers.cpp -o reduce-kmers

reduce-qmers: reduce-qmers.cpp
	$(CC) $(CFLAGS) reduce-qmers.cpp -o reduce-qmers

trim: trim.cpp Read.o bithash.o edit.o libgzstream.a
	$(CC) $(CFLAGS) trim.cpp Read.o bithash.o edit.o -o trim $(LDFLAGS)

build_bithash: build_bithash.cpp bithash.o
	$(CC) $(CFLAGS) build_bithash.cpp bithash.o -o build_bithash

correct_stats: stats.cpp
	$(CC) $(CFLAGS) stats.cpp -o correct_stats

Read.o: Read.cpp Read.h bithash.o
	$(CC) $(CFLAGS) -c Read.cpp

edit.o: edit.cpp edit.h
	$(CC) $(CFLAGS) -c edit.cpp

bithash.o: bithash.cpp bithash.h
	$(CC) $(CFLAGS) -c bithash.cpp

count.o: count.cpp
	$(CC) $(CFLAGS) -c count.cpp

gzstream.o: gzstream.C gzstream.h
	$(CC) $(CFLAGS) -c -o gzstream.o gzstream.C

libgzstream.a: gzstream.o
	ar cr libgzstream.a gzstream.o
