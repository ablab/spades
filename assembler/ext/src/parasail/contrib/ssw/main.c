/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 1.2.2
 *  Last revision by Mengyao Zhao on 2017-05-30.
 */
#include "config.h"

/* getopt needs _POSIX_C_SOURCE 2 */
#define _POSIX_C_SOURCE 2

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#if defined(_MSC_VER)
#include "wingetopt/src/getopt.h"
#else
#include <unistd.h>
#endif
#include "ssw.h"
#include "kseq.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

#if defined(_MSC_VER)
#include <io.h>
#define READ_FUNCTION _read
#else
#define READ_FUNCTION read
#endif

#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

static void reverse_comple(const char* seq, char* rc) {
	int32_t end = strlen(seq), start = 0;
	static const int8_t rc_table[128] = {
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4, 4,
		4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4, 4,
		4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};
	rc[end] = '\0';
	-- end;
	while (LIKELY(start < end)) {
		rc[start] = (char)rc_table[(int8_t)seq[end]];
		rc[end] = (char)rc_table[(int8_t)seq[start]];
		++ start;
		-- end;
	}
	if (start == end) rc[start] = (char)rc_table[(int8_t)seq[start]];
}

static void ssw_write (s_align* a,
			const kseq_t* ref_seq,
			const kseq_t* read,
			const char* read_seq,	// strand == 0: original read; strand == 1: reverse complement read
			const int8_t* ref_num,
			const int8_t* read_num,
			const int8_t* table,
			int8_t strand,	// 0: forward aligned ; 1: reverse complement aligned
			int8_t sam) {	// 0: Blast like output; 1: Sam format output

	int32_t mismatch;
	if (sam == 0) {	// Blast like output
		fprintf(stdout, "target_name: %s\nquery_name: %s\noptimal_alignment_score: %d\t", ref_seq->name.s, read->name.s, a->score1);
		if (a->score2 > 0) fprintf(stdout, "suboptimal_alignment_score: %d\t", a->score2);
		if (strand == 0) fprintf(stdout, "strand: +\t");
		else fprintf(stdout, "strand: -\t");
		if (a->ref_begin1 + 1) fprintf(stdout, "target_begin: %d\t", a->ref_begin1 + 1);
		fprintf(stdout, "target_end: %d\t", a->ref_end1 + 1);
		if (a->read_begin1 + 1) fprintf(stdout, "query_begin: %d\t", a->read_begin1 + 1);
		fprintf(stdout, "query_end: %d\n\n", a->read_end1 + 1);
		if (a->cigar) {
			int32_t c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
			uint32_t i;
			while (e < a->cigarLen || left > 0) {
				int32_t count = 0;
				int32_t q = qb;
				int32_t p = pb;
				fprintf(stdout, "Target: %8d    ", q + 1);
				for (c = e; c < a->cigarLen; ++c) {
					char letter = cigar_int_to_op(a->cigar[c]);
					uint32_t length = cigar_int_to_len(a->cigar[c]);
					uint32_t l = (count == 0 && left > 0) ? left: length;
					for (i = 0; i < l; ++i) {
						if (letter == 'I') fprintf(stdout, "-");
						else {
							fprintf(stdout, "%c", *(ref_seq->seq.s + q));
							++ q;
						}
						++ count;
						if (count == 60) goto step2;
					}
				}
step2:
				fprintf(stdout, "    %d\n                    ", q);
				q = qb;
				count = 0;
				for (c = e; c < a->cigarLen; ++c) {
					char letter = cigar_int_to_op(a->cigar[c]);
					uint32_t length = cigar_int_to_len(a->cigar[c]);
					uint32_t l = (count == 0 && left > 0) ? left: length;
					for (i = 0; i < l; ++i){
						if (letter == 'M') {
							if (table[(int)*(ref_seq->seq.s + q)] == table[(int)*(read_seq + p)])fprintf(stdout, "|");
							else fprintf(stdout, "*");
							++q;
							++p;
						} else {
							fprintf(stdout, " ");
							if (letter == 'I') ++p;
							else ++q;
						}
						++ count;
						if (count == 60) {
							qb = q;
							goto step3;
						}
					}
				}
step3:
				p = pb;
				fprintf(stdout, "\nQuery:  %8d    ", p + 1);
				count = 0;
				for (c = e; c < a->cigarLen; ++c) {
					char letter = cigar_int_to_op(a->cigar[c]);
					uint32_t length = cigar_int_to_len(a->cigar[c]);
					uint32_t l = (count == 0 && left > 0) ? left: length;
					for (i = 0; i < l; ++i) {
						if (letter == 'D') fprintf(stdout, "-");
						else {
							fprintf(stdout, "%c", *(read_seq + p));
							++p;
						}
						++ count;
						if (count == 60) {
							pb = p;
							left = l - i - 1;
							e = (left == 0) ? (c + 1) : c;
							goto end;
						}
					}
				}
				e = c;
				left = 0;
end:
				fprintf(stdout, "    %d\n\n", p);
			}
		}
	}else {	// Sam format output
		fprintf(stdout, "%s\t", read->name.s);
		if (a->score1 == 0) fprintf(stdout, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
		else {
			int32_t c, p;
			uint32_t mapq = -4.343 * log(1 - (double)abs(a->score1 - a->score2)/(double)a->score1);
			mapq = (uint32_t) (mapq + 4.99);
			mapq = mapq < 254 ? mapq : 254;
			if (strand) fprintf(stdout, "16\t");
			else fprintf(stdout, "0\t");
			fprintf(stdout, "%s\t%d\t%d\t", ref_seq->name.s, a->ref_begin1 + 1, mapq);
			mismatch = mark_mismatch(a->ref_begin1, a->read_begin1, a->read_end1, ref_num, read_num, read->seq.l, &a->cigar, &a->cigarLen);
			for (c = 0; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				fprintf(stdout, "%lu%c", (unsigned long)length, letter);
			}
			fprintf(stdout, "\t*\t0\t0\t");
			fprintf(stdout, "%s", read_seq);
			fprintf(stdout, "\t");
			if (read->qual.s && strand) {
				for (p = read->qual.l - 1; p >= 0; --p) fprintf(stdout, "%c", read->qual.s[p]);
			}else if (read->qual.s) fprintf (stdout, "%s", read->qual.s);
			else fprintf(stdout, "*");
			fprintf(stdout, "\tAS:i:%d", a->score1);
			fprintf(stdout,"\tNM:i:%d\t", mismatch);
			if (a->score2 > 0) fprintf(stdout, "ZS:i:%d\n", a->score2);
			else fprintf(stdout, "\n");
		}
	}
}

int main (int argc, char * const argv[]) {
	clock_t start, end;
	float cpu_time;
	gzFile read_fp, ref_fp;
	kseq_t *read_seq, *ref_seq;
	int32_t l, m, k, match = 2, mismatch = 2, gap_open = 3, gap_extension = 1, path = 0, reverse = 0, n = 5, sam = 0, protein = 0, header = 0, s1 = 67108864, s2 = 128, filter = 0;
	int8_t* mata = (int8_t*)calloc(25, sizeof(int8_t));
	const int8_t* mat = mata;
	char mat_name[16];
	mat_name[0] = '\0';
	int8_t* ref_num = (int8_t*)malloc(s1);
	int8_t* num = (int8_t*)malloc(s2), *num_rc = 0;
	char* read_rc = 0;

	static const int8_t mat50[] = {
	//  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     	5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,	// A
       -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,	// R
       -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,	// N
       -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,	// D
       -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,	// C
       -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,	// Q
       -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,	// E
     	0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,	// G
       -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,	// H
       -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,	// I
       -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,	// L
       -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,	// K
       -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,	// M
       -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,	// F
       -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,	// P
     	1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,	// S
    	0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, 	// T
       -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, 	// W
       -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, 	// Y
     	0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, 	// V
       -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, 	// B
       -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, 	// Z
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, 	// X
       -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1 	// *
	};

	/* This table is used to transform amino acid letters into numbers. */
	int8_t aa_table[128] = {
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
		23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
		14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
	};

	/* This table is used to transform nucleotide letters into numbers. */
	int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};

	int8_t* table = nt_table;

	// Parse command line.
	while ((l = getopt(argc, argv, "m:x:o:e:a:f:pcrsh")) >= 0) {
		switch (l) {
			case 'm': match = atoi(optarg); break;
			case 'x': mismatch = atoi(optarg); break;
			case 'o': gap_open = atoi(optarg); break;
			case 'e': gap_extension = atoi(optarg); break;
			case 'a': strcpy(mat_name, optarg); break;
			case 'f': filter = atoi(optarg); break;
			case 'p': protein = 1; break;
			case 'c': path = 1; break;
			case 'r': reverse = 1; break;
			case 's': sam = 1; break;
			case 'h': header = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: ssw_test [options] ... <target.fasta> <query.fasta>(or <query.fastq>)\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "\t-m N\tN is a positive integer for weight match in genome sequence alignment. [default: 2]\n");
		fprintf(stderr, "\t-x N\tN is a positive integer. -N will be used as weight mismatch in genome sequence alignment. [default: 2]\n");
		fprintf(stderr, "\t-o N\tN is a positive integer. -N will be used as the weight for the gap opening. [default: 3]\n");
		fprintf(stderr, "\t-e N\tN is a positive integer. -N will be used as the weight for the gap extension. [default: 1]\n");
		fprintf(stderr, "\t-p\tDo protein sequence alignment. Without this option, the ssw_test will do genome sequence alignment.\n");
		fprintf(stderr, "\t-a FILE\tFILE is either the Blosum or Pam weight matrix. [default: Blosum50]\n");
		fprintf(stderr, "\t-c\tReturn the alignment path.\n");
		fprintf(stderr, "\t-f N\tN is a positive integer. Only output the alignments with the Smith-Waterman score >= N.\n");
		fprintf(stderr, "\t-r\tThe best alignment will be picked between the original read alignment and the reverse complement read alignment.\n");
		fprintf(stderr, "\t-s\tOutput in SAM format. [default: no header]\n");
		fprintf(stderr, "\t-h\tIf -s is used, include header in SAM output.\n\n");
		return 1;
	}

	// initialize scoring matrix for genome sequences
	for (l = k = 0; LIKELY(l < 4); ++l) {
		for (m = 0; LIKELY(m < 4); ++m) mata[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
		mata[k++] = 0; // ambiguous base
	}
	for (m = 0; LIKELY(m < 5); ++m) mata[k++] = 0;

	if (protein == 1 && (! strcmp(mat_name, "\0"))) {
		n = 24;
		table = aa_table;
		mat = mat50;
	} else if (strcmp(mat_name, "\0")) {

	// Parse score matrix.
		FILE *f_mat = fopen(mat_name, "r");
		char line[128];
		mata = (int8_t*)realloc(mata, 1024 * sizeof(int8_t));
		k = 0;
		m = 0;
		while (fgets(line, 128, f_mat)) {
			if (line[0] == '*' || (line[0] >= 'A' && line[0] <= 'Z')) {
				if (line[0] >= 'A' && line[0] <= 'Z') aa_table[(int)line[0]] = aa_table[(int)line[0] + 32] = m;
				char str[4], *s = str;
				str[0] = '\0';
				l = 1;
				while (line[l]) {
					if ((line[l] >= '0' && line[l] <= '9') || line[l] == '-') *s++ = line[l];
					else if (str[0] != '\0') {
						*s = '\0';
						mata[k++] = (int8_t)atoi(str);
						s = str;
						str[0] = '\0';
					}
					++l;
				}
				if (str[0] != '\0') {
					*s = '\0';
					mata[k++] = (int8_t)atoi(str);
					s = str;
					str[0] = '\0';
				}
				++m;
			}
		}
		if (k == 0) {
			fprintf(stderr, "Problem of reading the weight matrix file.\n");
			return 1;
		}
		fclose(f_mat);
		n = m;
		table = aa_table;
		mat = mata;
	}

	read_fp = gzopen(argv[optind + 1], "r");

    if (! read_fp) {
        fprintf (stderr, "gzopen of '%s' failed.\n", argv[optind + 1]);
            exit (EXIT_FAILURE);
    }

	read_seq = kseq_init(read_fp);
	if (sam && header && path) {
		fprintf(stdout, "@HD\tVN:1.4\tSO:queryname\n");
		ref_fp = gzopen(argv[optind], "r");
		ref_seq = kseq_init(ref_fp);
		while ((l = kseq_read(ref_seq)) >= 0) fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n", ref_seq->name.s, (int32_t)ref_seq->seq.l);
		kseq_destroy(ref_seq);
		gzclose(ref_fp);
	} else if (sam && !path) {
		fprintf(stderr, "SAM format output is only available together with option -c.\n");
		sam = 0;
	}

	// alignment
	if (reverse == 1 && n == 5) {
		read_rc = (char*)malloc(s2);
		num_rc = (int8_t*)malloc(s2);
	}
	start = clock();
	while (kseq_read(read_seq) >= 0) {
		s_profile* p, *p_rc = 0;
		int32_t readLen = read_seq->seq.l;
		int32_t maskLen = readLen / 2;

		while (readLen >= s2) {
			++s2;
			kroundup32(s2);
			num = (int8_t*)realloc(num, s2);
			if (reverse == 1 && n == 5) {
				read_rc = (char*)realloc(read_rc, s2);
				num_rc = (int8_t*)realloc(num_rc, s2);
			}
		}
		for (m = 0; m < readLen; ++m) num[m] = table[(int)read_seq->seq.s[m]];
		p = ssw_init(num, readLen, mat, n, 2);
		if (reverse == 1 && n == 5) {
			reverse_comple(read_seq->seq.s, read_rc);
			for (m = 0; m < readLen; ++m) num_rc[m] = table[(int)read_rc[m]];
			p_rc = ssw_init(num_rc, readLen, mat, n, 2);
		}else if (reverse == 1 && n == 24) {
			fprintf (stderr, "Reverse complement alignment is not available for protein sequences. \n");
			return 1;
		}

		ref_fp = gzopen(argv[optind], "r");
		ref_seq = kseq_init(ref_fp);
		while (kseq_read(ref_seq) >= 0) {
			s_align* result, *result_rc = 0;
			int32_t refLen = ref_seq->seq.l;
			int8_t flag = 0;
			while (refLen > s1) {
				++s1;
				kroundup32(s1);
				ref_num = (int8_t*)realloc(ref_num, s1);
			}
			for (m = 0; m < refLen; ++m) ref_num[m] = table[(int)ref_seq->seq.s[m]];
			if (path == 1) flag = 2;
			result = ssw_align (p, ref_num, refLen, gap_open, gap_extension, flag, filter, 0, maskLen);
			if (reverse == 1 && protein == 0)
				result_rc = ssw_align(p_rc, ref_num, refLen, gap_open, gap_extension, flag, filter, 0, maskLen);
			if (result_rc && result_rc->score1 > result->score1 && result_rc->score1 >= filter) {
				if (sam) ssw_write (result_rc, ref_seq, read_seq, read_rc, ref_num, num_rc, table, 1, 1);
				else ssw_write (result_rc, ref_seq, read_seq, read_rc, ref_num, num_rc, table, 1, 0);
			}else if (result && result->score1 >= filter){
				if (sam) ssw_write(result, ref_seq, read_seq, read_seq->seq.s, ref_num, num, table, 0, 1);
				else ssw_write(result, ref_seq, read_seq, read_seq->seq.s, ref_num, num, table, 0, 0);
			} else if (! result) return 1;
			if (result_rc) align_destroy(result_rc);
			align_destroy(result);
		}

		if(p_rc) init_destroy(p_rc);
		init_destroy(p);
		kseq_destroy(ref_seq);
		gzclose(ref_fp);
	}
	end = clock();
	cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	fprintf(stderr, "CPU time: %f seconds\n", cpu_time);

	if (num_rc) {
		free(num_rc);
		free(read_rc);
	}
	kseq_destroy(read_seq);
	gzclose(read_fp);
	free(num);
	free(ref_num);
 	free(mata);
	return 0;
}
