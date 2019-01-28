## Comparison with other aligners

**Special thanks to Martin Šošić for the original version of this
document in the OpAl software package. See
https://github.com/Martinsos/opal/blob/master/aligner_comparison.md for
the original. **

The following presents results of speed comparisons of parasail with
other aligners: SSW, ssearch (FASTA), OpAl, and SWIPE.  Since SSW and
SWIPE only do Smith-Waterman alignment, we compared all software only
for SW.

Aligners were tested by quering sequences against UniProtKB/Swiss-Prot
database (containing 547964 sequences).  The database can be obtained
from www.uniprot.org/downloads -> UniProtKB/Swiss-Prot.  

Specific sequences can also be obtained from www.uniprot.org by
searching for them by name (Search tab).

All aligners were tested with following parameters:
* number of threads = 1
* gap opening = 3
* gap extension = 1
* score matrix = BLOSUM50

Only scores were calculated (not alignments). Time spent to read query
and database sequences was not measured. In the case of SSW, output was
programmatically suppressed to avoid its effect on alignment timing.

Links to the software used:
* SSW: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
* OpAl: https://github.com/Martinsos/opal
* SWIPE: https://github.com/torognes/swipe
* ssearch: http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.7b.tar.gz

How aligners were called:
* SSW: `./ssw_test -p uniprot_sprot.fasta <query_file>`
* OpAl: `./opal_aligner -s <query_file> uniprot_sprot.fasta`
* SWIPE: `./swipe -a 1 -p 1 -G 3 -E 1 -M BLOSUM50 -b 0 -i <query_file> -d uniprot_sprot`
  * NOTE: database had to be preprocessed for SWIPE using _makeblastdb_
* ssearch: `./ssearch36 -d 0 -T 1 -p -f -3 -g -1 -s BL50 <query_file> uniprot_sprot.fasta`
* parasail: `parasail_aligner -x -t 1 -a <function> -o 3 -e 1 -m blosum50 -f uniprot_sprot -q <query_file>`
  * SSE4.1, 16-bit: -a sw_striped_profile_sse41_128_16
  * SSE4.1, 8-bit: -a sw_striped_profile_sse41_128_8
  * SSE4.1, 8-bit with saturation check: -a sw_striped_profile_sse41_128_sat
  * AVX2, 16-bit: -a sw_striped_profile_avx2_256_16
  * AVX2, 8-bit: -a sw_striped_profile_avx2_256_8
  * AVX2, 8-bit with saturation check: -a sw_striped_profile_avx2_256_sat

The following tables show how much time it took for different sequences to be
aligned against the UniProtKB/Swiss-Prot database. All times are in seconds. The times are an average of three runs.

The following tests were performed on a MacBook Pro i5 CPU @ 2.53GHz with 8GB
RAM (SSE4.1 support). The compiler was Apple LLVM version 6.0 (clang-600.0.57.

|                                      |O74807  |P19930  |Q3ZAI3  |P18080|
|--------------------------------------|--------|--------|--------|------|
|query length|110|195|390|513|
|SSW (SSE2)|13.2|27.1|41.7|52.4|
|opal (SSE4.1)|15.2|21.5|35.9|44.6|
|SWIPE (SSSE3)|7.6|13.3|24.7|32.0|
|ssearch36 (SSE2)|12.9|20.4|29.6|38.1|
|parasail (SSE4.1) 16-bit|11.3|15.9|25.8|36.9|
|parasail (SSE4.1) sat|10.1|21.2|32.1|43.1|

The parasail saturation-checking function, like SSW, starts by trying the 8-bit precision implementation. If saturation/overflow is detected, the 8-bit function aborts early and the 16-bit implementation will then run. On this particular platform, running the 16-bit version alone is the fastest pairwise implementation. SWIPE remains fastest for this database search problem, though parasail is not far behind.

![](images/perf_mac.png)

The following tests were performed on an Intel Haswell E5-2670 v3 CPU running
at 2.3 Ghz with 64 GB 2133 Mhz DDR4 memory. The compiler used was Intel ICC
15.0.1 using level three optimization (-O3).

|                                      |O74807  |P19930  |Q3ZAI3  |P18080|
|--------------------------------------|--------|--------|--------|------|
|query length|110.0|195.0|390.0|513.0|
|SSW (SSE2)|12.0|25.6|41.5|48.7|
|opal (SSE4.1)|17.8|24.2|38.9|48.2|
|opal (AVX2)|12.2|15.4|22.9|28.4|
|SWIPE (SSSE3)|9.6|16.5|30.9|39.9|
|ssearch36 (SSE2)|11.7|21.1|30.6|37.2|
|parasail (SSE41) 16-bit|10.8|15.9|27.4|33.8|
|parasail (SSE41) sat|8.9|25.2|44.4|55.2|
|parasail (SSE41) satabort|9.2|20.7|33.8|41.5|
|parasail (SSE41) satcont|9.1|15.7|26.7|33.3|
|parasail (AVX2) 16-bit|9.3|12.0|17.4|21.5|
|parasail (AVX2) sat|9.5|20.6|29.4|35.2|
|parasail (AVX2) satabort|9.3|16.8|22.7|27.5|
|parasail (AVX2) satcont|9.6|13.5|19.5|24.0|

There are a number of alternatives to evaluate, namely what to do once the 8-bit calculation overflows. The 8-bit calculation can be allowed to finish ("sat"), abort as soon as overflow occurs ("satabort"), or abort early while copying enough state to restart from where it left off ("satcont"). Each approach has different performance characteristics.

![](images/perf_haswell.png)

