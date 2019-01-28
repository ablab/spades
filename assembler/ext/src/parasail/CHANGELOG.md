# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

This project follows the [Gitflow Workflow model](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow).

## [Unreleased]
The Unreleased section will be empty for tagged releases. Unreleased functionality appears in the develop branch.

## [2.3] - 2018-10-23
### Added
- Support for MSYS and mingw builds.
- parasail_aligner outputs version information in verbose mode.
- Functions `parasail_result_get_traceback` and `parasail_traceback_free` to retrieve traceback C strings. See [README.md] for details.

## [2.2] - 2018-07-10
### Added
- ARM NEON vectorized functions via the [simde] project. Also thanks to [philres] for the hardware donation. 

### Fixed
- parasail_aligner works again with profile-based functions

### Closed Issues
- Support for ARM64? [\#59]

## [2.1.5] - 2018-06-21
### Fixed
- Compilation now works for arm platforms. Only non-vectorized functions at this time.

## [2.1.4] - 2018-05-29
### Fixed
- Alignments were always case-insensitive. Now cigars and tracebacks are also case-insensitive.

### Closed Issues
- matches between uppercase and lowercase symbols are treated as mismatches [\#57]

## [2.1.3] - 2018-03-28
### Changed
- Instruction sets can be selectively disabled during configure. The
  default behavior is the same as before, automatically testing for
  each. Using the --enable variant of the new configure options will
  cause configure to fail if the proper C compiler flags cannot be found
  for the compiler.
```
	  --disable-sse2          disable SSE2 support (default=auto)
	  --disable-sse4_1        disable SSE4.1 support (default=auto)
	  --disable-avx2          disable AVX2 support (default=auto)
	  --disable-avx512        disable AVX512 support (default=auto)
	  --disable-altivec       disable Altivec support (default=auto)
```

### Closed Pull Requests
- Fix automagic detection of intrinsics [\#56] ([SoapZA])

## [2.1.2] - 2018-03-20
### Fixed
- `parasail_matrix_from_file()` was failing to read input file

### Closed Issues
- Read results in non-parseable SAM alignment [\#48]
- DNA Substitution Matrices [\#24]

## [2.1.1] - 2018-03-05
### Added
- `parasail_traceback_generic_extra()` to specify index width, FILE stream
- `parasail_free_unaligned()` to free memory that wasn't allocated using
  `parasail_memalign()`, e.g., from `parasail_cigar_decode()`

### Changed
- parasail_aligner with tracebacks can now redirect to a file using `-g`.
  Default is still stdout.

### Fixed
- parasail_aligner would seg fault at the end if not producing trace output
- `parasail_traceback_generic()`
  - sequence name buffers no longer overrun
  - alignment indexes can now be longer than 7 digits
- Alignment routines that store data in a large array, such as any
  returning the DP table or traceback, now use a 64-bit offset into the
  array allowing for the alignment of longer sequences.
- parasail_traceback_generic() now properly truncates local alignment
  output (see \#55)
- parasail_aligner SAM output uses soft clipping only in local alignments

### Closed Issues
- parasail_aligner need not show any residues beyond the aligned
  segments when running a local alignment (SW) [\#55]
- Question: does "sat" function cause the solution width to go higher
  than 16? [\#54]
- Crashes of parasail_aligner running under window cause an error
  message box to be displayed [\#53]
- Traceback output does not get sent to file specified as argument to -g
  option [\#52]
- readme.md erroneously refers to "sse4" rather than "sse41" [\#51]
- test_isa informs me that avx2 is not available, but it should be [\#50]
- EMBOSS and SSW style tracebacks can put sequence and match lines
  out-of-register [\#49]
- Read results in non-parseable SAM alignment [\#48]
- Reliable segmentation fault with all traceback alignments [\#47]

## [2.1] - 2018-01-15
### Added
- parasail_aligner -b batch_size to help reduce memory overhead
- parasail_aligner can take an input file on stdin
- meson build system

### Changed
- Reduce memory used by all trace routines
- parasail_sequences_from_file(filename) can read from "stdin"

### Removed
- Generated autoconf, automake, libtool files
- parasail_aligner Intel Cilk support
- KNC ISA

### Closed Issues
- Large Memory Consumption with Traceback [\#44]
- Speedup and memory reduction of backtracing alignment [\#43]
- Add meson build system [\#39]

### Merged Pull Requests
- Add meson [\#45] ([SoapZA])

## [2.0.6] - 2018-01-11
### Fixed
- Semi-global trace functions were reporting the wrong end location

## [2.0.5] - 2018-01-05
### Fixed
- Trace functions properly align memory
- Intel compiler caused bug in 8- and 16-bit vector scan functions

### Closed Issues
- Segfault in sw_trace_striped_avx2_256_8 [\#46]

## [2.0.4] - 2017-11-30
### Fixed
- CMake add_subdirectory() of parasail project works again. Thanks to
  [armintoepfer] for the bug report.
- Preprocessor symbol clash for cigar tracebacks. [\#40]
- Patch Makefile.in to avoid automake bug during 'make check'.

### Closed Issues
- Ambiguous define with htslib [\#40]

## [2.0.3] - 2017-11-3
### Merged Pull Requests
- Fix installed includes [\#38] ([rkern])

## [2.0.2] - 2017-10-17
### Added
- [manylinux] release builds.

### Fixed
- Traceback/cigar now works for non-striped alignment functions. A significant bug caused incorrect cigar strings and tracebacks for any alignment routine besides 'striped'.  Thanks to [huxihao] for the issue report on [parasail-python].

### Closed Issues
- strdup implicitly declared in parser.c [\#37]
- make check fails without zlib [\#36]

## [2.0.1] - 2017-09-29
### Fixed
- SSW emulation seg fault when using score_size flag.

## [2.0] - 2017-09-26
### Added
- Alignment trace functions for generating SAM CIGAR output.
- SAM CIGAR encode, decode, and accessor functions.
- Support for AltiVec/POWER ISA.
- [SSW] emulation functions.
- `parasail_result_t` attribute accessor functions.
- `parasail_traceback_generic()` function for printing tracebacks to stdout.
- Revamped sequence parsing based on [kseq.h].
  - `parasail_sequence_t` and `parasail_sequences_t` objects, used with new `parasail_sequences_from_file()`
  - Support for FASTA, FASTQ.
  - Optional support for compressed input files if libz is found during the build process.
- parasail_aligner
  - Verbose mode `-v`. This re-enables output common to v1.x.
  - Output format `-O {SAMH,SAM,EMBOSS,SSW}`.
    Requires one format argument as well as the use of a trace-enabled alignment function.

### Changed
- Reduced size of `parasail_result_t` object. Users should treat the
  result as an opaque pointer and instead use the new attribute accessor
  functions.

### Deprecated
- The 'block' vectorized functions should not be used and will be removed.

### Removed
- parasail_aligner 'packed' input files.

### Fixed
- parasail_aligner now understands the stop codon.

### Closed Issues
- It appears that the `-s` option to specify SAM is not handled by parasail_aligner [\#35]
- Stop Codon Error [\#34]
- Understanding Result [\#32]
- provide functions that return the full traceback [\#12]
- smaller memory footprint of parasail_result_t [\#11]

## [1.3.1] - 2017-09-21
### Fixed
- parasail_aligner option '-s' works again.

## [1.3] - 2017-09-01
### Changed
- Added parasail_aligner option '-G' for output compatible with [GrappoloTK].
- Changed CMake option BUILD_SHARED_LIBS default to ON.
- Added automatic deployment of CI artifacts to releases.

## [1.2] - 2017-01-28

### Changed
- Added alignment function
  - parasail_nw_banded (note, different interface than the other alignment functions)
- Added matrices
  - nuc44
  - dnafull
- Added matrix functions
  - parasail_matrix_from_file
  - parasail_matrix_copy
  - parasail_matrix_set_value
- Added parasail_aligner options
  - -m matrix -- can be a built-in matrix name or a filename to be parsed
  - 'packed' input files
- parasail_matrix_t attribute `int need_free` is now `int *user_matrix`

### Closed Issues
- Needleman-Wunsch with affine gap penalties [\#26]

### Closed Pull Requests
- fix [\#29] -- wontfix

## [1.1.2] - 2016-12-07

### Fixed
- autoconf build; libparasail now correctly depends on libm where needed
- CMake build; do not incorrectly force libparasail to depend on libpthread

## [1.1.1] - 2016-11-30

### Fixed
- libparasail now correctly links when pow() not in system C library

### Merged Pull Requests
- Allow injection via cmake submodule [\#27] ([armintoepfer])

## [1.1] - 2016-08-12

### Changed
- Stats functions are now affine, not linear.
- Semi-global and global alignments now use a more negative value to
  represent negative infinity instead of half the value of the smallest
  representable integer for the given bit width.
- end_query and end_ref reported for all routines.

### Fixed
- Stats functions are now affine, not linear.

### Closed Issues
- provide Java bindings [\#22]
- stats functions should be affine, not linear [\#10]
- parasail results off by one error [\#4]

## [1.0.3] - 2016-03-25

### Changed
- Added TravisCI support for autotools Linux and OSX builds.
- Added AppVeyor support for CMake Windows builds.
- PARASAIL_API and PARASAIL_LOCAL removed from all parasail functions.
- CMake build 
  - Added BUILD_SHARED_LIBS option.
  - Added parasail.def for MSVC DLL creation.
  - Set CMAKE_POSITION_INDEPENDENT_CODE to ON if BUILD_SHARED_LIBS is ON.
  - /arch:AVX is the correct flag for MSVC, not /arch:AVX2.

### Fixed
- parasail_free() was not being used to free ISA-specific sequence profiles. Caused MSVC 64-bit library to crash.
- CMake shared library build was basically not functional on any platform. It now works.

## [1.0.2] - 2016-03-17

### Changed
- 32-bit builds replace missing functionailty.
  - SSE2 _mm_set1_epi64x, _mm_set_epi64x
  - AVX2 _mm256_set1_epi64x, _mm256_set_epi64x

### Removed
- Python bindings and pygen.py generateor were removed. Now a stand-alone project [parasail-python].

### Fixed
- Multi-arch build for OSX now correctly detects SSE4.1 and AVX2 after fixing [\#20].

### Closed Issues
- -O3 optimization causes incorrect results on OSX clang for _mm256_blendv_epi8 [\#21]
- epi64 instructions not available on 32-bit platforms [\#20]
- python ctypes interface instead of cython [\#19] **wontfix** -- moved to [parasail-python] project.
- create python wheel for pip install [\#18] **wontfix** -- moved to [parasail-python] project.
- Adding example for python [\#18] **wontfix** -- moved to [parasail-python] project.

## [1.0.1] - 2016-03-01

### Changed
- Many improvements and bug fixes to the CMake build.
  - Needed to bump CMAKE_MINIMUM_REQUIRED to VERSION 3.1 to fix static linking.
  - Visual Studio, OSX, and Linux have been verified to work.
- Windows platform natively supported.
- If an instruction set, e.g., AVX2 is not detected, then the functions are stubbed out and return NULL and set errno to ENOSYS.
- restrict keyword is conditionally preprocessed away if it's not supported by the compiler (e.g., C++, C89).  parasail internally still uses restrict if there is a suitable extension (e.g., __restrict) but this change allows greater flexibility for external libraries and applications.
- parasail_aligner application now uses long instead of int for indexing. This supports larger input datasets.

### Fixed
- Changed C++ style comments to C style to support MSVC build.
- Corrected mixed declarations and code to support MSVC build.
- Fixed various warnings from gcc -Wall -Wextra, clang, icc. MSVC build still produces many warnings.

### Closed Issues
- incorrect default SSE41_CFLAGS for gcc 4.4.7 [\#17]
- test_isa should also report what the compiler supported [\#15]
- update README et al. for new citation [\#13]
- Adding flag to disable/enable binaries in CMakeLists.txt [\#9]
- Profile thread safety? [\#8]
- Can't get parasail\_aligner to use \> 1 thread [\#7]
- Missing \#include \<string.h\> in tests [\#6]
- Documentation [\#5]
- AVX2: no such instruction [\#1](https://github.com/jeffdaily/parasail/issues/1)

## [1.0.0] - 2015-09-16
First stable, production-ready version of parasail.

[Unreleased]: https://github.com/jeffdaily/parasail/compare/v2.3...develop
[2.3]:   https://github.com/jeffdaily/parasail/compare/v2.2...v2.3
[2.2]:   https://github.com/jeffdaily/parasail/compare/v2.1.5...v2.2
[2.1.5]: https://github.com/jeffdaily/parasail/compare/v2.1.4...v2.1.5
[2.1.4]: https://github.com/jeffdaily/parasail/compare/v2.1.3...v2.1.4
[2.1.3]: https://github.com/jeffdaily/parasail/compare/v2.1.2...v2.1.3
[2.1.2]: https://github.com/jeffdaily/parasail/compare/v2.1.1...v2.1.2
[2.1.1]: https://github.com/jeffdaily/parasail/compare/v2.1...v2.1.1
[2.1]:   https://github.com/jeffdaily/parasail/compare/v2.0.6...v2.1
[2.0.6]: https://github.com/jeffdaily/parasail/compare/v2.0.5...v2.0.6
[2.0.5]: https://github.com/jeffdaily/parasail/compare/v2.0.4...v2.0.5
[2.0.4]: https://github.com/jeffdaily/parasail/compare/v2.0.3...v2.0.4
[2.0.3]: https://github.com/jeffdaily/parasail/compare/v2.0.2...v2.0.3
[2.0.2]: https://github.com/jeffdaily/parasail/compare/v2.0.1...v2.0.2
[2.0.1]: https://github.com/jeffdaily/parasail/compare/v2.0...v2.0.1
[2.0]:   https://github.com/jeffdaily/parasail/compare/v1.3.1...v2.0
[1.3.1]: https://github.com/jeffdaily/parasail/compare/v1.3...v1.3.1
[1.3]:   https://github.com/jeffdaily/parasail/compare/v1.2...v1.3
[1.2]:   https://github.com/jeffdaily/parasail/compare/v1.1.2...v1.2
[1.1.2]: https://github.com/jeffdaily/parasail/compare/v1.1.1...v1.1.2
[1.1.1]: https://github.com/jeffdaily/parasail/compare/v1.1...v1.1.1
[1.1]:   https://github.com/jeffdaily/parasail/compare/v1.0.3...v1.1
[1.0.3]: https://github.com/jeffdaily/parasail/compare/v1.0.2...v1.0.3
[1.0.2]: https://github.com/jeffdaily/parasail/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/jeffdaily/parasail/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/jeffdaily/parasail/releases/tag/v1.0.0

[\#63]: https://github.com/jeffdaily/parasail/issues/63
[\#62]: https://github.com/jeffdaily/parasail/issues/62
[\#61]: https://github.com/jeffdaily/parasail/issues/61
[\#60]: https://github.com/jeffdaily/parasail/issues/60
[\#59]: https://github.com/jeffdaily/parasail/issues/59
[\#58]: https://github.com/jeffdaily/parasail/issues/58
[\#57]: https://github.com/jeffdaily/parasail/issues/57
[\#56]: https://github.com/jeffdaily/parasail/pull/56
[\#55]: https://github.com/jeffdaily/parasail/issues/55
[\#54]: https://github.com/jeffdaily/parasail/issues/54
[\#53]: https://github.com/jeffdaily/parasail/issues/53
[\#52]: https://github.com/jeffdaily/parasail/issues/52
[\#51]: https://github.com/jeffdaily/parasail/issues/51
[\#50]: https://github.com/jeffdaily/parasail/issues/50
[\#49]: https://github.com/jeffdaily/parasail/issues/49
[\#48]: https://github.com/jeffdaily/parasail/issues/48
[\#47]: https://github.com/jeffdaily/parasail/issues/47
[\#46]: https://github.com/jeffdaily/parasail/issues/46
[\#45]: https://github.com/jeffdaily/parasail/pull/45
[\#44]: https://github.com/jeffdaily/parasail/issues/44
[\#43]: https://github.com/jeffdaily/parasail/issues/43
[\#42]: https://github.com/jeffdaily/parasail/issues/42
[\#41]: https://github.com/jeffdaily/parasail/issues/41
[\#40]: https://github.com/jeffdaily/parasail/issues/40
[\#39]: https://github.com/jeffdaily/parasail/issues/39
[\#38]: https://github.com/jeffdaily/parasail/issues/38
[\#37]: https://github.com/jeffdaily/parasail/issues/37
[\#36]: https://github.com/jeffdaily/parasail/issues/36
[\#35]: https://github.com/jeffdaily/parasail/issues/35
[\#34]: https://github.com/jeffdaily/parasail/issues/34
[\#33]: https://github.com/jeffdaily/parasail/issues/33
[\#32]: https://github.com/jeffdaily/parasail/issues/32
[\#31]: https://github.com/jeffdaily/parasail/issues/31
[\#30]: https://github.com/jeffdaily/parasail/issues/30
[\#29]: https://github.com/jeffdaily/parasail/pull/29
[\#28]: https://github.com/jeffdaily/parasail/issues/28
[\#27]: https://github.com/jeffdaily/parasail/pull/27
[\#26]: https://github.com/jeffdaily/parasail/issues/26
[\#25]: https://github.com/jeffdaily/parasail/pull/25
[\#24]: https://github.com/jeffdaily/parasail/issues/24
[\#23]: https://github.com/jeffdaily/parasail/issues/23
[\#22]: https://github.com/jeffdaily/parasail/issues/22
[\#21]: https://github.com/jeffdaily/parasail/issues/21
[\#20]: https://github.com/jeffdaily/parasail/issues/20
[\#19]: https://github.com/jeffdaily/parasail/issues/19
[\#18]: https://github.com/jeffdaily/parasail/issues/18
[\#17]: https://github.com/jeffdaily/parasail/issues/17
[\#16]: https://github.com/jeffdaily/parasail/issues/16
[\#15]: https://github.com/jeffdaily/parasail/issues/15
[\#14]: https://github.com/jeffdaily/parasail/issues/14
[\#13]: https://github.com/jeffdaily/parasail/issues/13
[\#12]: https://github.com/jeffdaily/parasail/issues/12
[\#11]: https://github.com/jeffdaily/parasail/issues/11
[\#10]: https://github.com/jeffdaily/parasail/issues/10
[\#9]: https://github.com/jeffdaily/parasail/issues/9
[\#8]: https://github.com/jeffdaily/parasail/issues/8
[\#7]: https://github.com/jeffdaily/parasail/issues/7
[\#6]: https://github.com/jeffdaily/parasail/issues/6
[\#5]: https://github.com/jeffdaily/parasail/issues/5
[\#4]: https://github.com/jeffdaily/parasail/issues/4
[\#3]: https://github.com/jeffdaily/parasail/issues/3
[\#2]: https://github.com/jeffdaily/parasail/issues/2
[\#1]: https://github.com/jeffdaily/parasail/issues/1

[README.md]: README.md
[philres]: https://github.com/philres
[simde]: https://github.com/nemequ/simde
[manylinux]: https://github.com/pypa/manylinux
[parasail-python]: https://github.com/jeffdaily/parasail-python
[huxihao]: https://github.com/huxihao
[armintoepfer]: https://github.com/armintoepfer
[rkern]: https://github.com/rkern
[SoapZA]: https://github.com/SoapZA
[GrappoloTK]: https://github.com/luhowardmark/GrappoloTK
[SSW]: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
[kseq.h]: http://lh3lh3.users.sourceforge.net/kseq.shtml
