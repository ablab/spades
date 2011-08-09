AUTOMAKE_OPTIONS = foreign no-dependencies

SUBDIRS = io_lib progs tests

man_MANS = \
	man/man1/srf_index_hash.1 \
	man/man1/srf_info.1 \
	man/man1/srf_list.1 \
	man/man1/srf2fasta.1 \
	man/man1/srf2fastq.1 \
	man/man3/read_scf_header.3 \
	man/man3/read_allocate.3 \
	man/man3/ExperimentFile.3 \
	man/man3/write_reading.3 \
	man/man3/write_scf.3 \
	man/man3/read_deallocate.3 \
	man/man3/fwrite_reading.3 \
	man/man3/read2scf.3 \
	man/man3/scf2read.3 \
	man/man3/fwrite_scf.3 \
	man/man3/write_scf_header.3 \
	man/man3/read_reading.3 \
	man/man3/read_scf.3 \
	man/man3/exp2read.3 \
	man/man3/fread_reading.3 \
	man/man3/read2exp.3 \
	man/man3/fread_scf.3 \
	man/man4/Read.4

pkginclude_HEADERS = \
	io_lib/Read.h \
	io_lib/scf_extras.h \
	io_lib/translate.h \
	io_lib/abi.h \
	io_lib/fpoint.h \
	io_lib/seqIOABI.h \
	io_lib/alf.h \
	io_lib/seqIOCTF.h \
	io_lib/expFileIO.h \
	io_lib/plain.h \
	io_lib/scf.h \
	io_lib/sff.h \
	io_lib/array.h \
	io_lib/compress.h \
	io_lib/error.h \
	io_lib/mach-io.h \
	io_lib/misc.h \
	io_lib/open_trace_file.h \
	io_lib/tar_format.h \
	io_lib/traceType.h \
	io_lib/xalloc.h \
	io_lib/mFILE.h \
	io_lib/stdio_hack.h \
	io_lib/vlen.h \
	io_lib/hash_table.h \
	io_lib/jenkins_lookup3.h \
	io_lib/os.h \
	io_lib/compression.h \
	io_lib/ztr.h \
	io_lib/deflate_interlaced.h \
	io_lib/srf.h \
	io_lib/pooled_alloc.h \
	io_lib_config.h

bin_SCRIPTS = io_lib-config

EXTRA_DIST = README COPYRIGHT ChangeLog CHANGES man options.mk bootstrap \
	docs/ZTR_format docs/Hash_File_Format io_lib-config.in io_lib/os.h.in

dist-hook:
	rm -rf `find $(distdir)/man -name .svn`
