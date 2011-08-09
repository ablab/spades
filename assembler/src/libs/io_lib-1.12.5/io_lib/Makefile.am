lib_LTLIBRARIES = libstaden-read.la

libstaden_read_la_SOURCES = \
	Read.c \
	Read.h \
	scf_extras.c \
	scf_extras.h \
	translate.c \
	translate.h \
	compression.c \
	compression.h \
	ztr.c \
	ztr.h \
	ztr_translate.c \
	deflate_interlaced.c \
	deflate_interlaced.h \
	abi.h \
	fpoint.c \
	fpoint.h \
	seqIOABI.c \
	seqIOABI.h \
	seqIOALF.c \
	alf.h \
	ctfCompress.c \
	seqIOCTF.c \
	seqIOCTF.h \
	expFileIO.c \
	expFileIO.h \
	plain.h \
	seqIOPlain.c \
	misc_scf.c \
	read_scf.c \
	scf.h \
	write_scf.c \
	sff.c \
	sff.h \
	srf.c \
	srf.h \
	array.c \
	array.h \
	compress.c \
	compress.h \
	error.c \
	error.h \
	files.c \
	find.c \
	mach-io.c \
	mach-io.h \
	misc.h \
	open_trace_file.c \
	open_trace_file.h \
	read_alloc.c \
	strings.c \
	tar_format.h \
	traceType.c \
	traceType.h \
	xalloc.c \
	xalloc.h \
	vlen.c \
	vlen.h \
	hash_table.c \
	hash_table.h \
	jenkins_lookup3.c \
	jenkins_lookup3.h \
	mFILE.c \
	mFILE.h \
	stdio_hack.h \
	pooled_alloc.c \
	pooled_alloc.h

libstaden_read_la_CPPFLAGS = @LIBCURL_CPPFLAGS@

libstaden_read_la_LDFLAGS = -version-info @VERS_CURRENT@:@VERS_REVISION@:@VERS_AGE@ 
libstaden_read_la_LIBADD = @LIBZ@ @LIBCURL@

INCLUDES = -I${top_srcdir}
