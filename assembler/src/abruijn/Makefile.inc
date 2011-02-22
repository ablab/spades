
# Before including this Makefile.inc the following vars can be set:
#    sources   - *.c, *.cpp, *.y, *.l
#    subdirs   - subdirs to enter 
#    prog      - an executable that will be compiled
#    cflags    - additional CFLAGS
#    release   - version of prog={release, profile, debug}
# The following MUST be set:
#    root      - directory where this Makefile.inc lies

############################################################################

CC=g++
cflags +=-Wall
CFLAGS=-I$(root) $(cflags) -DVERSION=$(VERSION_STR)  $(ADDCFLAGS) 


DEPS=.depend
CDEPSFLAGS=-MM -MG
AR=ar qc
LD=ld -r 
LIBNAME=lib.o


##############################################################################

.PHONY: $(subdirs)

ifndef prog
lib:=$(LIBNAME)
endif

# We will remake *.c, *.cpp, *.y, *.l
objects:=$(sources:%.c=%.o)
objects:=$(objects:%.cpp=%.o)

# Libraries are taken from the subdirs
subdirlibs:=$(subdirs:%=%/$(LIBNAME))

# C/C++ sources
csources=$(filter %.c %.cpp,$(sources))

default: $(release)

debug:
	$(MAKE) all ADDCFLAGS="-DDEBUG -g"

profile:
	$(MAKE) all ADDCFLAGS="-pg -g -DPROFILE"
	
release:
	$(MAKE) all ADDCFLAGS="-O3"

all: depend $(subdirs) $(lib) $(prog)

$(lib): $(objects) $(subdirlibs)
	$(LD) -o $(lib) $(objects) $(subdirlibs)

$(prog): $(objects) $(subdirlibs) 
	$(CC) $(CFLAGS) $+ -o $@

$(subdirlibs): 

$(subdirs):
	$(MAKE) -C $@ all

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

%.c %.h: %.y
	$(BISON) $(BISONFLAGS) -o $@ $<

%.c: %.fl
	$(FLEX) $(FLEXFLAGS) -o$@ $<

clean:
	rm -f $(objects) $(lib) $(prog) $(DEPS) $(ytrash) $(fltrash) $(ggotrash)
	rm -f *.html
	for i in $(subdirs); do $(MAKE) -C $$i clean; done


ifneq "$(csources)" ""

ifneq ("$(wildcard .depend)", "")
  include $(DEPS)
endif

depend: $(DEPS)
$(DEPS): $(csources) $(HEADERS)
	@echo "Building dependencies..."
	@$(CC) $(CFLAGS) $(csources) $(CDEPSFLAGS) >$(DEPS)

else

depend:

endif
