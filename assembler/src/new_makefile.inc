PROJECT := 
DEPENDS := 

builddir=../../build/

CPPS := $(foreach project,$(DEPENDS),$(wildcard project/*.cpp))

OBJS := $(addprefix $(builddir)/$(CPPS:.cpp=.o))
LD_FLAGS := ...
CC_FLAGS := ...

all:$(EXEC)
    @echo "vvvvvvvvvv MAKE SUCCESSFUL vvvvvvvvvv"

PG=#-pg
CC=g++
CFLAGS=-c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=c++0x $(PG) $(OMP) -I./../common -I./.. $(MYCFLAGS) #-DNDEBUG
LDFLAGS=-O3 -L/usr/local/lib -L/usr/local/apr/lib -lz -llog4cxx -lapr-1 -laprutil-1 $(PG) $(MYLDFLAGS)
DEPSFL=-MM -MG


# link
$(EXEC): $(OBJS)
    $(CC) $(OBJS) $(LDFLAGS) -o $(builddir)/$(PROJECT)/$@
#   @cp $(builddir)/$@ .

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)


%.o: ../$(subst $(builddir),,%.cpp) $(dir %.cpp)/.dir
	g++ $CC_FLAGS -c -o $@ $<
    @echo $(builddir)\\ > $*.d
	$(CC) $(DEPSFL) $(CFLAGS) $< >> $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $*.d

%.dir:
    mkdir -p $*
    touch $@

.PRECIOUS: %.dir

# remove compilation products
clean:
    rm -rf ./$(builddir)

