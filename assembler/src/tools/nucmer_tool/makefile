SRC=nucmer_coords_tool.cpp nucmer_tool_test.cpp
DIR=./
PROG=nucmer_tool_test

MYCFLAGS=

builddir=$(DIR)/build

EXEC=$(builddir)/$(PROG)

OBJS := $(addprefix $(builddir)/, $(notdir $(SRC:.cpp=.o)))

all:$(EXEC)
	@echo "vvvvvvvvvv MAKE SUCCESSFUL vvvvvvvvvv"

PG=#-pg
CC=g++
CFLAGS=-c -O3 -g3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=gnu++0x $(PG) $(OMP) -I./../include -I./../../ext/include -I./. $(MYCFLAGS) #-DNDEBUG
LDFLAGS=-O3 -L/usr/local/lib -lz $(PG) $(MYLDFLAGS)
DEPSFL=-MM -MG


# link
$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -o $(EXEC)
#	@cp $(builddir)/$@ .
    
# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# compile and generate dependency info;
# more complicated dependency computation, so all prereqs listed
# will also become command-less, prereq-less targets
#   sed:    strip the target (everything before colon)
#   sed:    remove any continuation backslashes
#   fmt -1: list words one per line
#   sed:    strip leading spaces
#   sed:    add trailing colons
$(builddir)/%.o: %.cpp $(builddir)/.d
	$(CC) $(CFLAGS) $< -o $@
	@echo $(builddir) > $(builddir)/$*.d.tmp
	@tr "\n" "/" <$(builddir)/$*.d.tmp >$(builddir)/$*.d
	$(CC) $(DEPSFL) $(CFLAGS) $< >> $(builddir)/$*.d
	@cp -f $(builddir)/$*.d $(builddir)/$*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $(builddir)/$*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $(builddir)/$*.d
	@rm -f $(builddir)/$*.d.tmp

   
$(builddir)/.d:
	mkdir -p $(builddir)
	touch $@        

.PRECIOUS: $(builddir)/.d

# remove compilation products
clean:
	rm -f ./$(builddir)/*


