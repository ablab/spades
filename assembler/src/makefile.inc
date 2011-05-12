builddir=../../build/$(DIR)

EXEC=$(builddir)/$(PROG)

OBJS := $(addprefix $(builddir)/, $(notdir $(SRC:.cpp=.o)))

all:$(EXEC)
	@echo "vvvvvvvvvv MAKE SUCCESSFUL vvvvvvvvvv"

PG=#-pg
CC=g++
CFLAGS=-c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=gnu++0x $(PG) $(OMP) -I./../common -I./.. $(MYCFLAGS) #-DNDEBUG
LDFLAGS=-O3 -L/usr/local/lib -L/usr/local/apr/lib -lz -llog4cxx -lapr-1 -laprutil-1 $(PG) $(MYLDFLAGS)
DEPSFL=-MM -MG


# link
$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -o $(builddir)/$@
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

$(builddir)/%.o: ../common/%.cpp $(builddir)/.d
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


