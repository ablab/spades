builddir=../../build/$(DIR)

EXEC=$(builddir)/$(PROG)

OBJS := $(addprefix $(builddir)/, $(notdir $(SRC:.cpp=.o)))

all:$(EXEC)
	@echo "vvvvvvvvvv MAKE SUCCESSFUL vvvvvvvvvv"

CC=g++
CFLAGS=-c -O3 -Wall -pedantic -Wno-long-long -Wno-deprecated -std=gnu++0x -I../../src/libs/cute -I../../src/libs -I../../src/common -I../../src/common/read -I../../src/common/sequence -I../../src -I../../src/omni -I../../src/$(DIR) -I../../src/$(DIR)
LDFLAGS=-O3 -L/usr/local/lib -L/usr/local/apr/lib -lz -llog4cxx -lapr-1 -laprutil-1
DEPSFL=-MM -MG


# link
$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LDFLAGS) -o $(builddir)/$@
#   @cp $(builddir)/$@ .
    
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

$(builddir)/%.o: ../../src/$(DIR)/%.cpp $(builddir)/.d
	$(CC) $(CFLAGS) $< -o $@
	@echo $(builddir) > $(builddir)/$*.d.tmp
	@tr "\n" "/" <$(builddir)/$*.d.tmp >$(builddir)/$*.d
	$(CC) $(DEPSFL) $(CFLAGS) $< >> $(builddir)/$*.d
	@cp -f $(builddir)/$*.d $(builddir)/$*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $(builddir)/$*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $(builddir)/$*.d
	@rm -f $(builddir)/$*.d.tmp

$(builddir)/%.o: ../../src/common/%.cpp $(builddir)/.d
	$(CC) $(CFLAGS) $< -o $@
	@echo $(builddir) > $(builddir)/$*.d.tmp
	@tr "\n" "/" <$(builddir)/$*.d.tmp >$(builddir)/$*.d
	$(CC) $(DEPSFL) $(CFLAGS) $< >> $(builddir)/$*.d
	@cp -f $(builddir)/$*.d $(builddir)/$*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $(builddir)/$*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $(builddir)/$*.d
	@rm -f $(builddir)/$*.d.tmp
    
    
$(builddir)/%.o: ../../src/common/read/%.cpp $(builddir)/.d
	$(CC) $(CFLAGS) $< -o $@
	@echo $(builddir) > $(builddir)/$*.d.tmp
	@tr "\n" "/" <$(builddir)/$*.d.tmp >$(builddir)/$*.d
	$(CC) $(DEPSFL) $(CFLAGS) $< >> $(builddir)/$*.d
	@cp -f $(builddir)/$*.d $(builddir)/$*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $(builddir)/$*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $(builddir)/$*.d
	@rm -f $(builddir)/$*.d.tmp
   
   
$(builddir)/%.o: ../../src/common/sequence/%.cpp $(builddir)/.d
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
