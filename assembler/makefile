CC=g++
CFLAGS=-c -O3 -Wall -pedantic -Wno-long-long
LDFLAGS=-lz

SOURCE_DIR=src
BUILD_DIR=build
SOURCES=parser.cpp assembler.cpp
FILES=$(addprefix $(SOURCE_DIR)/,$(SOURCES))
OBJECTS=$(addprefix $(BUILD_DIR)/,$(SOURCES:.cpp=.o))
EXECUTABLE=assembler

all: $(FILES) $(EXECUTABLE)
#	rm -rf $(BUILD_DIR)/*.o
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(BUILD_DIR)/$@

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(BUILD_DIR)/.d 
	$(CC) -c $(CFLAGS) $< -o $@

$(BUILD_DIR)/.d:
	mkdir -p $@
	touch $@

clean:
	rm -rf ./$(BUILD_DIR)/*

.PRECIOUS: $(BUILD_DIR)/.d

