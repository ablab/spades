CC=g++
CFLAGS=-c -Wall -pedantic -Wno-long-long
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
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(BUILD_DIR)/$@

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(BUILD_DIR) 
	$(CC) -c $(CFLAGS) $< -o $@

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUID_DIR)/*

