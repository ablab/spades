CC=g++
CFLAGS=-c -O3 -Wall -pedantic -Wno-long-long
LDFLAGS=-lz

SOURCE_DIR=src
BUILD_DIR=build
<<<<<<< HEAD
SOURCES=parser.cpp seq.cpp paireddebruijn/constructHashTable.cpp
HEADERS=$(addprefix $(SOURCE_DIR)/,parser.hpp seq.hpp condensed_graph.hpp paireddebruijn/common.h)
FILES=$(addprefix $(SOURCE_DIR)/,$(SOURCES))
OBJECTS=$(addprefix $(BUILD_DIR)/,$(SOURCES:.cpp=.o))
EXECUTABLE=paireddebruijn/constructHashTable

all: $(FILES) $(EXECUTABLE)

=======
SOURCES=parser.cpp seq.cpp assembler.cpp condensed_graph.cpp
HEADERS=$(addprefix $(SOURCE_DIR)/,parser.hpp seq.hpp condensed_graph.hpp)
FILES=$(addprefix $(SOURCE_DIR)/,$(SOURCES))
OBJECTS=$(addprefix $(BUILD_DIR)/,$(SOURCES:.cpp=.o))
EXECUTABLE=assembler

all: $(FILES) $(EXECUTABLE)
>>>>>>> 0cea8de4319242c968b78fa4ae6e29097e9b96f3
#	rm -rf $(BUILD_DIR)/*.o
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(BUILD_DIR)/$@

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(BUILD_DIR)/.d $(HEADERS) 
	$(CC) -c $(CFLAGS) $< -o $@

$(BUILD_DIR)/.d:
	mkdir -p $(BUILD_DIR)
	touch $@		

clean:
	rm -rf ./$(BUILD_DIR)/*

.PRECIOUS: $(BUILD_DIR)/.d

