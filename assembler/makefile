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
	rm -rf *.o
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(BUILD_DIR)/$@

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp 
	$(CC) -c $(CFLAGS) $< -o $@

#.cpp.o: 
#	$(CC) $(CFLAGS) $(SOURCE_DIR)/$< -o $(BUILD_DIR)/$@
	
clean:
	rm -rf *.o $(EXECUTABLE)

