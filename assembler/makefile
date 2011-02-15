CC=g++
CFLAGS=-c -Wall -pedantic -Wno-long-long
LDFLAGS=

SOURCES=assembler.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=assembler

SOURCES-EC=assembler-ec.cpp
OBJECTS-EC=$(SOURCES-EC:.cpp=.o)
EXECUTABLE-EC=assembler-ec

all: $(SOURCES) $(EXECUTABLE) $(SOURCES-EC) $(EXECUTABLE-EC)
	rm -rf *.o
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(EXECUTABLE-EC): $(OBJECTS-EC) 
	$(CC) $(LDFLAGS) $(OBJECTS-EC) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -rf *.o $(EXECUTABLE) $(EXECUTABLE-EC)

