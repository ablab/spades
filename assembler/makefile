.PHONY: a d p at dt pt ct

a: 
	make -C src/abruijn
	
d:
	make -C src/debruijn

p:
	make -C src/paireddebruijn	

at: 
	make -C test/abruijn
	
dt:
	make -C test/debruijn

pt:
	make -C test/paireddebruijn	

ct:
	make -C test/common
	
clean:
	rm -rf ./build/*

all: a d p