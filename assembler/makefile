.PHONY: a d p at dt pt ct

a: 
	make -C src/abruijn

d:
	make -C src/debruijn

p:
	make -C src/paireddebruijn	

b:
	make -C src/bayes

h:
	make -C src/hammer

at: 
	make -C test/abruijn

dt:
	make -C test/debruijn

pt:
	make -C test/paireddebruijn	

ct:
	make -C test/common

ht:
	make -C test/hammer

clean:
	rm -rf ./build/*

test: ct at dt ht

#temporary fix
all: a b d h test
