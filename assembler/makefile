.PHONY: dd rd

all: rd

dd:
	make -C build/debug/debruijn

rd:
	make -C build/release/debruijn

ddt:
	make -C build/debug/test/debruijn

rdt:
	make -C build/release/test/debruijn

clean:
	make -C build/debug clean
	make -C build/release clean

#temporary fix
#all: make -C build/release all
