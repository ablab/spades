.PHONY: dd rd

rd:
	make -C build/release/debruijn

dd:
	make -C build/debug/debruijn

ddt:
	make -C build/debug/test/debruijn

dit:
	make -C build/debug/test/include_test

dt:
	make ddt
	make dit

clean:
	make -C build/debug clean
	make -C build/release clean

all:
	make -C build/release all
	make -C build/debug all

#temporary fix
#all: make -C build/release all
