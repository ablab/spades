.PHONY: dd rd

rd:
	make -C build/release/debruijn

dd:
	make -C build/debug/debruijn

rdt:
	make -C build/release/test/debruijn

ddt:
	make -C build/debug/test/debruijn

rdtls:
	make -C build/release/test/debruijn_tools

ddtls:
	make -C build/debug/test/debruijn_tools

dit:
	make -C build/debug/test/include_test

rh:
	make -C build/release/hammer hammer

dh:
	make -C build/debug/hammer hammer

rv:
	make -C build/release/online_vis

dv:
	make -C build/debug/online_vis

rc:
	make -C build/release/cap

dc:
	make -C build/debug/cap

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
