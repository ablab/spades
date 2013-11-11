.PHONY: dd rd

rd:
	$(MAKE) -C build/release/debruijn

dd:
	$(MAKE) -C build/debug/debruijn

rdt:
	$(MAKE) -C build/release/test/debruijn

ddt:
	$(MAKE) -C build/debug/test/debruijn

rdtls:
	$(MAKE) -C build/release/test/debruijn_tools

ddtls:
	$(MAKE) -C build/debug/test/debruijn_tools

dit:
	$(MAKE) -C build/debug/test/include_test

rit:
	$(MAKE) -C build/release/test/include_test

rh:
	$(MAKE) -C build/release/hammer hammer

dh:
	$(MAKE) -C build/debug/hammer hammer

rv:
	$(MAKE) -C build/release/online_vis

dv:
	$(MAKE) -C build/debug/online_vis

rc:
	$(MAKE) -C build/release/cap cap

dc:
	$(MAKE) -C build/debug/cap cap

rctls:
	$(MAKE) -C build/release/cap cap-tools

dctls:
	$(MAKE) -C build/debug/cap cap-tools


rct:
	$(MAKE) -C build/release/test/cap

dct:
	$(MAKE) -C build/debug/test/cap

dt:
	$(MAKE) ddt
	$(MAKE) dit

clean:
	$(MAKE) -C build/debug clean
	$(MAKE) -C build/release clean

all:
	$(MAKE) -C build/release all
	$(MAKE) -C build/debug all

#temporary fix
#all: $(MAKE) -C build/release all
