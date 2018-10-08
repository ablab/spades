.PHONY: dd rd

rd:
	$(MAKE) -C build/release/projects/spades

dd:
	$(MAKE) -C build/debug/projects/spades

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
	$(MAKE) -C build/release/projects/hammer hammer

dh:
	$(MAKE) -C build/debug/projects/hammer hammer

rv:
	$(MAKE) -C build/release/projects/online_vis

dv:
	$(MAKE) -C build/debug/projects/online_vis

rp:	
	$(MAKE) -C build/release/projects/corrector

dp:      
	$(MAKE) -C build/debug/projects/corrector

rs:
	$(MAKE) -C build/release/projects/scaffold_correction

ds:
	$(MAKE) -C build/debug/projects/scaffold_correction

rts:
	$(MAKE) -C build/release/projects/truseq_analysis

dts:
	$(MAKE) -C build/debug/projects/truseq_analysis

dt:
	$(MAKE) ddt
	$(MAKE) dit
clean:
	$(MAKE) -C build/debug/projects/cclean
	$(MAKE) -C build/release/projects/cclean

all:
	$(MAKE) -C build/release all
	$(MAKE) -C build/debug all
