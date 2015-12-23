.PHONY: dd rd

rd:
	$(MAKE) -C build/release/debruijn

dd:
	$(MAKE) -C build/debug/debruijn

rds:
	$(MAKE) -C build/release/dipspades

dds:
	$(MAKE) -C build/debug/dipspades

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

rp:	
	$(MAKE) -C build/release/corrector

dp:      
	$(MAKE) -C build/debug/corrector

rs:
	$(MAKE) -C build/release/scaffold_correction

ds:
	$(MAKE) -C build/debug/scaffold_correction

rts:
	$(MAKE) -C build/release/truseq_analysis

dts:
	$(MAKE) -C build/debug/truseq_analysis

dt:
	$(MAKE) ddt
	$(MAKE) dit

all:
	$(MAKE) -C build/release all
	$(MAKE) -C build/debug all

#temporary fix
#all: $(MAKE) -C build/release all
