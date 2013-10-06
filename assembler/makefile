.PHONY: dd rd

rd:
	make -C build/release/debruijn ${MAKEFLAGS}

dd:
	make -C build/debug/debruijn ${MAKEFLAGS}

rdt:
	make -C build/release/test/debruijn ${MAKEFLAGS}

ddt:
	make -C build/debug/test/debruijn ${MAKEFLAGS}

rdtls:
	make -C build/release/test/debruijn_tools ${MAKEFLAGS}

ddtls:
	make -C build/debug/test/debruijn_tools ${MAKEFLAGS}

dit:
	make -C build/debug/test/include_test ${MAKEFLAGS}

rh:
	make -C build/release/hammer hammer ${MAKEFLAGS}

dh:
	make -C build/debug/hammer hammer ${MAKEFLAGS}

rv:
	make -C build/release/online_vis ${MAKEFLAGS}

dv:
	make -C build/debug/online_vis ${MAKEFLAGS}

rc:
	make -C build/release/cap cap ${MAKEFLAGS}

dc:
	make -C build/debug/cap cap ${MAKEFLAGS}

rctls:
	make -C build/release/cap cap-tools ${MAKEFLAGS}

dctls:
	make -C build/debug/cap cap-tools ${MAKEFLAGS}


rct:
	make -C build/release/test/cap ${MAKEFLAGS}

dct:
	make -C build/debug/test/cap ${MAKEFLAGS}

dt:
	make ddt ${MAKEFLAGS}
	make dit ${MAKEFLAGS}

clean:
	make -C build/debug clean
	make -C build/release clean

all:
	make -C build/release all ${MAKEFLAGS}
	make -C build/debug all ${MAKEFLAGS}

#temporary fix
#all: make -C build/release all
