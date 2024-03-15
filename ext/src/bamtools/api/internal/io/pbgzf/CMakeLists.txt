# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal/io/pbgzf
# ==========================

set ( InternalPBGZFDir "${InternalIODir}/pbgzf" )

#--------------------------
# platform-independent PBGZF
#--------------------------
set ( CommonPBGZFSources
		${InternalPBGZFDir}/bgzf.c
		${InternalPBGZFDir}/knetfile.c
		${InternalPBGZFDir}/block.c
		${InternalPBGZFDir}/consumer.c
		${InternalPBGZFDir}/pbgzf.c
		${InternalPBGZFDir}/queue.c
		${InternalPBGZFDir}/reader.c
		${InternalPBGZFDir}/util.c
		${InternalPBGZFDir}/writer.c
)

#---------------------------
# make build-specific list
#---------------------------
set ( InternalPBGZFSources 
        ${CommonPBGZFSources} 

        PARENT_SCOPE # <-- leave this last
)

