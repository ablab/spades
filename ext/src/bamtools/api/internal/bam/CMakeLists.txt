# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal/bam
# ==========================

set( InternalBamDir "${InternalDir}/bam" )

set( InternalBamSources
         ${InternalBamDir}/BamHeader_p.cpp
         ${InternalBamDir}/BamMultiReader_p.cpp
         ${InternalBamDir}/BamRandomAccessController_p.cpp
         ${InternalBamDir}/BamReader_p.cpp
         ${InternalBamDir}/BamWriter_p.cpp

         PARENT_SCOPE # <-- leave this last
   )

