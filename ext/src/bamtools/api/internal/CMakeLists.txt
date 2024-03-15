# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal
# ==========================

set( InternalDir "internal" )

add_subdirectory( bam )
add_subdirectory( index )
add_subdirectory( io )
add_subdirectory( sam )
add_subdirectory( utils )

set( InternalSources
         ${InternalBamSources}
         ${InternalIndexSources}
         ${InternalIOSources}
         ${InternalSamSources}
         ${InternalUtilsSources}

         PARENT_SCOPE # <-- leave this last
   )

