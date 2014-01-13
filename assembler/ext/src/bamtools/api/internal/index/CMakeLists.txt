# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal/index
# ==========================

set( InternalIndexDir "${InternalDir}/index" )

set( InternalIndexSources
        ${InternalIndexDir}/BamIndexFactory_p.cpp
        ${InternalIndexDir}/BamStandardIndex_p.cpp
        ${InternalIndexDir}/BamToolsIndex_p.cpp

        PARENT_SCOPE # <-- leave this last
)

