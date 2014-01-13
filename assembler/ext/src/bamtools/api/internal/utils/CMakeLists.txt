# ==========================
# BamTools CMakeLists.txt
# (c) 2011 Derek Barnett
#
# src/api/internal/utils
# ==========================

set( InternalUtilsDir "${InternalDir}/utils" )

set( InternalUtilsSources
        ${InternalUtilsDir}/BamException_p.cpp

        PARENT_SCOPE # <-- leave this last
)

