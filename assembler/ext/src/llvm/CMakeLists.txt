project(llvm-support)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
include(CheckSymbolExists)

check_include_files(execinfo.h  HAVE_EXECINFO_H)
check_include_files(signal.h    HAVE_SIGNAL_H)
check_include_files(sys/stat.h  HAVE_SYS_STAT_H)
check_include_files(cxxabi.h    HAVE_CXXABI_H)
check_include_files(dlfcn.h     HAVE_DLFCN_H)
check_include_files(mach/mach.h HAVE_MACH_MACH_H)
check_include_files(link.h      HAVE_LINK_H)

check_symbol_exists(backtrace "execinfo.h" HAVE_BACKTRACE)

find_library(DL_LIB NAMES "dl")

add_definitions(-D__STDC_CONSTANT_MACROS)
add_definitions(-D__STDC_LIMIT_MACROS)

# FIXME: Signal handler return type, currently hardcoded to 'void'
set(RETSIGTYPE void)

configure_file(${EXT_DIR}/include/llvm/Config.h.in
               ${SPADES_BUILT_INCLUDE_DIR}/llvm/Config.h)

file(GLOB sources "[a-zA-Z]*.cpp" "[a-zA-Z]*.c")

add_library(llvm-support STATIC
  ${sources})
if (DL_LIB)
  target_link_libraries(llvm-support ${DL_LIB})
endif()
