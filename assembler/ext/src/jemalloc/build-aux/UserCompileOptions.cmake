# User-facing configuration options exposed through CMake interfaces

function(negate var setTo)
    if(${var})
        set(${setTo} 0 PARENT_SCOPE)
    else()
        set(${setTo} 1 PARENT_SCOPE)
  endif()
endfunction()

option(JEMALLOC_GENERATE_NAMESPACE_HDR "Generate private_namespace.h" OFF)
#168-AC_ARG_WITH([xslroot],
#169:  [AS_HELP_STRING([--with-xslroot=<path>], [XSL stylesheet root path])], [
#170-if test "x$with_xslroot" = "xno" ; then
#--
#279-AC_ARG_ENABLE([cxx],
#280:  [AS_HELP_STRING([--disable-cxx], [Disable C++ integration])],
option(JEMALLOC_CXX_DISABLE "Disable C++ integration" OFF)
#--
#417-AC_ARG_WITH([lg_vaddr],
#418:  [AS_HELP_STRING([--with-lg-vaddr=<lg-vaddr>], [Number of significant virtual address bits])],
#419-  [LG_VADDR="$with_lg_vaddr"], [LG_VADDR="detect"])
set(LG_VADDR "detect"
    CACHE STRING "Number of significant virtual address bits")
#--
#548-AC_ARG_WITH([version],
#549:  [AS_HELP_STRING([--with-version=<major>.<minor>.<bugfix>-<nrev>-g<gid>],
#550-   [Version string])],
set(JEMALLOC_VERSION ""
    CACHE STRING "Version string <major>.<minor>.<bugfix>-<nrev>-g<gid>")
#--
#855-AC_ARG_WITH([rpath],
#856:  [AS_HELP_STRING([--with-rpath=<rpath>], [Colon-separated rpath (ELF systems only)])],
#set(JEMALLOC_RPATH ""
#    CACHE STRING "Colon-separated rpath (ELF systems only)")
#--
#867-AC_ARG_ENABLE([autogen],
#868:  [AS_HELP_STRING([--enable-autogen], [Automatically regenerate configure output])],
#--
#885-AC_ARG_ENABLE([shared],
#886:  [AS_HELP_STRING([--enable-shared], [Build shared libaries])],
#--
#898-AC_ARG_ENABLE([static],
#899:  [AS_HELP_STRING([--enable-static], [Build static libaries])],
#900-if test "x$enable_static" = "xno" ; then
#--
#915-AC_ARG_WITH([mangling],
#916:  [AS_HELP_STRING([--with-mangling=<map>], [Mangle symbols in <map>])],
#917-  [mangling_map="$with_mangling"], [mangling_map=""])
#--
#920-AC_ARG_WITH([jemalloc_prefix],
#921:  [AS_HELP_STRING([--with-jemalloc-prefix=<prefix>], [Prefix to prepend to all public APIs])],
#922-  [JEMALLOC_PREFIX="$with_jemalloc_prefix"],
set(JEMALLOC_PREFIX "je_"
    CACHE STRING "Prefix to prepent to all public APIs")
#--
#939-AC_ARG_WITH([export],
#940:  [AS_HELP_STRING([--without-export], [disable exporting jemalloc public APIs])],
set(JEMALLOC_EXPORT ""
    CACHE STRING "Disable exporting jemalloc public APIs")
#--
#990-AC_ARG_WITH([private_namespace],
#991:  [AS_HELP_STRING([--with-private-namespace=<prefix>], [Prefix to prepend to all library-private APIs])],
#992-  [JEMALLOC_PRIVATE_NAMESPACE="${with_private_namespace}je_"],
set(JEMALLOC_PRIVATE_NAMESPACE "je_"
    CACHE STRING "Prefix to prepend to all library-private APIs")
set(private_namespace ${JEMALLOC_PRIVATE_NAMESPACE})
#--
#1000-AC_ARG_WITH([install_suffix],
#1001:  [AS_HELP_STRING([--with-install-suffix=<suffix>], [Suffix to append to all installed files])],
#1002-  [INSTALL_SUFFIX="$with_install_suffix"],
#--
#1009-AC_ARG_WITH([malloc_conf],
#1010:  [AS_HELP_STRING([--with-malloc-conf=<malloc_conf>], [config.malloc_conf options string])],
set(JEMALLOC_CONFIG_MALLOC_CONF "\"\""
    CACHE STRING "config.malloc_conf options string")
#--
#1092-AC_ARG_ENABLE([debug],
#1093:  [AS_HELP_STRING([--enable-debug],
#1094-                  [Build debugging code])],
#--
#1127-AC_ARG_ENABLE([stats],
#1128:  [AS_HELP_STRING([--disable-stats],
#1129-                  [Disable statistics calculation/reporting])],
option(JEMALLOC_STATS_DISABLE "Disable statistics calculation/reporting" OFF)
negate(JEMALLOC_STATS_DISABLE JEMALLOC_STATS)
#--
#1144-AC_ARG_ENABLE([experimental_smallocx],
#1145:  [AS_HELP_STRING([--enable-experimental-smallocx], [Enable experimental smallocx API])],
option(JEMALLOC_SMALLOCX_ENABLE "Enable experimental smallocx API" ON)
#--
#1160-AC_ARG_ENABLE([prof],
#1161:  [AS_HELP_STRING([--enable-prof], [Enable allocation profiling])],
#--
#1176-AC_ARG_ENABLE([prof-libunwind],
#1177:  [AS_HELP_STRING([--enable-prof-libunwind], [Use libunwind for backtracing])],
#--
#1186-AC_ARG_WITH([static_libunwind],
#1187:  [AS_HELP_STRING([--with-static-libunwind=<libunwind.a>],
#1188-  [Path to static libunwind library; use rather than dynamically linking])],
#--
#1213-AC_ARG_ENABLE([prof-libgcc],
#1214:  [AS_HELP_STRING([--disable-prof-libgcc],
#1215-  [Do not use libgcc for backtracing])],
#--
#1238-AC_ARG_ENABLE([prof-gcc],
#1239:  [AS_HELP_STRING([--disable-prof-gcc],
#1240-  [Do not use gcc intrinsics for backtracing])],
#--
#1301-AC_ARG_ENABLE([fill],
#1302:  [AS_HELP_STRING([--disable-fill], [Disable support for junk/zero filling])],
option(JEMALLOC_FILL_DISABLE "Disable support for junk/zero filling" OFF)
negate(JEMALLOC_FILL_DISABLE JEMALLOC_FILL)
#--
#1317-AC_ARG_ENABLE([utrace],
#1318:  [AS_HELP_STRING([--enable-utrace], [Enable utrace(2)-based tracing])],
option(JEMALLOC_UTRACE "Enable utrace(2)-based tracking" OFF)

#--
#1345-AC_ARG_ENABLE([xmalloc],
#1346:  [AS_HELP_STRING([--enable-xmalloc], [Support xmalloc option])],
option(JEMALLOC_XMALLOC "Support xmalloc option" OFF)
#--
#1361-AC_ARG_ENABLE([cache-oblivious],
#1362:  [AS_HELP_STRING([--disable-cache-oblivious],
#1363-                  [Disable support for cache-oblivious allocation alignment])],
option(JEMALLOC_CACHE_OBLIVIOUS_DISABLE "Disable support for cache-oblivious allocation alignment" OFF)
negate(JEMALLOC_CACHE_OBLIVIOUS_DISABLE JEMALLOC_CACHE_OBLIVIOUS)
##--
#1378-AC_ARG_ENABLE([log],
#1379:  [AS_HELP_STRING([--enable-log], [Support debug logging])],
#1380-[if test "x$enable_log" = "xno" ; then
option(JEMALLOC_LOG "Enable debug logging" OFF)
#--
#1394-AC_ARG_ENABLE([readlinkat],
#1395:  [AS_HELP_STRING([--enable-readlinkat], [Use readlinkat over readlink])],
option(JEMALLOC_READLINKAT "Use readlinkat instead of readlink" OFF)
#--
#1410-AC_ARG_ENABLE([extra-size-check],
#1411:  [AS_HELP_STRING([--enable-extra-size-check],
#1412-  [Perform additonal size related sanity checks])],
option(JEMALLOC_EXTRA_SIZE_CHECK "Perform additional size sanity checks" OFF)
#--
#1496-AC_ARG_WITH([lg_quantum],
#1497:  [AS_HELP_STRING([--with-lg-quantum=<lg-quantum>],
#1498-   [Base 2 log of minimum allocation alignment])],
#--
#1505-AC_ARG_WITH([lg_page],
#1506:  [AS_HELP_STRING([--with-lg-page=<lg-page>], [Base 2 log of system page size])],
#1507-  [LG_PAGE="$with_lg_page"], [LG_PAGE="detect"])
#--
#1559-AC_ARG_WITH([lg_hugepage],
#1560:  [AS_HELP_STRING([--with-lg-hugepage=<lg-hugepage>],
#1561-   [Base 2 log of system huge page size])],
#--
#1595-AC_ARG_ENABLE([libdl],
#1596:  [AS_HELP_STRING([--disable-libdl],
#1597-  [Do not use libdl])],
option(JEMALLOC_DLSYM_DISABLE "Do not use libdl" OFF)
#--
#1713-AC_ARG_ENABLE([syscall],
#1714:  [AS_HELP_STRING([--disable-syscall], [Disable use of syscall(2)])],
if(NOT APPLE)
    option(JEMALLOC_USE_SYSCALL_DISABLE "Disable use of syscall(2)" OFF)
    negate(JEMALLOC_USE_SYSCALL_DISABLE JEMALLOC_USE_SYSCALL)
endif()
#--
#1804-AC_ARG_ENABLE([lazy_lock],
#1805:  [AS_HELP_STRING([--enable-lazy-lock],
#1806-  [Enable lazy locking (only lock when multi-threaded)])],
#--
#2078-AC_ARG_ENABLE([zone-allocator],
#2079:  [AS_HELP_STRING([--disable-zone-allocator],
#2080-                  [Disable zone allocator for Darwin])],
if(APPLE)
option(JEMALLOC_ZONE_DISABLE "Disable zone allocator for Darwin" OFF)
negate(JEMALLOC_ZONE_DISASBLE JEMALLOC_ZONE)
endif()
#--
#2103-AC_ARG_ENABLE([initial-exec-tls],
#2104:  [AS_HELP_STRING([--disable-initial-exec-tls],
#2105-                  [Disable the initial-exec tls model])],
