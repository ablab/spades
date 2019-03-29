# Utilities.cmake
# Supporting functions to build Jemalloc

include(CheckTypeSize)
include(CheckCCompilerFlag)
include(CheckCSourceCompiles)

##############################################################################
# CheckTypeSizeValid
# Given C type 'type', populate lg(sizeof(type)) in '_RESULT_LG' after validing
# at least one numeric argument given after '_RESULT_LG' equals sizeof(type).
function(UtilCheckTypeSizeValid type _RESULT_LG)
    # The variable for 'check_type_size' must be unique because it's
    # globally cached across all invocations.
    check_type_size(${type} ${_RESULT_LG}foundSize)

    set(sizeOk False)

    # Iterate over vararg parameters after _RESULT_LG (the size checks)
    foreach(sizeCompare ${ARGN})
        if(${_RESULT_LG}foundSize EQUAL sizeCompare)
            set(sizeOk True)
            lg(${${_RESULT_LG}foundSize} LG_OF_FOUND_SIZE)
            set(${_RESULT_LG} ${LG_OF_FOUND_SIZE} PARENT_SCOPE)
            break()
        endif()
    endforeach()

    if (NOT sizeOk)
        message(FATAL_ERROR "Unsupported ${type} size: ${${_RESULT_LG}foundSize}")
    endif()
endfunction()

##############################################################################
# Power of two
# returns result in a VAR whose name is in RESULT_NAME
function(pow2 e RESULT_NAME)
    set(pow2_result 1)
    while(${e} GREATER 0)
        math(EXPR pow2_result "${pow2_result} + ${pow2_result}")
        math(EXPR e "${e} - 1")
    endwhile()
    set(${RESULT_NAME} ${pow2_result} PARENT_SCOPE)
endfunction()

##############################################################################
# Logarithm base 2
# returns result in a VAR whose name is in RESULT_NAME
function(lg x RESULT_NAME)
    set(lg_result 0)
    while (${x} GREATER 1)
        math(EXPR lg_result "${lg_result} + 1")
        math(EXPR x "${x} / 2")
    endwhile ()
    set(${RESULT_NAME} ${lg_result} PARENT_SCOPE)
endfunction()

##############################################################################
# Read one file and append it to another
function(AppendFileContents input output)
    file(READ ${input} buffer)
    file(APPEND ${output} "${buffer}")
endfunction()

##############################################################################
# Generate public symbols list
function(GeneratePublicSymbolsList public_sym_list mangling_map symbol_prefix output_file)
    set(JEMALLOC_RENAME_HDR "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc_rename.h")
    set(JEMALLOC_MANGLE_HDR "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc_mangle.h")
    set(JEMALLOC_MANGLE_JET_HDR "${PROJECT_BINARY_DIR}/include/jemalloc/jemalloc_mangle_jet.h")

    # If requested file already exists, don't rebuild it!
    # Note: this means if you update the ingput public_sym_list, the result
    # won't show up until you manually deelete ${output_file}
    if(EXISTS "${output_file}" AND
            ${output_file} IS_NEWER_THAN ${JEMALLOC_RENAME_HDR} AND
            ${output_file} IS_NEWER_THAN ${JEMALLOC_MANGLE_HDR} AND
            ${output_file} IS_NEWER_THAN ${JEMALLOC_MANGLE_JET_HDR})
        return()
    endif()

    file(REMOVE "${output_file}")

    # First remove from public symbols those that appear in the mangling map
    if(mangling_map)
        foreach(map_entry ${mangling_map})
            # Extract the symbol
            string(REGEX REPLACE "([^ \t]*):[^ \t]*" "\\1" sym ${map_entry})
            list(REMOVE_ITEM  public_sym_list ${sym})
            file(APPEND "${output_file}" "${map_entry}\n")
        endforeach()
    endif()  

    foreach(pub_sym ${public_sym_list})
        file(APPEND "${output_file}" "${pub_sym}:${symbol_prefix}${pub_sym}\n")
    endforeach()

    # Generate files depending on symbols list
    GenerateJemallocRename("${PUBLIC_SYM_FILE}"
                           "${JEMALLOC_RENAME_HDR}")
    GenerateJemallocMangle("${PUBLIC_SYM_FILE}" ${je_}
                           "${JEMALLOC_MANGLE_HDR}")

    # Needed for tests
    GenerateJemallocMangle("${PUBLIC_SYM_FILE}" ${JEMALLOC_PREFIX_JET} ${JEMALLOC_MANGLE_JET_HDR})
endfunction()

##############################################################################
# Decorate symbols with a prefix
#
# This is per jemalloc_mangle.sh script.
#
# IMHO, the script has a bug that is currently reflected here
# If the public symbol as alternatively named in a mangling map it is not
# reflected here. Instead, all symbols are #defined using the passed symbol_prefix
function(GenerateJemallocMangle public_sym_list symbol_prefix output_file)
    # Header
    file(WRITE "${output_file}"
        "/*\n * By default application code must explicitly refer to mangled symbol names,\n"
        " * so that it is possible to use jemalloc in conjunction with another allocator\n"
        " * in the same application.  Define JEMALLOC_MANGLE in order to cause automatic\n"
        " * name mangling that matches the API prefixing that happened as a result of\n"
        " * --with-mangling and/or --with-jemalloc-prefix configuration settings.\n"
        " */\n"
        "#ifdef JEMALLOC_MANGLE\n"
        "#  ifndef JEMALLOC_NO_DEMANGLE\n"
        "#    define JEMALLOC_NO_DEMANGLE\n"
        "#  endif\n"
        )

    file(STRINGS "${public_sym_list}" INPUT_STRINGS)

    foreach(line ${INPUT_STRINGS})
        string(REGEX REPLACE "([^ \t]*):[^ \t]*" "#  define \\1 ${symbol_prefix}\\1" output ${line})      
        file(APPEND "${output_file}" "${output}\n")
    endforeach()

    file(APPEND "${output_file}"
        "#endif\n\n"
        "/*\n"
        " * The ${symbol_prefix}* macros can be used as stable alternative names for the\n"
        " * public jemalloc API if JEMALLOC_NO_DEMANGLE is defined.  This is primarily\n"
        " * meant for use in jemalloc itself, but it can be used by application code to\n"
        " * provide isolation from the name mangling specified via --with-mangling\n"
        " * and/or --with-jemalloc-prefix.\n"
        " */\n"
        "#ifndef JEMALLOC_NO_DEMANGLE\n"
        )

    foreach(line ${INPUT_STRINGS})
        string(REGEX REPLACE "([^ \t]*):[^ \t]*" "#  undef ${symbol_prefix}\\1" output ${line})      
        file(APPEND "${output_file}" "${output}\n")
    endforeach()

    # Footer
    file(APPEND "${output_file}" "#endif\n")
endfunction()

##############################################################################
# Generate jemalloc_rename.h per jemalloc_rename.sh
function(GenerateJemallocRename public_sym_list_file file_path)
    if(EXISTS "${file_path}" AND
	   "${file_path}" IS_NEWER_THAN "${public_sym_list_file}")
        return()
    endif()

    # Header
    file(WRITE "${file_path}"
        "/*\n * Name mangling for public symbols is controlled by --with-mangling and\n"
        " * --with-jemalloc-prefix.  With" "default settings the je_" "prefix is stripped by\n"
        " * these macro definitions.\n"
        " */\n#ifndef JEMALLOC_NO_RENAME\n\n"
        )

    file(STRINGS "${public_sym_list_file}" INPUT_STRINGS)
    foreach(line ${INPUT_STRINGS})
        string(REGEX REPLACE "([^ \t]*):([^ \t]*)" "#define je_\\1 \\2" output ${line})
        file(APPEND "${file_path}" "${output}\n")
    endforeach()

    # Footer
    file(APPEND "${file_path}"
        "#endif\n"
        )
endfunction()

##############################################################################
# Create a jemalloc.h header by concatenating the following headers
# Mimic processing from jemalloc.sh
function(CreateJemallocHeader pubsyms header_list output_file)
    if(EXISTS "${output_file}")
    foreach(pub_hdr ${header_list})
        set(HDR_PATH "${PROJECT_BINARY_DIR}/include/jemalloc/${pub_hdr}")
        if(EXISTS "${HDR_PATH}" AND
                "${HDR_PATH}" IS_NEWER_THAN "${output_file}")
            # If header is newer than created header, re-create header.
            break()
        endif()

        # If we get here, then the generated pub_header is newer than all
        # source headers, so we don't need to regenerate it.
        return()
    endforeach()
    endif()

    message(STATUS "Generating API header ${output_file}")

    file(TO_NATIVE_PATH "${output_file}" ntv_output_file)

    # File Header
    file(WRITE "${ntv_output_file}"
        "#ifndef JEMALLOC_H_\n"
        "#define    JEMALLOC_H_\n"
        "#ifdef __cplusplus\n"
        "extern \"C\" {\n"
        "#endif\n\n"
        )

    foreach(pub_hdr ${header_list})
        if(False)
            message(STATUS "Copying ${pub_hdr} into public header")
        endif()
        set(HDR_PATH "${PROJECT_BINARY_DIR}/include/jemalloc/${pub_hdr}")
        file(TO_NATIVE_PATH "${HDR_PATH}" ntv_pub_hdr)
        AppendFileContents(${ntv_pub_hdr} ${ntv_output_file})
    endforeach()

    # Footer
    file(APPEND "${ntv_output_file}"
        "#ifdef __cplusplus\n"
        "}\n"
        "#endif\n"
        "#endif /* JEMALLOC_H_ */\n"
        )
endfunction()

##############################################################################
# Redefines public symbols prefxied with je_ via a macro
# Based on public_namespace.sh which echoes the result to a stdout
function(PublicNamespace public_sym_list_file output_file)
    if(EXISTS ${output_file})
        return()
    endif()

    file(STRINGS "${public_sym_list_file}" INPUT_STRINGS)
    foreach(line ${INPUT_STRINGS})
        string(REGEX REPLACE "([^ \t]*):[^ \t]*" "#define    je_\\1 JEMALLOC_N(\\1)" output ${line})
        file(APPEND ${output_file} "${output}\n")
    endforeach()
endfunction()

##############################################################################
# #undefs public je_prefixed symbols
# Based on public_unnamespace.sh which echoes the result to a stdout
function(PublicUnnamespace public_sym_list_file output_file)
    if(EXISTS ${output_file})
        return()
    endif()

    file(STRINGS "${public_sym_list_file}" INPUT_STRINGS)
    foreach(line ${INPUT_STRINGS})
        string(REGEX REPLACE "([^ \t]*):[^ \t]*" "#undef    je_\\1" output ${line})
        file(APPEND ${output_file} "${output}\n")
    endforeach()
endfunction()


##############################################################################
# Redefines a private symbol via a macro
# Based on private_namespace.sh
function(PrivateNamespace private_sym_list_file output_file)
    if(EXISTS ${output_file} AND
            ${output_file} IS_NEWER_THAN ${private_sym_list_file})
        return()
    endif()

    file(STRINGS ${private_sym_list_file} INPUT_STRINGS)
    foreach(line ${INPUT_STRINGS})
        file(APPEND ${output_file} "#define    ${line} JEMALLOC_N(${line})\n")
    endforeach()
endfunction()

##############################################################################
# Redefines a private symbol via a macro
# Based on private_namespace.sh
function(PrivateUnnamespace private_sym_list_file output_file)
    file(REMOVE "${output_file}")
    file(STRINGS ${private_sym_list_file} INPUT_STRINGS)
    foreach(line ${INPUT_STRINGS})
        file(APPEND ${output_file} "#undef ${line}\n")
    endforeach()
endfunction()


##############################################################################
# A function that configures a file_path and outputs
# end result into output_path
# ExpandDefine True/False if we want to process the file and expand
# lines that start with #undef DEFINE into what is defined in CMAKE
function(ConfigureFile file_path output_path ExpandDefine)
    if(EXISTS ${output_path} AND ${output_path} IS_NEWER_THAN ${file_path})
        # If generated file is newer than source file, don't re-generate
        # Note: this means if ONLY variables change in CMake configurations,
        # the files are ALSO not regenerated. You'll need to erase your build/
        # directory if you want CMake-level changes (new detected features, etc)
        # to show up in generated files.
        return()
    endif()

    file(REMOVE "${output_path}")

    # Convert autoconf .in files to .cmake files to generate proper .h files
    if(EXISTS ${file_path})
        get_filename_component(fname ${file_path} NAME)
        get_filename_component(oname ${output_path} NAME)
        if(NOT ${ExpandDefine})
            message(STATUS "Configuring ${fname} -> ${oname}")
            configure_file(${file_path} ${output_path} @ONLY)
        else()
            message(STATUS
                "Translating ${fname} -> ${fname}.cmake => ${oname}")

            # Quotes around the variables below are _necessary_ or CMake
            # will remove semicolons from files, which is very bad for us.
            file(READ ${file_path} CONFIG_TRANSLATE)
            string(REGEX REPLACE
                "#undef[ \t]*([^ \t\r\n]+)" "#cmakedefine \\1 @\\1@"
                CONFIG_TRANSLATED "${CONFIG_TRANSLATE}")
            file(WRITE ${output_path}.cmake "${CONFIG_TRANSLATED}")

            configure_file(${output_path}.cmake ${output_path} @ONLY)
        endif()
    else()
        message(FATAL_ERROR "${file_path} not found")
    endif()
endfunction()

##############################################################################
## Run Git and parse the output to populate version settings above
function(GetAndParseVersion)
    if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        execute_process(COMMAND ${GIT_EXECUTABLE}
            describe --long --abbrev=40
            HEAD
            WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
            OUTPUT_VARIABLE jemalloc_version)

        # Figure out version components    
        string (REPLACE "\n" "" jemalloc_version  ${jemalloc_version})
        set(jemalloc_version ${jemalloc_version} PARENT_SCOPE)
        #        message(STATUS "Version is ${jemalloc_version}")

        # replace in this order to get a valid cmake list
        string (REPLACE "-g" "-" T_VERSION ${jemalloc_version})
        string (REPLACE "-" "." T_VERSION  ${T_VERSION})
        string (REPLACE "." ";" T_VERSION  ${T_VERSION})

        list(LENGTH T_VERSION L_LEN)

        if(${L_LEN} GREATER 0)
            list(GET T_VERSION 0 jemalloc_version_major)
            set(jemalloc_version_major ${jemalloc_version_major} PARENT_SCOPE)
            #            message(STATUS "jemalloc_version_major: ${jemalloc_version_major}")
        endif()

        if(${L_LEN} GREATER 1)
            list(GET T_VERSION 1 jemalloc_version_minor)
            set(jemalloc_version_minor ${jemalloc_version_minor} PARENT_SCOPE)
            #            message(STATUS "jemalloc_version_minor: ${jemalloc_version_minor}")
        endif()

        if(${L_LEN} GREATER 2)
            list(GET T_VERSION 2 jemalloc_version_bugfix)
            set(jemalloc_version_bugfix ${jemalloc_version_bugfix} PARENT_SCOPE)
            #            message(STATUS "jemalloc_version_bugfix: ${jemalloc_version_bugfix}")
        endif()

        if(${L_LEN} GREATER 3)
            list(GET T_VERSION 3 jemalloc_version_nrev)
            set(jemalloc_version_nrev ${jemalloc_version_nrev} PARENT_SCOPE)
            #            message(STATUS "jemalloc_version_nrev: ${jemalloc_version_nrev}")
        endif()

        if(${L_LEN} GREATER 4)
            list(GET T_VERSION 4 jemalloc_version_gid)
            set(jemalloc_version_gid ${jemalloc_version_gid} PARENT_SCOPE)
            #            message(STATUS "jemalloc_version_gid: ${jemalloc_version_gid}")
        endif()
    endif()
endfunction()
