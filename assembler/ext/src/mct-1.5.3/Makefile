# This file is part of Miscellaneous Container Templates.
#
#             https://launchpad.net/libmct/
#
# Copyright (c) 2012 Paul Pogonyshev.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# TODO: Automagically add cpp<-hpp dependencies.  Not very important as it is possible to
#       build tests without this information.  If you want to alter sources and recompile
#       you are better off with SCons currently.


prefix = /usr/local


override config_files = mct/config.hpp tests/config.hpp benchmark/config.hpp


tests_CXXFLAGS = $(CXXFLAGS) -I. -DBOOST_TEST_DYN_LINK=1
ifneq "$(patsubst g++,,$(patsubst gcc,,$(patsubst clang,,$(CXX))))" "$(CXX)"
    ifeq "$(findstring -W,$(tests_CXXFLAGS))" ""
	tests_CXXFLAGS += -pedantic -Wno-long-long -Wall -Wunused-parameter
    endif
endif

tests_LDFLAGS = $(LDFLAGS) -lboost_unit_test_framework


MAIN_TEST_SOURCES = $(filter-out tests/external-validator.cpp,$(wildcard tests/*.cpp))
MAIN_TEST_OBJECTS = $(patsubst %.cpp,%.o,$(MAIN_TEST_SOURCES))


benchmark_CXXFLAGS = $(CXXFLAGS) -I.
ifeq "$(findstring -O,$(benchmark_CXXFLAGS))" ""
    benchmark_CXXFLAGS += -O2
endif


all: help

help:
	@echo "The recommended buildtool is SCons, Makefile is provided as a fallback with"
	@echo "limited functionality.  See file INSTALL for details."
	@echo
	@echo "  make install -- install to $(prefix)"
	@echo "  make install prefix=... -- install to given location"
	@echo "  make test -- build and run all available tests"
	@echo "  make benchmark -- build all available benchmarks"
	@echo "  make configure -- force reconfiguration (usually not needed)"
	@echo "  make ... CXX=... CXXFLAGS=... -- you can change typical parameters"
	@echo "  make clean -- remove all built files"
	@echo
	@echo "Note that building documentation HTML and distribution files is currently only"
	@echo "possible with SCons."


configure: ALWAYS-RUN
	-rm -f $(config_files)
	$(MAKE) -k $(config_files)


install: mct/config.hpp Makefile
	mkdir    -p $(prefix)/include/mct/impl
	install  -t $(prefix)/include/mct       mct/*.hpp
	install  -t $(prefix)/include/mct/impl  mct/impl/*.hpp


check: test
test: ALWAYS-RUN run-tests
	./run-tests --show_progress


tests/external-validator.o: tests/external-validator.cpp tests/config.hpp aux/tool-flags
	if grep "HAVE_BOOST_INTERPROCESS\\s*1" tests/config.hpp >/dev/null; then		\
	    $(CXX) $(tests_CXXFLAGS) -c $< -o $@;						\
	else											\
	    touch $@;										\
	fi

tests/external-validator: tests/external-validator.o aux/tool-flags
	if grep "HAVE_BOOST_INTERPROCESS\\s*1" tests/config.hpp >/dev/null; then		\
	    $(CXX) $(tests_LDFLAGS) tests/external-validator.o -o $@;				\
	else											\
	    touch $@;										\
	fi

tests/%.o: tests/%.cpp tests/config.hpp aux/tool-flags
	$(CXX) $(tests_CXXFLAGS) -c $< -o $@

run-tests: $(MAIN_TEST_OBJECTS) tests/external-validator aux/tool-flags
	ldflags="$(tests_LDFLAGS)";								\
	if grep "HAVE_BOOST_SERIALIZATION\\s*1" tests/config.hpp >/dev/null; then		\
	    ldflags="$$ldflags -lboost_serialization";						\
	fi;											\
	$(CXX) $$ldflags $(MAIN_TEST_OBJECTS) -o $@


benchmark: benchmark/set

benchmark/set: benchmark/set.o


benchmark/%.o: benchmark/%.cpp benchmark/config.hpp aux/tool-flags
	$(CXX) $(benchmark_CXXFLAGS) -c $< -o $@

benchmark/%: benchmark/%.o aux/tool-flags
	$(CXX) $(LDFLAGS) $< -o $@

# Disable direct compilation/linking in one step.
%: %.cpp


clean:
	-rm -f $(config_files)
	-rm -f run-tests $(MAIN_TEST_OBJECTS) tests/external-validator tests/external-validator.o
	-rm -f benchmark/set benchmark/set.o


# Without using (barely readable, yeah) functions, configuration scripts become simply
# insane.  The functions are the reason we require GNU Make.

override define die_if # condition, error_message
if test $(1); then										\
    echo "*** error ***: $(strip $(2))";							\
    exit 1;											\
fi
endef

override define configuring_parameters_header # what
echo "Configuring $(1) with:";									\
echo "  CXX:       $(CXX)";									\
echo "  CXXFLAGS:  $(CXXFLAGS)";								\
echo "  LDFLAGS:   $(LDFLAGS)";									\
echo
endef

override define yes_no_message # value, message
if test ! -z "$(strip $(2))"; then								\
    echo "  $(strip $(2)):" `echo $(1) | sed -e 's/1/yes/' -e 's/0/no/'`;			\
fi
endef

override define try_compile # result_variable, filename_base, extra_cxxflags, message
export $(1)=`$(CXX) $(CXXFLAGS) $(strip $(3)) -c $(2).cpp -o $(2).o 2>/dev/null			\
	     && echo 1 || echo 0`;								\
$(call yes_no_message,$$$(1),$(4))
endef

override define try_link # result_variable, filename_base, extra_cxxflags, extra_ldflags, message
$(call try_compile,$(1),$(2),$(3));								\
if test $$$(1) = 1; then									\
    export $(1)=`$(CXX) $(LDFLAGS) $(4) $(2).o -o $(2) 2>/dev/null && echo 1 || echo 0`;	\
fi;												\
$(call yes_no_message,$$$(1),$(5))
endef

override define check_cxx_type # result_variable, header, type, message
>aux/autogenerated_$(1).cpp;									\
if test ! -z "$(2)"; then									\
    echo "#include <$(2)>" >>aux/autogenerated_$(1).cpp;					\
fi;												\
echo "$(strip $(3))  x;" >>aux/autogenerated_$(1).cpp;						\
$(call try_compile,$(1),aux/autogenerated_$(1),,$(4))
endef

override define check_cxx_header # result_variable, header, message
echo "#include <$(2)>" >aux/autogenerated_$(1).cpp;						\
$(call try_compile,$(1),aux/autogenerated_$(1),,$(3))
endef

override define check_cxx_library # result_variable, header, library, 				\
				  # extra_cxxflags, extra_ldflags, message
echo "#include <$(2)>" > aux/autogenerated_$(1).cpp;						\
echo "int main () { }" >>aux/autogenerated_$(1).cpp;						\
$(call try_link,$(1),aux/autogenerated_$(1),$(4),-l$(3) $(5),$(6))
endef

override define check_provider # header_variable, namespace_variable, header, namespace,	\
			       # type_template, suppress_message_flag
if test -z "$$$(1)"; then									\
    __message__="$(strip $(4))::$(5) availability";						\
    if test ! -z "$(6)"; then									\
	__message__=;										\
    fi;												\
    $(call check_cxx_type,have_$(subst ::,_,$(strip $(4)))_$(5),$(strip $(3)),			\
			  $(strip $(4))::$(5) <int>,$$__message__);				\
    if test $$have_$(subst ::,_,$(4))_$(5) = 1; then						\
	export $(1)=$(strip $(3));								\
	export $(2)=$(strip $(4));								\
    fi;												\
fi
endef

override define renew_if_changed # filename_base, renew_if_any_newer
if $(foreach file,$(2),test ! $(file) -nt $(1) &&)						\
   cmp $(1) $(1).tmp >/dev/null 2>/dev/null; then						\
    rm $(1).tmp;										\
else												\
    mv $(1).tmp $(1);										\
fi
endef

override define write_configuration # filename_base, renew_if_any_newer
sed -e 's/"/\\"/g' -e 's/.*/echo "&"/' -e 's/@{/$${/g' $(1).in | sh >$(1).tmp;			\
$(call renew_if_changed,$(1),$(2));								\
echo
endef

override define tests_configuration_find_robust # concept
export have_unordered_$(1)=0;									\
export have_tr1_unordered_$(1)=0;								\
export have_boost_unordered_$(1)=0;								\
												\
$(call check_cxx_header,have_unordered_$(1),unordered_$(1), std::unordered_$(1) availability);	\
if test $$have_unordered_$(1) = 0; then								\
    $(call check_cxx_header,have_tr1_unordered_$(1),tr1/unordered_$(1),				\
			    std::tr1::unordered_$(1) availability);				\
    if test $$have_tr1_unordered_$(1) = 0; then							\
	$(call check_cxx_header,have_boost_unordered_$(1),boost/unordered_$(1).hpp,		\
				boost::unordered_$(1) availability);				\
	$(call die_if,$$have_boost_unordered_$(1) = 0,						\
		      no robust $(1) implementation to test against);				\
    fi;												\
fi
endef


mct/config.hpp: mct/config.hpp.in version aux/tool-flags
	@											\
	read version_string < version;								\
	export version_string;									\
	export major_version=`echo $$version_string | sed -e 's/\..\+//'`;			\
	export minor_version=`echo $$version_string | sed -e 's/.\+\.\(.\+\)\..\+/\1/'`;	\
	export patch_version=`echo $$version_string | sed -e 's/.\+\.//'`;			\
												\
	$(call configuring_parameters_header,MCT);						\
												\
	$(call check_provider,hash_header,hash_namespace,functional,std,hash);			\
	$(call check_provider,hash_header,hash_namespace,tr1/functional,std::tr1,hash);		\
	$(call check_provider,hash_header,hash_namespace,boost/functional/hash.hpp,boost,hash);	\
	$(call die_if,-z "$$hash_header",no known hash function provider found);		\
												\
	$(call check_provider,type_traits_header,type_traits_namespace,				\
			      type_traits,std,is_pod);						\
	$(call check_provider,type_traits_header,type_traits_namespace,				\
			      tr1/type_traits,std::tr1,is_pod);					\
	$(call check_provider,type_traits_header,type_traits_namespace,				\
			      boost/type_traits.hpp,boost,is_pod);				\
	export have_type_traits=`test -z "$$type_traits_header" && echo 0 || echo 1`;		\
												\
	$(call try_compile,cxx0x_supported,aux/check-cxx0x,,C++0x support);			\
	$(call check_cxx_type,have_long_long,,long long,Have 'long long' type);			\
												\
	$(call write_configuration,mct/config.hpp,$^)


# Note that we intentionally repeat part of main configuration.  Unlike with SCons I don't
# see a way to support old non-dynamic Boost.Test mode.
tests/config.hpp: tests/config.hpp.in mct/config.hpp aux/tool-flags
	@											\
	$(call configuring_parameters_header,tests);						\
												\
	$(call check_cxx_library,have_boost_unit_test_framework,boost/test/unit_test.hpp,	\
				 boost_unit_test_framework,-DBOOST_TEST_DYN_LINK=1,,Boost.Test);\
	$(call die_if,-z "$$have_boost_unit_test_framework",Boost.Test is required);		\
												\
	$(call check_cxx_header,have_boost_mpl,boost/mpl/vector.hpp,Boost.MPL);			\
	$(call die_if,-z "$$have_boost_mpl",Boost.MPL is required);				\
												\
	$(call tests_configuration_find_robust,set);						\
	$(call tests_configuration_find_robust,map);						\
												\
	$(call check_provider,hash_header,hash_namespace,functional,std,hash,-);		\
	$(call check_provider,hash_header,hash_namespace,tr1/functional,std::tr1,hash,-);	\
	$(call check_provider,hash_header,hash_namespace,					\
			      boost/functional/hash.hpp,boost,hash,-);				\
												\
	hash_namespace_enter=`echo $$hash_namespace						\
			      | sed -e 's/::/ /g' -e 's/\w\+/namespace & {/g'`;			\
	hash_namespace_leave=`echo $$hash_namespace | sed -e 's/::/ /g' -e 's/\w\+/}/g'`;	\
	export hash_namespace_enter hash_namespace_leave;					\
												\
	$(call check_cxx_type,have_long_long_hash_specialization,$$hash_header,			\
	    $$hash_namespace::hash <long long>,							\
	    Have 'long long'  $$hash_namespace::hash specialization);				\
												\
	$(call try_compile,have_boost_interprocess,aux/check-boost-interprocess-1-48,,		\
			   Boost.Interprocess >= 1.48);						\
	$(call check_cxx_library,have_boost_serialization,boost/serialization/serialization.hpp,\
				 boost_serialization,,,Boost.Serialization);			\
												\
	$(call write_configuration,tests/config.hpp,$^)


benchmark/config.hpp: benchmark/config.hpp.in mct/config.hpp aux/tool-flags
	@											\
	$(call configuring_parameters_header,benchmarks);					\
												\
	$(call check_cxx_header,have_unordered_set,unordered_set,				\
				std::unordered_set availability);				\
	$(call check_cxx_header,have_tr1_unordered_set,tr1/unordered_set,			\
				std::tr1::unordered_set availability);				\
	$(call check_cxx_header,have_boost_unordered_set,boost/unordered_set.hpp,		\
				boost::unordered_set availability);				\
	$(call check_cxx_header,have_google_dense_hash_set,google/dense_hash_set,		\
				google::dense_hash_set availability);				\
												\
	$(call write_configuration,benchmark/config.hpp,$^)


# This file is used to rerun configure if any of CXX, CXXFLAGS, LDFLAGS or the Makefile
# itself change.
aux/tool-flags: ALWAYS-RUN
	@											\
	echo  "$(CXX)"       > aux/tool-flags.tmp;						\
	echo  "$(CXXFLAGS)"  >>aux/tool-flags.tmp;						\
	echo  "$(LDFLAGS)"   >>aux/tool-flags.tmp;						\
	$(call renew_if_changed,aux/tool-flags,Makefile)


ALWAYS-RUN:
.PHONY: all help configure test benchmark install ALWAYS-RUN


# Local variables:
# mode: makefile-gmake
# fill-column: 90
# End:
