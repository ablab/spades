project(
  'parasail',
  ['c', 'cpp'],
  version : '2.0.4',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++11',
    'c_std=gnu11'],
  license : 'BSD-3',
  meson_version : '>= 0.42')

parasail_cc = meson.get_compiler('c')
parasail_cxx = meson.get_compiler('cpp')

################
# dependencies #
################

# zlib
parasail_zlib_dep = dependency('zlib', required : true)

# libm
parasail_libm_dep = parasail_cc.find_library('m', required : false)

# TODO:
# - Add proper OpenMP support
#   The current roadmap to integrating OpenMP into meson
#   isn't clear yet:
#     https://github.com/mesonbuild/meson/pull/1852
# - Due to its somewhat exotic nature, Altivec, AVX512
#   and other non-x86 SIMD implementations are not
#   yet integrated into Meson.

###########
# headers #
###########

subdir('parasail')

parasail_include_directories = include_directories(
  '.',
  'src',
  'contrib')

###########
# contrib #
###########

subdir('contrib')

#####################
# sources + library #
#####################

subdir('src')

########
# apps #
########

if not meson.is_subproject()
  subdir('apps')

#########
# tests #
#########

  subdir('tests')

###################
# dependency info #
###################

  import('pkgconfig').generate(
    libraries : parasail_lib,
    version : meson.project_version(),
    name : 'parasail',
    filebase : 'parasail-1',
    description : 'Pairwise Sequence Alignment Library')
endif

parasail_dep = declare_dependency(
  include_directories : parasail_include_directories,
  link_with : parasail_lib,
  dependencies : parasail_zlib_dep,
  version : meson.project_version())
