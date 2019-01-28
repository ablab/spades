###########
# contrib #
###########

parasail_contrib_sais_incdir = include_directories('sais-lite-lcp')
parasail_contrib_sais = declare_dependency(
  include_directories : parasail_contrib_sais_incdir,
  link_with : static_library(
    'sais',
    files(['sais-lite-lcp/sais.c']))
  )

parasail_contrib_ssw_incdir = include_directories('ssw')
parasail_contrib_ssw = declare_dependency(
  include_directories : parasail_contrib_ssw_incdir,
  link_with : static_library(
    'ssw',
    files(['ssw/ssw.c']),
    include_directories : parasail_include_directories))

parasail_contrib_getMemorySize = declare_dependency(
  link_with : static_library(
    'getMemorySize',
    files(['getMemorySize.c']))
  )
