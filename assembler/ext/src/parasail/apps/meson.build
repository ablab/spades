########
# apps #
########

executable(
  'parasail_aligner', [
     parasail_config_h,
     files('parasail_aligner.cpp')],
  install : true,
  include_directories : parasail_include_directories,
  link_with : parasail_lib,
  dependencies : [parasail_contrib_sais, parasail_contrib_getMemorySize])

executable(
  'parasail_stats', [
    parasail_config_h,
    files('parasail_stats.c')],
  install : true,
  include_directories : parasail_include_directories,
  link_with : parasail_lib)
