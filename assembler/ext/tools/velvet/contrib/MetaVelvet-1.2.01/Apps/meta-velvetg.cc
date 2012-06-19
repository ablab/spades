#include "meta-velvetg.hh"

int main( int argc, char* argv[] ){
  MetaVelvetG* mvg = new MetaVelvetG( argc, argv );

  cout << "[meta-velvetg] " << "Load meta-graph ..." << endl;
  MetaGraph* graph = mvg->loadGraph();
  cout << "[meta-velveth] " << "...done (load meta-graph)." << endl << endl;

  cout << "[meta-velvetg] " << "Estimate coverage parameters..." << endl;
  mvg->estimateCoverage( graph );
  cout << "[meta-velvetg] " << "...done (estimate coverage parameters)." << endl << endl;

  cout << "[meta-velvetg] " << "Remove low & high coverage nodes ..." << endl;
  mvg->removeNodes( graph );
  cout << "[meta-velvetg] " << "...done (remove low & high coverage nodes)." << endl << endl;

  cout << "[meta-velvetg] " << "Scaffolding based on paired-end information ..." << endl;
  mvg->scaffolding( graph );
  cout << "[meta-velvetg] " << "...done (scaffolding)." << endl << endl;

  cout << "[meta-velvetg] " << "Finalize assembly ... " << endl;
  mvg->finalize( graph );
  cout << "[meta-velvetg] " << "...done (finalize assembly)." << endl << endl;

  cout << "[meta-velvetg] " << "Output results ..." << endl;
  mvg->output( graph );
  cout << "[meta-velvetg] " << "...done (output results)." << endl << endl;

  cout << "[meta-velvetg] " << "All tasks were safely completed..." << endl << endl;
  delete mvg;
  return 0;
}

MetaGraph* MetaVelvetG::loadGraph() const {
  MetaGraph* graph = new MetaGraph( FileNameUtils::getSeqFileName(prefix), FileNameUtils::getInitialGraphFileName(prefix) );
  setInsertLengths( graph );
  setSplitJudge( graph );
  setSubgraphOptions( graph );
  return graph;
}

void MetaVelvetG::setInsertLengths( MetaGraph* metaGraph ) const {
  for( int cat=0 ; cat<CATEGORIES ; ++cat ){
    metaGraph->setInsertLength( cat, insertLength[cat], insertLengthSD[cat] );
    cout << "[meta-velvetg] " << "Category = 'short" << cat+1 << "'\t" 
	 << "Ave. = " << insertLength[cat] << ", SD = " << insertLengthSD[cat] << endl;
  }
  metaGraph->setInsertLength( CATEGORIES, longInsertLength, longInsertLengthSD );
  cout << "[meta-velvetg] " << "Category = 'long'" << "\t" 
       << "Ave. = " << longInsertLength << ", SD = " << longInsertLengthSD << endl;
}

void MetaVelvetG::setSplitJudge( MetaGraph* metaGraph ) const {
  SplitJudge* judge = new SplitJudge( repeatCoverageSD, minSplitLength, flagUseConnections, numValidConnections, numNoiseConnections );
  judge->setStatsWriter( FileNameUtils::getSplitStatsFileName(prefix), FileNameUtils::getSplitDetailFileName(prefix), flagReportSplitDetail );
  metaGraph->setSplitJudge( judge );
}

void MetaVelvetG::setSubgraphOptions( MetaGraph* metaGraph ) const {
  metaGraph->setSubgraphOptions( FileNameUtils::getSubgraphPrefix(prefix), flagReportSubgraph );
}

void MetaVelvetG::estimateCoverage( MetaGraph* metaGraph ){
  estimateExpectedCoverage( metaGraph );
  estimateExpectedCoverages( metaGraph );
  estimateCoverageCutoff( metaGraph );
}

void MetaVelvetG::estimateExpectedCoverages( MetaGraph* metaGraph ){
  cout << "[mate-velvetg] " << "Estimate expected coverages ... ";
  if( flagEstimateExpectedCoverages ){ 
    cout << "yes." << endl;
    metaGraph->saveStats( FileNameUtils::getHistoStatsFileName(prefix) );
    PeakDetectorParameters* params = getPeakDetectorParameters();
    PeakDetector* peakDetector = PeakDetectorFactory::instantiatePeakDetector( "simple", params );
    numCoveragePeaks = peakDetector->detectCoveragePeaks( metaGraph, expectedCoverages, coverageBoundaries );
    delete params; delete peakDetector;
    if( numCoveragePeaks <= 1 ){
      cout << endl
	   << "[meta-velvetg] " << "Warning: Can't find multiple coverage peaks." << endl
	   << "[meta-velvetg] " << "Trun on single coverage peak mode. " << endl << endl;
      numCoveragePeaks = 1;
      expectedCoverages[0] = expectedCoverage;
    }
  } else {
    cout << "no." << endl;
  }
  metaGraph->setExpectedCoverages( numCoveragePeaks, expectedCoverages );
  metaGraph->showExpectedCoverages();
}

PeakDetectorParameters* MetaVelvetG::getPeakDetectorParameters() const {
  return new PeakDetectorParameters( minPeakCoverage, maxPeakCoverage, histoBinWidth, histoSnRatio );
}

void MetaVelvetG::estimateExpectedCoverage( MetaGraph* metaGraph ){
  cout << "[mate-velvetg] " << "Estimate expected coverage ... ";
  if( flagEstimateExpectedCoverage ){ 
    cout << "yes.";
    expectedCoverage = metaGraph->estimateExpectedCoverage( prefix );
  } else {
    cout << "no.";
  }
  cout << "\t" << "Expected coverage = " << expectedCoverage << endl;
}

void MetaVelvetG::estimateCoverageCutoff( MetaGraph* metaGraph ){
  cout << "[meta-velvetg] " << "Estimate coverage cutoff ... ";
  if( flagEstimateCoverageCutoff ){
    cout << "yes.";
    coverageCutoff = expectedCoverages[numCoveragePeaks-1] / 2;
  } else {
    cout << "no.";
  }
  cout << "\t" << "Coverage cutoff = " << coverageCutoff << endl;
}

void MetaVelvetG::removeNodes( MetaGraph* metaGraph ){
  cout << "[meta-velvetg] " << "Min. coverage cutoff for short reads        = " << coverageCutoff << endl
       << "[meta-velvetg] " << "Min. coverage cutoff for long reads         = " << longCoverageCutoff << endl
       << "[meta-velvetg] " << "Max. coverage cutoff for short & long reads = " << maxCoverageCutoff << endl
       << "[meta-velvetg] " << "Min. contig length                          = " << minContigLength << endl;
  metaGraph->removeNodes( coverageCutoff, longCoverageCutoff, maxCoverageCutoff, 
			  minContigLength, flagExportFilteredNodes, 
			  FileNameUtils::getLowCoverageContigsFileName(prefix), 
			  FileNameUtils::getHighCoverageContigsFileName(prefix) );
}

void MetaVelvetG::scaffolding( MetaGraph* metaGraph ){
  metaGraph->scaffolding( flagMatePair, flagScaffolding, maxChimeraRate, flagDiscardChimera );
}

void MetaVelvetG::finalize( MetaGraph* metaGraph ){
  metaGraph->finalize( coverageCutoff, longCoverageCutoff );
}

void MetaVelvetG::output( MetaGraph* metaGraph ){
  metaGraph->save( FileNameUtils::getScaffoldFileName(prefix), minContigLength, coverageMask );
  metaGraph->saveStats( FileNameUtils::getStatsFileName(prefix) );
  metaGraph->save( FileNameUtils::getLastGraphFileName(prefix) );
  if( flagExportAlignments ){ metaGraph->saveAlignments( FileNameUtils::getAlignmentFileName(prefix), minContigLength, FileNameUtils::getSeqFileName(prefix) ); }
  if( flagExportAmos ){ metaGraph->saveAmos( FileNameUtils::getAmosFileName(prefix), minContigLength ); }
  if( flagExportUnusedReads ){ metaGraph->saveUnusedReads( prefix, minContigLength ); }
}

MetaVelvetG::MetaVelvetG( int argc, char* argv[] ){
  setDefaultParameters();
  cout << endl<< endl
       << "[meta-velvetg] " << "Check command line options..." << endl;
  if( !checkParameters( argc, argv ) or !setParameters( argc, argv ) or !checkParameters() ){
    cout << "\t" << "Your command line options are not appropriate." << endl
	 << "\t" << "Please read usage of meta-velvetg carefully." << endl
	 << endl
	 << "\t" << "Thanks, " << endl
	 << endl << endl;
    printUsage();
    exit(-1);
  }
  cout << "\t" << "OK. Your command line options seem to be good." << endl
       << endl << endl;
  guessNonMandatoryParameters();
}

void MetaVelvetG::setDefaultParameters(){
  // velvet parameters
  flagEstimateCoverageCutoff   = true;
  flagEstimateExpectedCoverage = true;
  flagReadTracking             = true;
  flagScaffolding              = true;
  flagExportUnusedReads        = false;
  flagExportFilteredNodes      = false;
  flagExportAmos               = false;
  flagExportAlignments         = false;
  accelerationBits   = 24;
  minContigLength    = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  coverageMask       = 1;
  coverageCutoff     = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  maxCoverageCutoff  = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  longCoverageCutoff = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  expectedCoverage   = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  longInsertLength   = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  longInsertLengthSD = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  for ( int cat = 0; cat < CATEGORIES; ++cat) {
    insertLength[cat]   = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
    insertLengthSD[cat] = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
    flagMatePair[cat]   = false;
  }
  
  // graph-splitting parameters
  flagDiscardChimera  = false;
  maxChimeraRate      = META_VELEVTG_DEFAULT_MAX_CHIMERA_RATE;
  repeatCoverageSD    = META_VELVETG_DEFAULT_REPEAT_COVERAGE_SD;
  minSplitLength      = META_VELVETG_DEFAULT_MIN_SPLIT_LENGTH;
  flagUseConnections  = true;
  numValidConnections = META_VELVETG_DEFAULT_NUM_VALID_CONNECTIONS;
  numNoiseConnections = META_VELVETG_DEFAULT_NUM_NOISE_CONNECTIONS;
  flagReportSplitDetail = false;
  flagReportSubgraph  = false;
  
  // peak parameters
  flagEstimateExpectedCoverages = true;
  numCoveragePeaks   = 1;
  for( int i=0 ; i<META_VELVET_G_MAX_NUM_COVERAGE_PEAKS ; ++i ){
    expectedCoverages[i] = META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED;
  }
  minPeakCoverage = META_VELVET_G_DEFAULT_MIN_PEAK_COVERAGE;
  maxPeakCoverage = META_VELVET_G_DEFAULT_MAX_PEAK_COVERAGE;
  histoBinWidth   = META_VELVET_G_DEFAULT_HISTO_BIN_WIDTH;
  histoSnRatio    = META_VELVET_G_DEFAULT_HISTO_SN_RATIO;
}

bool MetaVelvetG::checkParameters() const {
  return checkMetaHistoParameters();
}

bool MetaVelvetG::checkMetaHistoParameters() const {
  PeakDetectorParameters* param = getPeakDetectorParameters();
  bool checkResult = param->checkParameters();
  delete param;
  return checkResult;
}

bool MetaVelvetG::checkParameters( int argc, char* argv[] ) const {
  if( checkArgc(argc) and checkPrefix(argc, argv) ){
    return true;
  }
  return false;
}

bool MetaVelvetG::checkArgc( int argc ) const {
  if( argc < 2 ){
    return false;
  }
  return true;
}

bool MetaVelvetG::checkPrefix( int argc, char* argv[] ) const {
  if( strncmp( argv[1], "-", 1 ) == 0 ){
    cerr << "[meta-velvetg] Error: the first argument should be a directory path (should not be start with '-')" << endl;
    return false;
  }
  return checkMandatoryFiles( (string)argv[1] );
}

bool MetaVelvetG::checkMandatoryFiles( string prefix ) const {
  if( !FileUtils::isValidFile(FileNameUtils::getSeqFileName(prefix)) ){
    cerr << "[meta-velvetg] Error: There is no 'Sequences' file in the directory: " << prefix << endl
	 << "\t" << "You should run 'velveth' before running 'meta-velvetg'" << endl << endl;
    return false;
  }
  if( !FileUtils::isValidFile(FileNameUtils::getRoadmapFileName(prefix)) ){
    cerr << "[meta-velvetg] Error: There is no 'Roadmaps' file in the directory: " << prefix << endl
	 << "\t" << "You should run 'velveth' before running 'meta-velvetg'" << endl << endl;
    return false;
  }
  if( !FileUtils::isValidFile(FileNameUtils::getInitialGraphFileName(prefix)) ){
    cerr << "[meta-velvetg] Error: There is no 'Graph2' file in the directory: " << prefix << endl
	 << "\t" << "You should run 'velvetg' with '-read_trkg' opton before runing 'meta-velvetg'" << endl << endl;
    return false;
  }
  return true;
}

void MetaVelvetG::printUsage() const {
  cerr << endl << endl
       << "************************************************************" << endl
       << "Program: ----------------------------------------------" << endl
       << "  " << "meta-velvetg - contiging and scaffolding program for metagenomics NGS data" << endl
       << "\t" << "Version = " << META_VELVET_VERSION << endl
       << "\t" << "MAX_CATEGORIES = " << CATEGORIES << endl
       << "\t" << "MAX_KMER_LENGTH = " << MAXKMERLENGTH << endl
       << endl << endl
       << "Usage: ------------------------------------------------" << endl
       << "  " << "meta-velvetg directory [options]" << endl
       << "\t" << "directory                   " << "\t: directory name for output files" << endl
       << endl
       << "  " << "Graph-splitting options (metagenome-specific):" << endl
       << "\t" << "-discard_chimera <yes|no>    " << "\t: discard chimera sub-graph (default: no)" << endl
       << "\t" << "-max_chimera_rate <double>   " << "\t: maximum allowable chimera rate (default: 0.0)" << endl
       << "\t" << "-repeat_cov_sd               " << "\t: standard deviation of repeat node coverages (default: 0.1)" << endl
       << "\t" << "-min_split_length <int>      " << "\t: minimum node length required for repeat resolution (default: 0)" << endl
       << "\t" << "-valid_connections <int>     " << "\t: minimum allowable number of consistent paired-end connections (default: 1)" << endl
       << "\t" << "-noise_connections <int>     " << "\t: maximum allowable number of inconsistent paired-end connections (default: 0)" << endl
       << "\t" << "-use_connections <yes|no>    " << "\t: use paired-end connections for graph splitting (default: yes)" << endl
       << "\t" << "-report_split_detail <yes|no>" << "\t: report sequences around repeat nodes (default: no)" << endl
       << "\t" << "-report_subgraph <yes|no>    " << "\t: report node sequences for each subgraph (default: no)" << endl
       << endl
       << "  " << "Peak detection options (metagenome-specific):" << endl
       << "\t" << "-exp_covs <string|auto>      " << "\t: expected coverages for each species in microbiome (default: auto)" << endl
       << "\t" << "                             " << "\t    ex) -exp_covs 214_122_70_43_25_13.5" << endl
       << "\t" << "                             " << "\t    coverage values should be sorted in a descending order" << endl
       << "\t" << "-min_peak_cov <double>       " << "\t: minimum peak coverage (default: " << META_VELVET_G_DEFAULT_MIN_PEAK_COVERAGE << ")" << endl
       << "\t" << "-max_peak_cov <double>       " << "\t: maximum peak coverage (default: " << META_VELVET_G_DEFAULT_MAX_PEAK_COVERAGE << ")" << endl
       << "\t" << "-histo_bin_width <double>    " << "\t: bin width of peak coverage histogram (default: " << META_VELVET_G_DEFAULT_HISTO_BIN_WIDTH << ")" << endl
       << "\t" << "-histo_sn_ratio <double>     " << "\t: signal-noise ratio to remove peak noises (default: " << META_VELVET_G_DEFAULT_HISTO_SN_RATIO << ")" << endl
       << endl
       << "  " << "Contiging options: (common to single-genome)" << endl
       << "\t" << "-cov_cutoff <double|auto>    " << "\t: removal of low coverage nodes AFTER tour bus or allow the system to infer it (default: auto)" << endl
       << "\t" << "-max_coverage <double>       " << "\t: removal of high coverage nodes AFTER tour bus (default: no removal)" << endl
       << "\t" << "-long_cov_cutoff <double>    " << "\t: removal of nodes with low long-read coverage AFTER tour bus (default: no removal)" << endl
    //<< "\t" << "-read_trkg <yes|no>          " << "\t: tracking of short read positions in assembly (default: yes)" << endl
       << "\t" << "-max_branch_length <int>     " << "\t: maximum length in base pair of bubble (default: 100)" << endl
       << "\t" << "-max_divergence <double>     " << "\t: maximum divergence rate between two branches in a bubble (default: 0.2)" << endl
       << "\t" << "-max_gap_count <int>         " << "\t: maximum number of gaps allowed in the alignment of the two branches of a bubble (default: 3)" << endl
       << "\t" << "-min_contig_lgth <int>       " << "\t: minimum contig length exported to contigs.fa file (default: hash length * 2)" << endl
       << endl
       << "  " << "Scaffolding options: (common to single-genome)" << endl
       << "\t" << "-scaffolding <yes|no>        " << "\t: scaffolding of contigs used paired end information (default: on)" << endl
       << "\t" << "-exp_cov <double|auto>       " << "\t: expected coverage of unique regions or allow the system to infer it (default: auto)" << endl
       << "\t" << "-ins_length <double>         " << "\t: expected distance between two paired end reads for the category (default: no read pairing)" << endl
       << "\t" << "-ins_length_sd <double>      " << "\t: standard deviation of insert length for the category (default: insert length / 10)" << endl
       << "\t" << "-ins_length2 <double>        " << "\t: expected distance between two paired end reads for the category (default: no read pairing)" << endl
       << "\t" << "-ins_length2_sd <double>     " << "\t: standard deviation of insert length for the category (default: insert length / 10)" << endl
       << "\t" << "-ins_length_long <double>    " << "\t: expected distance between two long paired-end reads (default: no read pairing)" << endl
       << "\t" << "-ins_length_long_sd <double> " << "\t: standard deviation of insert length for the category" << endl
       << "\t" << "-min_pair_count <int>        " << "\t: minimum number of paired end connections to justify the scaffolding of two long contigs (default: 5)" << endl
       << "\t" << "-long_mult_cutoff <int>      " << "\t: minimum number of long reads required to merge contigs (default: 2)" << endl
       << "\t" << "-shortMatePaired <yes|no>    " << "\t: for mate-pair libraries, indicate that the library might be contaminated with paired-end reads (default no)" << endl
       << endl
       << "  " << "Output options: (common to single-genome)" << endl
       << "\t" << "-amos_file <yes|no>          " << "\t: export assembly to AMOS file (default: no export)" << endl
       << "\t" << "-coverage_mask <int>         " << "\t: minimum coverage required for confident regions of contigs (default: 1)" << endl
       << "\t" << "-unused_reads <yes|no>       " << "\t: export unused reads in UnusedReads.fa file (default: no)" << endl
       << "\t" << "-alignments <yes|no>         " << "\t: export a summary of contig alignment to the reference sequences (default: no)" << endl
       << "\t" << "-exportFiltered <yes|no>     " << "\t: export the long nodes which were eliminated by the coverage filters (default: no)" << endl
       << "\t" << "-paired_exp_fraction <double>" << "\t: remove all the paired end connections which less than the specified fraction of the expected count (default: 0.1)" << endl
       << endl << endl
       << "Output: -----------------------------------------------" << endl
       << "\t" << "directory/meta-velvetg.contigs.fa       \t: fasta file of contigs longer than twice hash length" << endl
       << "\t" << "directory/meta-velvetg.LastGraph        \t: special formatted file with all the information on the final graph" << endl
       << "\t" << "directory/meta-velvetg.Graph2-stats.txt \t: stats file (tab-delimited) useful for optimizing coverage peak values" << endl
       << "\t" << "directory/meta-velvetg.split-stats.txt  \t: stats file (tab-delimited) useful for optimizing graph-splitting parameters" << endl
       << endl 
       << "************************************************************" << endl
       << endl << endl;
}

bool MetaVelvetG::setParameters( int argc, char* argv[] ){
  prefix = argv[1];
  for ( int arg_index=2 ; arg_index<argc ; arg_index++ ) {
    char* arg = argv[arg_index++];
    if ( arg_index >= argc ) {
      cerr << "[meta-velvetg] Error: Unusual number of arguments!" << endl;
      return false;
    }
    if( !setParameter( arg,  argv[arg_index] ) ){
      return false;
    }
  }
  return true;
}

bool MetaVelvetG::setParameter( char* name, char* val ){
  if( strcmp(name, "-cov_cutoff")           == 0 ){ return setCoverageCutoff( val ); }
  if( strcmp(name, "-long_cov_cutoff")      == 0 ){ return setLongCoverageCutoff( val ); }
  if( strcmp(name, "-exp_cov")              == 0 ){ return setExpectedCoverage( val ); }
  if( strcmp(name, "-ins_length")           == 0 ){ return setInsertLength( 0, val ); }
  if( strcmp(name, "-ins_length_sd")        == 0 ){ return setInsertLengthSD( 0, val ); }
  if( strcmp(name, "-ins_length_long")      == 0 ){ return setLongInsertLength( val ); }
  if( strcmp(name, "-ins_length_long_sd")   == 0 ){ return setLongInsertLengthSD( val ); }
  if( strncmp(name, "-ins_length", 11)      == 0 && strchr(name, 'd') == NULL ){int cat = 0; sscanf(name, "-ins_length%d", &cat); return setInsertLength( cat-1, val );}
  if( strncmp(name, "-ins_length", 11)      == 0 ){ int cat = 0; sscanf(name, "-ins_length%d_sd", &cat); return setInsertLengthSD( cat-1, val ); }
  if( strcmp(name, "-read_trkg")            == 0 ){ return setFlagReadTracking( val ); }
  if( strcmp(name, "-scaffolding")          == 0 ){ return setFlagScaffolding( val ); }
  if( strcmp(name, "-exportFiltered")       == 0 ){ return setFlagExportFilteredNodes( val ); }
  if( strcmp(name, "-amos_file")            == 0 ){ return setFlagExportAmos( val ); }
  if( strcmp(name, "-alignments")           == 0 ){ return setFlagExportAlignments( val ); }
  if( strcmp(name, "-min_contig_lgth")      == 0 ){ return setMinContigLength( val ); }
  if( strcmp(name, "-coverage_mask")        == 0 ){ return setCoverageMask( val ); }
  if( strcmp(name, "-accel_bits")           == 0 ){ return setAccelerationBits( val ); }
  if( strcmp(name, "-max_branch_length")    == 0 ){ return setMaxBranchLength( val ); }
  if( strcmp(name, "-max_divergence")       == 0 ){ return _setMaxDivergence( val ); }
  if( strcmp(name, "-max_gap_count")        == 0 ){ return setMaxGapCount( val ); }
  if( strcmp(name, "-min_pair_count")       == 0 ){ return setMinPairCount( val ); }
  if( strcmp(name, "-max_coverage")         == 0 ){ return setMaxCoverage( val ); }
  if( strcmp(name, "-long_mult_cutoff")     == 0 ){ return setLongMultiCutoff( val ); }
  if( strcmp(name, "-paired_exp_fraction")  == 0 ){ return _setPairedExpFraction( val ); }
  if( strcmp(name, "-unused_reads")         == 0 ){ return setFlagExportUnusedReads(val); }
  if( strcmp(name, "-shortMatePaired")      == 0 ){ return setFlagMatePair( 0, val ); }
  if( strncmp(name, "-shortMatePaired", 16) == 0 ){ int cat = 0; sscanf(name, "-shortMatePaired%i", &cat); return setFlagMatePair( cat-1, val ); }
  if( strcmp(name, "-discard_chimera")      == 0 ){ return setFlagDiscardChimera( val ); }
  if( strcmp(name, "-max_chimera_rate")     == 0 ){ return setMaxChimeraRate( val ); }
  if( strcmp(name, "-exp_covs")             == 0 ){ return setExpectedCoverages( val ); }
  if( strcmp(name, "-manual_exp_cov_multi") == 0 ){ return setExpectedCoverages( val ); }
  if( strcmp(name, "-repeat_cov_sd")        == 0 ){ return setRepeatCoverageSD( val ); }
  if( strcmp(name, "-min_split_length")     == 0 ){ return setMinSplitLength( val ); }
  if( strcmp(name, "-min_peak_cov")         == 0 ){ return setMinPeakCoverage( val ); }
  if( strcmp(name, "-max_peak_cov")         == 0 ){ return setMaxPeakCoverage( val ); }
  if( strcmp(name, "-histo_bin_width")      == 0 ){ return setHistoBinWidth( val ); }
  if( strcmp(name, "-histo_sn_ratio")       == 0 ){ return setHistoSnRatio( val ); }
  if( strcmp(name, "-valid_connections")    == 0 ){ return setNumValidConnections( val ); }
  if( strcmp(name, "-noise_connections")    == 0 ){ return setNumNoiseConnections( val ); }
  if( strcmp(name, "-use_connections")      == 0 ){ return setFlagUseConnections( val ); }
  if( strcmp(name, "-report_split_detail")  == 0 ){ return setFlagReportSplitDetail( val ); }
  if( strcmp(name, "-report_subgraph")      == 0 ){ return setFlagReportSubgraph( val ); }
  if( strcmp(name, "--help")                == 0 ){ return false; }
  cerr << "[meta-velvetg] Unknown option: " << name << endl;
  return false;
}

bool MetaVelvetG::setCoverageCutoff( char* val ){
  if ( strcmp(val, "auto") == 0 ){
    flagEstimateCoverageCutoff = true;
  } else {
    flagEstimateCoverageCutoff = false;
    coverageCutoff = atof( val );
  }
  return true;
}

bool MetaVelvetG::setLongCoverageCutoff( char* val ){
  longCoverageCutoff = atof( val );
  return true;
}

bool MetaVelvetG::setExpectedCoverage( char* val ){
  if( strcmp(val, "auto") == 0 ){
    flagEstimateExpectedCoverage = true;
  } else {
    flagEstimateExpectedCoverage = false;
    expectedCoverage = atof( val );
  }
  return true;
}

bool MetaVelvetG::setInsertLength( int cat, char* val ){
  if( !checkCategory( cat ) ){ return false; }
  insertLength[cat] = atof( val );
  if( insertLength[cat] <= 0 ){
    cerr << "[meta-velvetg] Error: Invalid insert length for category " << cat << ": " << insertLength[cat] << endl;
    return false;
  }
  return true;
}

bool MetaVelvetG::setInsertLengthSD( int cat, char* val ){
  if( !checkCategory( cat ) ){ return false; }
  insertLengthSD[cat] = atof( val );
  if( insertLengthSD[cat] <= 0 ){
    cerr << "[meta-velvetg] Error: Invalid insert length SD for category " << cat << ": " << insertLengthSD[cat] << endl;
    return false;
  }
  return true;
}

bool MetaVelvetG::setLongInsertLength( char* val ){
  longInsertLength = atof( val );
  if( longInsertLength <= 0 ){
    cerr << "[meta-velvetg] Error: Invalid insert length for long reads: " << longInsertLength << endl;
    return false;
  }
  return true;
}

bool MetaVelvetG::setLongInsertLengthSD( char* val ){
  longInsertLengthSD = atof( val );
  if( longInsertLengthSD <= 0 ){
    cerr << "[meta-velvetg] Error: Invalid insert length SD for long reads: " << longInsertLengthSD << endl;
    return false;
  }
  return true;
}

bool MetaVelvetG::setFlagReadTracking( char* val ){
  flagReadTracking = (strcmp(val, "yes") == 0);
  return true;
}

bool MetaVelvetG::setFlagScaffolding( char* val ){
  flagScaffolding = (strcmp(val, "yes") == 0);
  return true;
}

bool MetaVelvetG::setFlagExportFilteredNodes( char* val ){
  flagExportFilteredNodes = (strcmp(val, "yes") == 0);
  return true;
}

bool MetaVelvetG::setFlagExportAmos( char* val ){
  flagExportAmos = (strcmp(val, "yes") == 0);
  return true;
}

bool MetaVelvetG::setFlagExportAlignments( char* val ){
  flagExportAlignments = (strcmp(val, "yes") == 0);
  return true;
}

bool MetaVelvetG::setMinContigLength( char* val ){
  minContigLength = atoi( val );
  return true;
}

bool MetaVelvetG::setCoverageMask( char* val ){
  coverageMask = atof( val );
  return true;
}

bool MetaVelvetG::setAccelerationBits( char* val ){
  accelerationBits = atoi( val );
  if( accelerationBits < 0 ){
    cerr << "[meta-velvetg] Error: Illegal acceleration parameter: " << accelerationBits << endl;
    return false;
  }
  return true;
}

bool MetaVelvetG::setMaxBranchLength( char* val ){
  setMaxReadLength( atoi(val) );
  setLocalMaxReadLength( atoi(val) );
  return true;
}

bool MetaVelvetG::_setMaxDivergence( char* val ){
  setMaxDivergence( atof(val) );
  setLocalMaxDivergence( atof(val) );
  return true;
}

bool MetaVelvetG::setMaxGapCount( char* val ){
  setMaxGaps( atoi(val) );
  setLocalMaxGaps( atoi(val) );
  return true;
}

bool MetaVelvetG::setMinPairCount( char* val ){
  setUnreliableConnectionCutoff( atoi(val) );
  return true;
}

bool MetaVelvetG::setMaxCoverage( char* val ){
  maxCoverageCutoff = atof( val );
  return true;
}

bool MetaVelvetG::setLongMultiCutoff( char * val ){
  setMultiplicityCutoff( atoi(val) );
  return true;
}

bool MetaVelvetG::_setPairedExpFraction( char * val ){
  setPairedExpFraction( atof(val) );
  return true;
}

bool MetaVelvetG::setFlagExportUnusedReads( char* val ){
  if( strcmp(val, "yes") == 0 ){
    flagExportUnusedReads = true;
    flagReadTracking = true;
  }
  return true;
}

bool MetaVelvetG::setFlagMatePair( int cat, char* val ){
  checkCategory( cat );
  flagMatePair[cat] = ( strcmp(val, "yes") == 0 );
  return true;
}

bool MetaVelvetG::setFlagDiscardChimera( char* val ){
  flagDiscardChimera = ( strcmp(val, "yes") == 0 );
  return true;
}

bool MetaVelvetG::setMaxChimeraRate( char* val ){
  maxChimeraRate = atof( val );
  return true;
}

bool MetaVelvetG::setExpectedCoverages( char* val ){
  flagEstimateExpectedCoverages = false;
  char* tokenPointer = val;
  char* tmp[META_VELVET_G_MAX_NUM_COVERAGE_PEAKS];
  for( numCoveragePeaks=0; numCoveragePeaks<META_VELVET_G_MAX_NUM_COVERAGE_PEAKS ; ++numCoveragePeaks ){
    if( (tmp[numCoveragePeaks] = strtok(tokenPointer, META_VELVET_G_COVERAGE_DELIM_CHAR)) != NULL ){
      expectedCoverages[numCoveragePeaks] = atof( tmp[numCoveragePeaks] );
    } else {
      break;
    }
    tokenPointer = NULL;
  }
  return checkCoverageOrder();
}

bool MetaVelvetG::setRepeatCoverageSD( char* val ){
  repeatCoverageSD = atof(val);
  return true;
}

bool MetaVelvetG::setMinSplitLength( char* val ){
  minSplitLength = atoi(val);
  return true;
}

bool MetaVelvetG::setMinPeakCoverage( char* val ){
  minPeakCoverage = atof(val);
  return true;
}

bool MetaVelvetG::setMaxPeakCoverage( char* val ){
  maxPeakCoverage = atof(val);
  return true;
}

bool MetaVelvetG::setHistoBinWidth( char* val ){
  histoBinWidth = atof(val);
  return true;
}

bool MetaVelvetG::setHistoSnRatio( char* val ){
  histoSnRatio = atof(val);
  return true;
}

bool MetaVelvetG::setNumValidConnections( char* val ){
  numValidConnections = atoi(val);
  return true;
}

bool MetaVelvetG::setNumNoiseConnections( char* val ){
  numNoiseConnections = atoi(val);
  return true;
}

bool MetaVelvetG::setFlagUseConnections( char* val ){
  flagUseConnections = ( strcmp(val, "yes") == 0 );
  return true;
}

bool MetaVelvetG::setFlagReportSplitDetail( char* val ){
  flagReportSplitDetail = ( strcmp(val, "yes") == 0 );
  return true;
}

bool MetaVelvetG::setFlagReportSubgraph( char* val ){
  flagReportSubgraph = ( strcmp(val, "yes") == 0 );
  return true;
}

bool MetaVelvetG::checkCategory( int cat ) const {
  if( (cat<0) or (cat>CATEGORIES) ){
    cerr << "[meta-velvetg] Error: Invalid category: " << cat << endl;
    return false;
  }
  return true;
}

bool MetaVelvetG::checkCoverageOrder() const {
  for( int i=1 ; i<numCoveragePeaks ; ++i ){
    if( expectedCoverages[i-1] <= expectedCoverages[i] ){
      cerr << "[meta-velvetg] Error: expected coverages should be sorted in a descending order: " 
	   << expectedCoverages[i-1] << " <-> " << expectedCoverages[i] << endl;
      return false;
    }
  }
  return true;
}

void MetaVelvetG::guessNonMandatoryParameters(){
  guessInsertLengthSD();
  guessCoverageFlags();
}

void MetaVelvetG::guessInsertLengthSD(){
  for( int cat=0 ; cat<CATEGORIES ; ++cat){
    if( isGuessInsertLengthSD( insertLength[cat], insertLengthSD[cat] ) )
      insertLengthSD[cat] = guessInsertLengthSD( insertLength[cat] );
  }
  if( isGuessInsertLengthSD( longInsertLength, longInsertLengthSD ) )
    longInsertLengthSD = guessInsertLengthSD( longInsertLength );
}

bool MetaVelvetG::isGuessInsertLengthSD( double ave, double sd ){
  return (ave != META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED) and (sd == META_VELVET_G_PARAMETER_VALUE_NOT_SPECIFIED);
}

double MetaVelvetG::guessInsertLengthSD( double ave ){
  return ave * META_VELVET_G_SD_RATIO_FOR_GUESS;
}

void MetaVelvetG::guessCoverageFlags(){
  if( expectedCoverage < 0 )
    flagEstimateExpectedCoverage = true;
  if( coverageCutoff < 0 )
    flagEstimateCoverageCutoff = true;
}


