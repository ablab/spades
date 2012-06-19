#include "meta-velveth.hh"

int main( int argc, char* argv[] ){
  MetaVelvetH* mvh = new MetaVelvetH( argc, argv );
  cout << "[meta-velveth] " << "Import sequences..." << endl;
  mvh->importSequences( argc, argv );
  cout << "[meta-velveth] " << "...done (import sequences)." << endl << endl;

  cout << "[meta-velveth] " << "Set roadmaps..." << endl;
  mvh->setRoadmaps();
  cout << "[meta-velveth] " << "...done (set roadmaps)." << endl << endl;
  
  cout << "[meta-velveth] " << "All tasks were safely completed..." << endl << endl;
  delete mvh;
  return 0;
}

MetaVelvetH::MetaVelvetH( int argc, char* argv[] ){
  setDefaultParameters();
  if( !checkParameters( argc, argv ) ){
    printUsage();
    exit(-1);
  }
  setParameters( argc, argv );
}

void MetaVelvetH::setDefaultParameters(){
  flagDoubleStrand = META_VELVET_H_FLAG_DOUBLE_STRAND;
  flagNoHash = META_VELVET_H_FLAG_NO_HASH;
}

void MetaVelvetH::importSequences( int argc, char* argv[] ){
  VelvetUtils::importSequences( FileNameUtils::getSeqFileName(prefix), argc, argv, flagDoubleStrand, flagNoHash );
}

void MetaVelvetH::setRoadmaps(){
  cout << "[meta-velveth] " << "  flag double strand ... " << (int)flagDoubleStrand << endl;
  VelvetUtils::setRoadmaps( FileNameUtils::getSeqFileName(prefix), FileNameUtils::getRoadmapFileName(prefix), hashLength, flagDoubleStrand, flagNoHash );
}

bool MetaVelvetH::checkParameters( int argc, char* argv[] ) const {
  cout << endl<< endl
       << "[meta-velveth] " << "Check command line options..." << endl;
  if( checkArgc(argc) and checkHashLength(argc, argv) ){
    cout << "\t" << "OK. Your command line options seem to be good." << endl
	 << endl << endl;
    return true;
  }
  cout << "\t" << "Your command line options are not appropriate." << endl
       << "\t" << "Please read usage of meta-velveth carefully." << endl
       << endl
       << "\t" << "Thanks, " << endl
       << endl << endl;
  return false;
}

bool MetaVelvetH::checkArgc( int argc ) const {
  if (argc < 4){
    return false;
  }
  return true;
}

bool MetaVelvetH::checkHashLength( int argc, char* argv[] ) const {
  return checkMultipleHashLength( argc, argv )
    and checkMaxHashLength( argc, argv )
    and checkMinHashLength( argc, argv )
    and checkEvenHashLength( argc, argv );
}

bool MetaVelvetH::checkMultipleHashLength( int argc, char* argv[] ) const {
  if ( strstr(argv[2],"," ) ){
    cerr << "ERROR: meta-velveth cannot use multiple k-mers" << endl;
    return false;
  }
  return true;
}

bool MetaVelvetH::checkMaxHashLength( int argc, char* argv[] ) const {
  int hashLength = atoi(argv[2]);
  if( hashLength >= MAXKMERLENGTH ){
    cerr << "ERROR: meta-velveth can't handle k-mers as long as " << MAXKMERLENGTH << endl
	 << "\t" << "please re-compile with setting 'MAXKMERLENGTH' parameter" << endl;
    return false;
  }
  return true;
}

bool MetaVelvetH::checkMinHashLength( int argc, char* argv[] ) const {
  int hashLength = atoi(argv[2]);
  if( hashLength <= 0 ){
    cerr << "ERROR: Invalid hash length: " << hashLength << endl;
    return false;
  }
  return true;
}

bool MetaVelvetH::checkEvenHashLength( int argc, char* argv[] ) const {
  int hashLength = atoi(argv[2]);
  if( hashLength % 2 == 0 ){
    cerr << "ERROR: meta-velveth can't use even length k-mers, such as " << hashLength << endl;
    return false;
  }
  return true;
}

void MetaVelvetH::printUsage() const {
  cerr << endl << endl
       << "************************************************************" << endl
       << "Program: ----------------------------------------------" << endl
       << "  " << "meta-velveth - simple hashing program" << endl
       << "\t" << "Version = " << META_VELVET_VERSION << endl
       << "\t" << "MAX_CATEGORIES = " << CATEGORIES << endl
       << "\t" << "MAX_KMER_LENGTH = " << MAXKMERLENGTH << endl
       << endl << endl
       << "Usage: ------------------------------------------------" << endl
       << "  " << "meta-velveth directory hash_length {[-file_format][-read_type] filename1 [filename2 ...]} {...} [options]" << endl
       << "\t" << "directory\t: directory name for output files" << endl
       << "\t" << "hash_length\t: an odd integer <= " << MAXKMERLENGTH << endl
       << "\t" << "filename\t: path to sequence file or - for standard input" << endl
       << endl
       << "  " << "File format options:" << endl
       << "\t" << "-fasta  -fastq  -raw  -fasta.gz  -fastq.gz  -raw.gz  -sam  -bam" << endl
       << endl
       << "  " << "Read type options:" << endl
       << "\t" << "-short\t-shortPaired" << endl;
  for( int cat=2 ; cat <= (int)CATEGORIES ; ++cat ){
    cerr << "\t" << "-short" << cat << "\t-shortPaired" << cat << endl;
  }
  cerr << "\t" << "-long\t-longPaired" << endl
       << "\t" << "-reference" << endl
       << endl
       << "  " << "Other options:" << endl
       << "\t" << "-strand_specific\t: for strand specific transcriptome sequencing data (default: off)" << endl
       << "\t" << "-reuse_Sequences\t: reuse Sequences file (or link) already in directory (no need to provide original filenames in this case (default: off)" << endl
       << "\t" << "-noHash\t\t\t: simply prepare Sequences file, do not hash reads or prepare Roadmaps file (default: off)" << endl
       << endl << endl
       << "Command line examples: --------------------------------" << endl
       << "  " << "Ex.1 - Short single end reads:" << endl
       << "\t" << "meta-velveth TestDir 29 -short -fastq s_1_sequence.txt" << endl
       << endl
       << "  " << "Ex.2 - Paired-end short reads (remember to interleave paired reads):" << endl
       << "\t" << "meta-velveth TestDir 31 -shortPaired -fasta interleaved.fna" << endl
       << endl
       << "  " << "Ex.3 - Two channels and some long reads:" << endl
       << "\t" << "meta-velveth TestDir 43 -short -fastq unmapped.fna -longPaired -fasta SangerReads.fasta" << endl
       << endl
       << "  " << "Ex.4 - Three channels:" << endl
       << "\t" << "meta-velveth TestDir 35 -shortPaired -fasta pe_lib1.fasta -shortPaired2 pe_lib2.fasta -short3 se_lib1.fa" << endl
       << endl << endl
       << "Output: -----------------------------------------------" << endl
       << "\t" << "directory/Sequences" << endl
       << "\t" << "directory/Roadmaps" << endl
       << "\t" << "[Both files are picked up by graph, so please leave them there]" << endl
       << endl 
       << "************************************************************" << endl
       << endl << endl;
}

void MetaVelvetH::setParameters( int argc, char* argv[] ){
  prefix = argv[1];
  hashLength = atoi(argv[2]);
  string mkdirCmd = "mkdir -p " + prefix;
  system( mkdirCmd.c_str() );
}


