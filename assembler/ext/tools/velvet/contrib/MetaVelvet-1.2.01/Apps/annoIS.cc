#include "annoIS.hh"

int main( int argc, char* argv[] ){
  AnnoIS* annois = new AnnoIS( argc, argv );

  cout << "[annoIS] " << "Load IS node info..." << endl;
  vector<ISNode*> isNodes = annois->loadISNodes();
  cout << "[annoIS] " << isNodes.size() << " IS nodes were found."  << endl;
  cout << "[annoIS] " << "...done (load IS node info)." << endl << endl;

  cout << "[annoIS] " << "Load last graph of velvet..." << endl;
  ISGraph* graph = annois->loadGraph();
  cout << "[annoIS] " << "...done (load last graph of velvet)." << endl << endl;

  cout << "[annoIS] " << "Annotate IS positions in last graph..." << endl;
  annois->annotateIS( graph, isNodes );
  cout << "[annoIS] " << "...done (annotate IS positions in last graph)." << endl << endl;

  cout << "[annoIS] " << "Output adjacent reliable node IDs..." << endl
       << "[annoIS] " << "\t" << "filename = " << annois->getARNodeFileName() << endl;
  annois->saveARNodes( isNodes );
  cout << "[annoIS] " << "...done (output adjacent reliable node IDs)." << endl << endl;

  delete annois;
  return 0;
}

vector<ISNode*> AnnoIS::loadISNodes() const {
  return ISNodeIO::load(nodeFileName);
}

ISGraph* AnnoIS::loadGraph() const {
  return new ISGraph( FileNameUtils::getSeqFileName(prefix), FileNameUtils::getVelvetLastGraphFileName(prefix) );
}

void AnnoIS::annotateIS( ISGraph* graph, const vector<ISNode*>& isNodes ) const {
  graph->annotateIS( isNodes, reliableNodeLength, flankingLength );
}

void AnnoIS::saveARNodes( const vector<ISNode*>& isNodes ) const {
  ISNodeIO::saveARNodes( this->getARNodeFileName(), isNodes );
}

AnnoIS::AnnoIS( int argc, char* argv[] ){
  setDefaultParameters();
  cout << endl << endl
       << "[annoIS] " << "Check command line options..." << endl;
  if( !checkParameters( argc, argv ) or !setParameters( argc, argv ) or !checkParameters() ){
    cout << "\t" << "Your command line options are not appropriate." << endl
	 << "\t" << "Please read usage of annoIS carefully." << endl
	 << endl
	 << "\t" << "Thanks, " << endl
	 << endl << endl;
    printUsage();
    exit(-1);
  }
  cout << "\t" << "OK. Your command line options seem to be good." << endl
       << endl << endl;
}

void AnnoIS::setDefaultParameters(){
  reliableNodeLength = ANNOIS_DEFAULT_RELIABLE_NODE_LENGTH;
  flankingLength     = ANNOIS_DEFAULT_FLANKING_LENGTH;
  flagUseMV          = false;
}

bool AnnoIS::checkParameters() const {
  return true;
}

bool AnnoIS::checkParameters( int argc, char* argv[] ) const {
  if( checkArgc(argc) and checkPrefix(argc, argv) ){
    return true;
  }
  return false;
}

bool AnnoIS::checkArgc( int argc ) const {
  if( argc < 2 ){
    return false;
  }
  return true;
}

bool AnnoIS::checkPrefix( int argc, char* argv[] ) const {
  if( strncmp( argv[1], "-", 1 ) == 0 ){
    cerr << "[annoIS] Error: the first argument should be a directory path (should not be start with '-')" << endl;
    return false;
  }
  return checkMandatoryFiles( (string)argv[1] );
}

bool AnnoIS::checkMandatoryFiles( string prefix ) const {
  if( !FileUtils::isValidFile(FileNameUtils::getSeqFileName(prefix)) ){
    cerr << "[annoIS] Error: Lack a mandatory file: " << FileNameUtils::getSeqFileName(prefix) << endl
	 << "\t" << "You should run 'velveth' before running 'annoIS'" << endl << endl;
    return false;
  }
  if( !FileUtils::isValidFile(this->getGraphFileName( prefix )) ){
    cerr << "[annoIS] Error: Lack a mandatory file: " << this->getGraphFileName( prefix ) << endl
	 << "\t" << "You should run 'velvetg' or 'meta-velvetg' before runing 'annoIS'" << endl << endl;
    return false;
  }
  return true;
}

bool AnnoIS::checkCategory( int cat ) const {
  if( (cat<0) or (cat>CATEGORIES) ){
    cerr << "[annoIS] Error: Invalid category: " << cat << endl;
    return false;
  }
  return true;
}

void AnnoIS::printUsage() const {
  cerr << endl << endl
       << "************************************************************" << endl
       << "Program: ----------------------------------------------" << endl
       << "  " << "annoIS - annotating genomic positions harboring insertion sequences IS" << endl
       << "\t" << "Version = " << ANNOIS_VERSION << endl
       << "\t" << "MAX_CATEGORIES = " << CATEGORIES << endl
       << "\t" << "MAX_KMER_LENGTH = " << MAXKMERLENGTH << endl
       << endl << endl
       << "Usage: ------------------------------------------------" << endl
       << "  " << "annois directory nodefile [options]" << endl
       << "\t" << "directory                   " << "\t: directory name for output files" << endl
       << "\t" << "nodefile                    " << "\t: node info file" << endl
       << endl
       << "  " << "annoIS options:" << endl
       << "\t" << "-reliable_length <int>      " << "\t: minimum reliable node length (default: " << ANNOIS_DEFAULT_RELIABLE_NODE_LENGTH << ")" << endl
       << "\t" << "-flanking_length <int>      " << "\t: length fo flanking sequences (default: " << ANNOIS_DEFAULT_FLANKING_LENGTH << ")" << endl
       << "\t" << "-meta <yes|no>              " << "\t: Use meta-velvet results (default: no = velvet results)" << endl
       << endl << endl
       << "Output: -----------------------------------------------" << endl
       << "\t" << "directory/*.AR-nodes.txt    " << "\t: List of Adjacent Reliable (AR) node IDs" << endl
       << "************************************************************" << endl
       << endl << endl;
}

bool AnnoIS::setParameters( int argc, char* argv[] ){
  prefix = argv[1];
  nodeFileName = argv[2];
  arNodeFileName = argv[3];
  for ( int arg_index=4 ; arg_index<argc ; arg_index++ ) {
    char* arg = argv[arg_index++];
    if ( arg_index >= argc ) {
      cerr << "[annoIS] Error: Unusual number of arguments!" << endl;
      return false;
    }
    if( !setParameter( arg,  argv[arg_index] ) ){
      return false;
    }
  }
  return true;
}

bool AnnoIS::setParameter( char* name, char* val ){
  if( strcmp(name, "-reliable_length")      == 0 ){ return setReliableNodeLength(val); }
  if( strcmp(name, "-flanking_length")      == 0 ){ return setFlankingLength(val); }
  if( strcmp(name, "-meta")                 == 0 ){ return setFlagUseMV(val); }
  if( strcmp(name, "--help")                == 0 ){ return false; }
  cerr << "[annoIS] Unknown option: " << name << endl;
  return false;
}

bool AnnoIS::setReliableNodeLength( char* val ){
  reliableNodeLength = atoi( val );
  return true;
}

bool AnnoIS::setFlankingLength( char* val ){
  flankingLength = atoi( val );
  return true;
}

bool AnnoIS::setFlagUseMV( char* val ){
  flagUseMV = (strcmp(val, "yes") == 0);
  return true;
}

string AnnoIS::getGraphFileName( const string& prefix ) const {
  if( flagUseMV ){
    return FileNameUtils::getLastGraphFileName(prefix);
  }
  return FileNameUtils::getVelvetLastGraphFileName(prefix);
}


