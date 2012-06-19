#include "VelvetUtils.hh"

void VelvetUtils::importSequences( const string& seqFileName, int argc, char* argv[], boolean& flagDoubleStrand, boolean& flagNoHash ){
  parseDataAndReadFiles( VUtils::getCharFileName(seqFileName), argc - 2, &(argv[2]), &flagDoubleStrand, &flagNoHash );
}

void VelvetUtils::setRoadmaps( const string& seqFileName, const string& roadmapFileName, int hashLength, boolean flagDoubleStrand, boolean noHash ){
  resetWordFilter( hashLength );
  SplayTable* splayTable = newSplayTable( hashLength, (boolean)true );
  ReadSet* allSequences = importReadSet( VUtils::getCharFileName(seqFileName) );
  inputSequenceArrayIntoSplayTableAndArchive( allSequences, splayTable, VUtils::getCharFileName(roadmapFileName), VUtils::getCharFileName(seqFileName) );
  destroySplayTable( splayTable );
}

void VelvetUtils::constructPregraph( const string& seqFileName, const string& roadmapFileName, const string& pregraphFileName ){
  RoadMapArray* roadmaps = importRoadMapArray( VUtils::getCharFileName(roadmapFileName) );
  PreGraph* preGraph = newPreGraph_pg( roadmaps, VUtils::getCharFileName(seqFileName) );
  concatenatePreGraph_pg( preGraph );
  clipTips_pg( preGraph );
  exportPreGraph_pg( VUtils::getCharFileName(pregraphFileName), preGraph );
  destroyPreGraph_pg( preGraph );
}


