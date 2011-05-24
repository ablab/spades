#include "cute.h"
#include "strobe_read.hpp"
#include "graphBuilder.hpp"

using namespace abruijn;
using namespace std;

void TestSpellingGenome1() {
    string const filename1 = "./data/input/MG1655-K12_emul1.fastq.gz";
    string const filename2 = "./data/input/MG1655-K12_emul1.fastq.gz";
    size_t const reads_cut = -1;
    int const mode = 1;
    size_t const take = 1;

    string const ref_genome_filename = "./data/input/MG1655-K12.fasta.gz";
    size_t const ref_genome_cut = -1;

	string file_names[2] = { filename1, filename2 };
	StrobeReader<2, Read, ireadstream> sr(file_names);
	PairedReader<ireadstream> paired_stream(sr, 220);
	SimpleReaderWrapper<PairedReader<ireadstream> > srw(paired_stream);
	CuttingReader<SimpleReaderWrapper<PairedReader<ireadstream> > > cr(srw, reads_cut);

	abruijn::GraphBuilderMaster<CuttingReader<SimpleReaderWrapper<PairedReader<ireadstream>>>> gbm(cr, take, mode);
	gbm.build();

	bool ok = gbm.SpellGenomeThroughGraph ( ref_genome_filename, ref_genome_cut );

	ASSERT ( ok );

}

cute::suite SpellingGenomeSuite(){
	cute::suite s;
	s.push_back(CUTE(TestSpellingGenome1));
	return s;
}
