#ifndef READSREFORMATTER_HPP_
#define READSREFORMATTER_HPP_

#include "ireadstream.hpp"
#include "read.hpp"

/*
 * @params
 * file1, file2: files in fastq.gz
 * outputFile: resulting file, with reads separated by \n, paired reads goes one after another
 *
 */
int forgetQualityPairedData(string file1, string file2, string outputFile) {
	ireadstream stream1(file1, false);
	ireadstream stream2(file2, true);

	FILE* outFile = fopen(outputFile.c_str(), "w");
	Read r1;
	Read r2;
	while (!(stream1.eof()) && !(stream2.eof())) {
		stream1 >> r1;
		stream2 >> r2;
		if (r1.isValid() && r2.isValid()) {
			Sequence s1 = r1.getSequence();
			Sequence s2 = r2.getSequence();

			fprintf(outFile, "%s %s\n", s1.str().c_str(), s2.str().c_str());
		}
	}
	fclose(outFile);
}
#endif// READSREFORMATTER_HPP_
