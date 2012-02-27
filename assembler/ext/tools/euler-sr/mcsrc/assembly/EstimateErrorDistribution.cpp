/***************************************************************************
 * Title:          EstimateErrorDistribution.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "utils.h"
#include "Spectrum.h"
#include "SimpleStats.h"
#include <cmath>
#include <math.h>
#include <map>
#include "IntegralTupleStatic.h"

#include "compatibility.h"


double ExpPdf(double x,double lambda);
double Factorial(ssize_t f);
double GaussPdf(double mu, double simasq, double x);
double PoissonPdf(ssize_t k, double lambda);
ssize_t   FindDistribIntersection(double lambda, double mu, double sigmasq);

void PrintUsage() {
		std::cout << "usage: estErrorDist spectrumFile [-numSamples] " << std::endl
							<< " [-ig g]           intial guess for separating exponential from gaussian" << std::endl
							<< " [-delta d]        convergence parameter d " << std::endl
				      << " [-printAll]       Print all estimated parameters: exponential error, and mean" << std::endl
							<< "                   and var. of gaussian correct." << std::endl;
		std::cout << " [-verbose]        Print lots of output" << std::endl;
		std::cout << " [-printCoverage]  Print the estimated coverage." << std::endl;
		std::cout << " [-printHistogram] Print a histogram of word counts." << std::endl;
		std::cout << " [-binary]         Read spectrum as binary" << std::endl;
		std::cout << " [-maxIts n]        Maximum number of iterations" << std::endl;
}
int main(int argc, char* argv[]) {

	std::string spectrumFileName;
	if (argc < 2 ) {
		PrintUsage();
		exit(1);
	}

	spectrumFileName = argv[1];

	ssize_t numSamples = 100000;

	int argi = 2;
	ssize_t expMax = 2;
	double delta = 0.00001;
	ssize_t verbose = 0;
  ssize_t printAll = 0;
	ssize_t printHistogram = 0;
	ssize_t printCoverage = 0;
	ssize_t printMaxError = 1;
	ssize_t maxIts = 100;
	ssize_t binarySpectrum = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-numSamples")== 0) {
			numSamples = atosz(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxIts")== 0) {
			maxIts = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-ig")== 0) {
			expMax = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-delta") == 0) {
			delta = atof(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-verbose") == 0) {
			verbose = 1;
		}
    else if (strcmp(argv[argi], "-printAll") == 0) {
      printAll = 1;
    }
		else if (strcmp(argv[argi], "-printHistogram") == 0) {
			printHistogram = 1;
		}
		else if (strcmp(argv[argi], "-printCoverage") == 0) {
			printCoverage = 1;
			printMaxError = 0;
		}
		else if (strcmp(argv[argi], "-binary") == 0) {
			binarySpectrum = 1;
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
		++argi;
	}


	std::string reportFileName = FormReportName(spectrumFileName);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);


	std::ifstream spectrumIn;
	if (binarySpectrum == 0) 
		openck(spectrumFileName, spectrumIn, std::ios::in, report);
	else 
		openck(spectrumFileName, spectrumIn, std::ios::in| std::ios::binary, report);
	

	ssize_t numTuples;
	ssize_t tupleSize_SSZT;
	if (binarySpectrum == 0) {
		spectrumIn >> tupleSize_SSZT;
		spectrumIn >> numTuples;
	}
	else {
		spectrumIn.read((char*) &tupleSize_SSZT, sizeof(ssize_t) );
		spectrumIn.read((char*) &numTuples, sizeof(ssize_t) );
	}
	CountedIntegralTuple::SetTupleSize(tupleSize_SSZT);

	std::string tuple;
	//UNUSED// ssize_t mult;
	std::vector<ssize_t> spectMult;
	ssize_t i;
	ssize_t numGenerated = 0;
	std::map<ssize_t, ssize_t> histogram;

	// Simply read the entire spectrum.
	if (numTuples < numSamples) {
		//		spectrumIn.close();
		spectMult.resize(numTuples);
		if (binarySpectrum == 0) {
			i = 0;
			while(spectrumIn) {
				spectrumIn >> tuple >> spectMult[i];
				if (histogram.find(spectMult[i]) == histogram.end()) {
					histogram[spectMult[i]] = 1;
				}
				i++;
			}
		}
		else {
			CountedIntegralTuple countedTuple;
			for (i = 0; i < numTuples; i++) { 
				spectrumIn.read((char*) &countedTuple, sizeof(CountedIntegralTuple));
				spectMult[i] = countedTuple.count;
				if (histogram.find(spectMult[i]) == histogram.end()) {
					histogram[spectMult[i]] = 1;
				}
			}
		}

		numGenerated = numTuples;
	}
	else {
		//		SeedRandom();
		std::vector<ssize_t> sample;
		sample.resize(numSamples);
		//UNUSED// ssize_t minIndex, maxIndex;
		double pSample;

		// Generate *roughly* numSamples indices.
		pSample = double(numSamples) / numTuples;
		for (i = 0; i < numTuples and numGenerated < numSamples; i++) {
			if (Uniform() < pSample) { 
				sample[numGenerated] = i;
				numGenerated++;
			}
		}

		ssize_t curTuple = 0;
		std::string line;
		//UNUSED// ssize_t mult;
		spectMult.resize(numGenerated);
		if (binarySpectrum == 0) {
			for (i = 0; i < numTuples and curTuple < numGenerated; i++) {
				while (sample[curTuple] > i) {
					std::getline(spectrumIn, line);
					i++;
				}
				spectrumIn >> tuple >> spectMult[curTuple];
				if (histogram.find(spectMult[curTuple]) == histogram.end()) {
					histogram[spectMult[curTuple]] = 1;
				}
				else {
					histogram[spectMult[curTuple]]++;
				}
				curTuple++;
			}
		}
		else {
			CountedIntegralTuple countedTuple;
			for (i = 0; i < numTuples and curTuple < numGenerated; i++) {
				while (sample[curTuple] > i) {
					spectrumIn.read((char*) &countedTuple, sizeof(CountedIntegralTuple));
					i++;
				}
				//spectrumIn >> tuple >> spectMult[curTuple];
				spectrumIn.read((char*) &countedTuple, sizeof(CountedIntegralTuple));
				spectMult[curTuple] = countedTuple.count;
				if (histogram.find(countedTuple.count) == histogram.end()) {
					histogram[countedTuple.count] = 1;
				}
				else {
					histogram[countedTuple.count]++;
				}
				curTuple++;
			}
		}
	}

	std::map<ssize_t,ssize_t>::iterator mapit;
	//UNUSED// ssize_t begin = 99999999;
	ssize_t best  = -1;
	ssize_t bestIndex = -1;
	//UNUSED// ssize_t second = -1;
	ssize_t prev = -1;
	ssize_t hump = 0;
	if (histogram.begin() != histogram.end()) {
		best = histogram.begin()->second;
		bestIndex = 1;
		prev = best;
	}
	if (printHistogram) {
		for (mapit = histogram.begin(); mapit != histogram.end(); ++mapit) {
			std::cout << mapit->first << " " << mapit->second << " " << best << std::endl;
			if (!hump and prev < mapit->second) {
				hump = 1;
			}
			if (! hump and best > mapit->second) {
				best = mapit->second;
				bestIndex = mapit->first;
			}
			prev = mapit->second;
		}
		std::cout << "best: " << bestIndex << std::endl;
	}
	//	exit(0);
	
	
	// Initial guesses.  The exponential process does not generate 
	// anything above 3, count the frequencies and means of 
	// exponential and gaussian process.

	
	ssize_t numExp, numGauss;
	ssize_t sumExp = 0; 
	ssize_t sumGauss = 0;
	double sumSq = 0;
	numExp = 0; numGauss = 0;
	for (i = 0; i < numGenerated; i++) {
		if (spectMult[i] <= expMax) {
			numExp++;
			sumExp += spectMult[i];
		}
		else {
			numGauss++;
			sumGauss += spectMult[i];
			sumSq += (double) spectMult[i] * (double) spectMult[i];
		}
	}
	//UNUSED// ssize_t window[3];
	std::map<ssize_t,ssize_t>::iterator histEnd, histIt;
	histEnd = histogram.end();
	//UNUSED// ssize_t h = 0;
	if (verbose) {
		for (histIt = histogram.begin(); histIt != histEnd; histIt++) {
			std::cout << (*histIt).second << " ";
		}
		std::cout << std::endl;
	}

	if (numGauss == 0) {
		std::cout << 1 << std::endl;
		exit(0);
	}
	if (numExp == 0) {
		std::cout << 1 << std::endl;
		exit(0);
	}


	double expMean = double(sumExp) /  numExp;
	double gaussMean = double(sumGauss) / numGauss;
	double expLambda;
	double gaussVar, gaussStddev;
	//	std::cout << "expMean: " << expMean << std::endl;
	expLambda = 1/expMean;

	gaussVar =  sumSq / numGauss - gaussMean * gaussMean ;
	gaussStddev = sqrt(gaussVar);
	double mixGauss = double(numGauss) / numGenerated;
	double mixExp   = 1.0 - mixGauss;

	ssize_t curIndex = 0;
	
	for (i = 0; i < numGenerated; i++) {
		if (spectMult[i] < (4*gaussStddev + gaussMean)) {
			spectMult[curIndex] = spectMult[i];
			curIndex++;
		}
	}

	numGenerated = curIndex;
	
	std::vector<double> pExp;
	pExp.resize(numGenerated);

	// Initialize the probablity of being 
	// generated by an error process.
	//UNUSED// ssize_t numExcluded = 0;
	for (i = 0; i < numGenerated; i++) {
		if (spectMult[i] <= expMax) 
			pExp[i] = 1.0;
		else
			pExp[i] = 0.0;
	}

	if (verbose) {
		std::cout << "initialized to: ";
		std::cout << "0 " << expLambda << " " << gaussMean << " " << gaussVar 
							<< " " << mixExp << std::endl;
	}
	ssize_t converged = 0;

	double expLambdaHat, gaussVarHat, gaussMeanHat;

	double expResp;
	ssize_t iter = 1;
	double ex, ga;
	
	while (! converged ) {
		//
		// Compute expectations.
		//
		for (i = 0; i < numGenerated; i++) { 
			//			ex = PoissonPdf(spectMult[i], expLambda);
			ex = ExpPdf(spectMult[i],expLambda);
			ga = GaussPdf(gaussMean, gaussVar, spectMult[i]);
			if (pExp[i] >= 0) {
				if (((ga == 0) or 
						 (ga + ex == 0) or 
						 ((mixExp *ex + (1-mixExp) *ga) == 0)) and
						spectMult[i] > gaussMean) {
					// prevent some overflow conditions.
					pExp[i] = 0;
				}
				else {
					/*
					pExp[i] = (mixExp) * PoissonPdf(spectMult[i], expLambda)/ 
						((mixExp) * PoissonPdf(spectMult[i], expLambda) + 
						 (1 - mixExp) * GaussPdf(gaussMean, gaussVar, spectMult[i]));
					*/
					pExp[i] = (mixExp) * ExpPdf(spectMult[i], expLambda)/ 
						((mixExp) * ExpPdf(spectMult[i], expLambda) + 
						 (1 - mixExp) * GaussPdf(gaussMean, gaussVar, spectMult[i]));
					
					if (isNan(pExp[i])) {
						std::cout << "Unable to compute expectation for " << spectMult[i] <<  " " 
											<< mixExp << " " << expLambda<< std::endl;
						std::cout << ExpPdf(spectMult[i], expLambda) << " " 
											<<  GaussPdf(gaussMean, gaussVar, spectMult[i]) << " " << ex + ga << std::endl;
						std::cout << "Make sure you are running on the entire spectrum, including" << std::endl
											<< "frequency-1 tuples, since the exponential distribution uses these."
											<< std::endl;
						exit(1);
					}
				}
			}
			/*			std::cout << spectMult[i] << " " 
								<< ex << " " << ga << " " << " " << mixExp 
								<< " " <<  pExp[i] << std::endl;
			*/
		}

		//
		// Compute maximization.
		//
		expResp = 0;
		for (i = 0; i < numGenerated; i++) {
			if (pExp[i] >= 0) 
				expResp += pExp[i];
		}

		// exponential parameter

		expMean = 0;
		for (i = 0; i < numGenerated; i++) {
			expMean += pExp[i]*spectMult[i];
		}
		expLambdaHat = expResp / expMean;
		
		// gaussian parameter mu

		double sumGauss = 0;
		for (i = 0; i < numGenerated; i++) {
			sumGauss += (1-pExp[i]) * spectMult[i];
		}
		gaussMeanHat = sumGauss / (numGenerated - expResp);

		// gaussian parameter sigma squared

		double gaussSumSq = 0;
		for (i = 0; i < numGenerated; i++) {
			gaussSumSq += (1-pExp[i]) * (spectMult[i] - gaussMean)* (spectMult[i]-gaussMean);
		}
		gaussVarHat = gaussSumSq / (numGenerated - expResp);

		mixExp = expResp / numGenerated;
		if (verbose) {
			std::cout << iter << " " 
								<< expLambdaHat << " "
								<< gaussMeanHat << " "
								<< gaussVarHat << " "
								<< mixExp << std::endl;
		}
		if (iter == maxIts) 
			converged = 1;
		if (fabs(expLambda - expLambdaHat) < delta or
				fabs(gaussMean - gaussMeanHat) < delta or
				fabs(gaussVar - gaussVarHat) < delta) {
			converged = 1;
		}
		expLambda = expLambdaHat;
		gaussMean = gaussMeanHat;
		gaussVar  = gaussVarHat;
		++iter;
	}


	ssize_t threshIsect;
	threshIsect = FindDistribIntersection(expLambda, gaussMean, gaussVar);
	if (printMaxError) {
   	std::cout << threshIsect << std::endl;
  }
	else if (printCoverage) {
		std::cout << gaussMean << std::endl;
	}
	if (printAll) {
		std::cout << expLambda << " " << gaussMean << " " << gaussVar << " " << mixExp << std::endl;
	}
//	std::ofstream pExpOut;
//	pExpOut.open("pExpOut.txt");
//	for(i = 0;i < numGenerated; i++) {
//		pExpOut << pExp[i] << " " << spectMult[i] << std::endl;
//	}
//	pExpOut.close();

	EndReport(report);
	report.close();
	return 0;
}


double Factorial(ssize_t f) {
	if (f<2) {
		return 1.0;
	}

	double res = f;
	f--;
	while (f > 1) {
		res *= f;
		f--;
	}
	return res;
}

double PoissonPdf(ssize_t k, double lambda) {
	if (k < 11)
		return exp(k*log(lambda))*exp(-lambda)/Factorial(k);
	else
		return 0.0;
}



double ExpPdf(double x,double lambda) {
	return lambda * exp(-x*lambda);
}

double GaussPdf(double mu, double sigmasq, double x) {
	//UNUSED// double expval = exp(-((x - mu)*(x-mu))/(2*sigmasq));
	//	assert(expval != 0.0);
	return 1/sqrt(sigmasq*2*M_PI) * exp(-((x - mu)*(x-mu))/(2*sigmasq));
}

ssize_t FindDistribIntersection(double lambda, double mu, double sigmasq) {

	double a, b, c;

	a = 1;
	b = -2*mu-2*lambda*sigmasq;
	c = 2*sigmasq*log(lambda) + sigmasq*log(2*M_PI*sigmasq) + mu*mu;

	if (b*b - 4*a*c  < 0) {
		return -1;
	}

	else 
		return (ssize_t) floor( (-b - sqrt(b*b-4*a*c))/2*a);
}

