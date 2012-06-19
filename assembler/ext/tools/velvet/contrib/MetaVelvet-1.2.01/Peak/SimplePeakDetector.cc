#include "SimplePeakDetector.hh"

void SimplePeakDetector::initializeCounts(){
  for( double binVal=minCoverage ; binVal<=maxCoverage ; binVal+=binWidth ){
    setCount( binVal, 0.0 );
  }
}

void SimplePeakDetector::setCounts( MetaGraph* graph ){
  for( int i=1 ; i<graph->getNumNodes() ; ++i ){
    Node* node = graph->getNode(i);
    if( node==NULL || getNodeLength(node) <= 0 ){
      continue;
    }
    double binVal = binWidth * (int)floor( VUtils::getNodeDensity(node) / binWidth );
    if( binVal <= maxCoverage ){
      addCount( binVal, getNodeLength(node) );
    }
  }
}

int SimplePeakDetector::detectCoveragePeaks( MetaGraph* graph, double* coveragePeaks, double* coverageBoundaries ){
  setCounts( graph );
  double firstValley = searchFirstValley();
  double noiseCutoff = getNoiseCutoff();
  vector<double> peaks;
  cout << "[MetaHisto] " << "Noise cutoff coverage = " << noiseCutoff << endl;
  for( double curCoverage=(firstValley+binWidth) ; curCoverage<=(maxCoverage-binWidth) ; curCoverage+=binWidth ){
    double prevCount = getCount(curCoverage-binWidth);
    double curCount  = getCount(curCoverage);
    double nextCount = getCount(curCoverage+binWidth);
    if( curCount>prevCount && curCount>nextCount ){
      if( curCount >= noiseCutoff ){
	peaks.push_back( curCoverage );
	cout << "[MetaHisto] " << "Find " << peaks.size() << "-th coverage peak: " << curCoverage << " (frequency count = " << curCount << ")" << endl;
      }
    }
  }
  int numPeaks = peaks.size();
  for( int i=numPeaks-1 ; i>=0 ; --i ){
    coveragePeaks[i] = peaks.at(numPeaks-i-1);
  }
  return numPeaks;
}

double SimplePeakDetector::getNoiseCutoff() const {
  double firstValley = searchFirstValley();
  double peakBin = firstValley;
  double maxCount = getCount(firstValley);
  for( double curCoverage=(firstValley+binWidth) ; curCoverage<=maxCoverage ; curCoverage+=binWidth ){
    double curCount = getCount(curCoverage);
    if( curCount > maxCount ){
      maxCount = curCount;
      peakBin  = curCoverage;
    }
  }
  cout << "[MetaHisto] " << "First valley = " << firstValley << endl;
  cout << "[MetaHisto] " << "Largest peak coverage = " << peakBin << " (frequency count = " << maxCount << ")" << endl;
  return maxCount / snRatio;
}

double SimplePeakDetector::searchFirstValley() const {
  for( double curCoverage=(minCoverage+binWidth) ; curCoverage<=(maxCoverage-binWidth) ; curCoverage+=binWidth ){
    double prevCount = getCount(curCoverage-binWidth);
    double curCount  = getCount(curCoverage);
    double nextCount = getCount(curCoverage+binWidth);
    if( curCount<prevCount && curCount>nextCount ){
      return curCoverage;
    }
    prevCount = curCount;
  }
  return minCoverage;
}


