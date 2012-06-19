#include "PeakDetectorParameters.hh"

bool PeakDetectorParameters::checkParameters() const {
  if( maxCoverage < 0 ){
    cout << "[PeakDetectorParameters] " << "Error: Invalid max. peak coverage: " << maxCoverage << endl
	 << "[PeakDetectorParameters] " << "max. peak coverage should be non-negative." << endl;
    return false;
  }
  if( minCoverage < 0 ){
    cout << "[PeakDetectorParameters] " << "Error: Invalid min. peak coverage: " << minCoverage << endl
	 << "[PeakDetectorParameters] " << "min. peak coverage should be non-negative." << endl;
    return false;
  }
  if( maxCoverage <= minCoverage ){
    cout << "[PeakDetectorParameters] " << "Error: Invalid min. & max. peak coverages (max: " << maxCoverage << ", min: " << minCoverage << ")" << endl
	 << "[PeakDetectorParameters] " << "max. peak coverage should be greater than min. peak coverage." << endl;
    return false;
  }
  if( binWidth <= 0 ){
    cout << "[PeakDetectorParameters] " << "Error: Invalid histogram bin width: " << binWidth << endl
	 << "[PeakDetectorParameters] " << "bin width should be greater than 0." << endl;
    return false;
  }
  if( snRatio <= 1 ){
    cout << "[PeakDetectorParameters] " << "Error: Invalid histogram signal-noise ratio: " << snRatio << endl
	 << "[PeakDetectorParameters] " << "sn ratio should be greater than 1." << endl;
    return false;
  }
  return true;
}


