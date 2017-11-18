#ifndef ANGULARCORRELATION_H
#define ANGULARCORRELATION_H

#include "JetEvent/Jet.h"
#include "JetEvent/JetConstituentIterator.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "UserAnalysis/Numeric.h" // Common numeric stuff

#include <map> // For ACF exclusively

// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

// http://arxiv.org/pdf/1104.1646
class AngularCorrelation {
public:
  AngularCorrelation(const Jet*& j, double R, double dR = 0.05):
  _jet(j), _R(R), _dR(dR), _sumK(0.), _sumTheta(0.), _muJ2(0.) {
    if(_jet->size() <= 1) {
      _sumK = _sumTheta = _muJ2 = NaN;
    }
    compute();
  }

  const Jet* jet() const { return _jet; } 
  double deltaG() const { return _sumTheta != 0.? _R * Divide(_sumK, _sumTheta): NaN; } // PUBLISH_JET_SHAPE
  double G() const { return Divide(_sumTheta, _muJ2); } // PUBLISH_JET_SHAPE
  double muJ2() const { return _muJ2; } // PUBLISH_JET_SHAPE
  double mStar2() const { return sqrt(CLHEP::pi * _dR * _dR) * Divide(deltaG() * G(), _R) * muJ2(); } // PUBLISH_JET_SHAPE
  double mStar() const { return sqrt(mStar2()); } // PUBLISH_JET_SHAPE
  
  void compute() {
    if(_jet->size() <= 1) return;

    JetConstituentIterator iFirst = JetConstituentIterator::first(_jet);
    JetConstituentIterator iLast = JetConstituentIterator::last(_jet);
    JetConstituentIterator jFirst = JetConstituentIterator::first(_jet);
    JetConstituentIterator jLast = JetConstituentIterator::last(_jet);

    for(int i = 0; iFirst != iLast; ++iFirst, i++) {
      for(int j = 0; jFirst != jLast; ++jFirst, j++) {
	if(iFirst == jFirst) continue;
	double et_i = iFirst.et();
	double et_j = jFirst.et();
	  
	double rapidity_i = iFirst.rapidity();
	double rapidity_j = jFirst.rapidity();
	double phi_i = iFirst.phi();
	double phi_j = jFirst.phi();

	double dRapidity = fabs(rapidity_i - rapidity_j);
	double dPhi = _deltaPhi(phi_i, phi_j);

	double deltaR2 = (dRapidity * dRapidity) + (dPhi * dPhi);
	double deltaR = sqrt(deltaR2);
	double K = Kernel(_R - deltaR, _dR);
	double Theta = HeavisideTheta(_R - deltaR);

	_sumK += et_i * et_j * deltaR2 * K;
	_sumTheta += et_i * et_j * deltaR2 * Theta;
	_muJ2 +=  et_i * et_j * deltaR2;
      } // j loop
    } // i loop
  }

private:
  DeltaPhi _deltaPhi;
  const Jet* _jet;
  double _R;
  double _dR;
  double _sumK, _sumTheta, _muJ2;
};

class AngularStructure {
public:
  AngularStructure(const Jet*& j, double dR = 0.05, double h0 = 4.):
  _jet(j), _Rmax(JetMisc(j).R()), _dR(dR), _h0(h0) {
    _RtoGmap.clear();
    _RtoDeltaGmap.clear();
    _RtoMuJ2map.clear();
    _DeltaGtoRmap.clear();
    _RtoMmap.clear();
    compute();
  }
  
  void compute() {
    this->scan();
    this->findPeaks();
  }

  //private:
  void scan() {
    double G, deltaG, muJ2, mStar;
    for(double r = _dR; r <= _Rmax; r += _dR) {
      AngularCorrelation ac(_jet, r, _dR);
      G = ac.G();
      deltaG = ac.deltaG();
      muJ2 = ac.muJ2();
      mStar = ac.mStar();
      if(!isNaN(G)) _RtoGmap.insert(std::make_pair(r, G));
      if(!isNaN(deltaG)) _RtoDeltaGmap.insert(std::make_pair(r, deltaG));
      if(!isNaN(muJ2)) _RtoMuJ2map.insert(std::make_pair(r, muJ2));
      if(!isNaN(deltaG)) _DeltaGtoRmap.insert(std::make_pair(deltaG, r));
      if(!isNaN(mStar)) _RtoMmap.insert(std::make_pair(r, mStar));
    }
  }

  bool isLocalMaxima(std::multimap<double, double>::reverse_iterator it) {
    Map::iterator iter = _RtoDeltaGmap.find(it->second); // Find iterator for this R
    double deltaGnow = iter->second; // deltaG of candidate peak
    double deltaGprev = NaN;
    double deltaGnext = NaN;
    if(iter != _RtoDeltaGmap.begin()) { // beginning of the map?
      std::advance(iter, -1); // Go back 1
      deltaGprev = iter->second;
      std::advance(iter, 1); // Reset
    }
    std::advance(iter, 1); // Go forward 1
    if(iter != _RtoDeltaGmap.end()) deltaGnext = iter->second; // end of the map?
    std::advance(iter, -1); // Reset
    return !(deltaGnow < deltaGprev) && !(deltaGnow < deltaGnext);    
  }

  double G(double R) const {
    double result = NaN;
    int count = _RtoGmap.count(R);
    if(count > 0) result = _RtoGmap.find(R)->second;
    return result;
  }
  
  int numPeaks() const { // PUBLISH_JET_SHAPE
    return _peaks.size();
  }
  
  double RiStar(unsigned int i) { // PUBLISH_JET_SHAPE
    double result = NaN;
    if((_peaks.size() == 0) || (_peaks.size() > i)) return NaN;
    Map::iterator iter = _peaks.begin();
    std::advance(iter, i-1);
    if(iter != _peaks.end()) result = iter->first;
    return result;
  }
  
  double MiStar(unsigned int i) { // PUBLISH_JET_SHAPE
    double result = NaN;
    if((_peaks.size() == 0) || (_peaks.size() > i)) return NaN;
    Map::iterator iter = _peaks.begin();
    std::advance(iter, i-1);
    if(iter != _peaks.end()) result = _RtoMmap.find(iter->first)->second; 
    return result;
  }

  void findPeaks() {
    if(_RtoDeltaGmap.size() == 0) return; // no candidate peaks

    MultiMap::reverse_iterator iHighest = _DeltaGtoRmap.rbegin();
    MultiMap::reverse_iterator iLowest = _DeltaGtoRmap.rend();
    MultiMap::reverse_iterator jHighest = _DeltaGtoRmap.rbegin();
    MultiMap::reverse_iterator jLowest = _DeltaGtoRmap.rend();

    for(jHighest = _DeltaGtoRmap.rbegin(); jHighest != jLowest; ++jHighest) { // Select a candidate peak to examine
      if(!isLocalMaxima(jHighest)) continue; // this is not a local maxima
      _descents.clear();
      for(iHighest = _DeltaGtoRmap.rbegin(); iHighest != jHighest; ++iHighest) { // Select a larger candidate peak
	if(iHighest == jHighest) continue; // Don't let the two candidate peaks be the same
	if(!isLocalMaxima(iHighest)) continue; // this is not a local maxima
	double iR = iHighest->second; // R position of (larger) candidate peak
	double jR = jHighest->second; // R position of (smaller) candidate peak
	Map::iterator iterR1 = _RtoDeltaGmap.find(std::min(iR, jR)); // Sort R1 < R2
	Map::iterator iterR2 = _RtoDeltaGmap.find(std::max(iR, jR)); // Sort R1 < R2

	Map::iterator min; // Iterator pointing to a minimum in the range R1 to R2
	double minDeltaG; // Its deltaG

	min = std::min_element(iterR1, ++iterR2, sorter); // Find smallest deltaG between two points
	minDeltaG = min->second;

	double descent = fabs(minDeltaG - jHighest->first);
	// Make sure that minimum is not at R1 or R2
	if(descent > 0.) _descents.insert(std::make_pair(descent, jHighest->second)); // Keep a record of vertical decents
      } // end inner loop over higher peaks
      MultiMap::iterator minDescentIter = _descents.begin(); // Minimum at beginning of sorted map
      double minDescent = NaN;
      double RminDescent = NaN;
      if(!_descents.empty()) { // Get values of minimum vertical descent, provided such exist
	minDescent = minDescentIter->first;
	RminDescent = minDescentIter->second;
      }
      // If prominence > h0 or local maxima is highest peak
      if(minDescent > _h0 || jHighest == _DeltaGtoRmap.rbegin()) _peaks.insert(std::make_pair(jHighest->second, jHighest->first));
    } // end outer loop of examined peaks
    // We've found all possible peaks
  }

  typedef std::map<double, double> Map;
  struct SorterOp {
    bool operator()(Map::value_type &i1, Map::value_type &i2) {
      return i1.second < i2.second; }
  } sorter;

  typedef std::multimap<double, double> MultiMap;

  const Jet* _jet;
  double _Rmax, _dR, _h0;
  Map _RtoGmap;
  Map _RtoDeltaGmap;
  Map _RtoMuJ2map;
  MultiMap _DeltaGtoRmap;
  Map _RtoMmap;
  MultiMap _descents;
  Map _peaks;
};

#endif
