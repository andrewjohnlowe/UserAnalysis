#ifndef EVENTSHAPES_H
#define EVENTSHAPES_H

#include "JetEvent/Jet.h"
#include "JetEvent/JetCollection.h"
#include "CLHEP/Vector/ThreeVector.h" // Common to most variables
#include "CLHEP/Vector/TwoVector.h" // For DeltaS
#include "CLHEP/Vector/LorentzVector.h" // For event shapes only
#include "CLHEP/Units/PhysicalConstants.h" // Used for angularities (tau)

#include "UserAnalysis/Numeric.h" // Common numeric stuff
#include "UserAnalysis/JetFlags.h" // Jet reco flags

#include <vector>
#include <algorithm> // for std::sort

// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

using namespace Numeric;

class EventShape { // XX = EJ, Jj
public:
  EventShape(const JetCollection* theJetCollection, const Jet* parentJet = 0):
  _Njets(numJets(theJetCollection)),
  _dijets(_Njets >= 2) { DEBUG(0);
    if(theJetCollection) { DEBUG(0);
      JetCollectionHelper::jetcollection_t theJets(theJetCollection->begin(), theJetCollection->end()); DEBUG(0);
      JetCollectionHelper::jetcollection_t::iterator firstJet = theJets.begin(); DEBUG(0);
      JetCollectionHelper::jetcollection_t::iterator lastJet  = theJets.end(); DEBUG(0);
      _jets.clear(); DEBUG(0);
      for(unsigned int i = 0; firstJet != lastJet; ++firstJet, i++) { DEBUG(0);
	if((*firstJet) == 0) continue; DEBUG(0);
	const Jet* aJet = (*firstJet)->clone(true, true); DEBUG(0); // Implicit new here
	if(aJet == 0) continue; DEBUG(0);
	_jets.push_back(new CLHEP::HepLorentzVector(boost(aJet->hlv(), parentJet))); DEBUG(0);
	delete aJet; // We're done with this jet; do everything onwards with HLVs
      }
      std::sort(_jets.begin(), _jets.end(), sortByEtDown); DEBUG(0); // Now we can be super-sure the jets are sorted after boosting
      _firstJet = _Njets >= 1? _jets[0]: 0; DEBUG(0);
      _secondJet = _Njets >= 2? _jets[1]: 0; DEBUG(0); 
      _thirdJet = _Njets >= 3? _jets[2]: 0; DEBUG(0);
      _fourthJet = _Njets >= 4? _jets[3]: 0; DEBUG(0); 
      
      RadialParameter R;
      double radialParameter = 0.1 * R(theJetCollection);
      _deltaR = radialParameter + 0.2; // Something to try later maybe; ensures jets don't overlap much
      //QuartetSortByMprime f;
      
      _jetQuartet = std::make_pair(std::make_pair(_firstJet, _secondJet), std::make_pair(_thirdJet, _fourthJet)); DEBUG(0);
      _jetQuartet2 = std::make_pair(std::make_pair(_firstJet, _secondJet), std::make_pair(_thirdJet, _fourthJet)); DEBUG(0);
      _jetQuartet3 = std::make_pair(std::make_pair(_firstJet, _secondJet), std::make_pair(_thirdJet, _fourthJet)); DEBUG(0);
      _jetQuartet4 = std::make_pair(std::make_pair(_firstJet, _secondJet), std::make_pair(_thirdJet, _fourthJet)); DEBUG(0);
      _jetQuartet5 = std::make_pair(std::make_pair(_firstJet, _secondJet), std::make_pair(_thirdJet, _fourthJet)); DEBUG(0);
      _jetQuartet6 = std::make_pair(std::make_pair(_firstJet, _secondJet), std::make_pair(_thirdJet, _fourthJet)); DEBUG(0);
      
      _4jetCombinations.clear(); DEBUG(0);
      _4jetCombinations2.clear(); DEBUG(0);
      _4jetCombinations3.clear(); DEBUG(0);
       
      if(_Njets >= 4) { DEBUG(0);
	JetPair j1j2 = std::make_pair(_firstJet, _secondJet); DEBUG(0);
	JetPair j3j4 = std::make_pair(_thirdJet, _fourthJet); DEBUG(0);
	JetPair j1j3 = std::make_pair(_firstJet, _thirdJet); DEBUG(0);
	JetPair j2j4 = std::make_pair(_secondJet, _fourthJet); DEBUG(0);
	JetPair j1j4 = std::make_pair(_firstJet, _fourthJet); DEBUG(0);
	JetPair j2j3 = std::make_pair(_secondJet, _thirdJet); DEBUG(0);
	_4jetCombinations.push_back(std::make_pair(rapidityOrder(j1j2), rapidityOrder(j3j4))); DEBUG(0);
	_4jetCombinations.push_back(std::make_pair(rapidityOrder(j1j3), rapidityOrder(j2j4))); DEBUG(0);
	_4jetCombinations.push_back(std::make_pair(rapidityOrder(j1j4), rapidityOrder(j2j3))); DEBUG(0);
	  
	//if(f.dR(j1j2) > deltaR && f.dR(j3j4) > deltaR)
	_4jetCombinations2.push_back(std::make_pair(rapidityOrder(j1j2), rapidityOrder(j3j4))); DEBUG(0);
	//if(f.dR(j1j3) > deltaR && f.dR(j2j4) > deltaR)
	_4jetCombinations2.push_back(std::make_pair(rapidityOrder(j1j3), rapidityOrder(j2j4))); DEBUG(0);
	//if(f.dR(j1j4) > deltaR && f.dR(j2j3) > deltaR)
	_4jetCombinations2.push_back(std::make_pair(rapidityOrder(j1j4), rapidityOrder(j2j3))); DEBUG(0);

	_jetQuartet = *std::min_element(_4jetCombinations.begin(), _4jetCombinations.end(), dRsort); DEBUG(0);
	_jetQuartet2 = *std::min_element(_4jetCombinations.begin(), _4jetCombinations.end(), qSort); DEBUG(0);
	_jetQuartet3 = *std::min_element(_4jetCombinations.begin(), _4jetCombinations.end(), scdfSort); DEBUG(0);
	_jetQuartet4 = *std::min_element(_4jetCombinations.begin(), _4jetCombinations.end(), sd0Sort); DEBUG(0);
	if(!_4jetCombinations2.empty()) {
	  _jetQuartet5 = *std::min_element(_4jetCombinations2.begin(), _4jetCombinations2.end(), mprimeSort); DEBUG(0);
	}
      }
      if(_Njets >= 2) { DEBUG(_Njets);
	_2jetCombinations.clear(); DEBUG(0);
	for(unsigned int i = 0; i != _jets.size(); i++) {
	  for(unsigned int j = 0; j != _jets.size(); j++) {
	    if(i >= j) continue; DEBUG(i); DEBUG(j);
	    CLHEP::HepLorentzVector* first = _jets[i]; DEBUG(0);
	    CLHEP::HepLorentzVector* second = _jets[j]; DEBUG(0);
	    _2jetCombinations.push_back(std::make_pair(first, second)); DEBUG(first); DEBUG(second); DEBUG(_2jetCombinations.size());

	    for(unsigned int k = 0; k != _jets.size(); k++) {
	      if(j >= k) continue; DEBUG(i); DEBUG(j);
	      for(unsigned int l = 0; l != _jets.size(); l++) {
		if(k >= l) continue; DEBUG(i); DEBUG(j);
		JetPair JiJj = std::make_pair(_jets[i], _jets[j]);
		JetPair JkJl = std::make_pair(_jets[k], _jets[l]);
		//if(f.dR(JiJj) > deltaR && f.dR(JkJl) > deltaR) 
		_4jetCombinations3.push_back(std::make_pair(rapidityOrder(JiJj), rapidityOrder(JkJl))); DEBUG(0);
	      } // l loop
	    } // k loop
	  } // j loop
	} // i loop
	_jetPair = *std::min_element(_2jetCombinations.begin(), _2jetCombinations.end(), deltaMsort); DEBUG(_jetPair.first);
	if(!_4jetCombinations3.empty()) {
	  _jetQuartet6 = *std::min_element(_4jetCombinations3.begin(), _4jetCombinations3.end(), mprimeSort); DEBUG(0);
	}
      } // if njets >= 2
    } // if theJetCollection
  } // CTOR end
  
  ~EventShape() {
    for(unsigned int i = 0; i != _jets.size(); i++) { DEBUG(0);
      delete _jets[i]; DEBUG(0);
      _jets[i] = 0; DEBUG(0);
    }
    _jets.clear(); DEBUG(0);
  }
  
  int/* jStr_*/numJets(const JetCollection* theJetCollection) const { DEBUG(0);
    return theJetCollection? theJetCollection->size(): 0;
  }

  bool isDijets() const { return _dijets; } 

  double/* jStr_evtXX_*/etaJ1() const {
    double result = NaN;
    if(_Njets >= 1) {
      result = fabs(_firstJet->eta());
    }
    return result;
  }

  double/* jStr_evtXX_*/etaJ2() const {
    double result = NaN;
    if(_Njets >= 2) {
      result = fabs(_secondJet->eta());
    }
    return result;
  }

  double/* jStr_evtXX_*/etaJ3() const {
    double result = NaN;
    if(_Njets >= 3) {
      result = fabs(_thirdJet->eta());
    }
    return result;
  }

  double/* jStr_evtXX_*/etaJ4() const {
    double result = NaN;
    if(_Njets >= 4) {
      result = fabs(_fourthJet->eta());
    }
    return result;
  }

  double/* jStr_evtXX_*/pTJ1() const {
    double result = NaN;
    if(_Njets >= 1) {
      result = fabs(_firstJet->perp());
    }
    return result;
  }

  double/* jStr_evtXX_*/pTJ2() const {
    double result = NaN;
    if(_Njets >= 2) {
      result = fabs(_secondJet->perp());
    }
    return result;
  }

  double/* jStr_evtXX_*/pTJ3() const {
    double result = NaN;
    if(_Njets >= 3) {
      result = fabs(_thirdJet->perp());
    }
    return result;
  }

  double/* jStr_evtXX_*/pTJ4() const {
    double result = NaN;
    if(_Njets >= 4) {
      result = fabs(_fourthJet->perp());
    }
    return result;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=10&resId=0&materialId=slides&confId=150186
  double/* jStr_evtXX_*/fE1() const {
    double result = NaN;
    if(_Njets >= 1) {
      double e1 = _firstJet->e();
      result = Divide(e1, jetSumE());
    }
    return result;
  }
  
  double/* jStr_evtXX_*/fE2() const {
    double result = NaN;
    if(_Njets >= 2) {
      double e2 = _secondJet->e();
      result = Divide(e2, jetSumE());
    }
    return result;
  }

  double/* jStr_evtXX_*/fE3() const {
    double result = NaN;
    if(_Njets >= 3) {
      double e3 = _thirdJet->e();
      result = Divide(e3, jetSumE());
    }
    return result;
  }

  double/* jStr_evtXX_*/fET1() const {
    double result = NaN;
    if(_Njets >= 1) {
      double et1 = _firstJet->et();
      result = Divide(et1, jetSumET());
    }
    return result;
  }
  
  double/* jStr_evtXX_*/fET2() const {
    double result = NaN;
    if(_Njets >= 2) {
      double et2 = _secondJet->et();
      result = Divide(et2, jetSumET());
    }
    return result;
  }

  double/* jStr_evtXX_*/fET3() const {
    double result = NaN;
    if(_Njets >= 3) {
      double et3 = _thirdJet->et();
      result = Divide(et3, jetSumET());
    }
    return result;
  }

  // Dijet contransverse mass
  double/* jStr_evtXX_*/mct() const {
    double result = NaN;
    if(_dijets) {
      double et1 = _firstJet->et();
      double et2 = _secondJet->et();
      double pt1 = _firstJet->perp();
      double pt2 = _secondJet->perp();
      double et = (et1 + et2);
      double pt = (pt1 - pt2);
      result = Sqrt((et * et) - (pt * pt));
    }
    return result;
  }

  // http://arxiv.org/pdf/hep-ph/0403297v1
  double/* jStr_evtXX_*/y1y2() const {
    double result = NaN;
    if(_dijets) {
      result = _firstJet->rapidity() * _secondJet->rapidity();
    }
    return result;
  }

  // http://www.actaphys.uj.edu.pl/vol36/pdf/v36p0321.pdf
  // https://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=176611
  double/* jStr_evtXX_*/q() const {
    double result = NaN;
    if(_dijets) {
      double d = (*_firstJet - *_secondJet).mag2();
      result = Sqrt(-d);
    }
    return result;
  }
  
  // Dijet invariant mass
  double/* jStr_evtXX_*/mjj() const {
    double result = NaN;
    if(_dijets) {
      result = (*_firstJet + *_secondJet).m();
    }
    return result;
  }

  double/* jStr_evtXX_*/mTjj() const {
    double result = NaN;
    if(_dijets) {
      result = (*_firstJet + *_secondJet).mt();
    }
    return result;
  }

  // Trijet invariant mass
  // https://twindico.hep.anl.gov/indico/getFile.py/access?contribId=3&sessionId=0&resId=1&materialId=slides&confId=650
  double/* jStr_evtXX_*/mjjj() const {
    double result = NaN;
    if(_Njets >= 3) {
      result = (*_firstJet + *_secondJet + *_thirdJet).m();
    }
    return result;
  }

  double/* jStr_evtXX_*/mTjjj() const {
    double result = NaN;
    if(_Njets >= 3) {
      result = (*_firstJet + *_secondJet + *_thirdJet).mt();
    }
    return result;
  }

  double/* jStr_evtXX_*/mjjjj() const {
    double result = NaN;
    if(_Njets >= 4) {
      result = (*_firstJet + *_secondJet + *_thirdJet + *_fourthJet).m();
    }
    return result;
  }

  double/* jStr_evtXX_*/mTjjjj() const {
    double result = NaN;
    if(_Njets >= 4) {
      result = (*_firstJet + *_secondJet + *_thirdJet + *_fourthJet).mt();
    }
    return result;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=164974
  double/* jStr_evtXX_*/zetaPlus() const {
    double result = NaN;
    if(_Njets >= 3) {
      double eta1 = _firstJet->eta();
      double eta2 = _secondJet->eta();
      double eta3 = _thirdJet->eta();
      result = fabs(eta3) + fabs(eta1 - eta2);
    }
    return result;
  }

  // http://arxiv.org/abs/1104.1175
  double/* jStr_evtXX_*/zetaMinus() const {
    double result = NaN;
    if(_Njets >= 3) {
      double eta1 = _firstJet->eta();
      double eta2 = _secondJet->eta();
      double eta3 = _thirdJet->eta();
      result = fabs(eta3) - fabs(eta1 - eta2);
    }
    return result;
  }

  // http://dare.uva.nl/document/208893
  double/* jStr_evtXX_*/B12() const {
    double result = NaN;
    if(_dijets) {
      double m = (*_firstJet + *_secondJet).m();
      result = Divide(m, jetSumM());
    }
    return result;
  }

  // http://dare.uva.nl/document/208893
  double/* jStr_evtXX_*/BT12() const {
    double result = NaN;
    if(_dijets) {
      double mt = (*_firstJet + *_secondJet).mt();
      result = Divide(mt, jetSumMT());
    }
    return result;
  }

  // dijet searches for supersymmetry at the large hadron collider
  // Randall, Tucker-Smith, PRL 101, 221803 (2008)
  double/* jStr_evtXX_*/Alpha() const {
    double result = NaN;
    if(_dijets) {
      double et2 = _secondJet->et();
      result = Divide(et2, mjj());
    }
    return result;
  }

  double/* jStr_evtXX_*/AlphaT() const {
    double result = NaN;
    if(_dijets) {
      double et2 = _secondJet->et();
      result = Divide(et2, mTjj());
    }
    return result;
  }

  // arXiv:1006.1650v1 [hep-ph] (Unburied Higgs)
  double/* jStr_evtXX_*/massDemocracy() const {
    double result = NaN;
    if(_dijets) {
      double m1 = _firstJet->m();
      double m2 = _secondJet->m();
      result = std::min(Divide(m1, m2), Divide(m2, m1));
    }
    return result;
  }

  // arXiv:1006.1650v1 [hep-ph] (Unburied Higgs) // changed pt --> et
  double/* jStr_evtXX_*/betaflow(double etmin = 0.) const {
    double result = NaN;
    if(_Njets >= 3) {
      double et1 = _firstJet->et();
      double et2 = _secondJet->et();
      double et3 = _thirdJet->et();
      double Q = Divide(et3, (et1 + et2));
      result = etmin == 0.? Q: // the default
	et3 < etmin? 0.: Q; // suppress very soft radiation; suitable etmin: 1 GeV, 5 GeV (paper) 
    }
    return result;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=6&resId=0&materialId=slides&confId=160270
  double/* jStr_evtXX_*/y23(double etmin = 0.) const {
    double b = betaflow(etmin);
    return b * b;
  }

  double/* jStr_evtXX_*/lny23(double etmin = 0.) const {
    double b = betaflow(etmin);
    return Ln(b * b);
  }

  // Angle between jets
  double/* jStr_evtXX_*/theta() const {
    double result = NaN;
    if(_dijets) {
      double mjj2 = mjj() * mjj();
      result = Arccos(1. - Divide(mjj2, 2. * _firstJet->e() * _firstJet->e()));
    }
    return result;
  }

  // Asymmetry
  // arXiv:1011.6182v2 [hep-ex]
  double/* jStr_evtXX_*/asym() const {
    double result = NaN;
    if(_dijets) {
      double et1 = _firstJet->et();
      double et2 = _secondJet->et();
      result = Divide(et1 - et2, et1 + et2);
    }
    return result;
  }

  // http://cdsweb.cern.ch/record/1326634
  double/* jStr_evtXX_*/yB() const {
    double result = NaN;
    if(_dijets) {
      double y1 = _firstJet->rapidity();
      double y2 = _secondJet->rapidity();
      result = 0.5*(y1 + y2);
    }
    return result;
  }

  // http://cdsweb.cern.ch/record/1326634
  double/* jStr_evtXX_*/yStar() const {
    double result = NaN;
    if(_dijets) {
      double y1 = _firstJet->rapidity();
      double y2 = _secondJet->rapidity();
      result = 0.5*(y1 - y2);
    }
    return result;
  }

  // http://cdsweb.cern.ch/record/1326634
  double/* jStr_evtXX_*/thetaStar() const {
    return _dijets? Arccos(tanh(yStar())): NaN;
  }

  // http://cdsweb.cern.ch/record/1326634
  double/* jStr_evtXX_*/chi() const {
    double result = NaN;
    if(_dijets) {
      double y1 = _firstJet->rapidity();
      double y2 = _secondJet->rapidity();
      result = exp(fabs(y1 - y2));
    }
    return result;
  }
 
  // Angular deltas
  // http://arxiv.org/pdf/hep-ph/0105325v2
  double/* jStr_evtXX_*/deltaPhiJJ() const {
    double result = NaN;
    if(_dijets) {
      double phi1 = _firstJet->phi();
      double phi2 = _secondJet->phi();
      DeltaPhi deltaPhi;
      result = deltaPhi(phi1, phi2);
    }
    return result;
  }

  double/* jStr_evtXX_*/deltaThetaJJ() const {
    double result = NaN;
    if(_dijets) {
      DeltaPhi deltaPhi;
      double theta1 = _firstJet->theta();
      double theta2 = _secondJet->theta();
      result = deltaPhi(theta1, theta2);
    }
    return result;
  }

  double/* jStr_evtXX_*/deltaEtaJJ() const {
    double result = NaN;
    if(_dijets) {
      double eta1 = _firstJet->eta();
      double eta2 = _secondJet->eta();
      result = eta1 - eta2;
    }
    return result;
  }

  double/* jStr_evtXX_*/deltaRapidityJJ() const {
    double result = NaN;
    if(_dijets) {
      double rapidity1 = _firstJet->rapidity();
      double rapidity2 = _secondJet->rapidity();
      result = rapidity1 - rapidity2;
    }
    return result;
  }

  // Angular delta R
  double/* jStr_evtXX_*/deltaRJJ() const {
    double result = NaN;
    if(_dijets) {
      double dEta = deltaEtaJJ();
      double dPhi = deltaPhiJJ();
      result = Sqrt((dEta * dEta) + (dPhi * dPhi));
    }
    return result;
  }

  double/* jStr_evtXX_*/deltaRJJY() const {
    double result = NaN;
    if(_dijets) {
      double dRapidity = deltaRapidityJJ();
      double dPhi = deltaPhiJJ();
      result = Sqrt((dRapidity * dRapidity) + (dPhi * dPhi));
    }
    return result;
  }

  // Angular sigmas
  double/* jStr_evtXX_*/sigmaPhiJJ() const {
    double result = NaN;
    if(_dijets) {
      PhiCorr phiCorr;
      double phi1 = phiCorr(_firstJet->phi());
      double phi2 = phiCorr(_secondJet->phi());
      result = phiCorr(phi1 + phi2);
    }
    return result;
  }

  double/* jStr_evtXX_*/sigmaThetaJJ() const {
    double result = NaN;
    if(_dijets) {
      PhiCorr phiCorr;
      double theta1 = phiCorr(_firstJet->theta());
      double theta2 = phiCorr(_secondJet->theta());
      result = phiCorr(theta1 + theta2);
    }
    return result;
  }

  double/* jStr_evtXX_*/sigmaEtaJJ() const {
    double result = NaN;
    if(_dijets) {
      double eta1 = _firstJet->eta();
      double eta2 = _secondJet->eta();
      result = eta1 + eta2;
    }
    return result;
  }

  double/* jStr_evtXX_*/sigmaRapidityJJ() const {
    double result = NaN;
    if(_dijets) {
      double rapidity1 = _firstJet->rapidity();
      double rapidity2 = _secondJet->rapidity();
      result = rapidity1 + rapidity2;
    }
    return result;
  }

  double/* jStr_evtXX_*/sigmaPtJJ() const {
    double result = NaN;
    if(_dijets) {
      double pt1 = _firstJet->perp();
      double pt2 = _secondJet->perp();
      result = pt1 + pt2;
    }
    return result;
  }

  double/* jStr_evtXX_*/sigmaEtJJ() const {
    double result = NaN;
    if(_dijets) {
      double et1 = _firstJet->et();
      double et2 = _secondJet->et();
      result = et1 + et2;
    }
    return result;
  }

  double/* jStr_evtXX_*/sigmaEt12() const {
    double result = NaN;
    if(_dijets) {
      double et1 = _firstJet->et();
      double et2 = _secondJet->et();
      result = et1 + et2;
    }
    return result;
  }

  double/* jStr_evtXX_*/sigmaEt34() const {
    double result = NaN;
    if(_Njets >= 4) {
      double et3 = _thirdJet->et();
      double et4 = _fourthJet->et();
      result = et3 + et4;
    }
    return result;
  }

  // http://dare.uva.nl/document/208893
  // PHYSICAL REVIEW D 82, 032002 (2010)
  // Measurement of the tt cross section using high-multiplicity jet events
  double/* jStr_evtXX_*/A234() const {
    double result = NaN;
    if(_Njets >= 4) {
      double et2 = _secondJet->et();
      double et3 = _thirdJet->et();
      double et4 = _fourthJet->et();
      result = Divide((et2 + et3) - et4, et2 + et3 + et4);
    }
    return result;
  }

  // Angular asymmetries
  double/* jStr_evtXX_*/asymPhiJJ() const {
    return _dijets? Divide(deltaPhiJJ(), sigmaPhiJJ()): NaN;
  }

  double/* jStr_evtXX_*/asymThetaJJ() const {
    return _dijets? Divide(deltaThetaJJ(), sigmaThetaJJ()): NaN;
  }

  double/* jStr_evtXX_*/asymEtaJJ() const {
    return _dijets? Divide(deltaEtaJJ(), sigmaEtaJJ()): NaN;
  }

  double/* jStr_evtXX_*/asymRapidityJJ() const {
    return _dijets? Divide(deltaRapidityJJ(), sigmaRapidityJJ()): NaN;
  }

  // arXiv:1010.3698v1 [hep-ph]
  double/* jStr_evtXX_*/acoplanarity() const {
    double result = NaN;
    if(_dijets) {
      PhiCorr phiCorr;
      DeltaPhi deltaPhi;
      double t1 = deltaPhi(CLHEP::pi, fabs(deltaPhiJJ()));
      double t2 = deltaPhi(CLHEP::pi, sigmaThetaJJ());
      result = phiCorr(fabs(t1) + fabs(t2));
    }
    return result;
  }

  double/* jStr_evtXX_*/twist() const {
    return _dijets? Arctan2(fabs(deltaPhiJJ()), fabs(deltaEtaJJ())): NaN;
  }

  double/* jStr_evtXX_*/twistY() const {
    return _dijets? Arctan2(fabs(deltaPhiJJ()), fabs(deltaRapidityJJ())): NaN;
  }

  // Sum of E over all jets
  double/* jStr_evtXX_*/jetSumE() const { DEBUG(0); // HE
    if(_jets.empty()) return NaN; DEBUG(0); 
    return std::for_each(_jets.begin(), _jets.end(), AddE()).sum;
  }
  
  // Sum of ET over all jets
  double/* jStr_evtXX_*/jetSumET() const { DEBUG(0); // HET
    if(_jets.empty()) return NaN; DEBUG(0);
    return std::for_each(_jets.begin(), _jets.end(), AddET()).sum;
  }

  // Sum of PT over all jets
  double/* jStr_evtXX_*/jetSumPT() const { DEBUG(0); // HPT
    if(_jets.empty()) return NaN; DEBUG(0);
    return std::for_each(_jets.begin(), _jets.end(), AddPT()).sum;
  }

  // Sum of M over all jets
  double/* jStr_evtXX_*/jetSumM() const { DEBUG(0); // HM
    if(_jets.empty()) return NaN; DEBUG(0);
    return std::for_each(_jets.begin(), _jets.end(), AddM()).sum;
  }

  // Sum of PX over all jets
  double jetSumPX() const { DEBUG(0); // HPX
    if(_jets.empty()) return NaN; DEBUG(0);
    return std::for_each(_jets.begin(), _jets.end(), AddPX()).sum;
  }

  // Sum of PY over all jets
  double jetSumPY() const { DEBUG(0); // HPY
    if(_jets.empty()) return NaN; DEBUG(0);
    return std::for_each(_jets.begin(), _jets.end(), AddPY()).sum;
  }

  double/* jStr_evtXX_*/jetSumMT() const { DEBUG(0);
    double sumET = jetSumET(); DEBUG(0);
    double sumPx = jetSumPX(); DEBUG(0);
    double sumPy = jetSumPY(); DEBUG(0);
    double sumET2 = sumET * sumET; DEBUG(0);
    double sumPx2 = sumPx * sumPx; DEBUG(0);
    double sumPy2 = sumPy * sumPy; DEBUG(0);
    double Mt = Sqrt(sumET2 - (sumPx2 + sumPy2)); DEBUG(0);
    return Mt;
  }

  // http://dare.uva.nl/document/208893
  double/* jStr_evtXX_*/HTprime() const { DEBUG(0);
    double result = NaN;
    if(_dijets) { DEBUG(0);
      result = jetSumET() - sigmaEt12(); DEBUG(0);
    }
    return result;
  }

  double/* jStr_evtXX_*/centrality() const { DEBUG(0);
    return Divide(jetSumET(), jetSumE());
  }

  double/* jStr_evtXX_*/centralityP() const { DEBUG(0);
    return Divide(jetSumPT(), jetSumE());
  }
  
  // Splittings
  double/* jStr_evtXX_*/zminJ1J2() const {
    double result = NaN;
    if(_dijets) {
      double e1 = _firstJet->e();
      double e2 = _secondJet->e();
      result = Divide(std::min(e1, e2), e1 + e2);
    }
    return result;
  }
  
  double/* jStr_evtXX_*/zmaxJ1J2() const {
    double result = NaN;
    if(_dijets) {
      double e1 = _firstJet->e();
      double e2 = _secondJet->e();
      result = Divide(std::max(e1, e2), e1 + e2);
    }
    return result;
  }

  // Run over ALL jets 
  double/* jStr_evtXX_*/zminAllJets() const { DEBUG(0);
    if(_jets.empty()) return NaN; DEBUG(0);
    double jetSum = jetSumE(); DEBUG(0);
    return Divide((*std::min_element(_jets.begin(), _jets.end(), jetEnergySort))->e(), jetSum);
  }
  
  // Run over ALL jets 
  double/* jStr_evtXX_*/zmaxAllJets() const { DEBUG(0);
    if(_jets.empty()) return NaN; DEBUG(0);
    double jetSum = jetSumE(); DEBUG(0);
    return Divide((*std::max_element(_jets.begin(), _jets.end(), jetEnergySort))->e(), jetSum);
  }

  // Angular splittings
  double/* jStr_evtXX_*/zminJ1J2Phi() const {
    double result = NaN;
    if(_dijets) {
      PhiCorr phiCorr;
      double phi1 = phiCorr(_firstJet->phi());
      double phi2 = phiCorr(_secondJet->phi());
      result = Divide(std::min(phi1, phi2), phiCorr(phi1 + phi2));
    }
    return result;
  }

  double/* jStr_evtXX_*/zminJ1J2Theta() const {
    double result = NaN;
    if(_dijets) {
      PhiCorr phiCorr;
      double theta1 = phiCorr(_firstJet->theta());
      double theta2 = phiCorr(_secondJet->theta());
      result = Divide(std::min(theta1, theta2), phiCorr(theta1 + theta2));
    }
    return result;
  }

  double/* jStr_evtXX_*/zminJ1J2Eta() const {
    double result = NaN;
    if(_dijets) {
      double eta1 = _firstJet->eta();
      double eta2 = _secondJet->eta();
      result = Divide(std::min(eta1, eta2), eta1 + eta2);
    }
    return result;
  }

  double/* jStr_evtXX_*/zminJ1J2Rapidity() const {
    double result = NaN;
    if(_dijets) {
      double rapidity1 = _firstJet->rapidity();
      double rapidity2 = _secondJet->rapidity();
      result = Divide(std::min(rapidity1, rapidity2), rapidity1 + rapidity2);
    }
    return result;
  }

  // arXiv:1010.3698v1 [hep-ph]
  double/* jStr_evtXX_*/cosHelicityJ1() {
    return _dijets? cosHelicity(_firstJet): NaN;
  }

  double/* jStr_evtXX_*/helicityJ1() {
    return _dijets? helicity(_firstJet): NaN;
  }

  double/* jStr_evtXX_*/azilicityJ1() {
    return _dijets? azilicity(_firstJet): NaN;
  }

  double/* jStr_evtXX_*/cosHelicityJ2() {
    return _dijets? cosHelicity(_secondJet): NaN;
  }

  double/* jStr_evtXX_*/helicityJ2() {
    return _dijets? helicity(_secondJet): NaN;
  }

  double/* jStr_evtXX_*/azilicityJ2() {
    return _dijets? azilicity(_secondJet): NaN;
  }

  double/* jStr_evtXX_*/cosThetaJ1() {
    return _dijets? cosTheta(_firstJet): NaN;
  }
 
  double/* jStr_evtXX_*/cosThetaJ2() {
    return _dijets? cosTheta(_secondJet): NaN;
  } 

  double/* jStr_evtXX_*/deltaRapidityXtoJ1CM() {
    return _dijets? deltaRapidityXtoJCM(_firstJet): NaN;
  }
 
  double/* jStr_evtXX_*/deltaRapidityXtoJ2CM() {
    return _dijets? deltaRapidityXtoJCM(_secondJet): NaN;
  }

  double/* jStr_evtXX_*/deltaRapidityXtoJ1() {
    return _dijets? deltaRapidityXtoJ(_firstJet): NaN;
  }
 
  double/* jStr_evtXX_*/deltaRapidityXtoJ2() {
    return _dijets? deltaRapidityXtoJ(_secondJet): NaN;
  }

  // http://arxiv.org/pdf/1110.2693v1
  double/* jStr_evtXX_*/DeltaM() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 2) { DEBUG(_jetPair.first); DEBUG(_jetPair.second);
      double m1 = _jetPair.first->m(); DEBUG(m1);
      double m2 = _jetPair.second->m(); DEBUG(m2);
      result = fabs(m1 - m2); DEBUG(result);
    }
    return result;
  }

  // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SusyXXbarFourJetSearch
  // https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=162608
  // http://arxiv.org/pdf/1110.2693v1
  double/* jStr_evtXX_*/M() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 2) { DEBUG(_jetPair.first); DEBUG(_jetPair.second);
      double m1 = _jetPair.first->m(); DEBUG(m1);
      double m2 = _jetPair.second->m(); DEBUG(m2);
      result = 0.5 * (m1 + m2); DEBUG(result);
    }
    return result;
  }

  // http://arxiv.org/pdf/1110.2693v1
  double/* jStr_evtXX_*/asymM() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 2) { DEBUG(0);
      result = Divide(0.5 * DeltaM(), M()); DEBUG(DeltaM()); DEBUG(M()); DEBUG(result);
    }
    return result;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=168793
  double/* jStr_evtXX_*/Delta12() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 2) { DEBUG(0);
      QuartetSortBySCDF q;
      result = q.pTvSum(_firstJet, _secondJet); DEBUG(result);
    }
    return result;
  }
  // https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=168793
  double/* jStr_evtXX_*/Delta12n() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 2) { DEBUG(0);
      QuartetSortBySCDF q;
      result = Divide(q.pTvSum(_firstJet, _secondJet), q.pTsum(_firstJet, _secondJet)); DEBUG(result);
    }
    return result;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=168793
  // http://prd.aps.org/pdf/PRD/v47/i11/p4857_1 
  // Phys. Rev. D 47, 4857, 1993
  double/* jStr_evtXX_*/DeltaSCDF() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySCDF q;
      result = fabs(q.DeltaS(_jetQuartet3)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/dPhiIJCDF() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySCDF q;
      result = fabs(q.dPhiIJ(_jetQuartet3)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/dPhiKLCDF() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySCDF q;
      result = fabs(q.dPhiKL(_jetQuartet3)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/PtIPtJCDF() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySCDF q;
      result = fabs(q.PtIPtJ(_jetQuartet3)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/PtKPtLCDF() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySCDF q;
      result = fabs(q.PtKPtL(_jetQuartet3)); DEBUG(result);
    }
    return result;
  }

  // BASED ON:
  // https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=168793
  // http://prd.aps.org/pdf/PRD/v47/i11/p4857_1 
  // Phys. Rev. D 47, 4857, 1993
  double/* jStr_evtXX_*/DeltaSD0() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySD0 q;
      result = fabs(q.DeltaS(_jetQuartet4)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/dPhiIJD0() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySD0 q;
      result = fabs(q.dPhiIJ(_jetQuartet4)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/dPhiKLD0() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySD0 q;
      result = fabs(q.dPhiKL(_jetQuartet4)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/PtIPtJD0() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySD0 q;
      result = fabs(q.PtIPtJ(_jetQuartet4)); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/PtKPtLD0() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySD0 q;
      result = fabs(q.PtKPtL(_jetQuartet4)); DEBUG(result);
    }
    return result;
  }

  // http://arxiv.org/pdf/1112.0003v1
  double/* jStr_evtXX_*/Q() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortByMass q;
      result = q.Q(_jetQuartet2); DEBUG(result);
    }
    return result;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=168793
  // https://indico.cern.ch/conferenceDisplay.py?confId=168793
  double/* jStr_evtXX_*/SCDF() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySCDF q;
      result = q.SCDF(_jetQuartet3); DEBUG(result);
    }
    return result;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=slides&confId=168793
  // https://indico.cern.ch/conferenceDisplay.py?confId=168793
  double/* jStr_evtXX_*/SD0() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      QuartetSortBySD0 q;
      result = q.SD0(_jetQuartet4); DEBUG(result);
    }
    return result;
  }

  // http://cdsweb.cern.ch/record/1416058/files/EXO-11-016-pas.pdf
  // https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=162608
  // https://indico.cern.ch/getFile.py/access?resId=0&materialId=slides&confId=174936
  // http://arxiv.org/pdf/1110.2693v1.pdf
  // http://www.springerlink.com/content/76478x2716nr466l/fulltext.pdf JHEP 09 (2011)074
  double/* jStr_evtXX_*/Mprime4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.Mprime(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/dMprime4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.dM(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/MprimeAvg4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.Mavg(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaMin4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaMin(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaMax4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaMax(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaPhiXX4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaPhiXX(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaYXX4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaYXX(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaRXX4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaRXX(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/TwistXX4() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.TwistXX(_jetQuartet5); DEBUG(result);
    }
    return result;
  }

  int/* jStr_evtXX_*/separated4() { return separated4(_deltaR); }
  int separated4(double dR) { DEBUG(0);
    int result = -1;
    if(_Njets >= 4 && !_4jetCombinations2.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.separated(_jetQuartet5, dR); DEBUG(result);
    }
    return result;
  }

  // http://arxiv.org/pdf/1110.2693v1.pdf
  double/* jStr_evtXX_*/Mprime() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.Mprime(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/dMprime() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.dM(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/MprimeAvg() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.Mavg(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaMin() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaMin(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaMax() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaMax(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaPhiXX() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaPhiXX(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaYXX() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaYXX(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/DeltaRXX() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.DeltaRXX(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  double/* jStr_evtXX_*/TwistXX() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.TwistXX(_jetQuartet6); DEBUG(result);
    }
    return result;
  }

  int/* jStr_evtXX_*/separated() { return separated(_deltaR); }
  int separated(double dR) { DEBUG(0);
    int result = -1;
    if(_Njets >= 4 && !_4jetCombinations3.empty()) { DEBUG(0);
      QuartetSortByMprime q;
      result = q.separated(_jetQuartet6, dR); DEBUG(result);
    }
    return result;
  }

  // http://arxiv.org/pdf/1112.0003v1
  double/* jStr_evtXX_*/M14() { DEBUG(0);
    double result = NaN;
    if(_Njets >= 4) { DEBUG(0);
      result = (*_firstJet + *_fourthJet).m();
    }
    return result;
  }

  double/* jStr_evtXX_*/cosTheta1() {
    double result = NaN;
    if(_Njets >= 4) {
      CLHEP::HepLorentzVector* a = _jetQuartet.first.first; 
      if(z1boostAxis().restMass2() <= m2min ||
	 a->restMass2() <= m2min ||
	 xboostAxis().restMass2() <= m2min ||
	 z1boostAxis().findBoostToCM().mag2() >= b2max) {
	result = -NaN;
      }
      else {
	CLHEP::Hep3Vector pa = a->boost(z1boostAxis().findBoostToCM()).v();
	CLHEP::Hep3Vector pX = xboostAxis().boost(z1boostAxis().findBoostToCM()).v();
	result = Divide(pa.dot(pX), Sqrt(pa.dot(pa) * pX.dot(pX)));
      }
    }
    return result;
  }

  double/* jStr_evtXX_*/cosTheta2() {
    double result = NaN;
    if(_Njets >= 4) {
      CLHEP::HepLorentzVector* c = _jetQuartet.second.first;
      if(z2boostAxis().restMass2() <= m2min || 
	 c->restMass2() <= m2min ||
	 xboostAxis().restMass2() <= m2min ||
	 z2boostAxis().findBoostToCM().mag2() >= b2max) {
	result = -NaN;
      }
      else {
	CLHEP::Hep3Vector pc = c->boost(z2boostAxis().findBoostToCM()).v();
	CLHEP::Hep3Vector pX = xboostAxis().boost(z2boostAxis().findBoostToCM()).v();
	result = Divide(pc.dot(pX), Sqrt(pc.dot(pc) * pX.dot(pX)));
      }
    }
    return result;
  }

  double/* jStr_evtXX_*/cosThetaStar1() {
    double result = NaN;
    if(_Njets >= 4) {
      if(xboostAxis().restMass2() <= m2min || xboostAxis().findBoostToCM().mag2() >= b2max) {
	result = -NaN;
      }
      else {
	CLHEP::Hep3Vector pZ1 = z1boostAxis().boost(xboostAxis().findBoostToCM()).v();
	result = Divide(pZ1.dot(eZPxFrame()), Sqrt(pZ1.dot(pZ1)));
      }
    }
    return result;
  }

  double/* jStr_evtXX_*/cosThetaStar2() {
    double result = NaN;
    if(_Njets >= 4) {
      if(xboostAxis().restMass2() <= m2min || xboostAxis().findBoostToCM().mag2() >= b2max) {
	result = -NaN;
      }
      else {
	CLHEP::Hep3Vector pZ2 = z2boostAxis().boost(xboostAxis().findBoostToCM()).v();
	result = Divide(pZ2.dot(eZPxFrame()), Sqrt(pZ2.dot(pZ2)));
      }
    }
    return result;
  }

  double/* jStr_evtXX_*/cosPhiTilde1() {
    double result = NaN;
    if(_Njets >= 4) {
      CLHEP::HepLorentzVector* a = _jetQuartet.first.first;
      CLHEP::HepLorentzVector* b = _jetQuartet.first.second;
      if(xboostAxis().restMass2() <= m2min || xboostAxis().findBoostToCM().mag2() >= b2max) {
	result = -NaN;
      }
      else {
	CLHEP::Hep3Vector pa = a->boost(xboostAxis().findBoostToCM()).v();
	CLHEP::Hep3Vector pb = b->boost(xboostAxis().findBoostToCM()).v();
	result = Divide((eZxFrame().cross(eZPxFrame())).dot(pa.cross(pb)), Sqrt((pa.cross(pb)).dot(pa.cross(pb))));
      }
    }
    return result;
  }

  double/* jStr_evtXX_*/cosPhiTilde2() {
    double result = NaN;
    if(_Njets >= 4) {
      CLHEP::HepLorentzVector* c = _jetQuartet.second.first;
      CLHEP::HepLorentzVector* d = _jetQuartet.second.second;
      if(xboostAxis().restMass2() <= m2min || xboostAxis().findBoostToCM().mag2() >= b2max) {
	result = -NaN;
      }
      else {
	CLHEP::Hep3Vector pc = c->boost(xboostAxis().findBoostToCM()).v();
	CLHEP::Hep3Vector pd = d->boost(xboostAxis().findBoostToCM()).v();
	result = Divide((eZxFrame().cross(eZPxFrame())).dot(pc.cross(pd)), Sqrt((pc.cross(pd)).dot(pc.cross(pd))));
      }
    }
    return result;
  }

  // http://arxiv.org/pdf/1001.3396
  // http://arxiv.org/pdf/1010.0676v2
  double/* jStr_evtXX_*/cosPhi() {
    double result = NaN;
    if(_Njets >= 4) {
      CLHEP::HepLorentzVector* a = _jetQuartet.first.first;
      CLHEP::HepLorentzVector* b = _jetQuartet.first.second;
      CLHEP::HepLorentzVector* c = _jetQuartet.second.first;
      CLHEP::HepLorentzVector* d = _jetQuartet.second.second;
      if(xboostAxis().restMass2() <= m2min || xboostAxis().findBoostToCM().mag2() >= b2max) {
	result = -NaN;
      }
      else {
	CLHEP::Hep3Vector pa = a->boost(xboostAxis().findBoostToCM()).v();
	CLHEP::Hep3Vector pb = b->boost(xboostAxis().findBoostToCM()).v();
	CLHEP::Hep3Vector pc = c->boost(xboostAxis().findBoostToCM()).v();
	CLHEP::Hep3Vector pd = d->boost(xboostAxis().findBoostToCM()).v();
	result = Divide((pa.cross(pb)).dot(pc.cross(pd)), Sqrt((pa.cross(pb)).dot(pa.cross(pb)) * (pc.cross(pd)).dot(pc.cross(pd))));
      }
    }
    return result;
  }

  double/* jStr_evtXX_*/PhiTilde1() { return Arccos(cosPhiTilde1()); }
  double/* jStr_evtXX_*/PhiTilde2() { return Arccos(cosPhiTilde2()); }
  double/* jStr_evtXX_*/Phi() { return Arccos(cosPhi()); }

private:
  CLHEP::HepLorentzVector& boost(CLHEP::HepLorentzVector i, const Jet* parentJet = 0) { DEBUG(0);
    if(parentJet != 0 && (parentJet->hlv().restMass2() <= m2min || parentJet->hlv().findBoostToCM().mag2() >= b2max)) { // !ZMxpvTachyonic
      i.setX(NaN);
      i.setY(NaN);
      i.setZ(NaN);
      i.setT(NaN);
      parentJet = 0; // signal that we are to not boost, pass back NaN HLV instead
    }
    return !parentJet? i: i.boost(parentJet->hlv().findBoostToCM());
  }
  
  const CLHEP::HepLorentzVector boostAxis() const {
    // Find boost axis of object that decayed into 2 jets
    return  *_firstJet + *_secondJet;
  }
  
  const CLHEP::Hep3Vector jetAxisCM(CLHEP::HepLorentzVector* lv) {
    // Reverse boost to translate dijet into CM frame, get back-to-back jets
    // Boost needed to get to center-of-mass  frame:
    // w.findBoostToCM() == - w.boostVector()
    // w.boost(w.findBoostToCM()) == w.rest4Vector()
    // return lv.boost(boostAxis().findBoostToCM()).v();
    return lv->boost(_firstJet->findBoostToCM(*_secondJet)).v();
  }

  // arXiv:1010.3698v1 [hep-ph]
  // Lots of unconst functions 'cos we're boosting (altering) jets that are data members
  double cosTheta(CLHEP::HepLorentzVector* j) {
    double result = NaN;
    if(_dijets) {
      CLHEP::Hep3Vector unboostedJet = jetAxisCM(j);
      result = unboostedJet.cosTheta();
    }
    return result;
  }

  // arXiv:1010.3698v1 [hep-ph]
  // Delta when jet translated to CM frame
  double deltaRapidityXtoJCM(CLHEP::HepLorentzVector* j) {
    double result = NaN;
    if(_dijets) {
      CLHEP::Hep3Vector unboostedJet = jetAxisCM(j);
      unboostedJet /= unboostedJet.mag();
      CLHEP::Hep3Vector boostaxis = boostAxis().v();
      boostaxis /= boostaxis.mag2()? boostaxis.mag(): 1.;
      result = (fabs(boostaxis.z()) >= 1. || fabs(unboostedJet.z()) >= 1.)? NaN:
	boostaxis.rapidity() - unboostedJet.rapidity();
    }
    return result;
  }

  // Don't translate to CM frame
  double deltaRapidityXtoJ(CLHEP::HepLorentzVector* j) {
    double result = NaN;
    if(_dijets) {
      CLHEP::Hep3Vector boostedJet = j->v();
      boostedJet /= boostedJet.mag();
      CLHEP::Hep3Vector boostaxis = boostAxis().v();
      boostaxis /= boostaxis.mag2()? boostaxis.mag(): 1.;
      result = (fabs(boostaxis.z()) >= 1. || fabs(boostedJet.z()) >= 1.)? NaN:
	boostaxis.rapidity() - boostedJet.rapidity();
    }
    return result;
  }
 
  // arXiv:1010.3698v1 [hep-ph]
  double cosHelicity(CLHEP::HepLorentzVector* j) {
    return _dijets? cos(helicity(j)): NaN;
  }

  // arXiv:1010.3698v1 [hep-ph]
  // Some guesswork with Hep3Vector functions here... looks OK
  double helicity(CLHEP::HepLorentzVector* j) {
    double result = NaN;
    CLHEP::HepLorentzVector boostaxis = boostAxis();
    if(_dijets && boostaxis.v().mag2()) {
      CLHEP::Hep3Vector unboostedJet = jetAxisCM(j);
      result = unboostedJet.polarAngle(boostaxis.v(), boostaxis.v());
    }
    return result;
  }

  // Some guesswork with Hep3Vector functions here... looks OK
  double azilicity(CLHEP::HepLorentzVector* j) {
    double result = NaN;
    CLHEP::HepLorentzVector boostaxis = boostAxis();
    if(_dijets && boostaxis.v().mag2()) {
      CLHEP::Hep3Vector unboostedJet = jetAxisCM(j);
      const CLHEP::Hep3Vector z = CLHEP::Hep3Vector(0., 0., 1.);
      result = z.perpPart(boostaxis.v()).mag2()? unboostedJet.azimAngle(z, boostaxis.v()): NaN;
      //if(result < 0) result += CLHEP::pi;
    }
    return result;
  }

  inline CLHEP::HepLorentzVector zboostAxis(CLHEP::HepLorentzVector* x, CLHEP::HepLorentzVector* y) {
    return (*x) + (*y);
  }
  
  inline CLHEP::HepLorentzVector xboostAxis(CLHEP::HepLorentzVector* a, CLHEP::HepLorentzVector* b, CLHEP::HepLorentzVector* c, CLHEP::HepLorentzVector* d) {
    return (*a) + (*b) + (*c) + (*d);
  }

  inline CLHEP::HepLorentzVector z1boostAxis() {
    CLHEP::HepLorentzVector* a = _jetQuartet.first.first;
    CLHEP::HepLorentzVector* b = _jetQuartet.first.second;
    return zboostAxis(a, b);
  }
  
  inline CLHEP::HepLorentzVector z2boostAxis() {
    CLHEP::HepLorentzVector* c = _jetQuartet.second.first;
    CLHEP::HepLorentzVector* d = _jetQuartet.second.second;
    return zboostAxis(c, d);
  }

  inline CLHEP::HepLorentzVector xboostAxis() {
    CLHEP::HepLorentzVector* a = _jetQuartet.first.first;
    CLHEP::HepLorentzVector* b = _jetQuartet.first.second;
    CLHEP::HepLorentzVector* c = _jetQuartet.second.first;
    CLHEP::HepLorentzVector* d = _jetQuartet.second.second;
    return (*a) + (*b) + (*c) + (*d);
  }

  inline CLHEP::Hep3Vector pZlabFrame(CLHEP::HepLorentzVector* x, CLHEP::HepLorentzVector* y) {
    return x->v() + y->v();
  }

  inline CLHEP::Hep3Vector pZ1labFrame() {
    CLHEP::HepLorentzVector* a = _jetQuartet.first.first;
    CLHEP::HepLorentzVector* b = _jetQuartet.first.second;
    return pZlabFrame(a, b);
  }

  inline CLHEP::Hep3Vector pZ2labFrame() {
    CLHEP::HepLorentzVector* c = _jetQuartet.second.first;
    CLHEP::HepLorentzVector* d = _jetQuartet.second.second;
    return pZlabFrame(c, d);
  }

  inline CLHEP::Hep3Vector pXlabFrame(CLHEP::HepLorentzVector* a, CLHEP::HepLorentzVector* b, CLHEP::HepLorentzVector* c, CLHEP::HepLorentzVector* d) {
    return a->v() + b->v() + c->v() + d->v();
  }

  inline CLHEP::Hep3Vector pXlabFrame() {
    CLHEP::HepLorentzVector* a = _jetQuartet.first.first;
    CLHEP::HepLorentzVector* b = _jetQuartet.first.second;
    CLHEP::HepLorentzVector* c = _jetQuartet.second.first;
    CLHEP::HepLorentzVector* d = _jetQuartet.second.second;
    return pXlabFrame(a, b, c, d);
  }

  inline CLHEP::Hep3Vector eZxFrame() {
    CLHEP::Hep3Vector result = CLHEP::Hep3Vector(NaN, NaN, NaN);
    CLHEP::HepLorentzVector z = CLHEP::HepLorentzVector(0., 0., 1.);
    if(xboostAxis().restMass2() <= m2min || xboostAxis().findBoostToCM().mag2() >= b2max) { 
      result = -result;
    }
    else {
      result = z.boost(xboostAxis().findBoostToCM()).v().unit();
    }
    return result;
  }

  inline CLHEP::Hep3Vector eZPxFrame() {
    CLHEP::Hep3Vector result = CLHEP::Hep3Vector(NaN, NaN, NaN);
    if(xboostAxis().restMass2() <= m2min || xboostAxis().findBoostToCM().mag2() >= b2max) { 
      result = -result;
    }
    else {
      CLHEP::HepLorentzVector Z = CLHEP::HepLorentzVector(z1boostAxis().v().unit());
      result = Z.boost(xboostAxis().findBoostToCM()).v().unit();
    }
    return result;
  }
  
  // Adder functor
  struct AddE: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddE(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += x->e(); }
  };

  struct AddET: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddET(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += x->et(); }
  };

  struct AddPT: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddPT(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += x->perp(); }
  };

  struct AddM: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddM(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += x->m(); }
  };

  struct AddPX: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddPX(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += x->px(); }
  };

  struct AddPY: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddPY(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += x->py(); }
  };

  // Sorter functor
  struct SortByEtDown {
    bool operator()(CLHEP::HepLorentzVector* a, CLHEP::HepLorentzVector* b) {
      return a->et() > b->et();
    }
  } sortByEtDown;
  
  // Sorter functor
  struct Jsort {
    bool operator()(CLHEP::HepLorentzVector* a, CLHEP::HepLorentzVector* b) { return a->e() < b->e(); }
  } jetEnergySort;

  typedef std::pair<CLHEP::HepLorentzVector*, CLHEP::HepLorentzVector*> JetPair;
  typedef std::pair<JetPair, JetPair> JetQuartet;

  struct PairSortByRapidity {
    bool operator()(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return a->rapidity() < b->rapidity();
    }
  } ySort;

  inline JetPair rapidityOrder(JetPair& ab) {
    CLHEP::HepLorentzVector* a = ab.first;
    CLHEP::HepLorentzVector* b = ab.second;
    return a->rapidity() < b->rapidity()? std::make_pair(a, b): std::make_pair(b, a);
  }

  struct PairSortByDeltaM {
    double dM(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return fabs(a->m() - b->m());
    }
    double dM(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first; DEBUG(a);
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(b);
      return dM(a, b);
    }
    bool operator()(JetPair& x, JetPair& y) { DEBUG(dM(x)); DEBUG(dM(y));
      return dM(x) < dM(y);
    }
  } deltaMsort;

  struct QuartetSortByMprime {
    double dR(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      double dEta = fabs(a->eta() - b->eta());
      double dPhi = DeltaPhi()(a->phi(), b->phi());
      return Sqrt((dEta * dEta) + (dPhi * dPhi));
    }
    double dR(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return dR(a, b);
    }
    double mass(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a + *b).m();
    }
    double mass(JetPair& ab) { // i.e., M1 or M2
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return mass(a, b);
    }
    double dM(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second;
      return fabs(mass(ab) - mass(cd)); // |M1 - M2|
    }
    double Mavg(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second;
      return 0.5 * (mass(ab) + mass(cd)); // |M1 + M2|/2
    }
    double Mprime(JetQuartet& abcd) {
      return Divide(dM(abcd), Mavg(abcd)); 
    }
    double sumPt(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a).perp() + (*b).perp();
    }
    double sumPt(JetPair& ab) { // i.e., M1 or M2
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return sumPt(a, b);
    }
    double Delta(JetPair& xy, JetQuartet& abcd) {
      return sumPt(xy) - Mavg(abcd); 
    }
    double Delta1(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      return Delta(ab, abcd); 
    }
    double Delta2(JetQuartet& abcd) {
      JetPair cd = abcd.second;
      return Delta(cd, abcd); 
    }
    double DeltaMin(JetQuartet& abcd) {
      return std::min(Delta1(abcd), Delta2(abcd)); 
    }
    double DeltaMax(JetQuartet& abcd) {
      return std::max(Delta1(abcd), Delta2(abcd)); 
    }
    double YX(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a + *b).rapidity();
    }
    double YX(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return YX(a, b);
    }
    double DeltaYXX(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second;
      return fabs(YX(ab) - YX(cd)); 
    }
    double PhiX(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a + *b).phi();
    }
    double PhiX(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return PhiX(a, b);
    }
    double DeltaPhiXX(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second;
      return DeltaPhi()(PhiX(ab), PhiX(cd));
    }
    double DeltaRXX(JetQuartet& abcd) {
      double dPhi = DeltaPhiXX(abcd);
      double dY = DeltaYXX(abcd);
      return Sqrt((dY * dY) + (dPhi * dPhi));
    }
    double TwistXX(JetQuartet& abcd) {
      double dPhi = DeltaPhiXX(abcd);
      double dY = DeltaYXX(abcd);
      return Arctan2(fabs(dPhi), fabs(dY));
    }
    int separated(JetQuartet& abcd, double R) {
      CLHEP::HepLorentzVector* a = abcd.first.first;
      CLHEP::HepLorentzVector* b = abcd.first.second;
      CLHEP::HepLorentzVector* c = abcd.second.first;
      CLHEP::HepLorentzVector* d = abcd.second.second;
      JetPair j1j2 = std::make_pair(a, b); DEBUG(0);
      JetPair j3j4 = std::make_pair(c, d); DEBUG(0);
      JetPair j1j3 = std::make_pair(a, c); DEBUG(0);
      JetPair j2j4 = std::make_pair(b, d); DEBUG(0);
      JetPair j1j4 = std::make_pair(a, d); DEBUG(0);
      JetPair j2j3 = std::make_pair(b, c); DEBUG(0);
      return
	(dR(j1j2) > R &&
	 dR(j3j4) > R &&
	 dR(j1j3) > R &&
	 dR(j2j4) > R &&
	 dR(j1j4) > R &&
	 dR(j2j3) > R)? 1: 0;
    }
    bool operator()(JetQuartet& x, JetQuartet& y) { 
      return Mprime(x) < Mprime(y);
    }
  } mprimeSort;

  inline JetPair massOrder(JetPair& ab) {
    CLHEP::HepLorentzVector* a = ab.first;
    CLHEP::HepLorentzVector* b = ab.second;
    return a->m() < b->m()? std::make_pair(a, b): std::make_pair(b, a);
  }

  struct QuartetSortByDeltaR {
    double dR(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      double dEta = fabs(a->eta() - b->eta());
      double dPhi = DeltaPhi()(a->phi(), b->phi());
      return Sqrt((dEta * dEta) + (dPhi * dPhi));
    }
    double dR(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return dR(a, b);
    }
    double dR(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second;
      return fabs(dR(ab)) + fabs(dR(cd));
    }
    bool operator()(JetQuartet& x, JetQuartet& y) {
      return dR(x) < dR(y);
    }
  } dRsort;

  struct PairSortByM {
    double mass(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a + *b).m();
    }
    double mass(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return mass(a, b);
    }
    bool operator()(JetPair& x, JetPair& y) {
      return mass(x) < mass(y);
    }
  } mSort;

  struct QuartetSortByMass {
    double mass(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a + *b).m();
    }
    double mass(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(mass(a,b));
      return mass(a, b);
    }
    double dM(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second; DEBUG(mass(ab)); DEBUG(mass(cd)); DEBUG(fabs(mass(ab) - mass(cd)));
      return fabs(mass(ab) - mass(cd));
    }
    double max(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second; DEBUG(mass(ab)); DEBUG(mass(cd)); DEBUG(fabs(std::max(mass(ab), mass(cd))));
      return fabs(std::max(mass(ab), mass(cd)));
    }
    double Q(JetQuartet& abcd) { DEBUG(dM(abcd)); DEBUG(max(abcd)); DEBUG(Divide(dM(abcd), max(abcd)));
      return Divide(dM(abcd), max(abcd));
    }
    bool operator()(JetQuartet& x, JetQuartet& y) { DEBUG(Q(x)); DEBUG(Q(y));
      return Q(x) < Q(y);
    }
    bool operator()(JetQuartet& x) { DEBUG(Q(x));
      return Q(x);
    }
  } qSort;

  struct QuartetSortBySCDF {
    double pTvSum(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      double px = (*a).x() + (*b).x();
      double py = (*a).y() + (*b).y();
      return sqrt((px * px) + (py * py));
    }
    double pTvSum(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(pTvSum(a,b));
      return pTvSum(a, b);
    }
    double pTsum(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a).perp() + (*b).perp();
    }
    double pTsumSqrt(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return sqrt(pTsum(a, b));
    }
    double pTsum(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(pTsum(a,b));
      return pTsum(a, b);
    }
    double pTsumSqrt(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(pTsumSqrt(a,b));
      return pTsumSqrt(a, b);
    }
    double SCDF(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second; DEBUG(pTvSum(ab)); DEBUG(pTvSum(cd)); DEBUG(fabs(pTvSum(ab) - pTvSum(cd)));
      double Qab = Divide(pTvSum(ab), pTsumSqrt(ab));
      double Qcd = Divide(pTvSum(cd), pTsumSqrt(cd));
      return sqrt(0.5* ((Qab * Qab) + (Qcd * Qcd)));
    }
    double SD0(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second; DEBUG(pTvSum(ab)); DEBUG(pTvSum(cd)); DEBUG(fabs(pTvSum(ab) - pTvSum(cd)));
      double Qab = Divide(pTvSum(ab), pTsum(ab));
      double Qcd = Divide(pTvSum(cd), pTsum(cd));
      return sqrt(0.5* ((Qab * Qab) + (Qcd * Qcd)));
    }
    double dPhiIJ(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      CLHEP::Hep2Vector pa = CLHEP::Hep2Vector((*a).x(), (*a).y());
      CLHEP::Hep2Vector pb = CLHEP::Hep2Vector((*b).x(), (*b).y());
      return DeltaPhi()(pa.phi(), pb.phi());
    }
    double dPhiKL(JetQuartet& abcd) {
      JetPair cd = abcd.second;
      CLHEP::HepLorentzVector* c = cd.first;
      CLHEP::HepLorentzVector* d = cd.second;
      CLHEP::Hep2Vector pc = CLHEP::Hep2Vector((*c).x(), (*c).y());
      CLHEP::Hep2Vector pd = CLHEP::Hep2Vector((*d).x(), (*d).y());
      return DeltaPhi()(pc.phi(), pd.phi());
    }
    double PtIPtJ(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return Divide(a->perp(), b->perp());
    }
    double PtKPtL(JetQuartet& abcd) {
      JetPair cd = abcd.second;
      CLHEP::HepLorentzVector* c = cd.first;
      CLHEP::HepLorentzVector* d = cd.second;
      return Divide(c->perp(), d->perp());
    }
    double DeltaS(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second;
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      CLHEP::HepLorentzVector* c = cd.first;
      CLHEP::HepLorentzVector* d = cd.second;
      CLHEP::Hep2Vector pa = CLHEP::Hep2Vector((*a).x(), (*a).y());
      CLHEP::Hep2Vector pb = CLHEP::Hep2Vector((*b).x(), (*b).y());
      CLHEP::Hep2Vector pc = CLHEP::Hep2Vector((*c).x(), (*c).y());
      CLHEP::Hep2Vector pd = CLHEP::Hep2Vector((*d).x(), (*d).y());
      double phiab = (pa + pb).phi();
      double phicd = (pc + pd).phi();
      return DeltaPhi()(phiab, phicd);
    }
    bool operator()(JetQuartet& x, JetQuartet& y) { DEBUG(SCDF(x)); DEBUG(SCDF(y));
      return SCDF(x) < SCDF(y);
    }
    bool operator()(JetQuartet& x) { DEBUG(SCDF(x));
      return SCDF(x);
    }
  } scdfSort;

  struct QuartetSortBySD0 {
    double pTvSum(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      double px = (*a).x() + (*b).x();
      double py = (*a).y() + (*b).y();
      return sqrt((px * px) + (py * py));
    }
    double pTvSum(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(pTvSum(a,b));
      return pTvSum(a, b);
    }
    double pTsum(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return (*a).perp() + (*b).perp();
    }
    double pTsumSqrt(CLHEP::HepLorentzVector*& a, CLHEP::HepLorentzVector*& b) {
      return sqrt(pTsum(a, b));
    }
    double pTsum(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(pTsum(a,b));
      return pTsum(a, b);
    }
    double pTsumSqrt(JetPair& ab) {
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second; DEBUG(pTsumSqrt(a,b));
      return pTsumSqrt(a, b);
    }
    double SCDF(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second; DEBUG(pTvSum(ab)); DEBUG(pTvSum(cd)); DEBUG(fabs(pTvSum(ab) - pTvSum(cd)));
      double Qab = Divide(pTvSum(ab), pTsumSqrt(ab));
      double Qcd = Divide(pTvSum(cd), pTsumSqrt(cd));
      return sqrt(0.5* ((Qab * Qab) + (Qcd * Qcd)));
    }
    double SD0(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second; DEBUG(pTvSum(ab)); DEBUG(pTvSum(cd)); DEBUG(fabs(pTvSum(ab) - pTvSum(cd)));
      double Qab = Divide(pTvSum(ab), pTsum(ab));
      double Qcd = Divide(pTvSum(cd), pTsum(cd));
      return sqrt(0.5* ((Qab * Qab) + (Qcd * Qcd)));
    }
    double dPhiIJ(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      CLHEP::Hep2Vector pa = CLHEP::Hep2Vector((*a).x(), (*a).y());
      CLHEP::Hep2Vector pb = CLHEP::Hep2Vector((*b).x(), (*b).y());
      return DeltaPhi()(pa.phi(), pb.phi());
    }
    double dPhiKL(JetQuartet& abcd) {
      JetPair cd = abcd.second;
      CLHEP::HepLorentzVector* c = cd.first;
      CLHEP::HepLorentzVector* d = cd.second;
      CLHEP::Hep2Vector pc = CLHEP::Hep2Vector((*c).x(), (*c).y());
      CLHEP::Hep2Vector pd = CLHEP::Hep2Vector((*d).x(), (*d).y());
      return DeltaPhi()(pc.phi(), pd.phi());
    }
    double PtIPtJ(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      return Divide(a->perp(), b->perp());
    }
    double PtKPtL(JetQuartet& abcd) {
      JetPair cd = abcd.second;
      CLHEP::HepLorentzVector* c = cd.first;
      CLHEP::HepLorentzVector* d = cd.second;
      return Divide(c->perp(), d->perp());
    }
    double DeltaS(JetQuartet& abcd) {
      JetPair ab = abcd.first;
      JetPair cd = abcd.second;
      CLHEP::HepLorentzVector* a = ab.first;
      CLHEP::HepLorentzVector* b = ab.second;
      CLHEP::HepLorentzVector* c = cd.first;
      CLHEP::HepLorentzVector* d = cd.second;
      CLHEP::Hep2Vector pa = CLHEP::Hep2Vector((*a).x(), (*a).y());
      CLHEP::Hep2Vector pb = CLHEP::Hep2Vector((*b).x(), (*b).y());
      CLHEP::Hep2Vector pc = CLHEP::Hep2Vector((*c).x(), (*c).y());
      CLHEP::Hep2Vector pd = CLHEP::Hep2Vector((*d).x(), (*d).y());
      double phiab = (pa + pb).phi();
      double phicd = (pc + pd).phi();
      return DeltaPhi()(phiab, phicd);
    }
    bool operator()(JetQuartet& x, JetQuartet& y) { DEBUG(SD0(x)); DEBUG(SD0(y));
      return SD0(x) < SD0(y);
    }
    bool operator()(JetQuartet& x) { DEBUG(SD0(x));
      return SD0(x);
    }
  } sd0Sort;

  double _deltaR;
  std::vector<CLHEP::HepLorentzVector*> _jets;
  const unsigned int _Njets;
  const bool _dijets;
  CLHEP::HepLorentzVector* _firstJet;
  CLHEP::HepLorentzVector* _secondJet;
  CLHEP::HepLorentzVector* _thirdJet;
  CLHEP::HepLorentzVector* _fourthJet;
  
  JetQuartet _jetQuartet, _jetQuartet2, _jetQuartet3, _jetQuartet4, _jetQuartet5, _jetQuartet6;
  std::vector<JetQuartet> _4jetCombinations, _4jetCombinations2, _4jetCombinations3;
  JetPair _jetPair;
  std::vector<JetPair> _2jetCombinations;
};

#endif

