#ifndef SPHERICITYSPHEROCITY_H
#define SPHERICITYSPHEROCITY_H

#include "JetEvent/Jet.h"
#include "JetEvent/JetConstituentIterator.h"
#include "JetUtils/JetCollectionHelper.h"
#include "CLHEP/Matrix/Matrix.h" // For SphericitySpherocity

#include "UserAnalysis/Hep2Matrix.h" // Linearised momentum tensor
#include "UserAnalysis/EigensolverSym3x3.h" // For SphericitySpherocity
#include "UserAnalysis/Numeric.h" // Common numeric stuff

#include <cmath> // For isnan(x)

// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

using namespace Numeric;

// arXiv:1010.3698v1 [hep-ph]
class SphericitySpherocity {
public:
  SphericitySpherocity(const JetCollection* theJetCollection, const Jet* parentJet = 0, bool useConstituents = false): 
  _theJetCollection(theJetCollection), _pmag(0), _pmag2(0), _pmagSum(0), _pmag2Sum(0), _pTsum(0), _pT2sum(0), _sumB(NaN),
  _sphericityLambda1(NaN), _sphericityLambda2(NaN), _sphericityLambda3(NaN),
  _spherocityLambda1(NaN), _spherocityLambda2(NaN), _spherocityLambda3(NaN) {
    _sphericityTensor = CLHEP::HepMatrix(3, 3, 0.);
    _spherocityTensor = CLHEP::HepMatrix(3, 3, 0.);

    _Mlin = Hep2Matrix(0.);
    
    if(_theJetCollection) {
      compute(parentJet, useConstituents);
      
      double lambda[3] = {NaN, NaN, NaN};
      if(normalisedEigenvalues(sphericityTensor(), lambda)) {
	_sphericityLambda1 = lambda[0];
	_sphericityLambda2 = lambda[1];
	_sphericityLambda3 = lambda[2];
      }
      
      if(normalisedEigenvalues(spherocityTensor(), lambda)) {
	_spherocityLambda1 = lambda[0];
	_spherocityLambda2 = lambda[1];
	_spherocityLambda3 = lambda[2];
      }
    }
  }
  
  // http://arxiv.org/pdf/0806.0023
  double/* jStr_sphXX_*/detSphericity() const { return sphericityTensor()? sphericityTensor()->determinant(): NaN; }
  double/* jStr_sphXX_*/detSpherocity() const { return spherocityTensor()? spherocityTensor()->determinant(): NaN; }

  double/* jStr_sphXX_*/sphericityLambda1() const { return _sphericityLambda1; }
  double/* jStr_sphXX_*/sphericityLambda2() const { return _sphericityLambda2; }
  double/* jStr_sphXX_*/sphericityLambda3() const { return _sphericityLambda3; }
  double/* jStr_sphXX_*/spherocityLambda1() const { return _spherocityLambda1; } 
  double/* jStr_sphXX_*/spherocityLambda2() const { return _spherocityLambda2; }
  double/* jStr_sphXX_*/spherocityLambda3() const { return _spherocityLambda3; }

  // CMS-AN-2012/004
  // http://www.phys.uniroma1.it/DipWeb/dottorato/DOTT_FISICA/MENU/03DOTTORANDI/Seminari24/Semin_24_ott2011/Pandolfi2011.pdf
  // http://castello.web.cern.ch/castello/cms-note-higgs2l2j.pdf
  double/* jStr_sphXX_*/pft() const {
    return sqrt(Divide(_pT2sum, _pTsum * _pTsum));;
  }

  // https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=166288
  // pt-sum that suppressed forward tracks
  // http://arxiv.org/pdf/1005.4060v2
  double/* jStr_sphXX_*/beamThrust() const {
    return _sumB;
  }

  // http://elpub.bib.uni-wuppertal.de/servlets/DerivateServlet/Derivate-1480/dc1006.pdf
  double/* jStr_sphXX_*/circularity() const {
    return 2.* Divide(std::min(_sphericityLambda1, _sphericityLambda2), _sphericityLambda1 + _sphericityLambda2);
  }
  
  double/* jStr_sphXX_*/sphericity() const {
    return (3./2.) * (_sphericityLambda2 + _sphericityLambda3);
  }
  
  double/* jStr_sphXX_*/spherocity() const {
    return (3./2.) * (_spherocityLambda2 + _spherocityLambda3);
  }
  
  double/* jStr_sphXX_*/aplanarity() const {
    return (3./2.) * _sphericityLambda3;
  }
  
  double/* jStr_sphXX_*/aplanority() const {
    return (3./2.) * _spherocityLambda3;
  }
  
  double/* jStr_sphXX_*/Y() const {
    return (sqrt(3.)/2.) * (_sphericityLambda2 - _sphericityLambda3);
  }
  
  double/* jStr_sphXX_*/planarity() const {
    return (_sphericityLambda2 - _sphericityLambda3);
  }

  double/* jStr_sphXX_*/planority() const {
    return (_spherocityLambda2 - _spherocityLambda3);
  }

  double/* jStr_sphXX_*/Dshape() const {
    return 27.* _spherocityLambda1 * _spherocityLambda2 * _spherocityLambda3;
  }

  // http://www-conf.kek.jp/dis06/transparencies/WG4/hfs-tasevsky.ppt
  // http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node213.html
  double/* jStr_sphXX_*/Cshape() const {
    return 3.* ((_spherocityLambda1 * _spherocityLambda2) + (_spherocityLambda2 * _spherocityLambda3) + (_spherocityLambda3 * _spherocityLambda1));
  }

  // Check; C = 1 - 2nd Fox-Wolfram moment
  double/* jStr_sphXX_*/H2() const { 
    return 1. - Cshape();
  }
  // http://arxiv.org/pdf/1001.4082
  double/* jStr_sphXX_*/Fshape() const { 
    double e1 = _Mlin.eigenvalues().first;
    double e2 = _Mlin.eigenvalues().second;
    return Divide(e2, e1);
  }
  double/* jStr_sphXX_*/G() const { 
    double f = Fshape();
    return Divide(4. * f, (1. + f) * (1. + f));
  }
  // http://arxiv.org/pdf/1001.4082 2D transverse sphericity
  double/* jStr_sphXX_*/ST2D() const { 
    double e1 = _Mlin.eigenvalues().first;
    double e2 = _Mlin.eigenvalues().second;
    return Divide(2. * e2, e1 + e2);
  }

  double/* jStr_sphXX_*/detMlin() const { 
    return _Mlin.determinant();
  }

private:
  CLHEP::HepLorentzVector& boost(CLHEP::HepLorentzVector i, const Jet* parentJet = 0) {
    if(parentJet != 0 && (parentJet->hlv().restMass2() <= m2min || parentJet->hlv().findBoostToCM().mag2() >= b2max)) { // !ZMxpvTachyonic
      i.setX(NaN);
      i.setY(NaN);
      i.setZ(NaN);
      i.setT(NaN);
      parentJet = 0; // signal that we are to not boost, pass back NaN HLV instead
    }
    return !parentJet? i: i.boost(parentJet->hlv().findBoostToCM());
  }

  void compute(const Jet* parentJet = 0, bool useConstituents = false) {
    if(!_theJetCollection || _theJetCollection->size() == 0) {
      return;
    }
    double px, py, pz;
    
    if(useConstituents && parentJet && parentJet->size()) { // loop over constituents in jet boosted into own rest frame
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(parentJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(parentJet);
      _sumB = 0.;
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	px = boost(firstConstituent.hlv(), parentJet).px();
	py = boost(firstConstituent.hlv(), parentJet).py();
	pz = boost(firstConstituent.hlv(), parentJet).pz();
	fillTensors(px, py, pz);
      }
    }
    else { // loop over jets, or subjets if looking at jet boosted into own rest frame
      JetCollectionHelper::jetcollection_t theJets(_theJetCollection->begin(), _theJetCollection->end());
      JetCollectionHelper::jetcollection_t::iterator first = theJets.begin();
      JetCollectionHelper::jetcollection_t::iterator last  = theJets.end();
      _sumB = 0.;
      for(; first != last; ++first) {
	if((*first) == 0) continue;
	const Jet* aJet = (*first)->clone(true, true);
	if(aJet == 0) continue;
	px = boost(aJet->hlv(), parentJet).px();
	py = boost(aJet->hlv(), parentJet).py();
	pz = boost(aJet->hlv(), parentJet).pz();
	fillTensors(px, py, pz);
	delete aJet;
      }
    }
    // Matricies are symmetric:
    _sphericityTensor[1][0] = _sphericityTensor[0][1];
    _sphericityTensor[2][0] = _sphericityTensor[0][2];
    _sphericityTensor[2][1] = _sphericityTensor[1][2];

    _spherocityTensor[1][0] = _spherocityTensor[0][1];
    _spherocityTensor[2][0] = _spherocityTensor[0][2];
    _spherocityTensor[2][1] = _spherocityTensor[1][2];

    _sphericityTensor /= _pmag2Sum? _pmag2Sum: NaN;
    _spherocityTensor /= _pmagSum? _pmagSum: NaN;
  }

  void fillTensors(double px, double py, double pz) {
    _pmag2 = (px * px) + (py * py) + (pz * pz);
    _pmag = Sqrt(_pmag2);
    _invpmag = Divide(1., _pmag);
    _pmagSum += _pmag;
    _pmag2Sum += _pmag2;

    _pT = sqrt((px * px) + (py * py));
    _pT2 = _pT * _pT;
    _pTsum += _pT;
    _pT2sum += _pT2;
    
    _sphericityTensor[0][0] += px * px;
    _sphericityTensor[0][1] += px * py;
    _sphericityTensor[0][2] += px * pz;
    _sphericityTensor[1][1] += py * py;
    _sphericityTensor[1][2] += py * pz;
    _sphericityTensor[2][2] += pz * pz;
    
    _spherocityTensor[0][0] += px * px * _invpmag;
    _spherocityTensor[0][1] += px * py * _invpmag;
    _spherocityTensor[0][2] += px * pz * _invpmag;
    _spherocityTensor[1][1] += py * py * _invpmag;
    _spherocityTensor[1][2] += py * pz * _invpmag;
    _spherocityTensor[2][2] += pz * pz * _invpmag;

    double pt = Sqrt((px * px) + (py * py));
    _Mlin[0][0] += Divide(px * px, pt);
    _Mlin[0][1] += Divide(px * py, pt);
    _Mlin[1][0] += Divide(px * py, pt);
    _Mlin[1][1] += Divide(py * py, pt);

    double eta = fabs(pseudoRapidity(_pmag, pz));
    _sumB += pt * exp(-1.* pt * eta);
  }

  const CLHEP::HepMatrix* sphericityTensor() const {
    return _theJetCollection? &_sphericityTensor: 0;
  }

  const CLHEP::HepMatrix* spherocityTensor() const {
    return _theJetCollection? &_spherocityTensor: 0;
  }

  bool normalisedEigenvalues(const CLHEP::HepMatrix* M, double* root) const {
    EigenSolverSym3x3 solve(*M);
    solve(root);
    double sum = NaN;
    if(!isnan(root[0]) && !isnan(root[1]) && !isnan(root[2])) sum = root[0] + root[1] + root[2];
    double norm = Divide(1., sum);
    root[0] = Multiply(root[0], norm);
    root[1] = Multiply(root[1], norm);
    root[2] = Multiply(root[2], norm);
    return sum != 0.;
  }

  const JetCollection* _theJetCollection;
  CLHEP::HepMatrix _sphericityTensor;
  CLHEP::HepMatrix _spherocityTensor;
  Hep2Matrix _Mlin;
  double _pmag;
  double _pmag2;
  double _pmagSum;
  double _pmag2Sum;
  double _pT;
  double _pT2;
  double _pTsum;
  double _pT2sum;
  double _sumB;
  double _invpmag;
  double _sphericityLambda1, _sphericityLambda2, _sphericityLambda3;
  double _spherocityLambda1, _spherocityLambda2, _spherocityLambda3;
};

// Older implementation
// arXiv:1010.3698v1 [hep-ph]
// class SphericitySpherocity {
// public:
//   SphericitySpherocity(const JetCollection* theJetCollection): 
//   _theJetCollection(theJetCollection), _pmag(0), _pmag2(0), _pmagSum(0), _pmag2Sum(0),
//   _sphericityLambda1(NaN), _sphericityLambda2(NaN), _sphericityLambda3(NaN),
//   _spherocityLambda1(NaN), _spherocityLambda2(NaN), _spherocityLambda3(NaN) {
//     _sphericityTensor = CLHEP::HepMatrix(3, 3, 0.);
//     _spherocityTensor = CLHEP::HepMatrix(3, 3, 0.);
    
//     if(_theJetCollection) {
//       compute();
      
//       double lambda[3] = {NaN, NaN, NaN};
//       if(normalisedEigenvalues(sphericityTensor(), lambda)) {
// 	_sphericityLambda1 = lambda[0];
// 	_sphericityLambda2 = lambda[1];
// 	_sphericityLambda3 = lambda[2];
//       }
      
//       if(normalisedEigenvalues(spherocityTensor(), lambda)) {
// 	_spherocityLambda1 = lambda[0];
// 	_spherocityLambda2 = lambda[1];
// 	_spherocityLambda3 = lambda[2];
//       }
//     }
//   }

//   // http://elpub.bib.uni-wuppertal.de/servlets/DerivateServlet/Derivate-1480/dc1006.pdf
//   double circularity() const {
//     return 2.* Divide(std::min(_sphericityLambda1, _sphericityLambda2), _sphericityLambda1 + _sphericityLambda2);
//   }

//   double sphericity() const {
//     return (3./2.) * (_sphericityLambda2 + _sphericityLambda3);
//   }

//   double spherocity() const {
//     return (3./2.) * (_spherocityLambda2 + _spherocityLambda3);
//   }

//   double aplanarity() const {
//     return (3./2.) * _sphericityLambda3;
//   }

//   double aplanority() const {
//     return (3./2.) * _spherocityLambda3;
//   }

//   double Y() const {
//     return (sqrt(3.)/2.) * (_sphericityLambda2 - _sphericityLambda3);
//   }

//   double DShape() const {
//     return 27.* _spherocityLambda1 * _spherocityLambda2 * _spherocityLambda3;
//   }
  
//   double Cshape() const {
//     return 3.* ((_spherocityLambda1 * _spherocityLambda2) + (_spherocityLambda2 * _spherocityLambda3) + (_spherocityLambda3 * _spherocityLambda1));
//   }

// private:
//   void compute() {
//     if(!_theJetCollection || _theJetCollection->size() == 0) {
//       return;
//     }
//     JetCollectionHelper::jetcollection_t theJets(_theJetCollection->begin(), _theJetCollection->end());
//     JetCollectionHelper::jetcollection_t::iterator first = theJets.begin();
//     JetCollectionHelper::jetcollection_t::iterator last  = theJets.end();
    
//     double px, py, pz;

//     for(; first != last; ++first) {
//       if((*first) == 0) continue;
//       const Jet* aJet = (*first)->clone(true, true);
//       if(aJet == 0) continue;

//       px = aJet->px();
//       py = aJet->py();
//       pz = aJet->pz();

//       _pmag2 = (px * px) + (py * py) + (pz * pz);
//       _pmag = Sqrt(_pmag2);
//       _invpmag = Divide(1., _pmag);
//       _pmagSum += _pmag;
//       _pmag2Sum += _pmag2;

//       _sphericityTensor[0][0] += px * px;
//       _sphericityTensor[0][1] += px * py;
//       _sphericityTensor[0][2] += px * pz;
//       _sphericityTensor[1][1] += py * py;
//       _sphericityTensor[1][2] += py * pz;
//       _sphericityTensor[2][2] += pz * pz;

//       _spherocityTensor[0][0] += px * px * _invpmag;
//       _spherocityTensor[0][1] += px * py * _invpmag;
//       _spherocityTensor[0][2] += px * pz * _invpmag;
//       _spherocityTensor[1][1] += py * py * _invpmag;
//       _spherocityTensor[1][2] += py * pz * _invpmag;
//       _spherocityTensor[2][2] += pz * pz * _invpmag;
//       delete aJet;
//     }
//     // Matricies are symmetric:
//     _sphericityTensor[1][0] = _sphericityTensor[0][1];
//     _sphericityTensor[2][0] = _sphericityTensor[0][2];
//     _sphericityTensor[2][1] = _sphericityTensor[1][2];

//     _spherocityTensor[1][0] = _spherocityTensor[0][1];
//     _spherocityTensor[2][0] = _spherocityTensor[0][2];
//     _spherocityTensor[2][1] = _spherocityTensor[1][2];

//     _sphericityTensor /= _pmag2Sum;
//     _spherocityTensor /= _pmagSum;
//   }

//   const CLHEP::HepMatrix* sphericityTensor() const {
//     return _theJetCollection? &_sphericityTensor: 0;
//   }

//   const CLHEP::HepMatrix* spherocityTensor() const {
//     return _theJetCollection? &_spherocityTensor: 0;
//   }

//   bool normalisedEigenvalues(const CLHEP::HepMatrix* M, double* root) const {
//     EigenSolverSym3x3 solve(*M);
//     solve(root);
//     double sum = NaN;
//     if(!isnan(root[0]) && !isnan(root[1]) && !isnan(root[2])) sum = root[0] + root[1] + root[2];
//     double norm = Divide(1., sum);
//     root[0] = Multiply(root[0], norm);
//     root[1] = Multiply(root[1], norm);
//     root[2] = Multiply(root[2], norm);
//     return sum != 0.;
//   }

//   const JetCollection* _theJetCollection;
//   CLHEP::HepMatrix _sphericityTensor;
//   CLHEP::HepMatrix _spherocityTensor;
//   double _pmag;
//   double _pmag2;
//   double _pmagSum;
//   double _pmag2Sum;
//   double _invpmag;
//   double _sphericityLambda1, _sphericityLambda2, _sphericityLambda3;
//   double _spherocityLambda1, _spherocityLambda2, _spherocityLambda3;
// };

#endif // SPHERICITYSPHEROCITY_H

