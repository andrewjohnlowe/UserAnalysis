#ifndef FOXWOLFRAM_H
#define FOXWOLFRAM_H

#include "JetEvent/Jet.h"
#include "JetEvent/JetConstituentIterator.h" // For Fox-Wolfram in CM frame exclusively
#include "CLHEP/Vector/ThreeVector.h" // Common to most variables
#include "CLHEP/Vector/LorentzVector.h" // For event shapes only

#include "UserAnalysis/Numeric.h" // Common numeric stuff

#include <vector>
#include <algorithm> // for sorting

// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

// #define CHECK_FPE							\
//   if(FPEtest()) {							\
//     std::cout << "in " << __FUNCTION__					\
// 	      << ", " << __FILE__ << ":" << __LINE__ << std::endl;	\
//     feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);		\
//   }

using namespace Numeric;

namespace SApolynomials {
  template<unsigned int n> inline double S(CLHEP::Hep3Vector& i, CLHEP::Hep3Vector& j, CLHEP::Hep3Vector& k) { 
    return
      n == 1u? S<1u>(i, j, k):
      n == 2u? S<2u>(i, j, k):
      n == 3u? S<3u>(i, j, k):
      n == 4u? S<4u>(i, j, k):
      NaN;
  }
  template<> inline double S<1u>(CLHEP::Hep3Vector& /*i*/, CLHEP::Hep3Vector& /*j*/, CLHEP::Hep3Vector& /*k*/) { return 1.; }
  template<> inline double S<2u>(CLHEP::Hep3Vector& i, CLHEP::Hep3Vector& j, CLHEP::Hep3Vector& k) { 
    double ij = i.dot(j);
    double jk = j.dot(k);
    double ki = k.dot(i);
    return ij + jk + ki;
  }
  template<> inline double S<3u>(CLHEP::Hep3Vector& i, CLHEP::Hep3Vector& j, CLHEP::Hep3Vector& k) { 
    double ij = i.dot(j);
    double jk = j.dot(k);
    double ki = k.dot(i);
    return (ij * jk) + (ki * ij) + (jk * ki);
  }
  template<> inline double S<4u>(CLHEP::Hep3Vector& i, CLHEP::Hep3Vector& j, CLHEP::Hep3Vector& k) { 
    double ij = i.dot(j);
    double jk = j.dot(k);
    double ki = k.dot(i);
    return (ij * ij) + (jk * jk) + (ki * ki);
  }

  inline double A(CLHEP::Hep3Vector& i, CLHEP::Hep3Vector& j, CLHEP::Hep3Vector& k) {
    double ij = (i).dot(j);
    double jk = (j).dot(k);
    double ki = (k).dot(i);
    return ((ki * ki * jk) + (ij * ij * ki) + (jk * jk * ij)) - ((ij * ij * jk) + (ki * ki * ij) + (jk * jk * ki));
  }
}

using namespace SApolynomials;

class FoxWolfram {
public:
  FoxWolfram(const JetCollection* theJetCollection = 0, const Jet* parentJet = 0, bool useConstituents = false, unsigned int maxNum4Vectors = 1000):
  _z(CLHEP::Hep3Vector(0., 0., 1.)),
  _E(NaN),
  _xiPlus(NaN), _xiMinus(NaN),
  _xPlus(NaN), _xMinus(NaN),
  _fourMomenta(0), _threeMomenta(0), _unitVectors(0) { DEBUG(0);
    if(useConstituents) { DEBUG(0);
      _fourMomenta = fourMomenta(parentJet);
      if(_fourMomenta && _fourMomenta->size() <= 1) { DEBUG(0);
	for(unsigned int i = 0; i != _fourMomenta->size(); i++) { DEBUG(0);
	  delete _fourMomenta->at(i); DEBUG(0);
	  _fourMomenta->at(i) = 0; DEBUG(0);
	}
	_fourMomenta->clear(); DEBUG(0);
	delete _fourMomenta; // guess I need to do this if I'm doing:
	_fourMomenta = 0; // weird crash in opt build if I don't do this
      }
    }
    else { DEBUG(0);
      _fourMomenta = fourMomenta(theJetCollection, parentJet);
    }
    if(_fourMomenta) { DEBUG(0);
      std::sort(_fourMomenta->begin(), _fourMomenta->end(), sortByMomentumDown);
      unsigned int N = maxNum4Vectors < _fourMomenta->size()? maxNum4Vectors: _fourMomenta->size();
      std::vector<CLHEP::HepLorentzVector*>::iterator first = _fourMomenta->begin() + N;
      std::vector<CLHEP::HepLorentzVector*>::iterator last = _fourMomenta->end();
      for(; first != last; ++first) { DEBUG(0);
      	delete *first; DEBUG(0);
      	*first = 0; DEBUG(0);
      }
      _fourMomenta->erase(_fourMomenta->begin() + N, _fourMomenta->end());
    }
    _threeMomenta = threeMomenta(_fourMomenta);
    _unitVectors = normalise(_threeMomenta);
    
    _sumCr.resize(9);
    _sumCi.resize(9);
    _momC.resize(9);
    _sumB.resize(9);
    _sumK.resize(9);
    _sumD.resize(9);
    _sumH.resize(9);
    _sumPi.resize(4);
    
    if(_unitVectors) { DEBUG(0);
      std::fill(_sumCr.begin(), _sumCr.end(), 0.0);
      std::fill(_sumCi.begin(), _sumCi.end(), 0.0);
      std::fill(_momC.begin(), _momC.end(), 0.0);
      std::fill(_sumB.begin(), _sumB.end(), 0.0);
      std::fill(_sumK.begin(), _sumK.end(), 0.0);
      std::fill(_sumD.begin(), _sumD.end(), 0.0);
      std::fill(_sumH.begin(), _sumH.end(), 0.0);
      std::fill(_sumPi.begin(), _sumPi.end(), 0.0);
      _sumIJ = 0.0;
      _E = jetSumE(_fourMomenta);
      _xiPlus = Divide(xiPlus(_fourMomenta), _E);
      _xiMinus = Divide(xiMinus(_fourMomenta), _E);
      _xPlus = Divide(xPlus(_fourMomenta), _E);
      _xMinus = Divide(xMinus(_fourMomenta), _E);
      compute();
    }
    else { DEBUG(0);
      std::fill(_sumCr.begin(), _sumCr.end(), NaN);
      std::fill(_sumCi.begin(), _sumCi.end(), NaN);
      std::fill(_momC.begin(), _momC.end(), NaN);
      std::fill(_sumB.begin(), _sumB.end(), NaN);
      std::fill(_sumK.begin(), _sumK.end(), NaN);
      std::fill(_sumD.begin(), _sumD.end(), NaN);
      std::fill(_sumH.begin(), _sumH.end(), NaN);
      std::fill(_sumPi.begin(), _sumPi.end(), NaN);
      _sumIJ = NaN;
      _E = NaN;
      _xiPlus = NaN;
      _xiMinus = NaN;
      _xPlus = NaN;
      _xMinus = NaN;
    }
  }
  
  ~FoxWolfram() { DEBUG(0);
    if(_fourMomenta) { DEBUG(0);
      for(unsigned int i = 0; i != _fourMomenta->size(); i++) { DEBUG(0);
      	delete _fourMomenta->at(i); DEBUG(0);
      	_fourMomenta->at(i) = 0; DEBUG(0);
      }
      _fourMomenta->clear(); DEBUG(0);
    }
    delete _fourMomenta; 
    delete _threeMomenta; 
    delete _unitVectors;
  }

  double/* jStr_fwmXX_*/Psi1() const { return Psi(); }

  double/* jStr_fwmXX_*/B0() const { return B(0); }
  double/* jStr_fwmXX_*/B1() const { return B(1); }
  double/* jStr_fwmXX_*/B2() const { return B(2); }
  double/* jStr_fwmXX_*/B3() const { return B(3); }
  double/* jStr_fwmXX_*/B4() const { return B(4); }
  double/* jStr_fwmXX_*/B5() const { return B(5); }
  double/* jStr_fwmXX_*/B6() const { return B(6); }
  double/* jStr_fwmXX_*/B7() const { return B(7); }
  double/* jStr_fwmXX_*/B8() const { return B(8); }

  double/* jStr_fwmXX_*/C0() const { return C(0); }
  double/* jStr_fwmXX_*/C1() const { return C(1); }
  double/* jStr_fwmXX_*/C2() const { return C(2); }
  double/* jStr_fwmXX_*/C3() const { return C(3); }
  double/* jStr_fwmXX_*/C4() const { return C(4); }
  double/* jStr_fwmXX_*/C5() const { return C(5); }
  double/* jStr_fwmXX_*/C6() const { return C(6); }
  double/* jStr_fwmXX_*/C7() const { return C(7); }
  double/* jStr_fwmXX_*/C8() const { return C(8); }

  double/* jStr_fwmXX_*/K0() const { return K(0); }
  double/* jStr_fwmXX_*/K1() const { return K(1); }
  double/* jStr_fwmXX_*/K2() const { return K(2); }
  double/* jStr_fwmXX_*/K3() const { return K(3); }
  double/* jStr_fwmXX_*/K4() const { return K(4); }
  double/* jStr_fwmXX_*/K5() const { return K(5); }
  double/* jStr_fwmXX_*/K6() const { return K(6); }
  double/* jStr_fwmXX_*/K7() const { return K(7); }
  double/* jStr_fwmXX_*/K8() const { return K(8); }

  double/* jStr_fwmXX_*/D0() const { return D(0); }
  double/* jStr_fwmXX_*/D1() const { return D(1); }
  double/* jStr_fwmXX_*/D2() const { return D(2); }
  double/* jStr_fwmXX_*/D3() const { return D(3); }
  double/* jStr_fwmXX_*/D4() const { return D(4); }
  double/* jStr_fwmXX_*/D5() const { return D(5); }
  double/* jStr_fwmXX_*/D6() const { return D(6); }
  double/* jStr_fwmXX_*/D7() const { return D(7); }
  double/* jStr_fwmXX_*/D8() const { return D(8); }

  double/* jStr_fwmXX_*/H0() const { return H(0); }
  double/* jStr_fwmXX_*/H1() const { return H(1); }
  double/* jStr_fwmXX_*/H2() const { return H(2); }
  double/* jStr_fwmXX_*/H3() const { return H(3); }
  double/* jStr_fwmXX_*/H4() const { return H(4); }
  double/* jStr_fwmXX_*/H5() const { return H(5); }
  double/* jStr_fwmXX_*/H6() const { return H(6); }
  double/* jStr_fwmXX_*/H7() const { return H(7); }
  double/* jStr_fwmXX_*/H8() const { return H(8); }

  double/* jStr_fwmXX_*/Q0() const { return Q(0); }
  double/* jStr_fwmXX_*/Q1() const { return Q(1); }
  double/* jStr_fwmXX_*/Q2() const { return Q(2); }
  double/* jStr_fwmXX_*/Q3() const { return Q(3); }
  double/* jStr_fwmXX_*/Q4() const { return Q(4); }
  double/* jStr_fwmXX_*/Q5() const { return Q(5); }
  double/* jStr_fwmXX_*/Q6() const { return Q(6); }
  double/* jStr_fwmXX_*/Q7() const { return Q(7); }
  double/* jStr_fwmXX_*/Q8() const { return Q(8); }

  double/* jStr_fwmXX_*/Pi1() const { return Pi(1); }
  double/* jStr_fwmXX_*/Pi2() const { return Pi(2); }
  double/* jStr_fwmXX_*/Pi3() const { return Pi(3); }
  double/* jStr_fwmXX_*/Pi4() const { return Pi(4); }

  double/* jStr_fwmXX_*/B10() const { return B(1, 0); }
  double/* jStr_fwmXX_*/B20() const { return B(2, 0); }
  double/* jStr_fwmXX_*/B30() const { return B(3, 0); }
  double/* jStr_fwmXX_*/B40() const { return B(4, 0); }
  double/* jStr_fwmXX_*/B50() const { return B(5, 0); }
  double/* jStr_fwmXX_*/B60() const { return B(6, 0); }
  double/* jStr_fwmXX_*/B70() const { return B(7, 0); }
  double/* jStr_fwmXX_*/B80() const { return B(8, 0); }

  double/* jStr_fwmXX_*/C10() const { return C(1, 0); }
  double/* jStr_fwmXX_*/C20() const { return C(2, 0); }
  double/* jStr_fwmXX_*/C30() const { return C(3, 0); }
  double/* jStr_fwmXX_*/C40() const { return C(4, 0); }
  double/* jStr_fwmXX_*/C50() const { return C(5, 0); }
  double/* jStr_fwmXX_*/C60() const { return C(6, 0); }
  double/* jStr_fwmXX_*/C70() const { return C(7, 0); }
  double/* jStr_fwmXX_*/C80() const { return C(8, 0); }

  double/* jStr_fwmXX_*/K10() const { return K(1, 0); }
  double/* jStr_fwmXX_*/K20() const { return K(2, 0); }
  double/* jStr_fwmXX_*/K30() const { return K(3, 0); }
  double/* jStr_fwmXX_*/K40() const { return K(4, 0); }
  double/* jStr_fwmXX_*/K50() const { return K(5, 0); }
  double/* jStr_fwmXX_*/K60() const { return K(6, 0); }
  double/* jStr_fwmXX_*/K70() const { return K(7, 0); }
  double/* jStr_fwmXX_*/K80() const { return K(8, 0); }

  double/* jStr_fwmXX_*/D10() const { return D(1, 0); }
  double/* jStr_fwmXX_*/D20() const { return D(2, 0); }
  double/* jStr_fwmXX_*/D30() const { return D(3, 0); }
  double/* jStr_fwmXX_*/D40() const { return D(4, 0); }
  double/* jStr_fwmXX_*/D50() const { return D(5, 0); }
  double/* jStr_fwmXX_*/D60() const { return D(6, 0); }
  double/* jStr_fwmXX_*/D70() const { return D(7, 0); }
  double/* jStr_fwmXX_*/D80() const { return D(8, 0); }

  double/* jStr_fwmXX_*/H10() const { return H(1, 0); }
  double/* jStr_fwmXX_*/H20() const { return H(2, 0); }
  double/* jStr_fwmXX_*/H30() const { return H(3, 0); }
  double/* jStr_fwmXX_*/H40() const { return H(4, 0); }
  double/* jStr_fwmXX_*/H50() const { return H(5, 0); }
  double/* jStr_fwmXX_*/H60() const { return H(6, 0); }
  double/* jStr_fwmXX_*/H70() const { return H(7, 0); }
  double/* jStr_fwmXX_*/H80() const { return H(8, 0); }

  double/* jStr_fwmXX_*/Q10() const { return Q(1, 0); }
  double/* jStr_fwmXX_*/Q20() const { return Q(2, 0); }
  double/* jStr_fwmXX_*/Q30() const { return Q(3, 0); }
  double/* jStr_fwmXX_*/Q40() const { return Q(4, 0); }
  double/* jStr_fwmXX_*/Q50() const { return Q(5, 0); }
  double/* jStr_fwmXX_*/Q60() const { return Q(6, 0); }
  double/* jStr_fwmXX_*/Q70() const { return Q(7, 0); }
  double/* jStr_fwmXX_*/Q80() const { return Q(8, 0); }

  // Idea borrowed from:
  // https://cdsweb.cern.ch/record/1421606/files/FWD-10-004-pas.pdf
  // Not a Fox-Wolfram moment; just borrowing some functionality as a short-cut
  // http://www-h1.desy.de/psfiles/papers/desy11-166.pdf
  double/* jStr_fwmXX_*/xiPlus() const { return _xiPlus; }
  double/* jStr_fwmXX_*/xiMinus() const { return _xiMinus; }

  // https://indico.cern.ch/getFile.py/access?contribId=2&resId=0&materialId=slides&confId=176612
  // Piggybacking on already available functionality; not FW moment
  double/* jStr_fwmXX_*/xPlus() const { return _xPlus; }
  double/* jStr_fwmXX_*/xMinus() const { return _xMinus; }

private:
  double B(unsigned int N) const { DEBUG(Divide(_sumB[N], _E)); return Divide(_sumB[N], _E); }
  double C(unsigned int N) const { 
    double Cr = Divide(_sumCr[N], _E); 
    double Ci = Divide(_sumCi[N], _E); 
    return (Cr * Cr) + (Ci * Ci);
  }
  double K(unsigned int N) const { DEBUG(Divide(_sumK[N], _E * _E));  return Divide(_sumK[N], _E * _E); }
  double D(unsigned int N) const { DEBUG(Divide(_sumD[N], _E)); return Divide(_sumD[N], _E); }
  double H(unsigned int N) const { DEBUG(Divide(_sumH[N], _E * _E)); return Divide(_sumH[N], _E * _E); }
  double Pi(unsigned int N) const { DEBUG(Divide(_sumPi[N - 1], _E * _E * _E)); return Divide(_sumPi[N - 1], _E * _E * _E); }
  double Psi() const { DEBUG(Divide(_sumPsi, _E * _E * _E)); return Divide(_sumPsi, _E * _E * _E); }

  // http://hep1.phys.ntu.edu.tw/~mwang/thesis_eric.pdf
  // http://prd.aps.org/pdf/PRD/v70/i1/e012001
  double Q(unsigned int N) const { return Divide(_sumH[N], _sumIJ); }

  double B(unsigned int l, unsigned int m) const { return Divide(B(l), B(m)); }
  double C(unsigned int l, unsigned int m) const { return Divide(C(l), C(m)); }
  double K(unsigned int l, unsigned int m) const { return Divide(K(l), K(m)); }
  double D(unsigned int l, unsigned int m) const { return Divide(D(l), D(m)); }
  double H(unsigned int l, unsigned int m) const { return Divide(H(l), H(m)); }
  double Q(unsigned int l, unsigned int m) const { return Divide(Q(l), Q(m)); }
  double Pi(unsigned int l, unsigned int m) const { return Divide(Pi(l), Pi(m)); }

  double cosTheta(CLHEP::Hep3Vector& a, CLHEP::Hep3Vector& b) {
    return a.dot(b);
  }

  CLHEP::HepLorentzVector& boost(CLHEP::HepLorentzVector hlv, const Jet*& jet) { DEBUG(parentJet);
    const Jet* parentJet = jet;
    if(parentJet != 0 && (parentJet->hlv().restMass2() <= m2min || parentJet->hlv().findBoostToCM().mag2() >= b2max)) { // !ZMxpvTachyonic
      hlv.setX(NaN);
      hlv.setY(NaN);
      hlv.setZ(NaN);
      hlv.setT(NaN);
      parentJet = 0; // signal that we are to not boost, pass back NaN HLV instead
      DEBUG(hlv);
    }
    return !parentJet? hlv: hlv.boost(parentJet->hlv().findBoostToCM());
  }

  double stp(CLHEP::Hep3Vector& i, CLHEP::Hep3Vector& j, CLHEP::Hep3Vector& k) const {
    return (i.cross(j)).dot(k);
  }

  // http://www.sciencedirect.com/science/article/pii/0370269379904441
  // Physics Letters B
  // Volume 82, Issue 1, 12 March 1979, Pages 134-138
  // doi:10.1016/0370-2693(79)90444-1 | How to Cite or Link Using DOI
  //   Permissions & Reprints

  // Tests for planar events in e+e− annihilation☆
  // also: http://ccdb4fs.kek.jp/cgi-bin/img/allpdf?197811220
  bool compute() { DEBUG(0);
    // fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    // feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    for(unsigned int iIndex = 0; iIndex != _threeMomenta->size(); iIndex++) {
      CLHEP::Hep3Vector i = _threeMomenta->at(iIndex); // CHECK_FPE;
      double p_i = i.mag(); // CHECK_FPE;
      double z = cosTheta(_unitVectors->at(iIndex), _z); // CHECK_FPE;

      _sumB[0] += p_i * P<0u>(z); // CHECK_FPE;
      _sumB[1] += p_i * P<1u>(z); // CHECK_FPE;
      _sumB[2] += p_i * P<2u>(z); // CHECK_FPE;
      _sumB[3] += p_i * P<3u>(z); // CHECK_FPE;
      _sumB[4] += p_i * P<4u>(z); // CHECK_FPE;
      _sumB[5] += p_i * P<5u>(z); // CHECK_FPE;
      _sumB[6] += p_i * P<6u>(z); // CHECK_FPE;
      _sumB[7] += p_i * P<7u>(z); // CHECK_FPE;
      _sumB[8] += p_i * P<8u>(z); // CHECK_FPE;

      double pt_i = i.perp(); // CHECK_FPE;
      double phi_i = PhiCorr()(i.phi()); // CHECK_FPE;

      double sn = 0.; // CHECK_FPE;
      double cs = 0.; // CHECK_FPE;

      _sumCr[0] += pt_i; _sumCi[0] += 0.; // CHECK_FPE;

      sincos(1. * phi_i, &sn, &cs); _sumCr[1] += pt_i * cs; _sumCi[1] += pt_i * sn; // CHECK_FPE; 
      sincos(2. * phi_i, &sn, &cs); _sumCr[2] += pt_i * cs; _sumCi[2] += pt_i * sn; // CHECK_FPE; 
      sincos(3. * phi_i, &sn, &cs); _sumCr[3] += pt_i * cs; _sumCi[3] += pt_i * sn; // CHECK_FPE; 
      sincos(4. * phi_i, &sn, &cs); _sumCr[4] += pt_i * cs; _sumCi[4] += pt_i * sn; // CHECK_FPE; 
      sincos(5. * phi_i, &sn, &cs); _sumCr[5] += pt_i * cs; _sumCi[5] += pt_i * sn; // CHECK_FPE; 
      sincos(6. * phi_i, &sn, &cs); _sumCr[6] += pt_i * cs; _sumCi[6] += pt_i * sn; // CHECK_FPE; 
      sincos(7. * phi_i, &sn, &cs); _sumCr[7] += pt_i * cs; _sumCi[7] += pt_i * sn; // CHECK_FPE; 
      sincos(8. * phi_i, &sn, &cs); _sumCr[8] += pt_i * cs; _sumCi[8] += pt_i * sn; // CHECK_FPE; 

      _sumD[0] += pt_i; // CHECK_FPE;
      _sumD[1] += pt_i * cos(phi_i); // CHECK_FPE;
      _sumD[2] += pt_i * cos(2. * phi_i); // CHECK_FPE;
      _sumD[3] += pt_i * cos(3. * phi_i); // CHECK_FPE;
      _sumD[4] += pt_i * cos(4. * phi_i); // CHECK_FPE;
      _sumD[5] += pt_i * cos(5. * phi_i); // CHECK_FPE;
      _sumD[6] += pt_i * cos(6. * phi_i); // CHECK_FPE;
      _sumD[7] += pt_i * cos(7. * phi_i); // CHECK_FPE;
      _sumD[8] += pt_i * cos(8. * phi_i); // CHECK_FPE;

      for(unsigned int jIndex = 0; jIndex != _threeMomenta->size(); jIndex++) {
	CLHEP::Hep3Vector j = _threeMomenta->at(jIndex); // CHECK_FPE;
	double p_j = j.mag(); // CHECK_FPE;
	double p_ij = p_i * p_j; // CHECK_FPE;
	double pt_j = j.perp(); // CHECK_FPE;
	double pt_ij = pt_i * pt_j; // CHECK_FPE;
	double phi_j = PhiCorr()(j.phi()); // CHECK_FPE;
	double deltaPhi = DeltaPhi()(phi_i, phi_j); // CHECK_FPE;
	double x = cosTheta(_unitVectors->at(iIndex), _unitVectors->at(jIndex)); // CHECK_FPE;
	
	_sumK[0] += pt_ij; // CHECK_FPE;
	_sumK[1] += pt_ij * cos(deltaPhi); // CHECK_FPE;
	_sumK[2] += pt_ij * cos(2. * deltaPhi); // CHECK_FPE;
	_sumK[3] += pt_ij * cos(3. * deltaPhi); // CHECK_FPE;
	_sumK[4] += pt_ij * cos(4. * deltaPhi); // CHECK_FPE;
	_sumK[5] += pt_ij * cos(5. * deltaPhi); // CHECK_FPE;
	_sumK[6] += pt_ij * cos(6. * deltaPhi); // CHECK_FPE;
	_sumK[7] += pt_ij * cos(7. * deltaPhi); // CHECK_FPE;
	_sumK[8] += pt_ij * cos(8. * deltaPhi); // CHECK_FPE;

	_sumH[0] += p_ij * P<0u>(x); // CHECK_FPE;
	_sumH[1] += p_ij * P<1u>(x); // CHECK_FPE;
	_sumH[2] += p_ij * P<2u>(x); // CHECK_FPE;
	_sumH[3] += p_ij * P<3u>(x); // CHECK_FPE;
	_sumH[4] += p_ij * P<4u>(x); // CHECK_FPE;
	_sumH[5] += p_ij * P<5u>(x); // CHECK_FPE;
	_sumH[6] += p_ij * P<6u>(x); // CHECK_FPE;
	_sumH[7] += p_ij * P<7u>(x); // CHECK_FPE;
	_sumH[8] += p_ij * P<8u>(x); // CHECK_FPE;
	
	_sumIJ += p_ij; // CHECK_FPE;

	for(unsigned int kIndex = 0; kIndex != _threeMomenta->size(); kIndex++) {
	  CLHEP::Hep3Vector k = _threeMomenta->at(kIndex); // CHECK_FPE;
	  double p_k = k.mag(); // CHECK_FPE;
	  double p_ijk = p_ij * p_k; // CHECK_FPE;
	  double ijk = stp(_unitVectors->at(iIndex), _unitVectors->at(jIndex), _unitVectors->at(kIndex)); // CHECK_FPE;
	  double p3ijk = p_ijk * ijk; // CHECK_FPE;
	  _sumPsi += p3ijk * A(_unitVectors->at(iIndex), _unitVectors->at(jIndex), _unitVectors->at(kIndex)); // CHECK_FPE;
	  double p3ijk2 = p3ijk * ijk; // CHECK_FPE;
	  
	  _sumPi[0] += p3ijk2 * S<1u>(_unitVectors->at(iIndex), _unitVectors->at(jIndex), _unitVectors->at(kIndex)); // CHECK_FPE;
	  _sumPi[1] += p3ijk2 * S<2u>(_unitVectors->at(iIndex), _unitVectors->at(jIndex), _unitVectors->at(kIndex)); // CHECK_FPE;
	  _sumPi[2] += p3ijk2 * S<3u>(_unitVectors->at(iIndex), _unitVectors->at(jIndex), _unitVectors->at(kIndex)); // CHECK_FPE;
	  _sumPi[3] += p3ijk2 * S<4u>(_unitVectors->at(iIndex), _unitVectors->at(jIndex), _unitVectors->at(kIndex)); // CHECK_FPE;
	} // k loop
      } // j loop
    } // i loop
    // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    // feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    return true;
  }

  struct SortByMomentumDown {
    bool operator()(CLHEP::HepLorentzVector* a, CLHEP::HepLorentzVector* b) { return a->v().mag() > b->v().mag(); }
    bool operator()(const CLHEP::Hep3Vector& a, const CLHEP::Hep3Vector& b) { return a.mag() > b.mag(); }
  } sortByMomentumDown;

  std::vector<CLHEP::HepLorentzVector*>* fourMomenta(const JetCollection* jetCollection, const Jet* parentJet = 0) { DEBUG(0);
    std::vector<CLHEP::HepLorentzVector*>* fourMomenta = 0;
    if(jetCollection) { DEBUG(0);
      JetCollectionHelper::jetcollection_t theJets(jetCollection->begin(), jetCollection->end()); DEBUG(0);
      JetCollectionHelper::jetcollection_t::iterator firstJet = theJets.begin(); DEBUG(0);
      JetCollectionHelper::jetcollection_t::iterator lastJet  = theJets.end(); DEBUG(0);
      fourMomenta = new std::vector<CLHEP::HepLorentzVector*>;
      fourMomenta->clear(); DEBUG(0);
      for(; firstJet != lastJet; ++firstJet) { DEBUG(0);
	if((*firstJet) == 0) continue; DEBUG(0);
	const Jet* aJet = (*firstJet)->clone(true, true); DEBUG(0); // Implicit new here
	if(aJet == 0) continue; DEBUG(0);
	fourMomenta->push_back(new CLHEP::HepLorentzVector(boost(aJet->hlv(), parentJet))); DEBUG(0);
	delete aJet; // We're done with this jet; do everything onwards with HLVs
      }
    }
    return fourMomenta;
  }
  
  std::vector<CLHEP::HepLorentzVector*>* fourMomenta(const Jet* jet) { DEBUG(0);
    std::vector<CLHEP::HepLorentzVector*>* fourMomenta = 0;
    if(jet) { DEBUG(0);
      fourMomenta = new std::vector<CLHEP::HepLorentzVector*>;
      fourMomenta->clear();
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(jet); DEBUG(0);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(jet); DEBUG(0);
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	fourMomenta->push_back(new CLHEP::HepLorentzVector(boost(firstConstituent.hlv(), jet))); DEBUG(0);	
      }
    }
    return fourMomenta;
  }
  
  std::vector<CLHEP::Hep3Vector>* threeMomenta(std::vector<CLHEP::HepLorentzVector*>* fourMomenta) { DEBUG(0);
    std::vector<CLHEP::Hep3Vector>* threeMomenta = 0;
    if(fourMomenta && fourMomenta->size()) { DEBUG(0);
      threeMomenta = new std::vector<CLHEP::Hep3Vector>(fourMomenta->size());
      threeMomenta->clear();
      for(unsigned int i = 0; i != fourMomenta->size(); i++) {
	if(fourMomenta->at(i) == 0) continue; 
	threeMomenta->push_back(fourMomenta->at(i)->v());
      }
    }
    return threeMomenta;
  }

  std::vector<CLHEP::Hep3Vector>* normalise(std::vector<CLHEP::Hep3Vector>* threeMomenta) { DEBUG(0);
    std::vector<CLHEP::Hep3Vector>* unitVectors = 0;
    if(threeMomenta && threeMomenta->size()) { DEBUG(0);
      unitVectors = new std::vector<CLHEP::Hep3Vector>(threeMomenta->size());
      unitVectors->clear();
      for(unsigned int i = 0; i != threeMomenta->size(); i++) {
	unitVectors->push_back(threeMomenta->at(i).unit());
      }
    }
    return unitVectors;
  }

  // Adder functor
  struct AddE: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddE(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += x->e(); }
  };

  double jetSumE(std::vector<CLHEP::HepLorentzVector*>* fourMomenta) const { DEBUG(0);
    if(fourMomenta == 0 || fourMomenta->empty()) return NaN;
    return std::for_each(fourMomenta->begin(), fourMomenta->end(), AddE()).sum;
  }
  
  struct AddXiPlusTerm: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddXiPlusTerm(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += (x->e() + x->pz()); }
  };
  
  struct AddXiMinusTerm: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddXiMinusTerm(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += (x->e() - x->pz()); }
  };
  
  double xiPlus(std::vector<CLHEP::HepLorentzVector*>* fourMomenta) const { DEBUG(0);
    if(fourMomenta == 0 || fourMomenta->empty()) return NaN;
    return std::for_each(fourMomenta->begin(), fourMomenta->end(), AddXiPlusTerm()).sum;
  }

  double xiMinus(std::vector<CLHEP::HepLorentzVector*>* fourMomenta) const { DEBUG(0);
    if(fourMomenta == 0 || fourMomenta->empty()) return NaN;
    return std::for_each(fourMomenta->begin(), fourMomenta->end(), AddXiMinusTerm()).sum;
  }

  struct AddXPlusTerm: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddXPlusTerm(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += (x->perp() * exp(fabs(x->pseudoRapidity()))); }
  };
  
  struct AddXMinusTerm: public std::unary_function<CLHEP::HepLorentzVector*, void> {
    AddXMinusTerm(): sum(0) {}
    double sum;
    result_type operator()(argument_type x) { sum += (x->perp() * exp(-fabs(x->pseudoRapidity()))); }
  };
  
  double xPlus(std::vector<CLHEP::HepLorentzVector*>* fourMomenta) const { DEBUG(0);
    if(fourMomenta == 0 || fourMomenta->empty()) return NaN;
    return std::for_each(fourMomenta->begin(), fourMomenta->end(), AddXPlusTerm()).sum;
  }

  double xMinus(std::vector<CLHEP::HepLorentzVector*>* fourMomenta) const { DEBUG(0);
    if(fourMomenta == 0 || fourMomenta->empty()) return NaN;
    return std::for_each(fourMomenta->begin(), fourMomenta->end(), AddXMinusTerm()).sum;
  }

  CLHEP::Hep3Vector _z;
  double _E;
  double _xiPlus, _xiMinus;
  double _xPlus, _xMinus;
  std::vector<CLHEP::HepLorentzVector*>* _fourMomenta;
  std::vector<CLHEP::Hep3Vector>* _threeMomenta;
  std::vector<CLHEP::Hep3Vector>* _unitVectors;
  std::vector<double> _sumCr;
  std::vector<double> _sumCi;
  std::vector<double> _momC;
  std::vector<double> _sumB;
  std::vector<double> _sumK;
  std::vector<double> _sumD;
  std::vector<double> _sumH;
  std::vector<double> _sumPi;
  double _sumIJ;
  double _sumPsi;
};

#endif
