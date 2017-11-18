#ifndef NUMERIC_H
#define NUMERIC_H

#include <fenv.h> // For FPE trapping <cfenv>
#include <iostream> // So FPE trapping can report to stdout
#include <cmath> // For FP classification
#include <sstream> // For string --> number
#include <string> // For string --> number
#include <limits> // For NaNs
#include "CLHEP/Units/PhysicalConstants.h" // For pi

using namespace std;

volatile double __volatiledouble__;

namespace Numeric { // NaNs, infinities, free functions with argument bounds checking...
  const double NaN = std::numeric_limits<double>::quiet_NaN();
  const double Infinity = NaN;// std::numeric_limits<double>::infinity();
  const double b2max = 1.;
  const double m2min = 0.;

  inline std::string FPclassify(double x) {
    int i = std::fpclassify(x);
    std::string s;
    switch (i) {
    case FP_INFINITE: s = "FP_INFINITE"; break;
    case FP_NAN: s = "FP_NAN"; break;
    case FP_NORMAL: s = "FP_NORMAL"; break;
    case FP_SUBNORMAL: s = "FP_SUBNORMAL"; break;
    case FP_ZERO: s = "FP_ZERO"; break;
    default: s = "UNDEFINED"; break;
    }
    return s;
  }

  inline std::string FPcategory(double x) {
    int ia = std::isfinite(x); std::string sa = ia? "isfinite ": "";
    int ib = std::isinf(x);    std::string sb = ib? "isinf ": "";
    int ic = std::isnan(x);    std::string sc = ic? "isnan ": "";
    int id = std::isnormal(x); std::string sd = id? "isnormal ": "";
    int ie = std::signbit(x);  std::string se = ie? "-": "";
    return sa+sb+sc+sd+se;
  }

  inline int FPEtest(void) {
    int set_excepts = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    if(set_excepts & FE_DIVBYZERO)
      std::cout << "Detected FE_DIVBYZERO" << std::endl;
    if(set_excepts & FE_INVALID)
      std::cout << "Detected FE_INVALID" << std::endl;
    if(set_excepts & FE_OVERFLOW)
      std::cout << "Detected FE_OVERFLOW" << std::endl;
    return set_excepts;
  }
  
  inline std::string FPEtype(int i) {
    std::string s;
    switch (i) {
    case NULL: s = ""; break;
    case FE_DIVBYZERO: s = "FE_DIVBYZERO"; break;
    case FE_INEXACT: s = "FE_INEXACT"; break;
    case FE_INVALID: s = "FE_INVALID"; break;
    case FE_OVERFLOW: s = "FE_OVERFLOW"; break;
    case FE_UNDERFLOW: s = "FE_UNDERFLOW"; break;
    default: s = "Undefined"; break;
    }
    return s;
  }

  template <typename numbertype>
  inline bool isEqual(numbertype const& a, numbertype const& b) {
    // See http://www.parashift.com/c++-faq-lite/newbie.html#faq-29.17,
    return std::abs(a - b) <= std::numeric_limits<numbertype>::epsilon() * std::abs(a);
    // see Knuth section 4.2.2 pages 217-218
  } 

  template<typename numbertype>
  inline bool hasNaN(const numbertype& /*N*/) {
    return numeric_limits<numbertype>::has_quiet_NaN;
  }
  
  template<typename numbertype>
  inline bool isInteger(const numbertype& /*N*/) {
    return numeric_limits<numbertype>::is_integer;
  }
  
  template<typename numbertype>
  inline bool hasInfinity(const numbertype& /*N*/) {
    return numeric_limits<numbertype>::has_infinity;
  }

  inline double sgn(const double& N) {
    return signbit(N)? -1.: 1.;
  }

  inline double Sqrt(double x) {
    return x >= 0. && !isinf(x)? sqrt(x): NaN;
  }
  
  inline double Ln(double x) {
    return x > 0. && !isinf(x)? log(x): 
      x == 0.? -Infinity: NaN;
  }
  
  inline double Arccos(double x) {
    return fabs(x) <= 1.? acos(x): NaN;
  }
  
  inline double Arcsin(double x) {
    return fabs(x) <= 1.? asin(x): NaN;
  }
  
  inline double Arctan2(double y, double x) {
    return std::fpclassify(y) != FP_ZERO && std::fpclassify(x) != FP_ZERO? atan2(y, x): NaN;
  }

  inline double Divide(double x, double y) {
    return (
  	    std::fpclassify(x) == FP_INFINITE ||
  	    std::fpclassify(x) == FP_NAN ||
  	    std::fpclassify(y) == FP_ZERO ||
  	    std::fpclassify(y) == FP_INFINITE ||
  	    std::fpclassify(y) == FP_NAN	    
  	    )? NaN: x / y;
  }
  
  // inline double Divide(double x, double y) {
  //   double result;
  //   // fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //   // feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  //   result = (
  // 	      std::fpclassify(x) == FP_INFINITE ||
  // 	      std::fpclassify(x) == FP_NAN ||
  // 	      std::fpclassify(y) == FP_ZERO ||
  // 	      std::fpclassify(y) == FP_INFINITE ||
  // 	      std::fpclassify(y) == FP_NAN	    
  // 	      )? NaN: x / y;
  //   __volatiledouble__ = result;
  //   if(FPEtest()) {
  //     std::cout << "in " << __FUNCTION__
  // 		<< ", " << __FILE__ << ":" << __LINE__ << ", x: "
  // 		<< std::setprecision(10) << x << std::setprecision(4)
  // 		<< " " << FPclassify(x) << " / y: "
  // 		<< std::setprecision(10) << y << std::setprecision(4)
  // 		<< " " << FPclassify(y) << std::endl;
  //     result = NaN;
  //   }
  //   __volatiledouble__ = result;
  //   // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  //   // feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
  //   return result;
  // }

  inline double Multiply(double x, double y) {
    return (
	    std::fpclassify(x) == FP_INFINITE ||
	    std::fpclassify(x) == FP_NAN ||
	    std::fpclassify(y) == FP_INFINITE ||
	    std::fpclassify(y) == FP_NAN	    
	    )? NaN: x * y;
  }
  
  inline double Power(double x, double y) {
    return
      fabs(x) == 0. && y < 0.? Infinity:
      fabs(x) == 1. && fabs(y) == Infinity? NaN:
      fabs(x) == 0. && fabs(y) == 0.? NaN:
      fabs(x) == Infinity && fabs(y) == 0.? NaN: pow(x, y);
  }

  inline double Asymmetry(double x, double y) {
    return Divide(x - y, x + y);
  }
  
  inline double pseudoRapidity(double mag, double z) {
    return mag == 0? 0.0:
      mag == z? Infinity:
      mag == -z? -Infinity:
      0.5 * log((mag + z) / (mag - z));
  }

  inline double HeavisideTheta(double x) {
    return x >= 0.? 1.0: 0.0;
  }

  inline double DiracDelta(double x) {
    return x == 0.? Infinity: 0.;
  }

  inline double Kernel(double x, double dR) {
    return exp(-(x * x)/(dR * dR)) / sqrt(CLHEP::pi * dR * dR);
  }

  // Legendre polynomials
  template<unsigned int l> inline double P(double x) { return P<l>(x); }
  template<> inline double P<0u>(double x) { return 1.; }
  template<> inline double P<1u>(double x) { return x; }
  template<> inline double P<2u>(double x) { return 0.5* ((3. * x * x) - 1.); }
  template<> inline double P<3u>(double x) { return 0.5* ((5. * x * x) - (3. * x)); }
  template<> inline double P<4u>(double x) { return 0.125* ((35.* x * x * x * x) - (30. * x * x) + 3.); }
  template<> inline double P<5u>(double x) { return 0.125* ((63.* x * x * x * x * x) - (70. * x * x * x) + (15.* x)); }
  template<> inline double P<6u>(double x) { return 0.0625* ((231.* x * x * x * x * x * x) - (315. * x * x * x * x) + (105.* x * x) - 5.); }
  template<> inline double P<7u>(double x) { return 0.0625* ((429.* x * x * x * x * x * x * x) - (693. * x * x * x * x * x) + (315.* x * x * x) - (35. * x)); }
  template<> inline double P<8u>(double x) { return 0.0078125* ((6435.* x * x * x * x * x * x * x * x) - (12012. * x * x * x * x * x * x) + (6930.* x * x * x * x) - (1260. * x * x) + 35.); }

  template<typename numbertype>
  inline numbertype str2num(std::string const& s) {
    std::istringstream i(s);
    numbertype x;
    if(!(i >> x)) x = std::numeric_limits<numbertype>::quiet_NaN();
    return x;
  }
  
  template<typename numbertype>
  inline std::string num2str(const numbertype& N) {
    std::stringstream out;
    out << N;
    return out.str();
  }

  struct PhiCorr { // Helper function
    template <typename T>
    T operator()(T phi) {
      if(phi < -CLHEP::pi) phi += CLHEP::twopi;
      if(phi > CLHEP::pi) phi -= CLHEP::twopi;
      return phi;
    }
    template <typename T>
    T operator()(T phi) const {
      if(phi < -CLHEP::pi) phi += CLHEP::twopi;
      if(phi > CLHEP::pi) phi -= CLHEP::twopi;
      return phi;
    }
  };

  class DeltaPhi { // Helper function
  public:
    template <typename T>
    T operator()(T a, T b) {
      return phiCorr(phiCorr(a) - phiCorr(b)); 
    }
    template <typename T>
    T operator()(T a, T b) const {
      return phiCorr(phiCorr(a) - phiCorr(b)); 
    }

  private:
    PhiCorr phiCorr;
  };

  class SigmaPhi { // Helper function
  public:
    template <typename T>
    T operator()(T a, T b) {
      return phiCorr(phiCorr(a) + phiCorr(b)); 
    }
    template <typename T>
    T operator()(T a, T b) const {
      return phiCorr(phiCorr(a) + phiCorr(b)); 
    }

  private:
    PhiCorr phiCorr;
  };

  inline bool isBHadron(int pdg) {
    int mpdg = abs(pdg);
    return (   ( 500 < mpdg && mpdg < 599 )   ||
	       ( 10500 < mpdg && mpdg < 10599 ) ||
	       (  5000 < mpdg && mpdg < 5999  ) ||
	       ( 20500 < mpdg && mpdg < 20599 ) );
  }
  
  inline bool isDHadron(int pdg) {
    int mpdg = abs(pdg);
    return (   ( 400 < mpdg && mpdg < 499 )   || 
	       ( 10400 < mpdg && mpdg < 10499 ) ||
	       (  4000 < mpdg && mpdg < 4999  ) ||
	       ( 20400 < mpdg && mpdg < 20499 ) );
  }
  
}

#endif
