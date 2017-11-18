#ifndef JETFLAGS_H
#define JETFLAGS_H

#include "JetEvent/Jet.h"
#include "JetEvent/JetCollection.h"
#include "JetTagInfo/TruthInfo.h" // For MC truth tagging only

#include "UserAnalysis/Numeric.h"
#include <string> // To determine jet inputs and algorithm
 
// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

using namespace Numeric;

struct JetAuthor {
  std::string/* jStr_authorJ_*/operator()(const Jet*& j) { return j->jetAuthor(); }
  std::string/* jStr_authorE_*/operator()(const JetCollection*& j) { return j->author(); }
  std::string tag(const Jet*& j) { return j->getCalibTag(1); }
};

struct JetAlg { // PUBLISH_JET_SHAPE, PUBLISH_EVENT_SHAPE
  int found(std::string& author) {
    size_t found;
    bool antikt = false;
    bool cam = false;
    bool kt = false;
    bool cone = false;
    bool siscone = false;
    
    found = author.find("AntiKt");
    if(found != std::string::npos) antikt = true;
    found = author.find("Cam");
    if(found != std::string::npos) cam = true;
    found = author.find("Kt");
    if(found != std::string::npos) kt = true;
    found = author.find("Cone");
    if(found != std::string::npos) cone = true;
    found = author.find("SISCone");
    if(found != std::string::npos) siscone = true;
    
    if(antikt) kt = false;
    if(siscone) cone = false;
    
    int alg = 0;
    if(antikt) alg = 1;
    if(cam) alg = 2;
    if(kt) alg = 3;
    if(cone) alg = 4;
    if(siscone) alg = 5;
    
    return alg;
  }

  int/* jStr_algJ_*/operator()(const Jet*& j) {
    std::string author = j->jetAuthor();
    return found(author);
  }

  int/* jStr_algE_*/operator()(const JetCollection*& j) {
    std::string author = j->author();
    return found(author);
  }
};

struct Input { // PUBLISH_JET_SHAPE, PUBLISH_EVENT_SHAPE
  int found(std::string& author) {
    size_t found;
    bool tower = false;
    bool topo = false;
    bool truth = false;
    bool topotower = false;
    bool lctopo = false;
    
    found = author.find("Tower");
    if(found != std::string::npos) tower = true;
    found = author.find("Topo");
    if(found != std::string::npos) topo = true;
    found = author.find("Truth");
    if(found != std::string::npos) truth = true;
    found = author.find("LCTopo");
    if(found != std::string::npos) lctopo = true;
    
    if(topo && tower) {
      topotower = true;
      topo = false;
      tower = false;
    }
    if(lctopo) topo = false;
    
    int input = 0;
    if(topotower) input = 1;
    if(tower) input = 2;
    if(topo) input = 3;
    if(lctopo) input = 4;
    if(truth) input = 5;
    
    return input;
  }

  int/* jStr_inputJ_*/operator()(const Jet*& j) {
    std::string author = j->jetAuthor();
    return found(author);
  }

  int/* jStr_inputE_*/operator()(const JetCollection*& j) {
    std::string author = j->author();
    return found(author);
  }
};
  
struct Bjet { // PUBLISH_JET_SHAPE, PUBLISH_EVENT_SHAPE
  int found(std::string& author) {
    size_t found;
    int bjet = 0;
    found = author.find("B");
    if(found != std::string::npos) bjet = 1;
    return bjet;
  }

  int/* jStr_bjetJ_*/operator()(const Jet*& j) {
    std::string author = j->jetAuthor();
    return found(author);
  }

  int/* jStr_bjetE_*/operator()(const JetCollection*& j) {
    std::string author = j->author();
    return found(author);
  }
};

// A kludge to allow filtering of specially-flagged jet collections; flag them with 'X'
struct Xjet { // PUBLISH_JET_SHAPE, PUBLISH_EVENT_SHAPE
  int found(std::string& author) {
    size_t found;
    int xjet = 0;
    found = author.find("X");
    if(found != std::string::npos) xjet = 1;
    return xjet;
  }

  int/* jStr_xjetJ_*/operator()(const Jet*& j) {
    std::string author = j->jetAuthor();
    return found(author);
  }

  int/* jStr_xjetE_*/operator()(const JetCollection*& j) {
    std::string author = j->author();
    return found(author);
  }
};

// A kludge for slimming jet collections
struct FilterJets { // PUBLISH_EVENT_SHAPE
  int/* jStr_bjetE_*/operator()(const JetCollection*& j, double etacut = 999., double pTcut = 0.) {
    JetCollection* jetColl = const_cast<JetCollection*>(j);
    for(JetCollection::iterator it = jetColl->begin(); it != jetColl->end();) {
      if((fabs((*it)->eta()) > etacut || (*it)->pt() < pTcut)) {
	it = jetColl->erase(it);
      }
      else {
	++it;
      }
    }
    j = const_cast<JetCollection*>(jetColl);
    return j->size();
  }
};

struct FilterBjets { // PUBLISH_EVENT_SHAPE
  int/* jStr_bjetE_*/operator()(const JetCollection*& j, double wcut) {
    JetCollection* jetColl = const_cast<JetCollection*>(j);
    for(JetCollection::iterator it = jetColl->begin(); it != jetColl->end();) {
      if((*it)->getFlavourTagWeight() <= wcut) {
	it = jetColl->erase(it);
      }
      else {
	++it;
      }
    }
    j = const_cast<JetCollection*>(jetColl);
    return j->size();
  }

  int/* jStr_bjetE_*/operator()(const JetCollection*& j) {
    JetCollection* jetColl = const_cast<JetCollection*>(j);
    for(JetCollection::iterator it = jetColl->begin(); it != jetColl->end();) {
      std::string label("N/A");
      const Analysis::TruthInfo* mcinfo = (*it)->tagInfo<Analysis::TruthInfo>("TruthInfo");
      if(mcinfo) {
  	label = mcinfo->jetTruthLabel();
      }
      if(label != "B") {
  	it = jetColl->erase(it);
      }
      else {
  	++it;
      }
    }
    j = const_cast<JetCollection*>(jetColl);
    return j->size();
  }
};

struct MyJets { // PUBLISH_JET_SHAPE, PUBLISH_EVENT_SHAPE
  int found(std::string& author) {
    size_t found;
    int myJets = 0;
    found = author.find("My");
    if(found != std::string::npos) myJets = 1;
    return myJets;
  }

  int/* jStr_myJetsJ_*/operator()(const Jet*& j) {
    std::string author = j->jetAuthor();
    return found(author);
  }

  int/* jStr_myJetsE_*/operator()(const JetCollection*& j) {
    std::string author = j->author();
    return found(author);
  }
};

struct RadialParameter { // PUBLISH_JET_SHAPE, PUBLISH_EVENT_SHAPE
  int found(std::string& author) {
    size_t found;
    std::string num = "";
    found = author.find_first_of("0123456789");
    while(found != std::string::npos) {
      char c = author[found];
      num.append(1, c);
      found = author.find_first_of("0123456789", found + 1);
    }
    return str2num<unsigned int>(num.c_str());
  }

  int/* jStr_radialParamJ_*/operator()(const Jet*& j) {
    std::string author = j->jetAuthor();
    return found(author);
  }
  
  int/* jStr_radialParamE_*/operator()(const JetCollection*& j) {
    std::string author = j->author();
    return found(author);
  }
};

#endif

