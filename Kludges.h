#ifndef KLUDGES_H
#define KLUDGES_H

#include "JetEvent/Jet.h"
#include "JetEvent/JetAssociationGeneric.h" // To get subjets

// // debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

// #include <tr1/memory> // Last jet turns the lights out, using shared_ptr
// #define TAKE_OWNERSHIP(x)						\
//   std::tr1::shared_ptr<const Jet> __##x = std::tr1::shared_ptr<const Jet>(x, Deleter());

class Deleter { // Garbage collection on behalf of leaky CoreSubJetsTool in JetMomentTools
public:
  void operator()(const JetCollection*& jetColl) const {
    for(unsigned int i = 0; i != jetColl->size(); i++) { DEBUG(i);
      const JetAssociationGeneric<JetCollection>* association = 
	jetColl->at(i)->getAssociation<JetAssociationGeneric<JetCollection> >("CoreSubJets"); DEBUG(association);
      if(association) {
	JetCollection* subjetColl = association->get(); DEBUG(subjetColl); DEBUG(subjetColl->size());
	if(subjetColl) subjetColl->clear(); DEBUG(subjetColl->size());
	delete subjetColl; DEBUG(subjetColl);
      }
    }
  }
};

#endif // KLUDGES_H

