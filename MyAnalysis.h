#ifndef MY_ANALYSIS_H
#define MY_ANALYSIS_H

#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ObjectVector.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "StoreGate/StoreGateSvc.h"
#include "GaudiKernel/ITHistSvc.h"
#include "AthenaBaseComps/AthAlgorithm.h"
#include "JetMomentTools/JetVertexAssociationTool.h"

#include <vector>
#include <string>

class VxContainer;
class JetCollection;
class TruthParticleContainer;

class MyAnalysis: public AthAlgorithm {

public:

  MyAnalysis(const std::string& name, ISvcLocator* pSvcLocator);
  ~MyAnalysis();

  virtual StatusCode beginRun();
  virtual StatusCode initialize();
  virtual StatusCode finalize();
  virtual StatusCode execute();
  virtual StatusCode initEvent();

private:

  // Handle  to JetVertexAssociationTool.
  ToolHandle<JetVertexAssociationTool> m_jvfTool;

  std::string m_jetcollname;
  std::vector<std::string> m_jetcollnames;
  double m_cutIP3DSV1;
  std::string m_vxCandidatesName;

  ITHistSvc* m_thistSvc;

  /** a handle on Store Gate for access to the Event Store */
  StoreGateSvc* m_storeGate;

  // to add event info to new ntuple (used to go by default in CollectionTree)
  StatusCode addEventInfo(); 

  /** get quark flavour of jets */
  int getQuarkJetFlavour(JetCollection::const_iterator jetItr);

  /** Tag jets */
  long long int QGtag(const TruthParticleContainer* mcpartTES, const Jet* j);
  long long int QGtag2(const TruthParticleContainer* mcpartTES, const Jet* j);

  // /** name of the AOD truth particle container to retrieve from StoreGate */
  std::string m_truthParticleContainerName;

  /** Athena-Aware Ntuple (AAN) variables - branches of the AAN TTree */

  // stuff for new ntuple
  // The standard AANT, CollectionTree, is bare bones
  TTree* m_tree_AS; 

  unsigned int m_eventNr;

  unsigned int    m_runNumber;
  unsigned int    m_eventNumber;
  unsigned int    m_eventTime;
  unsigned int    m_lumiBlock;
  unsigned int    m_bCID;
  double  m_eventWeight;
  
  unsigned long int _nparticles;
  unsigned long int _NFPE;

  int _jStr_numPVE;
  double _jStr_PVxE;
  double _jStr_PVyE;
  double _jStr_PVzE;
  double _jStr_PVrE;
  double _jStr_PVchiSqrE;
  double _jStr_PVnDoFE;
  double _jStr_PVfitE;
  int _jStr_PVnTrkE;
  int _jStr_PVtypeE;

  std::vector<int>* _jStr_numPVJ;
  std::vector<double>* _jStr_PVxJ;
  std::vector<double>* _jStr_PVyJ;
  std::vector<double>* _jStr_PVzJ;
  std::vector<double>* _jStr_PVrJ;
  std::vector<double>* _jStr_PVchiSqrJ;
  std::vector<double>* _jStr_PVnDoFJ;
  std::vector<double>* _jStr_PVfitJ;
  std::vector<int>* _jStr_PVnTrkJ;
  std::vector<int>* _jStr_PVtypeJ;

  int _jStr_numJets;

  // EventShape, whole event, jets as input:
  std::vector<double>* _jStr_evtEJ_etaJ1;
  std::vector<double>* _jStr_evtEJ_etaJ2;
  std::vector<double>* _jStr_evtEJ_etaJ3;
  std::vector<double>* _jStr_evtEJ_etaJ4;
  std::vector<double>* _jStr_evtEJ_pTJ1;
  std::vector<double>* _jStr_evtEJ_pTJ2;
  std::vector<double>* _jStr_evtEJ_pTJ3;
  std::vector<double>* _jStr_evtEJ_pTJ4;

  std::vector<double>* _jStr_evtEJ_fE1;
  std::vector<double>* _jStr_evtEJ_fE2;
  std::vector<double>* _jStr_evtEJ_fE3;
  std::vector<double>* _jStr_evtEJ_fET1;
  std::vector<double>* _jStr_evtEJ_fET2;
  std::vector<double>* _jStr_evtEJ_fET3;
  std::vector<double>* _jStr_evtEJ_mct;
  std::vector<double>* _jStr_evtEJ_q;
  std::vector<double>* _jStr_evtEJ_mjj;
  std::vector<double>* _jStr_evtEJ_mTjj;
  std::vector<double>* _jStr_evtEJ_mjjj;
  std::vector<double>* _jStr_evtEJ_mTjjj;
  std::vector<double>* _jStr_evtEJ_mjjjj;
  std::vector<double>* _jStr_evtEJ_mTjjjj;
  std::vector<double>* _jStr_evtEJ_zetaPlus;
  std::vector<double>* _jStr_evtEJ_zetaMinus;
  std::vector<double>* _jStr_evtEJ_B12;
  std::vector<double>* _jStr_evtEJ_BT12;
  std::vector<double>* _jStr_evtEJ_Alpha;
  std::vector<double>* _jStr_evtEJ_AlphaT;
  std::vector<double>* _jStr_evtEJ_massDemocracy;
  std::vector<double>* _jStr_evtEJ_betaflow;
  std::vector<double>* _jStr_evtEJ_betaflow_1GeV;
  std::vector<double>* _jStr_evtEJ_betaflow_5GeV;
  std::vector<double>* _jStr_evtEJ_y23;
  std::vector<double>* _jStr_evtEJ_y23_1GeV;
  std::vector<double>* _jStr_evtEJ_y23_5GeV;
  std::vector<double>* _jStr_evtEJ_lny23;
  std::vector<double>* _jStr_evtEJ_lny23_1GeV;
  std::vector<double>* _jStr_evtEJ_lny23_5GeV;
  std::vector<double>* _jStr_evtEJ_theta;
  std::vector<double>* _jStr_evtEJ_asym;
  std::vector<double>* _jStr_evtEJ_yB;
  std::vector<double>* _jStr_evtEJ_yStar;
  std::vector<double>* _jStr_evtEJ_thetaStar;
  std::vector<double>* _jStr_evtEJ_chi;
  std::vector<double>* _jStr_evtEJ_deltaPhiJJ;
  std::vector<double>* _jStr_evtEJ_deltaThetaJJ;
  std::vector<double>* _jStr_evtEJ_deltaEtaJJ;
  std::vector<double>* _jStr_evtEJ_deltaRapidityJJ;
  std::vector<double>* _jStr_evtEJ_deltaRJJ;
  std::vector<double>* _jStr_evtEJ_deltaRJJY;
  std::vector<double>* _jStr_evtEJ_sigmaPhiJJ;
  std::vector<double>* _jStr_evtEJ_sigmaThetaJJ;
  std::vector<double>* _jStr_evtEJ_sigmaEtaJJ;
  std::vector<double>* _jStr_evtEJ_sigmaRapidityJJ;
  std::vector<double>* _jStr_evtEJ_sigmaPtJJ;
  std::vector<double>* _jStr_evtEJ_sigmaEtJJ;
  std::vector<double>* _jStr_evtEJ_sigmaEt12;
  std::vector<double>* _jStr_evtEJ_sigmaEt34;
  std::vector<double>* _jStr_evtEJ_A234;
  std::vector<double>* _jStr_evtEJ_asymPhiJJ;
  std::vector<double>* _jStr_evtEJ_asymThetaJJ;
  std::vector<double>* _jStr_evtEJ_asymEtaJJ;
  std::vector<double>* _jStr_evtEJ_asymRapidityJJ;
  std::vector<double>* _jStr_evtEJ_acoplanarity;
  std::vector<double>* _jStr_evtEJ_twist;
  std::vector<double>* _jStr_evtEJ_twistY;
  std::vector<double>* _jStr_evtEJ_jetSumE;
  std::vector<double>* _jStr_evtEJ_jetSumET;
  std::vector<double>* _jStr_evtEJ_jetSumPT;
  std::vector<double>* _jStr_evtEJ_jetSumM;
  std::vector<double>* _jStr_evtEJ_jetSumMT;
  std::vector<double>* _jStr_evtEJ_HTprime;
  std::vector<double>* _jStr_evtEJ_centrality;
  std::vector<double>* _jStr_evtEJ_centralityP;
  std::vector<double>* _jStr_evtEJ_zminJ1J2;
  std::vector<double>* _jStr_evtEJ_zmaxJ1J2;
  std::vector<double>* _jStr_evtEJ_zminAllJets;
  std::vector<double>* _jStr_evtEJ_zmaxAllJets;
  std::vector<double>* _jStr_evtEJ_zminJ1J2Phi;
  std::vector<double>* _jStr_evtEJ_zminJ1J2Theta;
  std::vector<double>* _jStr_evtEJ_zminJ1J2Eta;
  std::vector<double>* _jStr_evtEJ_zminJ1J2Rapidity;
  std::vector<double>* _jStr_evtEJ_cosHelicityJ1;
  std::vector<double>* _jStr_evtEJ_helicityJ1;
  std::vector<double>* _jStr_evtEJ_azilicityJ1;
  std::vector<double>* _jStr_evtEJ_cosHelicityJ2;
  std::vector<double>* _jStr_evtEJ_helicityJ2;
  std::vector<double>* _jStr_evtEJ_azilicityJ2;
  std::vector<double>* _jStr_evtEJ_cosThetaJ1;
  std::vector<double>* _jStr_evtEJ_cosThetaJ2;
  std::vector<double>* _jStr_evtEJ_deltaRapidityXtoJ1CM;
  std::vector<double>* _jStr_evtEJ_deltaRapidityXtoJ2CM;
  std::vector<double>* _jStr_evtEJ_deltaRapidityXtoJ1;
  std::vector<double>* _jStr_evtEJ_deltaRapidityXtoJ2;
  std::vector<double>* _jStr_evtEJ_cosTheta1;
  std::vector<double>* _jStr_evtEJ_cosTheta2;
  std::vector<double>* _jStr_evtEJ_cosThetaStar1;
  std::vector<double>* _jStr_evtEJ_cosThetaStar2;
  std::vector<double>* _jStr_evtEJ_cosPhiTilde1;
  std::vector<double>* _jStr_evtEJ_cosPhiTilde2;
  std::vector<double>* _jStr_evtEJ_cosPhi;
  std::vector<double>* _jStr_evtEJ_PhiTilde1;
  std::vector<double>* _jStr_evtEJ_PhiTilde2;
  std::vector<double>* _jStr_evtEJ_Phi;

  std::vector<double>* _jStr_evtEJ_Mprime4;
  std::vector<double>* _jStr_evtEJ_dMprime4;
  std::vector<double>* _jStr_evtEJ_MprimeAvg4;
  std::vector<double>* _jStr_evtEJ_DeltaMin4;
  std::vector<double>* _jStr_evtEJ_DeltaMax4;
  std::vector<double>* _jStr_evtEJ_DeltaPhiXX4;
  std::vector<double>* _jStr_evtEJ_DeltaYXX4;
  std::vector<double>* _jStr_evtEJ_DeltaRXX4;
  std::vector<double>* _jStr_evtEJ_TwistXX4;
  std::vector<int>* _jStr_evtEJ_separated4;
  std::vector<double>* _jStr_evtEJ_Mprime;
  std::vector<double>* _jStr_evtEJ_dMprime;
  std::vector<double>* _jStr_evtEJ_MprimeAvg;
  std::vector<double>* _jStr_evtEJ_DeltaMin;
  std::vector<double>* _jStr_evtEJ_DeltaMax;
  std::vector<double>* _jStr_evtEJ_DeltaPhiXX;
  std::vector<double>* _jStr_evtEJ_DeltaYXX;
  std::vector<double>* _jStr_evtEJ_DeltaRXX;
  std::vector<double>* _jStr_evtEJ_TwistXX;
  std::vector<int>* _jStr_evtEJ_separated;

  std::vector<double>* _jStr_evtEJ_M;
  std::vector<double>* _jStr_evtEJ_DeltaM;
  std::vector<double>* _jStr_evtEJ_asymM;
  std::vector<double>* _jStr_evtEJ_Q;
  std::vector<double>* _jStr_evtEJ_SCDF;
  std::vector<double>* _jStr_evtEJ_dPhiIJCDF;
  std::vector<double>* _jStr_evtEJ_dPhiKLCDF;
  std::vector<double>* _jStr_evtEJ_PtIPtJCDF;
  std::vector<double>* _jStr_evtEJ_PtKPtLCDF;
  std::vector<double>* _jStr_evtEJ_SD0;
  std::vector<double>* _jStr_evtEJ_dPhiIJD0;
  std::vector<double>* _jStr_evtEJ_dPhiKLD0;
  std::vector<double>* _jStr_evtEJ_PtIPtJD0;
  std::vector<double>* _jStr_evtEJ_PtKPtLD0;
  std::vector<double>* _jStr_evtEJ_DeltaSCDF;
  std::vector<double>* _jStr_evtEJ_DeltaSD0;
  std::vector<double>* _jStr_evtEJ_Delta12;
  std::vector<double>* _jStr_evtEJ_Delta12n;
  std::vector<double>* _jStr_evtEJ_M14;
  std::vector<double>* _jStr_evtEJ_y1y2;

  // EventShape, mini event (jet in own CM frame), subjets as input:
  std::vector<double>* _jStr_evtJj_etaJ1;
  std::vector<double>* _jStr_evtJj_etaJ2;
  std::vector<double>* _jStr_evtJj_etaJ3;
  std::vector<double>* _jStr_evtJj_etaJ4;
  std::vector<double>* _jStr_evtJj_pTJ1;
  std::vector<double>* _jStr_evtJj_pTJ2;
  std::vector<double>* _jStr_evtJj_pTJ3;
  std::vector<double>* _jStr_evtJj_pTJ4;

  std::vector<double>* _jStr_evtJj_fE1;
  std::vector<double>* _jStr_evtJj_fE2;
  std::vector<double>* _jStr_evtJj_fE3;
  std::vector<double>* _jStr_evtJj_fET1;
  std::vector<double>* _jStr_evtJj_fET2;
  std::vector<double>* _jStr_evtJj_fET3;
  std::vector<double>* _jStr_evtJj_mct;
  std::vector<double>* _jStr_evtJj_q;
  std::vector<double>* _jStr_evtJj_mjj;
  std::vector<double>* _jStr_evtJj_mTjj;
  std::vector<double>* _jStr_evtJj_mjjj;
  std::vector<double>* _jStr_evtJj_mTjjj;
  std::vector<double>* _jStr_evtJj_mjjjj;
  std::vector<double>* _jStr_evtJj_mTjjjj;
  std::vector<double>* _jStr_evtJj_zetaPlus;
  std::vector<double>* _jStr_evtJj_zetaMinus;
  std::vector<double>* _jStr_evtJj_B12;
  std::vector<double>* _jStr_evtJj_BT12;
  std::vector<double>* _jStr_evtJj_Alpha;
  std::vector<double>* _jStr_evtJj_AlphaT;
  std::vector<double>* _jStr_evtJj_massDemocracy;
  std::vector<double>* _jStr_evtJj_betaflow;
  std::vector<double>* _jStr_evtJj_betaflow_1GeV;
  std::vector<double>* _jStr_evtJj_betaflow_5GeV;
  std::vector<double>* _jStr_evtJj_y23;
  std::vector<double>* _jStr_evtJj_y23_1GeV;
  std::vector<double>* _jStr_evtJj_y23_5GeV;
  std::vector<double>* _jStr_evtJj_lny23;
  std::vector<double>* _jStr_evtJj_lny23_1GeV;
  std::vector<double>* _jStr_evtJj_lny23_5GeV;
  std::vector<double>* _jStr_evtJj_theta;
  std::vector<double>* _jStr_evtJj_asym;
  std::vector<double>* _jStr_evtJj_yB;
  std::vector<double>* _jStr_evtJj_yStar;
  std::vector<double>* _jStr_evtJj_thetaStar;
  std::vector<double>* _jStr_evtJj_chi;
  std::vector<double>* _jStr_evtJj_deltaPhiJJ;
  std::vector<double>* _jStr_evtJj_deltaThetaJJ;
  std::vector<double>* _jStr_evtJj_deltaEtaJJ;
  std::vector<double>* _jStr_evtJj_deltaRapidityJJ;
  std::vector<double>* _jStr_evtJj_deltaRJJ;
  std::vector<double>* _jStr_evtJj_deltaRJJY;
  std::vector<double>* _jStr_evtJj_sigmaPhiJJ;
  std::vector<double>* _jStr_evtJj_sigmaThetaJJ;
  std::vector<double>* _jStr_evtJj_sigmaEtaJJ;
  std::vector<double>* _jStr_evtJj_sigmaRapidityJJ;
  std::vector<double>* _jStr_evtJj_sigmaPtJJ;
  std::vector<double>* _jStr_evtJj_sigmaEtJJ;
  std::vector<double>* _jStr_evtJj_sigmaEt12;
  std::vector<double>* _jStr_evtJj_sigmaEt34;
  std::vector<double>* _jStr_evtJj_A234;
  std::vector<double>* _jStr_evtJj_asymPhiJJ;
  std::vector<double>* _jStr_evtJj_asymThetaJJ;
  std::vector<double>* _jStr_evtJj_asymEtaJJ;
  std::vector<double>* _jStr_evtJj_asymRapidityJJ;
  std::vector<double>* _jStr_evtJj_acoplanarity;
  std::vector<double>* _jStr_evtJj_twist;
  std::vector<double>* _jStr_evtJj_twistY;
  std::vector<double>* _jStr_evtJj_jetSumE;
  std::vector<double>* _jStr_evtJj_jetSumET;
  std::vector<double>* _jStr_evtJj_jetSumPT;
  std::vector<double>* _jStr_evtJj_jetSumM;
  std::vector<double>* _jStr_evtJj_jetSumMT;
  std::vector<double>* _jStr_evtJj_HTprime;
  std::vector<double>* _jStr_evtJj_centrality;
  std::vector<double>* _jStr_evtJj_centralityP;
  std::vector<double>* _jStr_evtJj_zminJ1J2;
  std::vector<double>* _jStr_evtJj_zmaxJ1J2;
  std::vector<double>* _jStr_evtJj_zminAllJets;
  std::vector<double>* _jStr_evtJj_zmaxAllJets;
  std::vector<double>* _jStr_evtJj_zminJ1J2Phi;
  std::vector<double>* _jStr_evtJj_zminJ1J2Theta;
  std::vector<double>* _jStr_evtJj_zminJ1J2Eta;
  std::vector<double>* _jStr_evtJj_zminJ1J2Rapidity;
  std::vector<double>* _jStr_evtJj_cosHelicityJ1;
  std::vector<double>* _jStr_evtJj_helicityJ1;
  std::vector<double>* _jStr_evtJj_azilicityJ1;
  std::vector<double>* _jStr_evtJj_cosHelicityJ2;
  std::vector<double>* _jStr_evtJj_helicityJ2;
  std::vector<double>* _jStr_evtJj_azilicityJ2;
  std::vector<double>* _jStr_evtJj_cosThetaJ1;
  std::vector<double>* _jStr_evtJj_cosThetaJ2;
  std::vector<double>* _jStr_evtJj_deltaRapidityXtoJ1CM;
  std::vector<double>* _jStr_evtJj_deltaRapidityXtoJ2CM;
  std::vector<double>* _jStr_evtJj_deltaRapidityXtoJ1;
  std::vector<double>* _jStr_evtJj_deltaRapidityXtoJ2;
  std::vector<double>* _jStr_evtJj_cosTheta1;
  std::vector<double>* _jStr_evtJj_cosTheta2;
  std::vector<double>* _jStr_evtJj_cosThetaStar1;
  std::vector<double>* _jStr_evtJj_cosThetaStar2;
  std::vector<double>* _jStr_evtJj_cosPhiTilde1;
  std::vector<double>* _jStr_evtJj_cosPhiTilde2;
  std::vector<double>* _jStr_evtJj_cosPhi;
  std::vector<double>* _jStr_evtJj_PhiTilde1;
  std::vector<double>* _jStr_evtJj_PhiTilde2;
  std::vector<double>* _jStr_evtJj_Phi;

  std::vector<double>* _jStr_evtJj_Mprime4;
  std::vector<double>* _jStr_evtJj_dMprime4;
  std::vector<double>* _jStr_evtJj_MprimeAvg4;
  std::vector<double>* _jStr_evtJj_DeltaMin4;
  std::vector<double>* _jStr_evtJj_DeltaMax4;
  std::vector<double>* _jStr_evtJj_DeltaPhiXX4;
  std::vector<double>* _jStr_evtJj_DeltaYXX4;
  std::vector<double>* _jStr_evtJj_DeltaRXX4;
  std::vector<double>* _jStr_evtJj_TwistXX4;
  std::vector<int>* _jStr_evtJj_separated4;
  std::vector<double>* _jStr_evtJj_Mprime;
  std::vector<double>* _jStr_evtJj_dMprime;
  std::vector<double>* _jStr_evtJj_MprimeAvg;
  std::vector<double>* _jStr_evtJj_DeltaMin;
  std::vector<double>* _jStr_evtJj_DeltaMax;
  std::vector<double>* _jStr_evtJj_DeltaPhiXX;
  std::vector<double>* _jStr_evtJj_DeltaYXX;
  std::vector<double>* _jStr_evtJj_DeltaRXX;
  std::vector<double>* _jStr_evtJj_TwistXX;
  std::vector<int>* _jStr_evtJj_separated;

  std::vector<double>* _jStr_evtJj_M;
  std::vector<double>* _jStr_evtJj_DeltaM;
  std::vector<double>* _jStr_evtJj_asymM;
  std::vector<double>* _jStr_evtJj_Q;
  std::vector<double>* _jStr_evtJj_SCDF;
  std::vector<double>* _jStr_evtJj_dPhiIJCDF;
  std::vector<double>* _jStr_evtJj_dPhiKLCDF;
  std::vector<double>* _jStr_evtJj_PtIPtJCDF;
  std::vector<double>* _jStr_evtJj_PtKPtLCDF;
  std::vector<double>* _jStr_evtJj_SD0;
  std::vector<double>* _jStr_evtJj_dPhiIJD0;
  std::vector<double>* _jStr_evtJj_dPhiKLD0;
  std::vector<double>* _jStr_evtJj_PtIPtJD0;
  std::vector<double>* _jStr_evtJj_PtKPtLD0;
  std::vector<double>* _jStr_evtJj_DeltaSCDF;
  std::vector<double>* _jStr_evtJj_DeltaSD0;
  std::vector<double>* _jStr_evtJj_Delta12;
  std::vector<double>* _jStr_evtJj_Delta12n;
  std::vector<double>* _jStr_evtJj_M14;
  std::vector<double>* _jStr_evtJj_y1y2;

  std::vector<double>* _jStr_fwmEJ_xiPlus;
  std::vector<double>* _jStr_fwmEJ_xiMinus;

  std::vector<double>* _jStr_fwmEJ_xPlus;
  std::vector<double>* _jStr_fwmEJ_xMinus;

  // Fox-Wolfram, whole event, jets as input:
  std::vector<double>* _jStr_fwmEJ_Psi1;
  std::vector<double>* _jStr_fwmEJ_B0;
  std::vector<double>* _jStr_fwmEJ_B1;
  std::vector<double>* _jStr_fwmEJ_B2;
  std::vector<double>* _jStr_fwmEJ_B3;
  std::vector<double>* _jStr_fwmEJ_B4;
  std::vector<double>* _jStr_fwmEJ_B5;
  std::vector<double>* _jStr_fwmEJ_B6;
  std::vector<double>* _jStr_fwmEJ_B7;
  std::vector<double>* _jStr_fwmEJ_B8;
  std::vector<double>* _jStr_fwmEJ_C0;
  std::vector<double>* _jStr_fwmEJ_C1;
  std::vector<double>* _jStr_fwmEJ_C2;
  std::vector<double>* _jStr_fwmEJ_C3;
  std::vector<double>* _jStr_fwmEJ_C4;
  std::vector<double>* _jStr_fwmEJ_C5;
  std::vector<double>* _jStr_fwmEJ_C6;
  std::vector<double>* _jStr_fwmEJ_C7;
  std::vector<double>* _jStr_fwmEJ_C8;
  std::vector<double>* _jStr_fwmEJ_K0;
  std::vector<double>* _jStr_fwmEJ_K1;
  std::vector<double>* _jStr_fwmEJ_K2;
  std::vector<double>* _jStr_fwmEJ_K3;
  std::vector<double>* _jStr_fwmEJ_K4;
  std::vector<double>* _jStr_fwmEJ_K5;
  std::vector<double>* _jStr_fwmEJ_K6;
  std::vector<double>* _jStr_fwmEJ_K7;
  std::vector<double>* _jStr_fwmEJ_K8;
  std::vector<double>* _jStr_fwmEJ_D0;
  std::vector<double>* _jStr_fwmEJ_D1;
  std::vector<double>* _jStr_fwmEJ_D2;
  std::vector<double>* _jStr_fwmEJ_D3;
  std::vector<double>* _jStr_fwmEJ_D4;
  std::vector<double>* _jStr_fwmEJ_D5;
  std::vector<double>* _jStr_fwmEJ_D6;
  std::vector<double>* _jStr_fwmEJ_D7;
  std::vector<double>* _jStr_fwmEJ_D8;
  std::vector<double>* _jStr_fwmEJ_H0;
  std::vector<double>* _jStr_fwmEJ_H1;
  std::vector<double>* _jStr_fwmEJ_H2;
  std::vector<double>* _jStr_fwmEJ_H3;
  std::vector<double>* _jStr_fwmEJ_H4;
  std::vector<double>* _jStr_fwmEJ_H5;
  std::vector<double>* _jStr_fwmEJ_H6;
  std::vector<double>* _jStr_fwmEJ_H7;
  std::vector<double>* _jStr_fwmEJ_H8;
  std::vector<double>* _jStr_fwmEJ_Q0;
  std::vector<double>* _jStr_fwmEJ_Q1;
  std::vector<double>* _jStr_fwmEJ_Q2;
  std::vector<double>* _jStr_fwmEJ_Q3;
  std::vector<double>* _jStr_fwmEJ_Q4;
  std::vector<double>* _jStr_fwmEJ_Q5;
  std::vector<double>* _jStr_fwmEJ_Q6;
  std::vector<double>* _jStr_fwmEJ_Q7;
  std::vector<double>* _jStr_fwmEJ_Q8;
  std::vector<double>* _jStr_fwmEJ_Pi1;
  std::vector<double>* _jStr_fwmEJ_Pi2;
  std::vector<double>* _jStr_fwmEJ_Pi3;
  std::vector<double>* _jStr_fwmEJ_Pi4;
  std::vector<double>* _jStr_fwmEJ_B10;
  std::vector<double>* _jStr_fwmEJ_B20;
  std::vector<double>* _jStr_fwmEJ_B30;
  std::vector<double>* _jStr_fwmEJ_B40;
  std::vector<double>* _jStr_fwmEJ_B50;
  std::vector<double>* _jStr_fwmEJ_B60;
  std::vector<double>* _jStr_fwmEJ_B70;
  std::vector<double>* _jStr_fwmEJ_B80;
  std::vector<double>* _jStr_fwmEJ_C10;
  std::vector<double>* _jStr_fwmEJ_C20;
  std::vector<double>* _jStr_fwmEJ_C30;
  std::vector<double>* _jStr_fwmEJ_C40;
  std::vector<double>* _jStr_fwmEJ_C50;
  std::vector<double>* _jStr_fwmEJ_C60;
  std::vector<double>* _jStr_fwmEJ_C70;
  std::vector<double>* _jStr_fwmEJ_C80;
  std::vector<double>* _jStr_fwmEJ_K10;
  std::vector<double>* _jStr_fwmEJ_K20;
  std::vector<double>* _jStr_fwmEJ_K30;
  std::vector<double>* _jStr_fwmEJ_K40;
  std::vector<double>* _jStr_fwmEJ_K50;
  std::vector<double>* _jStr_fwmEJ_K60;
  std::vector<double>* _jStr_fwmEJ_K70;
  std::vector<double>* _jStr_fwmEJ_K80;
  std::vector<double>* _jStr_fwmEJ_D10;
  std::vector<double>* _jStr_fwmEJ_D20;
  std::vector<double>* _jStr_fwmEJ_D30;
  std::vector<double>* _jStr_fwmEJ_D40;
  std::vector<double>* _jStr_fwmEJ_D50;
  std::vector<double>* _jStr_fwmEJ_D60;
  std::vector<double>* _jStr_fwmEJ_D70;
  std::vector<double>* _jStr_fwmEJ_D80;
  std::vector<double>* _jStr_fwmEJ_H10;
  std::vector<double>* _jStr_fwmEJ_H20;
  std::vector<double>* _jStr_fwmEJ_H30;
  std::vector<double>* _jStr_fwmEJ_H40;
  std::vector<double>* _jStr_fwmEJ_H50;
  std::vector<double>* _jStr_fwmEJ_H60;
  std::vector<double>* _jStr_fwmEJ_H70;
  std::vector<double>* _jStr_fwmEJ_H80;
  std::vector<double>* _jStr_fwmEJ_Q10;
  std::vector<double>* _jStr_fwmEJ_Q20;
  std::vector<double>* _jStr_fwmEJ_Q30;
  std::vector<double>* _jStr_fwmEJ_Q40;
  std::vector<double>* _jStr_fwmEJ_Q50;
  std::vector<double>* _jStr_fwmEJ_Q60;
  std::vector<double>* _jStr_fwmEJ_Q70;
  std::vector<double>* _jStr_fwmEJ_Q80;

  // // Fox-Wolfram, whole event, jet constituents as input:
  // std::vector<double>* _jStr_fwmEi_Psi1;
  // std::vector<double>* _jStr_fwmEi_B0;
  // std::vector<double>* _jStr_fwmEi_B1;
  // std::vector<double>* _jStr_fwmEi_B2;
  // std::vector<double>* _jStr_fwmEi_B3;
  // std::vector<double>* _jStr_fwmEi_B4;
  // std::vector<double>* _jStr_fwmEi_B5;
  // std::vector<double>* _jStr_fwmEi_B6;
  // std::vector<double>* _jStr_fwmEi_B7;
  // std::vector<double>* _jStr_fwmEi_B8;
  // std::vector<double>* _jStr_fwmEi_C0;
  // std::vector<double>* _jStr_fwmEi_C1;
  // std::vector<double>* _jStr_fwmEi_C2;
  // std::vector<double>* _jStr_fwmEi_C3;
  // std::vector<double>* _jStr_fwmEi_C4;
  // std::vector<double>* _jStr_fwmEi_C5;
  // std::vector<double>* _jStr_fwmEi_C6;
  // std::vector<double>* _jStr_fwmEi_C7;
  // std::vector<double>* _jStr_fwmEi_C8;
  // std::vector<double>* _jStr_fwmEi_K0;
  // std::vector<double>* _jStr_fwmEi_K1;
  // std::vector<double>* _jStr_fwmEi_K2;
  // std::vector<double>* _jStr_fwmEi_K3;
  // std::vector<double>* _jStr_fwmEi_K4;
  // std::vector<double>* _jStr_fwmEi_K5;
  // std::vector<double>* _jStr_fwmEi_K6;
  // std::vector<double>* _jStr_fwmEi_K7;
  // std::vector<double>* _jStr_fwmEi_K8;
  // std::vector<double>* _jStr_fwmEi_D0;
  // std::vector<double>* _jStr_fwmEi_D1;
  // std::vector<double>* _jStr_fwmEi_D2;
  // std::vector<double>* _jStr_fwmEi_D3;
  // std::vector<double>* _jStr_fwmEi_D4;
  // std::vector<double>* _jStr_fwmEi_D5;
  // std::vector<double>* _jStr_fwmEi_D6;
  // std::vector<double>* _jStr_fwmEi_D7;
  // std::vector<double>* _jStr_fwmEi_D8;
  // std::vector<double>* _jStr_fwmEi_H0;
  // std::vector<double>* _jStr_fwmEi_H1;
  // std::vector<double>* _jStr_fwmEi_H2;
  // std::vector<double>* _jStr_fwmEi_H3;
  // std::vector<double>* _jStr_fwmEi_H4;
  // std::vector<double>* _jStr_fwmEi_H5;
  // std::vector<double>* _jStr_fwmEi_H6;
  // std::vector<double>* _jStr_fwmEi_H7;
  // std::vector<double>* _jStr_fwmEi_H8;
  // std::vector<double>* _jStr_fwmEi_Q0;
  // std::vector<double>* _jStr_fwmEi_Q1;
  // std::vector<double>* _jStr_fwmEi_Q2;
  // std::vector<double>* _jStr_fwmEi_Q3;
  // std::vector<double>* _jStr_fwmEi_Q4;
  // std::vector<double>* _jStr_fwmEi_Q5;
  // std::vector<double>* _jStr_fwmEi_Q6;
  // std::vector<double>* _jStr_fwmEi_Q7;
  // std::vector<double>* _jStr_fwmEi_Q8;
  // std::vector<double>* _jStr_fwmEi_Pi1;
  // std::vector<double>* _jStr_fwmEi_Pi2;
  // std::vector<double>* _jStr_fwmEi_Pi3;
  // std::vector<double>* _jStr_fwmEi_Pi4;
  // std::vector<double>* _jStr_fwmEi_B10;
  // std::vector<double>* _jStr_fwmEi_B20;
  // std::vector<double>* _jStr_fwmEi_B30;
  // std::vector<double>* _jStr_fwmEi_B40;
  // std::vector<double>* _jStr_fwmEi_B50;
  // std::vector<double>* _jStr_fwmEi_B60;
  // std::vector<double>* _jStr_fwmEi_B70;
  // std::vector<double>* _jStr_fwmEi_B80;
  // std::vector<double>* _jStr_fwmEi_C10;
  // std::vector<double>* _jStr_fwmEi_C20;
  // std::vector<double>* _jStr_fwmEi_C30;
  // std::vector<double>* _jStr_fwmEi_C40;
  // std::vector<double>* _jStr_fwmEi_C50;
  // std::vector<double>* _jStr_fwmEi_C60;
  // std::vector<double>* _jStr_fwmEi_C70;
  // std::vector<double>* _jStr_fwmEi_C80;
  // std::vector<double>* _jStr_fwmEi_K10;
  // std::vector<double>* _jStr_fwmEi_K20;
  // std::vector<double>* _jStr_fwmEi_K30;
  // std::vector<double>* _jStr_fwmEi_K40;
  // std::vector<double>* _jStr_fwmEi_K50;
  // std::vector<double>* _jStr_fwmEi_K60;
  // std::vector<double>* _jStr_fwmEi_K70;
  // std::vector<double>* _jStr_fwmEi_K80;
  // std::vector<double>* _jStr_fwmEi_D10;
  // std::vector<double>* _jStr_fwmEi_D20;
  // std::vector<double>* _jStr_fwmEi_D30;
  // std::vector<double>* _jStr_fwmEi_D40;
  // std::vector<double>* _jStr_fwmEi_D50;
  // std::vector<double>* _jStr_fwmEi_D60;
  // std::vector<double>* _jStr_fwmEi_D70;
  // std::vector<double>* _jStr_fwmEi_D80;
  // std::vector<double>* _jStr_fwmEi_H10;
  // std::vector<double>* _jStr_fwmEi_H20;
  // std::vector<double>* _jStr_fwmEi_H30;
  // std::vector<double>* _jStr_fwmEi_H40;
  // std::vector<double>* _jStr_fwmEi_H50;
  // std::vector<double>* _jStr_fwmEi_H60;
  // std::vector<double>* _jStr_fwmEi_H70;
  // std::vector<double>* _jStr_fwmEi_H80;
  // std::vector<double>* _jStr_fwmEi_Q10;
  // std::vector<double>* _jStr_fwmEi_Q20;
  // std::vector<double>* _jStr_fwmEi_Q30;
  // std::vector<double>* _jStr_fwmEi_Q40;
  // std::vector<double>* _jStr_fwmEi_Q50;
  // std::vector<double>* _jStr_fwmEi_Q60;
  // std::vector<double>* _jStr_fwmEi_Q70;
  // std::vector<double>* _jStr_fwmEi_Q80;

  std::vector<double>* _jStr_fwmJj_xiPlus;
  std::vector<double>* _jStr_fwmJj_xiMinus;

  std::vector<double>* _jStr_fwmJj_xPlus;
  std::vector<double>* _jStr_fwmJj_xMinus;

  // Fox-Wolfram, mini event, subjets as input:
  std::vector<double>* _jStr_fwmJj_Psi1;
  std::vector<double>* _jStr_fwmJj_B0;
  std::vector<double>* _jStr_fwmJj_B1;
  std::vector<double>* _jStr_fwmJj_B2;
  std::vector<double>* _jStr_fwmJj_B3;
  std::vector<double>* _jStr_fwmJj_B4;
  std::vector<double>* _jStr_fwmJj_B5;
  std::vector<double>* _jStr_fwmJj_B6;
  std::vector<double>* _jStr_fwmJj_B7;
  std::vector<double>* _jStr_fwmJj_B8;
  std::vector<double>* _jStr_fwmJj_C0;
  std::vector<double>* _jStr_fwmJj_C1;
  std::vector<double>* _jStr_fwmJj_C2;
  std::vector<double>* _jStr_fwmJj_C3;
  std::vector<double>* _jStr_fwmJj_C4;
  std::vector<double>* _jStr_fwmJj_C5;
  std::vector<double>* _jStr_fwmJj_C6;
  std::vector<double>* _jStr_fwmJj_C7;
  std::vector<double>* _jStr_fwmJj_C8;
  std::vector<double>* _jStr_fwmJj_K0;
  std::vector<double>* _jStr_fwmJj_K1;
  std::vector<double>* _jStr_fwmJj_K2;
  std::vector<double>* _jStr_fwmJj_K3;
  std::vector<double>* _jStr_fwmJj_K4;
  std::vector<double>* _jStr_fwmJj_K5;
  std::vector<double>* _jStr_fwmJj_K6;
  std::vector<double>* _jStr_fwmJj_K7;
  std::vector<double>* _jStr_fwmJj_K8;
  std::vector<double>* _jStr_fwmJj_D0;
  std::vector<double>* _jStr_fwmJj_D1;
  std::vector<double>* _jStr_fwmJj_D2;
  std::vector<double>* _jStr_fwmJj_D3;
  std::vector<double>* _jStr_fwmJj_D4;
  std::vector<double>* _jStr_fwmJj_D5;
  std::vector<double>* _jStr_fwmJj_D6;
  std::vector<double>* _jStr_fwmJj_D7;
  std::vector<double>* _jStr_fwmJj_D8;
  std::vector<double>* _jStr_fwmJj_H0;
  std::vector<double>* _jStr_fwmJj_H1;
  std::vector<double>* _jStr_fwmJj_H2;
  std::vector<double>* _jStr_fwmJj_H3;
  std::vector<double>* _jStr_fwmJj_H4;
  std::vector<double>* _jStr_fwmJj_H5;
  std::vector<double>* _jStr_fwmJj_H6;
  std::vector<double>* _jStr_fwmJj_H7;
  std::vector<double>* _jStr_fwmJj_H8;
  std::vector<double>* _jStr_fwmJj_Q0;
  std::vector<double>* _jStr_fwmJj_Q1;
  std::vector<double>* _jStr_fwmJj_Q2;
  std::vector<double>* _jStr_fwmJj_Q3;
  std::vector<double>* _jStr_fwmJj_Q4;
  std::vector<double>* _jStr_fwmJj_Q5;
  std::vector<double>* _jStr_fwmJj_Q6;
  std::vector<double>* _jStr_fwmJj_Q7;
  std::vector<double>* _jStr_fwmJj_Q8;
  std::vector<double>* _jStr_fwmJj_Pi1;
  std::vector<double>* _jStr_fwmJj_Pi2;
  std::vector<double>* _jStr_fwmJj_Pi3;
  std::vector<double>* _jStr_fwmJj_Pi4;
  std::vector<double>* _jStr_fwmJj_B10;
  std::vector<double>* _jStr_fwmJj_B20;
  std::vector<double>* _jStr_fwmJj_B30;
  std::vector<double>* _jStr_fwmJj_B40;
  std::vector<double>* _jStr_fwmJj_B50;
  std::vector<double>* _jStr_fwmJj_B60;
  std::vector<double>* _jStr_fwmJj_B70;
  std::vector<double>* _jStr_fwmJj_B80;
  std::vector<double>* _jStr_fwmJj_C10;
  std::vector<double>* _jStr_fwmJj_C20;
  std::vector<double>* _jStr_fwmJj_C30;
  std::vector<double>* _jStr_fwmJj_C40;
  std::vector<double>* _jStr_fwmJj_C50;
  std::vector<double>* _jStr_fwmJj_C60;
  std::vector<double>* _jStr_fwmJj_C70;
  std::vector<double>* _jStr_fwmJj_C80;
  std::vector<double>* _jStr_fwmJj_K10;
  std::vector<double>* _jStr_fwmJj_K20;
  std::vector<double>* _jStr_fwmJj_K30;
  std::vector<double>* _jStr_fwmJj_K40;
  std::vector<double>* _jStr_fwmJj_K50;
  std::vector<double>* _jStr_fwmJj_K60;
  std::vector<double>* _jStr_fwmJj_K70;
  std::vector<double>* _jStr_fwmJj_K80;
  std::vector<double>* _jStr_fwmJj_D10;
  std::vector<double>* _jStr_fwmJj_D20;
  std::vector<double>* _jStr_fwmJj_D30;
  std::vector<double>* _jStr_fwmJj_D40;
  std::vector<double>* _jStr_fwmJj_D50;
  std::vector<double>* _jStr_fwmJj_D60;
  std::vector<double>* _jStr_fwmJj_D70;
  std::vector<double>* _jStr_fwmJj_D80;
  std::vector<double>* _jStr_fwmJj_H10;
  std::vector<double>* _jStr_fwmJj_H20;
  std::vector<double>* _jStr_fwmJj_H30;
  std::vector<double>* _jStr_fwmJj_H40;
  std::vector<double>* _jStr_fwmJj_H50;
  std::vector<double>* _jStr_fwmJj_H60;
  std::vector<double>* _jStr_fwmJj_H70;
  std::vector<double>* _jStr_fwmJj_H80;
  std::vector<double>* _jStr_fwmJj_Q10;
  std::vector<double>* _jStr_fwmJj_Q20;
  std::vector<double>* _jStr_fwmJj_Q30;
  std::vector<double>* _jStr_fwmJj_Q40;
  std::vector<double>* _jStr_fwmJj_Q50;
  std::vector<double>* _jStr_fwmJj_Q60;
  std::vector<double>* _jStr_fwmJj_Q70;
  std::vector<double>* _jStr_fwmJj_Q80;

  std::vector<double>* _jStr_fwmJi_xiPlus;
  std::vector<double>* _jStr_fwmJi_xiMinus;

  std::vector<double>* _jStr_fwmJi_xPlus;
  std::vector<double>* _jStr_fwmJi_xMinus;

  // Fox-Wolfram, mini event, jet constituents as input:
  std::vector<double>* _jStr_fwmJi_Psi1;
  std::vector<double>* _jStr_fwmJi_B0;
  std::vector<double>* _jStr_fwmJi_B1;
  std::vector<double>* _jStr_fwmJi_B2;
  std::vector<double>* _jStr_fwmJi_B3;
  std::vector<double>* _jStr_fwmJi_B4;
  std::vector<double>* _jStr_fwmJi_B5;
  std::vector<double>* _jStr_fwmJi_B6;
  std::vector<double>* _jStr_fwmJi_B7;
  std::vector<double>* _jStr_fwmJi_B8;
  std::vector<double>* _jStr_fwmJi_C0;
  std::vector<double>* _jStr_fwmJi_C1;
  std::vector<double>* _jStr_fwmJi_C2;
  std::vector<double>* _jStr_fwmJi_C3;
  std::vector<double>* _jStr_fwmJi_C4;
  std::vector<double>* _jStr_fwmJi_C5;
  std::vector<double>* _jStr_fwmJi_C6;
  std::vector<double>* _jStr_fwmJi_C7;
  std::vector<double>* _jStr_fwmJi_C8;
  std::vector<double>* _jStr_fwmJi_K0;
  std::vector<double>* _jStr_fwmJi_K1;
  std::vector<double>* _jStr_fwmJi_K2;
  std::vector<double>* _jStr_fwmJi_K3;
  std::vector<double>* _jStr_fwmJi_K4;
  std::vector<double>* _jStr_fwmJi_K5;
  std::vector<double>* _jStr_fwmJi_K6;
  std::vector<double>* _jStr_fwmJi_K7;
  std::vector<double>* _jStr_fwmJi_K8;
  std::vector<double>* _jStr_fwmJi_D0;
  std::vector<double>* _jStr_fwmJi_D1;
  std::vector<double>* _jStr_fwmJi_D2;
  std::vector<double>* _jStr_fwmJi_D3;
  std::vector<double>* _jStr_fwmJi_D4;
  std::vector<double>* _jStr_fwmJi_D5;
  std::vector<double>* _jStr_fwmJi_D6;
  std::vector<double>* _jStr_fwmJi_D7;
  std::vector<double>* _jStr_fwmJi_D8;
  std::vector<double>* _jStr_fwmJi_H0;
  std::vector<double>* _jStr_fwmJi_H1;
  std::vector<double>* _jStr_fwmJi_H2;
  std::vector<double>* _jStr_fwmJi_H3;
  std::vector<double>* _jStr_fwmJi_H4;
  std::vector<double>* _jStr_fwmJi_H5;
  std::vector<double>* _jStr_fwmJi_H6;
  std::vector<double>* _jStr_fwmJi_H7;
  std::vector<double>* _jStr_fwmJi_H8;
  std::vector<double>* _jStr_fwmJi_Q0;
  std::vector<double>* _jStr_fwmJi_Q1;
  std::vector<double>* _jStr_fwmJi_Q2;
  std::vector<double>* _jStr_fwmJi_Q3;
  std::vector<double>* _jStr_fwmJi_Q4;
  std::vector<double>* _jStr_fwmJi_Q5;
  std::vector<double>* _jStr_fwmJi_Q6;
  std::vector<double>* _jStr_fwmJi_Q7;
  std::vector<double>* _jStr_fwmJi_Q8;
  std::vector<double>* _jStr_fwmJi_Pi1;
  std::vector<double>* _jStr_fwmJi_Pi2;
  std::vector<double>* _jStr_fwmJi_Pi3;
  std::vector<double>* _jStr_fwmJi_Pi4;
  std::vector<double>* _jStr_fwmJi_B10;
  std::vector<double>* _jStr_fwmJi_B20;
  std::vector<double>* _jStr_fwmJi_B30;
  std::vector<double>* _jStr_fwmJi_B40;
  std::vector<double>* _jStr_fwmJi_B50;
  std::vector<double>* _jStr_fwmJi_B60;
  std::vector<double>* _jStr_fwmJi_B70;
  std::vector<double>* _jStr_fwmJi_B80;
  std::vector<double>* _jStr_fwmJi_C10;
  std::vector<double>* _jStr_fwmJi_C20;
  std::vector<double>* _jStr_fwmJi_C30;
  std::vector<double>* _jStr_fwmJi_C40;
  std::vector<double>* _jStr_fwmJi_C50;
  std::vector<double>* _jStr_fwmJi_C60;
  std::vector<double>* _jStr_fwmJi_C70;
  std::vector<double>* _jStr_fwmJi_C80;
  std::vector<double>* _jStr_fwmJi_K10;
  std::vector<double>* _jStr_fwmJi_K20;
  std::vector<double>* _jStr_fwmJi_K30;
  std::vector<double>* _jStr_fwmJi_K40;
  std::vector<double>* _jStr_fwmJi_K50;
  std::vector<double>* _jStr_fwmJi_K60;
  std::vector<double>* _jStr_fwmJi_K70;
  std::vector<double>* _jStr_fwmJi_K80;
  std::vector<double>* _jStr_fwmJi_D10;
  std::vector<double>* _jStr_fwmJi_D20;
  std::vector<double>* _jStr_fwmJi_D30;
  std::vector<double>* _jStr_fwmJi_D40;
  std::vector<double>* _jStr_fwmJi_D50;
  std::vector<double>* _jStr_fwmJi_D60;
  std::vector<double>* _jStr_fwmJi_D70;
  std::vector<double>* _jStr_fwmJi_D80;
  std::vector<double>* _jStr_fwmJi_H10;
  std::vector<double>* _jStr_fwmJi_H20;
  std::vector<double>* _jStr_fwmJi_H30;
  std::vector<double>* _jStr_fwmJi_H40;
  std::vector<double>* _jStr_fwmJi_H50;
  std::vector<double>* _jStr_fwmJi_H60;
  std::vector<double>* _jStr_fwmJi_H70;
  std::vector<double>* _jStr_fwmJi_H80;
  std::vector<double>* _jStr_fwmJi_Q10;
  std::vector<double>* _jStr_fwmJi_Q20;
  std::vector<double>* _jStr_fwmJi_Q30;
  std::vector<double>* _jStr_fwmJi_Q40;
  std::vector<double>* _jStr_fwmJi_Q50;
  std::vector<double>* _jStr_fwmJi_Q60;
  std::vector<double>* _jStr_fwmJi_Q70;
  std::vector<double>* _jStr_fwmJi_Q80;

  // Flags:
  std::vector<std::string>* _jStr_authorE;
  std::vector<std::string>* _jStr_authorJ;
  std::vector<int>* _jStr_radialParamE;
  std::vector<int>* _jStr_radialParamJ;
  std::vector<int>* _jStr_algE;
  std::vector<int>* _jStr_algJ;
  std::vector<int>* _jStr_inputE;
  std::vector<int>* _jStr_inputJ;
  std::vector<int>* _jStr_bjetE;
  std::vector<int>* _jStr_bjetJ;
  std::vector<int>* _jStr_xjetE;
  std::vector<int>* _jStr_xjetJ;
  std::vector<int>* _jStr_myJetsE;
  std::vector<int>* _jStr_myJetsJ;
  std::vector<int>* _jStr_indexJ;

  // Misc jet properties:
  std::vector<int>* _jStr_jetMiscJ_numConstituents;
  std::vector<double>* _jStr_jetMiscJ_jetM;
  std::vector<double>* _jStr_jetMiscJ_jetMt;
  std::vector<double>* _jStr_jetMiscJ_jetE;
  std::vector<double>* _jStr_jetMiscJ_jetP;
  std::vector<double>* _jStr_jetMiscJ_jetEt;
  std::vector<double>* _jStr_jetMiscJ_jetPt;
  std::vector<double>* _jStr_jetMiscJ_jetPhi;
  std::vector<double>* _jStr_jetMiscJ_jetEta;
  std::vector<double>* _jStr_jetMiscJ_jetRapidity;
  std::vector<double>* _jStr_jetMiscJ_xJ;
  std::vector<double>* _jStr_jetMiscJ_gamma;
  std::vector<double>* _jStr_jetMiscJ_R;
  std::vector<double>* _jStr_jetMiscJ_Cbar;
  std::vector<int>* _jStr_jetMiscJ_numSubjets;

  // Jet pull vars:
  std::vector<double>* _jStr_pullJ_det;
  std::vector<double>* _jStr_pullJ_ratio;
  std::vector<double>* _jStr_pullJ_pullPf;
  std::vector<double>* _jStr_pullJ_angularEccentricity;
  std::vector<double>* _jStr_pullJ_orientation;
  std::vector<double>* _jStr_pullJ_girth;
  std::vector<double>* _jStr_pullJ_Cbar;
  std::vector<double>* _jStr_pullJ_g;
  std::vector<double>* _jStr_pullJ_e;
  std::vector<double>* _jStr_pullJ_B;
  std::vector<double>* _jStr_pullJ_logB;
  std::vector<double>* _jStr_pullJ_pullTheta;
  std::vector<double>* _jStr_pullJ_pullMag;

  // Unburied Higgs vars:
  std::vector<double>* _jStr_ubhJ_z;
  std::vector<double>* _jStr_ubhJ_z2;
  std::vector<double>* _jStr_ubhJ_a1;
  std::vector<double>* _jStr_ubhJ_a2;
  std::vector<double>* _jStr_ubhJ_a3;
  std::vector<double>* _jStr_ubhJ_meanpt;
  std::vector<double>* _jStr_ubhJ_meanet;
  std::vector<double>* _jStr_ubhJ_mbar;
  std::vector<double>* _jStr_ubhJ_massDemocracy;
  std::vector<double>* _jStr_ubhJ_fE1;
  std::vector<double>* _jStr_ubhJ_fE2;
  std::vector<double>* _jStr_ubhJ_fE3;
  std::vector<double>* _jStr_ubhJ_fET1;
  std::vector<double>* _jStr_ubhJ_fET2;
  std::vector<double>* _jStr_ubhJ_fET3;
  std::vector<double>* _jStr_ubhJ_Alpha;
  std::vector<double>* _jStr_ubhJ_AlphaT;
  std::vector<double>* _jStr_ubhJ_betaflow;
  std::vector<double>* _jStr_ubhJ_betaflow_1GeV;
  std::vector<double>* _jStr_ubhJ_betaflow_5GeV;
  std::vector<double>* _jStr_ubhJ_y23;
  std::vector<double>* _jStr_ubhJ_y23_1GeV;
  std::vector<double>* _jStr_ubhJ_y23_5GeV;
  std::vector<double>* _jStr_ubhJ_lny23;
  std::vector<double>* _jStr_ubhJ_lny23_1GeV;
  std::vector<double>* _jStr_ubhJ_lny23_5GeV;
  std::vector<double>* _jStr_ubhJ_subjetAsymmetry;
  std::vector<double>* _jStr_dipJ_dipolarity;

  // N-subjettiness vars:
  std::vector<double>* _jStr_nsubjnessJ_tau1;
  std::vector<double>* _jStr_nsubjnessJ_tau2;
  std::vector<double>* _jStr_nsubjnessJ_tau3;
  std::vector<double>* _jStr_nsubjnessJ_tau2tau1;
  std::vector<double>* _jStr_nsubjnessJ_tau3tau2;

  // psi jet shapes:
  std::vector<double>* _jStr_psiJ2_psi;
  std::vector<double>* _jStr_psiJ2_rho;
  std::vector<double>* _jStr_psiRatiosJ_psi1;
  std::vector<double>* _jStr_psiRatiosJ_psi2;
  std::vector<double>* _jStr_psiRatiosJ_psi3;
  std::vector<double>* _jStr_psiRatiosJ_psi7;
  std::vector<double>* _jStr_psiRatiosJ_psi717;
  std::vector<double>* _jStr_psiRatiosJ_psi127;
  std::vector<double>* _jStr_psiRatiosJ_psi37;

  // BCF jet shapes:
  std::vector<int>* _jStr_bcfJ_bcfVersion_v0a1;
  std::vector<int>* _jStr_bcfJ_a_v0a1;
  std::vector<double>* _jStr_bcfJ_bcfT_v0a1;
  std::vector<double>* _jStr_bcfJ_bcf_v0a1;
  std::vector<double>* _jStr_bcfJ_bcfAsymY_v0a1;
  std::vector<double>* _jStr_bcfJ_bcfAsymPhi_v0a1;
  std::vector<double>* _jStr_bcfJ_bcfAsymYPhi_v0a1;
  std::vector<double>* _jStr_bcfJ_bcfAsymYPhi2_v0a1;

  std::vector<int>* _jStr_bcfJ_bcfVersion_v0a2;
  std::vector<int>* _jStr_bcfJ_a_v0a2;
  std::vector<double>* _jStr_bcfJ_bcfT_v0a2;
  std::vector<double>* _jStr_bcfJ_bcf_v0a2;
  std::vector<double>* _jStr_bcfJ_bcfAsymY_v0a2;
  std::vector<double>* _jStr_bcfJ_bcfAsymPhi_v0a2;
  std::vector<double>* _jStr_bcfJ_bcfAsymYPhi_v0a2;
  std::vector<double>* _jStr_bcfJ_bcfAsymYPhi2_v0a2;

  // planar flow jet shapes:
  std::vector<double>* _jStr_pfJ_pf;
  std::vector<double>* _jStr_pfJ_detST;
  std::vector<double>* _jStr_pfJ_lambdaST;

  // splittings:
  std::vector<double>* _jStr_zminJ;
  std::vector<double>* _jStr_zmaxJ;

  // angularities
  std::vector<double>* _jStr_tauJ09_a;
  std::vector<double>* _jStr_tauJ09_tau;
  std::vector<double>* _jStr_tauJ20_a;
  std::vector<double>* _jStr_tauJ20_tau;
  std::vector<double>* _jStr_tauJ40_a;
  std::vector<double>* _jStr_tauJ40_tau;

  std::vector<double>* _jStr_radEi_bcfJ1_v0a1;
  std::vector<double>* _jStr_radEi_bcfJ2_v0a1;
  std::vector<double>* _jStr_radEi_bcfJ3_v0a1;
  std::vector<double>* _jStr_radEi_bcfTJ1_v0a1;
  std::vector<double>* _jStr_radEi_bcfTJ2_v0a1;
  std::vector<double>* _jStr_radEi_bcfTJ3_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYJ1_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYJ2_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYJ3_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymPhiJ1_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymPhiJ2_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymPhiJ3_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYPhiJ1_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYPhiJ2_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYPhiJ3_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYPhi2J1_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYPhi2J2_v0a1;
  std::vector<double>* _jStr_radEi_bcfAsymYPhi2J3_v0a1;

  std::vector<double>* _jStr_radEi_bcfJ1_v0a2;
  std::vector<double>* _jStr_radEi_bcfJ2_v0a2;
  std::vector<double>* _jStr_radEi_bcfJ3_v0a2;
  std::vector<double>* _jStr_radEi_bcfTJ1_v0a2;
  std::vector<double>* _jStr_radEi_bcfTJ2_v0a2;
  std::vector<double>* _jStr_radEi_bcfTJ3_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYJ1_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYJ2_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYJ3_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymPhiJ1_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymPhiJ2_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymPhiJ3_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYPhiJ1_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYPhiJ2_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYPhiJ3_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYPhi2J1_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYPhi2J2_v0a2;
  std::vector<double>* _jStr_radEi_bcfAsymYPhi2J3_v0a2;

  std::vector<double>* _jStr_radEi_BJ1;
  std::vector<double>* _jStr_radEi_BJ2;
  std::vector<double>* _jStr_radEi_BJ3;
  std::vector<double>* _jStr_radEi_girthJ1;
  std::vector<double>* _jStr_radEi_girthJ2;
  std::vector<double>* _jStr_radEi_girthJ3;
  std::vector<double>* _jStr_radEi_girth32;
  std::vector<double>* _jStr_radEi_girth21;
  std::vector<double>* _jStr_radEi_girthAsymJ1J2;
  std::vector<double>* _jStr_radEi_alpha1;
  std::vector<double>* _jStr_radEi_alpha2;
  std::vector<double>* _jStr_radEi_alpha;
  std::vector<double>* _jStr_radEi_beta1;
  std::vector<double>* _jStr_radEi_beta2;
  std::vector<double>* _jStr_radEi_beta;
  std::vector<double>* _jStr_radEi_thetaJ1J2;
  std::vector<double>* _jStr_radEi_dipolarityInfLine;
  std::vector<double>* _jStr_radEi_dipolarityLineSeg;
  std::vector<double>* _jStr_radEi_BCF1;
  std::vector<double>* _jStr_radEi_BCF2;
  std::vector<double>* _jStr_radEi_BCF3;
  std::vector<double>* _jStr_radEi_dipolarity;

  // whole event, jets as input:
  std::vector<double>* _jStr_sphEJ_detSphericity;
  std::vector<double>* _jStr_sphEJ_detSpherocity;
  std::vector<double>* _jStr_sphEJ_sphericityLambda1;
  std::vector<double>* _jStr_sphEJ_sphericityLambda2;
  std::vector<double>* _jStr_sphEJ_sphericityLambda3;
  std::vector<double>* _jStr_sphEJ_spherocityLambda1;
  std::vector<double>* _jStr_sphEJ_spherocityLambda2;
  std::vector<double>* _jStr_sphEJ_spherocityLambda3;
  std::vector<double>* _jStr_sphEJ_circularity;
  std::vector<double>* _jStr_sphEJ_sphericity;
  std::vector<double>* _jStr_sphEJ_spherocity;
  std::vector<double>* _jStr_sphEJ_aplanarity;
  std::vector<double>* _jStr_sphEJ_aplanority;
  std::vector<double>* _jStr_sphEJ_Y;
  std::vector<double>* _jStr_sphEJ_planarity;
  std::vector<double>* _jStr_sphEJ_planority;
  std::vector<double>* _jStr_sphEJ_Dshape;
  std::vector<double>* _jStr_sphEJ_Cshape;
  std::vector<double>* _jStr_sphEJ_H2;
  std::vector<double>* _jStr_sphEJ_Fshape;
  std::vector<double>* _jStr_sphEJ_beamThrust;
  std::vector<double>* _jStr_sphEJ_G;
  std::vector<double>* _jStr_sphEJ_ST2D;
  std::vector<double>* _jStr_sphEJ_detMlin;
  std::vector<double>* _jStr_sphEJ_pft;

  // whole event, jet constituents as input:
  std::vector<double>* _jStr_sphEi_detSphericity;
  std::vector<double>* _jStr_sphEi_detSpherocity;
  std::vector<double>* _jStr_sphEi_sphericityLambda1;
  std::vector<double>* _jStr_sphEi_sphericityLambda2;
  std::vector<double>* _jStr_sphEi_sphericityLambda3;
  std::vector<double>* _jStr_sphEi_spherocityLambda1;
  std::vector<double>* _jStr_sphEi_spherocityLambda2;
  std::vector<double>* _jStr_sphEi_spherocityLambda3;
  std::vector<double>* _jStr_sphEi_circularity;
  std::vector<double>* _jStr_sphEi_sphericity;
  std::vector<double>* _jStr_sphEi_spherocity;
  std::vector<double>* _jStr_sphEi_aplanarity;
  std::vector<double>* _jStr_sphEi_aplanority;
  std::vector<double>* _jStr_sphEi_Y;
  std::vector<double>* _jStr_sphEi_planarity;
  std::vector<double>* _jStr_sphEi_planority;
  std::vector<double>* _jStr_sphEi_Dshape;
  std::vector<double>* _jStr_sphEi_Cshape;
  std::vector<double>* _jStr_sphEi_H2;
  std::vector<double>* _jStr_sphEi_Fshape;
  std::vector<double>* _jStr_sphEi_beamThrust;
  std::vector<double>* _jStr_sphEi_G;
  std::vector<double>* _jStr_sphEi_ST2D;
  std::vector<double>* _jStr_sphEi_detMlin;
  std::vector<double>* _jStr_sphEi_pft;

  // mini event, subjets as input:
  std::vector<double>* _jStr_sphJj_detSphericity;
  std::vector<double>* _jStr_sphJj_detSpherocity;
  std::vector<double>* _jStr_sphJj_sphericityLambda1;
  std::vector<double>* _jStr_sphJj_sphericityLambda2;
  std::vector<double>* _jStr_sphJj_sphericityLambda3;
  std::vector<double>* _jStr_sphJj_spherocityLambda1;
  std::vector<double>* _jStr_sphJj_spherocityLambda2;
  std::vector<double>* _jStr_sphJj_spherocityLambda3;
  std::vector<double>* _jStr_sphJj_circularity;
  std::vector<double>* _jStr_sphJj_sphericity;
  std::vector<double>* _jStr_sphJj_spherocity;
  std::vector<double>* _jStr_sphJj_aplanarity;
  std::vector<double>* _jStr_sphJj_aplanority;
  std::vector<double>* _jStr_sphJj_Y;
  std::vector<double>* _jStr_sphJj_planarity;
  std::vector<double>* _jStr_sphJj_planority;
  std::vector<double>* _jStr_sphJj_Dshape;
  std::vector<double>* _jStr_sphJj_Cshape;
  std::vector<double>* _jStr_sphJj_H2;
  std::vector<double>* _jStr_sphJj_Fshape;
  std::vector<double>* _jStr_sphJj_beamThrust;
  std::vector<double>* _jStr_sphJj_G;
  std::vector<double>* _jStr_sphJj_ST2D;
  std::vector<double>* _jStr_sphJj_detMlin;
  std::vector<double>* _jStr_sphJj_pft;

  // mini event, jet constituents as input:
  std::vector<double>* _jStr_sphJi_detSphericity;
  std::vector<double>* _jStr_sphJi_detSpherocity;
  std::vector<double>* _jStr_sphJi_sphericityLambda1;
  std::vector<double>* _jStr_sphJi_sphericityLambda2;
  std::vector<double>* _jStr_sphJi_sphericityLambda3;
  std::vector<double>* _jStr_sphJi_spherocityLambda1;
  std::vector<double>* _jStr_sphJi_spherocityLambda2;
  std::vector<double>* _jStr_sphJi_spherocityLambda3;
  std::vector<double>* _jStr_sphJi_circularity;
  std::vector<double>* _jStr_sphJi_sphericity;
  std::vector<double>* _jStr_sphJi_spherocity;
  std::vector<double>* _jStr_sphJi_aplanarity;
  std::vector<double>* _jStr_sphJi_aplanority;
  std::vector<double>* _jStr_sphJi_Y;
  std::vector<double>* _jStr_sphJi_planarity;
  std::vector<double>* _jStr_sphJi_planority;
  std::vector<double>* _jStr_sphJi_Dshape;
  std::vector<double>* _jStr_sphJi_Cshape;
  std::vector<double>* _jStr_sphJi_H2;
  std::vector<double>* _jStr_sphJi_Fshape;
  std::vector<double>* _jStr_sphJi_beamThrust;
  std::vector<double>* _jStr_sphJi_G;
  std::vector<double>* _jStr_sphJi_ST2D;
  std::vector<double>* _jStr_sphJi_detMlin;
  std::vector<double>* _jStr_sphJi_pft;

  /*
    std::vector<int>* _radialParameterEvent;
    std::vector<std::string>* _authorEvent;

    std::vector<int>* _index;
    std::vector<double>* _jetE;
    std::vector<double>* _jetP;
    std::vector<double>* _jetEta;
    std::vector<double>* _jetPhi;

    // jet shape variables
    std::vector<double>* _angularEccentricity;
    std::vector<double>* _det;
    std::vector<double>* _detST;
    std::vector<double>* _e;
    std::vector<double>* _g;
    std::vector<double>* _B;
    std::vector<double>* _logB;
    std::vector<double>* _girth;
    std::vector<double>* _jetEt;
    std::vector<double>* _jetM;
    std::vector<double>* _lambdaST;
    std::vector<double>* _orientation;
    std::vector<double>* _pf;
    std::vector<double>* _psi;
    std::vector<double>* _pullMag;
    std::vector<double>* _pullPf;
    std::vector<double>* _pullTheta;
    std::vector<double>* _ratio;
    std::vector<double>* _rho;
    std::vector<double>* _tauNeg9;
    std::vector<double>* _tauNeg20;
    std::vector<int>* _numConstituents;
    std::vector<double>* _zmaxJ;
    std::vector<double>* _zminJ;
    std::vector<int>* _radialParameter;
    std::vector<std::string>* _author;
    std::vector<double>* _R;

    std::vector<int>* _numSubjets;
    std::vector<double>* _subjetMassDemocracy;
    std::vector<double>* _subjetBetaflow;
    std::vector<double>* _tau1;
    std::vector<double>*  _tau2;
    std::vector<double>* _tau3;
    std::vector<double>* _tau2tau1;
    std::vector<double>* _tau3tau2;
    std::vector<double>* _G2;
    std::vector<int>* _nP; 
    std::vector<double>* _Rstar1; 
    std::vector<double>* _Rstar2; 
    std::vector<double>* _Rstar3; 
    std::vector<double>* _Mstar1; 
    std::vector<double>* _Mstar2; 
    std::vector<double>* _Mstar3; 
    std::vector<double>* _jetDipolarity;
    std::vector<double>* _subjetAsymmetry;
    std::vector<double>* _subjetZ;
    std::vector<double>* _subjetZ2;
    std::vector<double>* _mbar;

    // event shape variables
    std::vector<double>* _Dshape;
    std::vector<double>* _Y;
    std::vector<double>* _acoplanarity;
    std::vector<double>* _BJ1;
    std::vector<double>* _BJ2;
    std::vector<double>* _BJ3;
    std::vector<double>* _girthJ1;
    std::vector<double>* _girthJ2;
    std::vector<double>* _girthJ3;
    std::vector<double>* _girthAsymJ1J2;
    std::vector<double>* _alpha;
    std::vector<double>* _alpha1;
    std::vector<double>* _alpha2;
    std::vector<double>* _aplanarity;
    std::vector<double>* _aplanority;
    std::vector<double>* _asym;
    std::vector<double>* _asymEtaJJ;
    std::vector<double>* _asymPhiJJ;
    std::vector<double>* _asymRapidityJJ;
    std::vector<double>* _asymThetaJJ;
    std::vector<double>* _azilicityJ1;
    std::vector<double>* _azilicityJ2;
    std::vector<double>* _beta;
    std::vector<double>* _beta1;
    std::vector<double>* _beta2;
    std::vector<double>* _betaflow;
    std::vector<double>* _chi;
    std::vector<double>* _cosHelicityJ1;
    std::vector<double>* _cosHelicityJ2;
    std::vector<double>* _cosThetaJ1;
    std::vector<double>* _cosThetaJ2;
    std::vector<double>* _deltaEtaJJ;
    std::vector<double>* _deltaPhiJJ;
    std::vector<double>* _deltaRJJ;
    std::vector<double>* _deltaRJJY;
    std::vector<double>* _deltaRapidityJJ;
    std::vector<double>* _deltaRapidityXtoJ1;
    std::vector<double>* _deltaRapidityXtoJ1CM;
    std::vector<double>* _deltaRapidityXtoJ2;
    std::vector<double>* _deltaRapidityXtoJ2CM;
    std::vector<double>* _deltaThetaJJ;
    std::vector<double>* _dipolarity;
    std::vector<double>* _dipolarityInfLine;
    std::vector<double>* _dipolarityLineSeg;
    std::vector<double>* _BCF1;
    std::vector<double>* _BCF2;
    std::vector<double>* _BCF3;

    std::vector<double>* _bcfT;
    std::vector<double>* _bcf;
    std::vector<double>* _bcfAsymY;
    std::vector<double>* _bcfAsymPhi;
    std::vector<double>* _bcfAsymYPhi;
    std::vector<double>* _bcfAsymYPhi2;
    std::vector<double>* _bcfJ1;
    std::vector<double>* _bcfJ2;
    std::vector<double>* _bcfJ3;
    std::vector<double>* _bcfTJ1;
    std::vector<double>* _bcfTJ2;
    std::vector<double>* _bcfTJ3;
    std::vector<double>* _bcfAsymYJ1;
    std::vector<double>* _bcfAsymYJ2;
    std::vector<double>* _bcfAsymYJ3;
    std::vector<double>* _bcfAsymPhiJ1;
    std::vector<double>* _bcfAsymPhiJ2;
    std::vector<double>* _bcfAsymPhiJ3;
    std::vector<double>* _bcfAsymYPhiJ1;
    std::vector<double>* _bcfAsymYPhiJ2;
    std::vector<double>* _bcfAsymYPhiJ3;
    std::vector<double>* _bcfAsymYPhi2J1;
    std::vector<double>* _bcfAsymYPhi2J2;
    std::vector<double>* _bcfAsymYPhi2J3;

    std::vector<double>* _bcf2T;
    std::vector<double>* _bcf2;
    std::vector<double>* _bcf2AsymY;
    std::vector<double>* _bcf2AsymPhi;
    std::vector<double>* _bcf2AsymYPhi;
    std::vector<double>* _bcf2AsymYPhi2;
    std::vector<double>* _bcf2J1;
    std::vector<double>* _bcf2J2;
    std::vector<double>* _bcf2J3;
    std::vector<double>* _bcf2TJ1;
    std::vector<double>* _bcf2TJ2;
    std::vector<double>* _bcf2TJ3;
    std::vector<double>* _bcf2AsymYJ1;
    std::vector<double>* _bcf2AsymYJ2;
    std::vector<double>* _bcf2AsymYJ3;
    std::vector<double>* _bcf2AsymPhiJ1;
    std::vector<double>* _bcf2AsymPhiJ2;
    std::vector<double>* _bcf2AsymPhiJ3;
    std::vector<double>* _bcf2AsymYPhiJ1;
    std::vector<double>* _bcf2AsymYPhiJ2;
    std::vector<double>* _bcf2AsymYPhiJ3;
    std::vector<double>* _bcf2AsymYPhi2J1;
    std::vector<double>* _bcf2AsymYPhi2J2;
    std::vector<double>* _bcf2AsymYPhi2J3;

    std::vector<double>* _gamma;
    std::vector<double>* _gamma1;
    std::vector<double>* _gamma2;
    std::vector<double>* _helicityJ1;
    std::vector<double>* _helicityJ2;
    std::vector<double>* _jetSumE;
    std::vector<double>* _massDemocracy;
    std::vector<double>* _mjj;
    std::vector<double>* _sigmaEtaJJ;
    std::vector<double>* _sigmaPhiJJ;
    std::vector<double>* _sigmaRapidityJJ;
    std::vector<double>* _sigmaThetaJJ;
    std::vector<double>* _sphericity;
    std::vector<double>* _spherocity;
    std::vector<double>* _tJetSep;
    std::vector<double>* _tJetSepInv;
    std::vector<double>* _theta;
    std::vector<double>* _thetaJ1J2;
    std::vector<double>* _thetaStar;
    std::vector<double>* _twist;
    std::vector<double>* _twistY;
    std::vector<double>* _yB;
    std::vector<double>* _yStar;
    std::vector<double>* _zmaxAllJets;
    std::vector<double>* _zmaxJ1J2;
    std::vector<double>* _zminAllJets;
    std::vector<double>* _zminJ1J2;
    std::vector<double>* _zminJ1J2Eta;
    std::vector<double>* _zminJ1J2Phi;
    std::vector<double>* _zminJ1J2Rapidity;
    std::vector<double>* _zminJ1J2Theta;
    int _numJets;

  */
  std::vector<double>* _YFlip12;
  std::vector<double>* _YFlip23;
  std::vector<double>* _EM_FRACTION_CLUSTER;
  std::vector<double>* _ELLIPTICAREA;
  std::vector<double>* _AMAREA;
  std::vector<double>* _HULL_LENGTH;
  std::vector<double>* _HULL_AREA;
  std::vector<double>* _EM_FRACTION_MCTRUTH;
  std::vector<double>* _LowEtConstituentsFrac;
  std::vector<double>* _JetEccentricity;
  std::vector<double>* _DRtoReco;
  std::vector<double>* _PtNearest;
  std::vector<double>* _WIDTH;
  std::vector<double>* _rbb;
  std::vector<double>* _rfilt;

  std::vector<int>* _nConst;
  std::vector<std::vector<double> >* _const_et;
  std::vector<std::vector<double> >* _const_rapidity;
  std::vector<std::vector<double> >* _const_phi;
  std::vector<std::vector<double> >* _const_Ret;
  std::vector<std::vector<double> >* _const_Drapidity;
  std::vector<std::vector<double> >* _const_Dphi;

  std::vector<double>* _jvf;
  std::vector<int>* _ntrk;
  std::vector<double>* _trkpt;
  
  std::vector<long long int>* _QGtag; // from MC truth
  std::vector<long long int>* _QGtag2; // from MC truth
  std::vector<int>* _jetTrueFlavour; // from MC truth

  // Key 	Nature of weight 	Tagger info
  // e.g. w_keyname    LR    tagger information 
  std::vector<double>* w_cmb; // DEFAULT TAGGER, IP3D+SV1
  std::vector<double>* w_TrackCounting2D; // - based on cuts on the first two most significant tracks
  std::vector<double>* w_JetProb; // probability based only on the IP resolution function of prompt tracks
  std::vector<double>* w_IP1D; // likelihood ratio (LR) based on the ImpactParameter significances of all tracks, only longitudinal
  std::vector<double>* w_IP2D; // LR based on the ImpactParameter significances of all tracks, only transverse
  std::vector<double>* w_IP3D; // LR based on the ImpactParameter significances of all tracks, transverse+longitudinal
  std::vector<double>* w_SV0; // distance/error Secondary Vertex Tagger based on BTagVrtSec vertex finder, normalized and signed distance (in 3D) between primary vertex and inclusive secondary vertex
  std::vector<double>* w_SV1; // LR Secondary Vertex Tagger based on BTagVrtSec vertex finder, 2d+1d likelihood
  std::vector<double>* w_SV2; // LR Secondary Vertex Tagger based on BTagVrtSec vertex finder, 3d likelihood
  std::vector<double>* w_BaselineTagger; // LR IP3D+ SV1, based on their sum of the weights (this is also the default you get with just getFlavourTagWeight() )
  std::vector<double>* w_JetFitterTag; // LR Secondary Vertex Tagger based on the JetFitter vertex finder, uses projective 1d likelihood with many different categories
  std::vector<double>* w_JetFitterCOMB; // LR IP3D+JetFitterTag, based on their sum of the weights
  std::vector<double>* w_JetFitterTagNN; // NN Secondary Vertex Tagger based on the JetFitter vertex finder, based on Neural Network trained with JetNet
  std::vector<double>* w_JetFitterCOMBNN; // NN IP3D+JetFitterTagNN, based on a separate Neural Network trained with JetNet
  std::vector<double>* w_SoftMuonTag; // LR soft muon Tagger
  std::vector<double>* w_SoftElectronTag; // LR soft electron Tagger 

};

#endif // MY_ANALYSIS_H

