#include <string>
#include <vector>
#include <list>
#include <map>

#include "CBNT_Utils/CBNT_AthenaAwareBase.h"

#include "GaudiKernel/ITHistSvc.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrkToolInterfaces/ITrackSelectorTool.h"
#include "TrkVertexAnalysisUtils/V0Tools.h"

#include "TrigDecisionTool/TrigDecisionTool.h"

class StoreGateSvc;

/////////////////////////////////////////////////////////////////////////////
class RealDataSkeleton : public CBNT_AthenaAwareBase {

public:
  RealDataSkeleton (const std::string& name, ISvcLocator* pSvcLocator);
  ~RealDataSkeleton() {};
 
   virtual StatusCode CBNT_initializeBeforeEventLoop();
   virtual StatusCode CBNT_initialize();
   virtual StatusCode CBNT_finalize();
   virtual StatusCode CBNT_execute();
   virtual StatusCode CBNT_clear(); 

private:

  // User quantites, set in jobOpts
  std::string m_trackCollName; // Reconstructed tracks
  std::string m_vertexCollName; // Primary  vertices
  std::string m_v0CollName; // V0 collection
  std::string m_muonCollectionKey; // Name of muon container

  // Event info
  std::vector<int>* m_eventNum;
  std::vector<int>* m_runNum; 
  std::vector<int>* m_lumiBlock;
  std::vector<int>* m_timeStamp;
	
  // Track selection
  ToolHandle < Trk::ITrackSelectorTool > m_trkSelector;
  // V0 Tools
  ToolHandle<Trk::V0Tools> m_V0Tools;

  ToolHandle<Trig::TrigDecisionTool> m_trigDec;

  bool m_doTrigger;
  std::vector<std::string> m_triggerChains;
  std::map<std::string,int> m_triggersPassed;

};

