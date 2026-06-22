// system include files
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Common/interface/ValidHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"

#include <TH2.h>

class SimTrackSimVertexDumperV2 : public edm::one::EDAnalyzer<> {
public:
  explicit SimTrackSimVertexDumperV2(const edm::ParameterSet&);
  ~SimTrackSimVertexDumperV2() override {}

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override {}
  void endJob() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   TH2F* h2_hits_xy;

private:
  edm::EDGetTokenT<edm::HepMCProduct> hepmcToken_;
  edm::EDGetTokenT<edm::SimTrackContainer> simTrackToken_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVertexToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer> HcalToken_;
  bool dumpHepMC_;
  std::map<int,std::string> _map_procType;
};

SimTrackSimVertexDumperV2::SimTrackSimVertexDumperV2(const edm::ParameterSet& iConfig)
    : hepmcToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("moduleLabelHepMC"))),
      simTrackToken_(consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("moduleLabelTk"))),
      simVertexToken_(consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("moduleLabelVtx"))),
      HcalToken_(consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("HcalHits"))),
      dumpHepMC_(iConfig.getUntrackedParameter<bool>("dumpHepMC")) {
        _map_procType[0]   = "Primary generator";
        _map_procType[91]  = "Transportation";
        _map_procType[92]  = "Coulpled transportation";
        _map_procType[1]   = "Coulomb scattering";
        _map_procType[2]   = "Ionisation";
        _map_procType[3]   = "Bremsstrahlung";
        _map_procType[4]   = "e+e- pair production";
        _map_procType[5]   = "Annihilation in 2 gamma";
        _map_procType[6]   = "Annihilation in mu+mu- ";
        _map_procType[7]   = "Annihilation in hadrons";
        _map_procType[10]  = "Multiple scattering";
        _map_procType[12]  = "Photoelectric";
        _map_procType[13]  = "Compton scattering";
        _map_procType[14]  = "Gamma conversion in e+e-";
        _map_procType[15]  = "Gamma conversion in mu+mu- ";
        _map_procType[23]  = "Synchrotron radiation";
        _map_procType[111] = "Hadron elastic scattering";
        _map_procType[121] = "Hadron inelastic";
        _map_procType[131] = "Neutron capture";
        _map_procType[141] = "Neutron fission";
        _map_procType[151] = "Hadron stopping at rest";
        _map_procType[201] = "Decay";
        _map_procType[202] = "Decay with spin";
        _map_procType[203] = "Pion decay with spin";
        _map_procType[204] = "Radioactive decay";
        _map_procType[205] = "Unknown decay";
        _map_procType[206] = "External decay";
        _map_procType[401] = "Step limiter";
        _map_procType[403] = "Neutron killer";
        _map_procType[301] = "GFlash";


        // The root tuple
        // usesResource(TFileService::kSharedResource);
        // edm::Service<TFileService> fs;
        // h2_hits_xy = new TH2F("h2_hits_xy","hitx",600, 0., 600.,300, 0., 300.);
      }

//
// member functions
//

// ------------ method called to produce the data  ------------
void SimTrackSimVertexDumperV2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace HepMC;

  std::vector<SimTrack> theSimTracks;
  std::vector<SimVertex> theSimVertexes;

  auto MCEvt = edm::makeValid(iEvent.getHandle(hepmcToken_));
  const HepMC::GenEvent* evt = MCEvt->GetEvent();

  auto SimTk = edm::makeValid(iEvent.getHandle(simTrackToken_));
  auto SimVtx = edm::makeValid(iEvent.getHandle(simVertexToken_));

  theSimTracks.insert(theSimTracks.end(), SimTk->begin(), SimTk->end());
  theSimVertexes.insert(theSimVertexes.end(), SimVtx->begin(), SimVtx->end());

 
  edm::LogPrint("DumpTkVtx") << "\n SimVertex / SimTrack structure dump \n";
  edm::LogPrint("DumpTkVtx") << " SimVertex in the event = " << theSimVertexes.size();
  edm::LogPrint("DumpTkVtx") << " SimTracks in the event = " << theSimTracks.size();
  edm::LogPrint("DumpTkVtx") << "\n";

  for (unsigned int isimvtx = 0; isimvtx < theSimVertexes.size(); isimvtx++) {
    edm::LogPrint("DumpTkVtx") << "SimVertex " << isimvtx << " = " << theSimVertexes[isimvtx] << " , proc = " << _map_procType[theSimVertexes[isimvtx].processType()] << "\n";

    for (unsigned int isimtk = 0; isimtk < theSimTracks.size(); isimtk++) {
      if (theSimTracks[isimtk].vertIndex() >= 0 && std::abs(theSimTracks[isimtk].vertIndex()) == (int)isimvtx) {
        edm::LogPrint("DumpTkVtx") << "  SimTrack " << isimtk << " = " << theSimTracks[isimtk]
                                   << " Track Id = " << theSimTracks[isimtk].trackId();

        // for debugging purposes
        if (dumpHepMC_) {
          if (theSimTracks[isimtk].genpartIndex() != -1) {
            HepMC::GenParticle* part = evt->barcode_to_particle(theSimTracks[isimtk].genpartIndex());
            if (part) {
              edm::LogPrint("DumpTkVtx") << "  ---> Corresponding to HepMC particle " << *part;
            } else {
              edm::LogPrint("DumpTkVtx") << " ---> Corresponding HepMC particle to barcode "
                                         << theSimTracks[isimtk].genpartIndex() << " not in selected event ";
            }
          }
        }
      }
    }
    edm::LogPrint("DumpTkVtx") << "\n";
  }

  for (std::vector<SimTrack>::iterator isimtk = theSimTracks.begin(); isimtk != theSimTracks.end(); ++isimtk) {
    if (isimtk->noVertex()) {
      edm::LogPrint("DumpTkVtx") << "SimTrack without an associated Vertex = " << *isimtk;
    }
  }

  edm::LogPrint("DumpTkVtx") << "\n";


  std::vector<std::pair<int, std::string> > theCaloComposition;

  int oldsize = 0;
  std::vector<PCaloHit> theCaloHits;
  auto HcalHits = iEvent.getHandle(HcalToken_);

  if (HcalHits.isValid()) {
    theCaloHits.insert(theCaloHits.end(), HcalHits->begin(), HcalHits->end());
    std::pair<int, std::string> label19(theCaloHits.size() - oldsize, "HcalHits");
    oldsize = theCaloHits.size();
    theCaloComposition.push_back(label19);
  }

  edm::LogPrint("SimHitCaloHitDumper") << "\n Calorimeter Hits in the event = " << theCaloHits.size();
  edm::LogPrint("SimHitCaloHitDumper") << "\n";

  int nhit = 0;
  for (std::vector<std::pair<int, std::string> >::iterator icoll = theCaloComposition.begin(); icoll != theCaloComposition.end(); ++icoll) {
    edm::LogPrint("SimHitCaloHitDumper") << "\n";
    edm::LogPrint("SimHitCaloHitDumper") << (*icoll).second << " hits in the event = " << (*icoll).first;
    edm::LogPrint("SimHitCaloHitDumper") << "\n";
    for (int ihit = 0; ihit < (*icoll).first; ++ihit) {
      edm::LogPrint("SimHitCaloHitDumper") << theCaloHits[nhit];
      nhit++;
    }
  }

  edm::LogPrint("DumpTkVtx") << "\n";

  // for (unsigned int isimvtx = 0; isimvtx < theSimVertexes.size(); isimvtx++) {
  //    h2_hits_xy->Fill(theSimVertexes[isimvtx])
  // }

  // int oldsize = 0;
  // std::vector<PCaloHit> theCaloHits;
  // auto HcalHits = iEvent.getHandle(HcalToken_);

  // if (HcalHits.isValid()) {
  //   theCaloHits.insert(theCaloHits.end(), HcalHits->begin(), HcalHits->end());
  //   std::pair<int, std::string> label19(theCaloHits.size() - oldsize, "HcalHits");
  //   oldsize = theCaloHits.size();
  //   theCaloComposition.push_back(label19);
  // }
  // int nhit = 0;
  // for (std::vector<std::pair<int, std::string> >::iterator icoll = theCaloComposition.begin(); icoll != theCaloComposition.end(); ++icoll) {
  //   for (int ihit = 0; ihit < (*icoll).first; ++ihit) {
  //     theCaloHits[nhit];
  //   }
  // }

  return;
}

void SimTrackSimVertexDumperV2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("moduleLabelHepMC", edm::InputTag("generatorSmeared"))
      ->setComment("Input generated HepMC event after vtx smearing");
  desc.add<edm::InputTag>("moduleLabelTk", edm::InputTag("g4SimHits"))
      ->setComment("Module for input SimTrack collection");
  desc.add<edm::InputTag>("moduleLabelVtx", edm::InputTag("g4SimHits"))
      ->setComment("Module for input SimVertex collection");
  desc.add<edm::InputTag>("HcalHits", edm::InputTag("g4SimHits"))
      ->setComment("Module for HcalHits collection");
  desc.addUntracked<bool>("dumpHepMC", false);
  descriptions.add("simTrackSimVertexDumper", desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(SimTrackSimVertexDumperV2);