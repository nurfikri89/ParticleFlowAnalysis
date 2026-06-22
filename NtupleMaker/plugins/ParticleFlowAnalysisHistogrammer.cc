#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "ParticleFlowAnalysis/NtupleMaker/plugins/ParticleFlowAnalysisHistogrammer.h"
#include "ParticleFlowAnalysis/NtupleMaker/include/SimHitInfo.h"
#include "ParticleFlowAnalysis/NtupleMaker/include/HBHETowerInfo.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Calibration/IsolatedParticles/interface/CaloPropagateTrack.h"

#include <TROOT.h>
#include <TVector3.h>
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;
using namespace edm;
using namespace reco;


ParticleFlowAnalysisHistogrammer::ParticleFlowAnalysisHistogrammer(const edm::ParameterSet& iConfig) {

  tokengenParticles_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles"));
  tokengenVertexXYZ_ = consumes<math::XYZPointF>(iConfig.getParameter<InputTag>("genVertex"));

  g4SimHitsPCaloHitsEBToken_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("g4SimHitsPCaloHitsEB"));
  g4SimHitsPCaloHitsEEToken_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("g4SimHitsPCaloHitsEE"));
  g4SimHitsPCaloHitsPSToken_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("g4SimHitsPCaloHitsPS"));
  g4SimHitsPCaloHitsHCALToken_ = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("g4SimHitsPCaloHitsHCAL"));
  // simCaloHitsHFToken_   = consumes<edm::PCaloHitContainer>(iConfig.getParameter<edm::InputTag>("simCaloHitsHF"));

  geometryToken_      = esConsumes<CaloGeometry, CaloGeometryRecord>();
  caloTopologyToken_  = esConsumes<CaloTopology, CaloTopologyRecord>();
  hcalTopologyToken_  = esConsumes<HcalTopology, HcalRecNumberingRecord>();
  hcalDDDrecToken_    = esConsumes<HcalDDDRecConstants, HcalRecNumberingRecord>();
  magneticFieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  pdtToken_           = esConsumes<HepPDT::ParticleDataTable, PDTRecord>();
  ecalPFRecHitThresholdsToken_= esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd>();
  hcalPFCutsToken_    = esConsumes<HcalPFCuts, HcalPFCutsRcd>(edm::ESInputTag("", "withTopo"));
  hcalDBToken_        = esConsumes<HcalDbService, HcalDbRecord>();
  hcalRespCorrsToken_ = esConsumes<HcalRespCorrs, HcalRespCorrsRcd>(edm::ESInputTag("", "withTopo"));

  input_simtrack_token_ = consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits"));
  input_simvertex_token_ = consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits"));

  //
  //
  //
  theParameterMap = new HcalSimParameterMap(iConfig);

  // The root tuple
  usesResource(TFileService::kSharedResource);

  bookTH2F("h2_genPartEta_vs_genPartE", "", 70, -3.5, 3.5, 1000, 0.f, 100.);

  bookProfile2D("hp2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy", "", 60, -30, 30, 7, 1, 8);

  bookTH2F("h2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);
  bookTH2F("h2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_depth", "", 60, -30, 30, 7, 1, 8);

  bookTH2F("h2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);
  bookTH2F("h2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_energy", "", 60, -30, 30, 1000, 0.f, 10.);

  bookTH2F("h2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);
  bookTH2F("h2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_energy_v2", "", 60, -30, 30, 300, 0.f, 1500.);

  // bookTH2F("h2_z_vs_y",    "", 125, 0., 250., 100, 0., 500.);
  // bookTH2F("h2_z_vs_y_v2", "", 125, 0., 250., 100, 0., 500.);

  // bookTH2F("h2_rho_vs_z",     "", 125, 0., 250., 100, 0., 500.);
  // bookTH2F("h2_rho_vs_z_v2",  "", 125, 0., 250., 100, 0., 500.);

  bookTH2F("h2_z_vs_y",    "", 1400, -700., 700., 1000, -300., 700.);
  bookTH2F("h2_z_vs_y_v2", "", 1400, -700., 700., 1000, -300., 700.);

  bookTH2F("h2_z_vs_rho",     "", 1400, -700., 700., 1000, 0., 700.);
  bookTH2F("h2_z_vs_rho_v2",  "", 1400, -700., 700., 1000, 0., 700.);

  // bookTH3F("h3_eventIdx_vs_z_vs_y",    "", 100, 1, 101, 125, 0., 250., 100, 0., 500.);
  // bookTH3F("h3_eventIdx_vs_z_vs_y_v2", "", 100, 1, 101, 125, 0., 250., 100, 0., 500.);

  eventCounter_=0;

  bookProfile2D("hp2_depth1_genPartEta_vs_hcalsimhit_ieta_vs_energy", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth2_genPartEta_vs_hcalsimhit_ieta_vs_energy", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth3_genPartEta_vs_hcalsimhit_ieta_vs_energy", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth4_genPartEta_vs_hcalsimhit_ieta_vs_energy", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth5_genPartEta_vs_hcalsimhit_ieta_vs_energy", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth6_genPartEta_vs_hcalsimhit_ieta_vs_energy", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth7_genPartEta_vs_hcalsimhit_ieta_vs_energy", "", 30, 0., 3., 60, -30, 30);

  bookProfile2D("hp2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);
  bookProfile2D("hp2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit", "", 60, -30, 30, 7, 1, 8);

  bookProfile2D("hp2_depth1_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth2_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth3_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth4_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth5_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth6_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit", "", 30, 0., 3., 60, -30, 30);
  bookProfile2D("hp2_depth7_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit", "", 30, 0., 3., 60, -30, 30);
}

ParticleFlowAnalysisHistogrammer::~ParticleFlowAnalysisHistogrammer(){}

void ParticleFlowAnalysisHistogrammer::beginRun(const edm::Run& run, const edm::EventSetup & es) {}
void ParticleFlowAnalysisHistogrammer::endRun(const edm::Run& run, const edm::EventSetup & es) {}

void ParticleFlowAnalysisHistogrammer::analyze(const Event& iEvent, const EventSetup& iSetup) {
  LogDebug("ParticleFlowAnalysisHistogrammer") <<"START event: "<<iEvent.id().event() <<" in run "<<iEvent.id().run()<<endl;


  run  = iEvent.id().run();
  evt  = iEvent.id().event();
  lumiBlock = iEvent.id().luminosityBlock();
  time = iEvent.time();

  orun = (size_t)run;
  oevt = (size_t)evt;
  olumiBlock = (size_t)lumiBlock;
  otime = (size_t)((iEvent.time().value())>>32);

  const HepPDT::ParticleDataTable* pdt = &iSetup.getData(pdtToken_);
  const MagneticField* bField = &iSetup.getData(magneticFieldToken_);
  const CaloGeometry* geometry = &iSetup.getData(geometryToken_);
  const CaloTopology* caloTopology = &iSetup.getData(caloTopologyToken_);
  const HcalTopology* theHBHETopology = &iSetup.getData(hcalTopologyToken_);
  const HcalDDDRecConstants* hcaloConst = &iSetup.getData(hcalDDDrecToken_);
  const EcalPFRecHitThresholds* ecalPFCuts = &iSetup.getData(ecalPFRecHitThresholdsToken_);
  const HcalPFCuts* hcalPFCuts = &iSetup.getData(hcalPFCutsToken_);
  const HcalRespCorrs* hcalRespCorrs = &iSetup.getData(hcalRespCorrsToken_);

  const HcalGeometry* cellGeometryHB = dynamic_cast<const HcalGeometry *>(geometry->getSubdetectorGeometry(DetId::Hcal, HcalSubdetector::HcalBarrel));
  const HcalGeometry* cellGeometryHE = dynamic_cast<const HcalGeometry *>(geometry->getSubdetectorGeometry(DetId::Hcal, HcalSubdetector::HcalEndcap));
  const HcalGeometry* cellGeometryHF = dynamic_cast<const HcalGeometry *>(geometry->getSubdetectorGeometry(DetId::Hcal, HcalSubdetector::HcalForward));

  //==========================================
  //
  // Get gen vertex
  //
  //==========================================
  Handle<math::XYZPointF> genVertexHandle;
  iEvent.getByToken(tokengenVertexXYZ_, genVertexHandle);
  auto genVertex = genVertexHandle.product();
  genVertex_x = genVertex->X();
  genVertex_y = genVertex->Y();
  genVertex_z = genVertex->Z();

  //===================================================================================================
  //
  // Since we are analysis single particle guns, the container should *always* have EXACTLY one entry.
  //
  //===================================================================================================
  Handle<GenParticleCollection> genParticlesHandle;
  iEvent.getByToken(tokengenParticles_, genParticlesHandle);
  auto genParticles = genParticlesHandle.product();

  if (genParticles->size() == 0){
    std::cout << "No genParts" << std::endl;
    return;
  }
  auto p = genParticles->at(0);
  genPart_pt     = p.pt();
  genPart_eta    = p.eta();
  genPart_phi    = p.phi();
  genPart_mass   = p.mass();
  genPart_energy = p.energy();
  genPart_p      = p.p();
  genPart_pdgId  = p.pdgId();
  genPart_charge = p.charge();

  m_Histos2D["h2_genPartEta_vs_genPartE"]->Fill(genPart_eta,genPart_energy);

  // const reco::GenParticle& genPart = (*(genPartIDs[0].trkItr));

  bool fillYZPlot=true;
  if (abs(genPart_eta) < 1.5){
    eventCounter_++;
    fillYZPlot=true;
  }

  //==========================================
  //
  // G4
  //
  //==========================================
  // std::cout<< "New Event"  << std::endl;

  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMCTruth
  edm::Handle<SimVertexContainer> simVertices;
  iEvent.getByToken(input_simvertex_token_, simVertices);
  const edm::SimVertexContainer* mSimVertexCollection = &*simVertices;
  // const SimVertex* vertex = &((*getSimVertexCollection())[track->vertIndex()]);

  edm::Handle<SimTrackContainer> simTracks;
  iEvent.getByToken(input_simtrack_token_, simTracks);
  const edm::SimTrackContainer*  mSimTrackCollection = &*simTracks;

  for (unsigned i = 0; i < mSimTrackCollection->size(); ++i) {
    auto track = mSimTrackCollection->at(i);
    // std::cout << i << " / " << mSimTrackCollection->size() << " | " << track.vertIndex() << " | genpartIndex = " << track.genpartIndex() << std::endl;

    if (fillYZPlot) {
      m_Histos2D["h2_z_vs_y_v2"]->Fill(track.getPositionAtBoundary().z(),track.getPositionAtBoundary().y());
    }

    if (fillYZPlot) {
      m_Histos2D["h2_z_vs_rho_v2"]->Fill(track.getPositionAtBoundary().z(),track.getPositionAtBoundary().Rho());
    }

    if ((eventCounter_ < 100) && fillYZPlot) {
      // m_Histos3D["h3_eventIdx_vs_z_vs_y_v2"]->Fill(eventCounter_,track.getPositionAtBoundary().z(),track.getPositionAtBoundary().y());
    }

    if (track.vertIndex() >= 0 && track.vertIndex() < (int) mSimVertexCollection->size()){
      // std::cout << " | " << mSimVertexCollection->at(track.vertIndex()) << " | procType = " << mSimVertexCollection->at(track.vertIndex()).processType() << std::endl;

      auto& thisVertex = mSimVertexCollection->at(track.vertIndex());
      if ((eventCounter_ < 100) && fillYZPlot) {
        // m_Histos3D["h3_eventIdx_vs_z_vs_y"]->Fill(eventCounter_,thisVertex.position().z(),thisVertex.position().y());
      }
      if (fillYZPlot) {
        m_Histos2D["h2_z_vs_y"]->Fill(thisVertex.position().z(),thisVertex.position().y());
      }

    if (fillYZPlot) {
      m_Histos2D["h2_z_vs_rho"]->Fill(thisVertex.position().z(),thisVertex.position().Rho());
    }


    }
  else{      std::cout << " | No Vertex"  << std::endl;
    }
  }

  // edm::Handle<edm::PCaloHitContainer> g4SimHitsPCaloHitsHCALHandle;
  // iEvent.getByToken(g4SimHitsPCaloHitsHCALToken_, g4SimHitsPCaloHitsHCALHandle);
  // auto g4SimHitsPCaloHitsHCAL = g4SimHitsPCaloHitsHCALHandle.product();
  // size_t n_g4SimHitsPCaloHitsHCAL = g4SimHitsPCaloHitsHCAL->size();

  // for (size_t idx = 0; idx < n_g4SimHitsPCaloHitsHCAL; idx++){
  //   const PCaloHit& caloHit = g4SimHitsPCaloHitsHCAL->at(idx);
  //   DetId newid = HcalHitRelabeller::relabel(caloHit.id(), hcaloConst);
  //   HcalDetId theHcalDetId(newid);
  //   for (unsigned i = 0; i < mSimTrackCollection->size(); ++i) {
  //     auto track = mSimTrackCollection->at(i);
  //     if ((int) track.trackId() == caloHit.geantTrackId()){
  //       if (track.vertIndex() >= 0 && track.vertIndex() < (int) mSimVertexCollection->size()){
  //         std::cout << " | " << mSimVertexCollection->at(track.vertIndex()) << " | procType = " << mSimVertexCollection->at(track.vertIndex()).processType() ;
  //         std::cout << ", caloHit.energyEM()  = " << caloHit.energyEM();
  //         std::cout << ", caloHit.energyHad() = " << caloHit.energyHad();
  //         std::cout << ", ieta = " << theHcalDetId.ieta();
  //         std::cout << ", iphi = " << theHcalDetId.iphi() << std::endl;
  //       }else{
  //         std::cout << " | No Vertex"  << std::endl;
  //       }
  //     }
  //   }
  // }

  //==========================================
  //
  // PCaloHit: HCAL
  //
  //==========================================
  std::map<unsigned int, SimHitInfo> map_simHitInfoHBHE;
  std::map<int, HBHETowerInfo> map_hbheToweInfoHBHE;

  for (int ieta=-29; ieta < 30;ieta++){
    map_hbheToweInfoHBHE[ieta] = HBHETowerInfo();
    map_hbheToweInfoHBHE[ieta].energy    = 0.f;
    map_hbheToweInfoHBHE[ieta].energyEM  = 0.f;
    map_hbheToweInfoHBHE[ieta].energyHAD = 0.f;
    map_hbheToweInfoHBHE[ieta].nCaloHits = 0;
    map_hbheToweInfoHBHE[ieta].energy_25ns    = 0.f;
    map_hbheToweInfoHBHE[ieta].energyEM_25ns  = 0.f;
    map_hbheToweInfoHBHE[ieta].energyHAD_25ns = 0.f;
    map_hbheToweInfoHBHE[ieta].nCaloHits_25ns = 0;
  }

  auto hcalCond = &iSetup.getData(hcalDBToken_);
  theParameterMap->setDbService(hcalCond);

  //
  //TEMP
  edm::Handle<edm::PCaloHitContainer> g4SimHitsPCaloHitsHCALHandle;
  iEvent.getByToken(g4SimHitsPCaloHitsHCALToken_, g4SimHitsPCaloHitsHCALHandle);
  auto g4SimHitsPCaloHitsHCAL = g4SimHitsPCaloHitsHCALHandle.product();
  size_t n_g4SimHitsPCaloHitsHCAL = g4SimHitsPCaloHitsHCAL->size();


  for (size_t idx = 0; idx < n_g4SimHitsPCaloHitsHCAL; idx++){
    const PCaloHit& caloHit = g4SimHitsPCaloHitsHCAL->at(idx);
   // std::cout << caloHit.geantTrackId() << std::endl;
    //   "sampFactor = " << samplingFactor <<
    //   " fCtoGeV = " << fCtoGeV <<
    //   " peToAnalog = " << photoelectronsToAnalog <<
    //   " simHitToPE = " << simHitToPhotoelectrons << std::endl;
    DetId newid = HcalHitRelabeller::relabel(caloHit.id(), hcaloConst);
    HcalDetId theHcalDetId(newid);
    // double eta = -9.;
    // double phi = -9.;
    // if (newid.subdetId() == HcalSubdetector::HcalBarrel){
    //   eta = cellGeometryHB->getPosition(theHcalDetId).eta();
    //   phi = cellGeometryHB->getPosition(theHcalDetId).phi();
    // }else if(newid.subdetId() == HcalSubdetector::HcalEndcap){
    //   eta = cellGeometryHE->getPosition(theHcalDetId).eta();
    //   phi = cellGeometryHE->getPosition(theHcalDetId).phi();
    // }
    //

    // std::cout << caloHit.geantTrackId() << " / mSimTrackCollection->size() = " << mSimTrackCollection->size() << std::endl;

    // TODO: Can try these examples to retrieve more parameters used for HCAL response simulation
    //
    // if (newid.subdetId() == HcalSubdetector::HcalBarrel || newid.subdetId() == HcalSubdetector::HcalEndcap){
    //   const HcalSimParameters& pars = dynamic_cast<const HcalSimParameters&>(theParameterMap->simParameters(theHcalDetId));
    //   float samplingFactor = pars.samplingFactor(theHcalDetId);
    //   float fCtoGeV = pars.fCtoGeV(theHcalDetId);
    //   float photoelectronsToAnalog = pars.photoelectronsToAnalog(theHcalDetId);
    //   float simHitToPhotoelectrons = pars.simHitToPhotoelectrons(theHcalDetId);
    //   std::cout <<
    //   "sampFactor = " << samplingFactor <<
    //   " fCtoGeV = " << fCtoGeV <<
    //   " peToAnalog = " << photoelectronsToAnalog <<
    //   " simHitToPE = " << simHitToPhotoelectrons << std::endl;
    // }
    if (newid.subdetId() == HcalSubdetector::HcalBarrel || newid.subdetId() == HcalSubdetector::HcalEndcap ){
      const HcalSimParameters& pars = dynamic_cast<const HcalSimParameters&>(theParameterMap->simParameters(theHcalDetId));
      float samplingFactor = pars.samplingFactor(theHcalDetId);

      auto it = map_simHitInfoHBHE.find(caloHit.id());

      // DetId not in map yet, so lets set SimHitInfo for this DetId.
      if (it == map_simHitInfoHBHE.end()) {
        map_simHitInfoHBHE[caloHit.id()] = SimHitInfo();
        map_simHitInfoHBHE[caloHit.id()].energy    = 0.f;
        map_simHitInfoHBHE[caloHit.id()].energyEM  = 0.f;
        map_simHitInfoHBHE[caloHit.id()].energyHAD = 0.f;
        map_simHitInfoHBHE[caloHit.id()].nCaloHits = 0;
        map_simHitInfoHBHE[caloHit.id()].energy_25ns    = 0.f;
        map_simHitInfoHBHE[caloHit.id()].energyEM_25ns  = 0.f;
        map_simHitInfoHBHE[caloHit.id()].energyHAD_25ns = 0.f;
        map_simHitInfoHBHE[caloHit.id()].nCaloHits_25ns = 0;
      }
      map_simHitInfoHBHE[caloHit.id()].energy    += caloHit.energy();
      map_simHitInfoHBHE[caloHit.id()].energyEM  += caloHit.energyEM();
      map_simHitInfoHBHE[caloHit.id()].energyHAD += caloHit.energyHad();
      map_simHitInfoHBHE[caloHit.id()].nCaloHits += 1;
      if (caloHit.time() <= 25){
        map_simHitInfoHBHE[caloHit.id()].energy_25ns    += caloHit.energy();
        map_simHitInfoHBHE[caloHit.id()].energyEM_25ns  += caloHit.energyEM();
        map_simHitInfoHBHE[caloHit.id()].energyHAD_25ns += caloHit.energyHad();
        map_simHitInfoHBHE[caloHit.id()].nCaloHits_25ns += 1;
      }

      int ieta = theHcalDetId.ieta();
      map_hbheToweInfoHBHE[ieta].energy    += samplingFactor * caloHit.energy();
      map_hbheToweInfoHBHE[ieta].energyEM  += samplingFactor * caloHit.energyEM();
      map_hbheToweInfoHBHE[ieta].energyHAD += samplingFactor * caloHit.energyHad();
      map_hbheToweInfoHBHE[ieta].nCaloHits += 1;
      if (caloHit.time() <= 25){
        map_hbheToweInfoHBHE[ieta].energy_25ns    += samplingFactor * caloHit.energy();
        map_hbheToweInfoHBHE[ieta].energyEM_25ns  += samplingFactor * caloHit.energyEM();
        map_hbheToweInfoHBHE[ieta].energyHAD_25ns += samplingFactor * caloHit.energyHad();
        map_hbheToweInfoHBHE[ieta].nCaloHits_25ns += 1;
      }
    }
  }

  for (int ieta=-29; ieta < 30;ieta++){
    // if (map_hbheToweInfoHBHE[ieta].energy > 0.1f){
      if      (genPart_eta >= 0.0f && genPart_eta < 0.2f) m_Histos2D["h2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.2f && genPart_eta < 0.4f) m_Histos2D["h2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.4f && genPart_eta < 0.6f) m_Histos2D["h2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.6f && genPart_eta < 0.8f) m_Histos2D["h2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.8f && genPart_eta < 1.0f) m_Histos2D["h2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.0f && genPart_eta < 1.2f) m_Histos2D["h2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.2f && genPart_eta < 1.4f) m_Histos2D["h2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.4f && genPart_eta < 1.6f) m_Histos2D["h2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.6f && genPart_eta < 1.8f) m_Histos2D["h2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.8f && genPart_eta < 2.0f) m_Histos2D["h2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.0f && genPart_eta < 2.2f) m_Histos2D["h2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.2f && genPart_eta < 2.4f) m_Histos2D["h2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.4f && genPart_eta < 2.6f) m_Histos2D["h2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.6f && genPart_eta < 2.8f) m_Histos2D["h2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.8f && genPart_eta < 3.0f) m_Histos2D["h2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_energy"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);

      if      (genPart_eta >= 0.0f && genPart_eta < 0.2f) m_Histos2D["h2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.2f && genPart_eta < 0.4f) m_Histos2D["h2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.4f && genPart_eta < 0.6f) m_Histos2D["h2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.6f && genPart_eta < 0.8f) m_Histos2D["h2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 0.8f && genPart_eta < 1.0f) m_Histos2D["h2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.0f && genPart_eta < 1.2f) m_Histos2D["h2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.2f && genPart_eta < 1.4f) m_Histos2D["h2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.4f && genPart_eta < 1.6f) m_Histos2D["h2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.6f && genPart_eta < 1.8f) m_Histos2D["h2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 1.8f && genPart_eta < 2.0f) m_Histos2D["h2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.0f && genPart_eta < 2.2f) m_Histos2D["h2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.2f && genPart_eta < 2.4f) m_Histos2D["h2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.4f && genPart_eta < 2.6f) m_Histos2D["h2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.6f && genPart_eta < 2.8f) m_Histos2D["h2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
      else if (genPart_eta >= 2.8f && genPart_eta < 3.0f) m_Histos2D["h2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_energy_v2"]->Fill(ieta,map_hbheToweInfoHBHE[ieta].energy);
    // }
  }

  for (std::map<unsigned int, SimHitInfo>::const_iterator it = map_simHitInfoHBHE.cbegin(); it != map_simHitInfoHBHE.cend(); ++it) {
    DetId newid = HcalHitRelabeller::relabel(it->first, hcaloConst);
    HcalDetId theHcalDetId(newid);

    // double eta = -9.;
    // double phi = -9.;
    // if (newid.subdetId() == HcalSubdetector::HcalBarrel){
    //   eta = cellGeometryHB->getPosition(theHcalDetId).eta();
    //   phi = cellGeometryHB->getPosition(theHcalDetId).phi();
    // }else if(newid.subdetId() == HcalSubdetector::HcalEndcap){
    //   eta = cellGeometryHE->getPosition(theHcalDetId).eta();
    //   phi = cellGeometryHE->getPosition(theHcalDetId).phi();
    // }
    const HcalSimParameters& pars = dynamic_cast<const HcalSimParameters&>(theParameterMap->simParameters(theHcalDetId));
    float samplingFactor = pars.samplingFactor(theHcalDetId);

    float energyThisCell = (it->second).energy * samplingFactor;
    float nCaloHits = (it->second).nCaloHits;
    int ieta = theHcalDetId.ieta();
    int depth = theHcalDetId.depth();

    if      (genPart_eta >= 0.0f && genPart_eta < 0.2f) m_Profiles2D["hp2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 0.2f && genPart_eta < 0.4f) m_Profiles2D["hp2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 0.4f && genPart_eta < 0.6f) m_Profiles2D["hp2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 0.6f && genPart_eta < 0.8f) m_Profiles2D["hp2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 0.8f && genPart_eta < 1.0f) m_Profiles2D["hp2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 1.0f && genPart_eta < 1.2f) m_Profiles2D["hp2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 1.2f && genPart_eta < 1.4f) m_Profiles2D["hp2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 1.4f && genPart_eta < 1.6f) m_Profiles2D["hp2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 1.6f && genPart_eta < 1.8f) m_Profiles2D["hp2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 1.8f && genPart_eta < 2.0f) m_Profiles2D["hp2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 2.0f && genPart_eta < 2.2f) m_Profiles2D["hp2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 2.2f && genPart_eta < 2.4f) m_Profiles2D["hp2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 2.4f && genPart_eta < 2.6f) m_Profiles2D["hp2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 2.6f && genPart_eta < 2.8f) m_Profiles2D["hp2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);
    else if (genPart_eta >= 2.8f && genPart_eta < 3.0f) m_Profiles2D["hp2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_energy"]->Fill(ieta,depth,energyThisCell);

    if(energyThisCell > 0.1){
      if      (genPart_eta >= 0.0f && genPart_eta < 0.2f) m_Histos2D["h2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 0.2f && genPart_eta < 0.4f) m_Histos2D["h2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 0.4f && genPart_eta < 0.6f) m_Histos2D["h2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 0.6f && genPart_eta < 0.8f) m_Histos2D["h2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 0.8f && genPart_eta < 1.0f) m_Histos2D["h2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 1.0f && genPart_eta < 1.2f) m_Histos2D["h2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 1.2f && genPart_eta < 1.4f) m_Histos2D["h2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 1.4f && genPart_eta < 1.6f) m_Histos2D["h2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 1.6f && genPart_eta < 1.8f) m_Histos2D["h2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 1.8f && genPart_eta < 2.0f) m_Histos2D["h2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 2.0f && genPart_eta < 2.2f) m_Histos2D["h2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 2.2f && genPart_eta < 2.4f) m_Histos2D["h2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 2.4f && genPart_eta < 2.6f) m_Histos2D["h2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 2.6f && genPart_eta < 2.8f) m_Histos2D["h2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
      else if (genPart_eta >= 2.8f && genPart_eta < 3.0f) m_Histos2D["h2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_depth"]->Fill(ieta,depth);
    }

    if      (depth==1) m_Profiles2D["hp2_depth1_genPartEta_vs_hcalsimhit_ieta_vs_energy"]->Fill(genPart_eta,ieta,energyThisCell);
    else if (depth==2) m_Profiles2D["hp2_depth2_genPartEta_vs_hcalsimhit_ieta_vs_energy"]->Fill(genPart_eta,ieta,energyThisCell);
    else if (depth==3) m_Profiles2D["hp2_depth3_genPartEta_vs_hcalsimhit_ieta_vs_energy"]->Fill(genPart_eta,ieta,energyThisCell);
    else if (depth==4) m_Profiles2D["hp2_depth4_genPartEta_vs_hcalsimhit_ieta_vs_energy"]->Fill(genPart_eta,ieta,energyThisCell);
    else if (depth==5) m_Profiles2D["hp2_depth5_genPartEta_vs_hcalsimhit_ieta_vs_energy"]->Fill(genPart_eta,ieta,energyThisCell);
    else if (depth==6) m_Profiles2D["hp2_depth6_genPartEta_vs_hcalsimhit_ieta_vs_energy"]->Fill(genPart_eta,ieta,energyThisCell);
    else if (depth==7) m_Profiles2D["hp2_depth7_genPartEta_vs_hcalsimhit_ieta_vs_energy"]->Fill(genPart_eta,ieta,energyThisCell);


    if      (genPart_eta >= 0.0f && genPart_eta < 0.2f) m_Profiles2D["hp2_genPartEta0p0To0p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 0.2f && genPart_eta < 0.4f) m_Profiles2D["hp2_genPartEta0p2To0p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 0.4f && genPart_eta < 0.6f) m_Profiles2D["hp2_genPartEta0p4To0p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 0.6f && genPart_eta < 0.8f) m_Profiles2D["hp2_genPartEta0p6To0p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 0.8f && genPart_eta < 1.0f) m_Profiles2D["hp2_genPartEta0p8To1p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 1.0f && genPart_eta < 1.2f) m_Profiles2D["hp2_genPartEta1p0To1p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 1.2f && genPart_eta < 1.4f) m_Profiles2D["hp2_genPartEta1p2To1p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 1.4f && genPart_eta < 1.6f) m_Profiles2D["hp2_genPartEta1p4To1p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 1.6f && genPart_eta < 1.8f) m_Profiles2D["hp2_genPartEta1p6To1p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 1.8f && genPart_eta < 2.0f) m_Profiles2D["hp2_genPartEta1p8To2p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 2.0f && genPart_eta < 2.2f) m_Profiles2D["hp2_genPartEta2p0To2p2_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 2.2f && genPart_eta < 2.4f) m_Profiles2D["hp2_genPartEta2p2To2p4_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 2.4f && genPart_eta < 2.6f) m_Profiles2D["hp2_genPartEta2p4To2p6_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 2.6f && genPart_eta < 2.8f) m_Profiles2D["hp2_genPartEta2p6To2p8_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);
    else if (genPart_eta >= 2.8f && genPart_eta < 3.0f) m_Profiles2D["hp2_genPartEta2p8To3p0_hcalsimhit_ieta_vs_hcalsimhit_depth_vs_ncalohit"]->Fill(ieta,depth,nCaloHits);

    if      (depth==1) m_Profiles2D["hp2_depth1_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit"]->Fill(genPart_eta,ieta,nCaloHits);
    else if (depth==2) m_Profiles2D["hp2_depth2_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit"]->Fill(genPart_eta,ieta,nCaloHits);
    else if (depth==3) m_Profiles2D["hp2_depth3_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit"]->Fill(genPart_eta,ieta,nCaloHits);
    else if (depth==4) m_Profiles2D["hp2_depth4_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit"]->Fill(genPart_eta,ieta,nCaloHits);
    else if (depth==5) m_Profiles2D["hp2_depth5_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit"]->Fill(genPart_eta,ieta,nCaloHits);
    else if (depth==6) m_Profiles2D["hp2_depth6_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit"]->Fill(genPart_eta,ieta,nCaloHits);
    else if (depth==7) m_Profiles2D["hp2_depth7_genPartEta_vs_hcalsimhit_ieta_vs_ncalohit"]->Fill(genPart_eta,ieta,nCaloHits);
  }

  //==========================================
  //
  // PCaloHit: EB
  //
  //==========================================
  /*
  std::map<unsigned int, SimHitInfo> map_simHitInfoEB;

  edm::Handle<edm::PCaloHitContainer> g4SimHitsPCaloHitsEBHandle;
  iEvent.getByToken(g4SimHitsPCaloHitsEBToken_, g4SimHitsPCaloHitsEBHandle);
  auto g4SimHitsPCaloHitsEB = g4SimHitsPCaloHitsEBHandle.product();
  size_t n_g4SimHitsPCaloHitsEB = g4SimHitsPCaloHitsEB->size();

  nG4SimHitPCaloHitEB=0;
  for (size_t idx = 0; idx < n_g4SimHitsPCaloHitsEB; idx++){
    const PCaloHit& caloHit = g4SimHitsPCaloHitsEB->at(idx);
    DetId detIdRaw(caloHit.id());
    EBDetId theECalEBDetId(detIdRaw);
    auto cellGeometry = geometry->getSubdetectorGeometry(theECalEBDetId)->getGeometry(theECalEBDetId);
    double eta = cellGeometry->getPosition().eta();
    double phi = cellGeometry->getPosition().phi();

    auto it = map_simHitInfoEB.find(caloHit.id());
    // DetId not in map yet, so lets set SimHitInfo for this DetId.
    if (it == map_simHitInfoEB.end()) {
      map_simHitInfoEB[caloHit.id()] = SimHitInfo();
      map_simHitInfoEB[caloHit.id()].energy    = 0.f;
      map_simHitInfoEB[caloHit.id()].energyEM  = 0.f;
      map_simHitInfoEB[caloHit.id()].energyHAD = 0.f;
      map_simHitInfoEB[caloHit.id()].nCaloHits = 0;
    }

    map_simHitInfoEB[caloHit.id()].energy    += caloHit.energy();
    map_simHitInfoEB[caloHit.id()].energyEM  += caloHit.energyEM();
    map_simHitInfoEB[caloHit.id()].energyHAD += caloHit.energyHad();
    map_simHitInfoEB[caloHit.id()].nCaloHits += 1;
  }
  */

  //==========================================
  //
  // PCaloHit: EE
  //
  //==========================================
  /*
  std::map<unsigned int, SimHitInfo> map_simHitInfoEE;

  edm::Handle<edm::PCaloHitContainer> g4SimHitsPCaloHitsEEHandle;
  iEvent.getByToken(g4SimHitsPCaloHitsEEToken_, g4SimHitsPCaloHitsEEHandle);
  auto g4SimHitsPCaloHitsEE = g4SimHitsPCaloHitsEEHandle.product();
  size_t n_g4SimHitsPCaloHitsEE = g4SimHitsPCaloHitsEE->size();

  for (size_t idx = 0; idx < n_g4SimHitsPCaloHitsEE; idx++){
    const PCaloHit& caloHit = g4SimHitsPCaloHitsEE->at(idx);
    DetId detIdRaw(caloHit.id());
    EEDetId theECalEEDetId(detIdRaw);
    auto cellGeometry = geometry->getSubdetectorGeometry(theECalEEDetId)->getGeometry(theECalEEDetId);
    double eta = cellGeometry->getPosition().eta();
    double phi = cellGeometry->getPosition().phi();

    auto it = map_simHitInfoEE.find(caloHit.id());
    // DetId not in map yet, so lets set SimHitInfo for this DetId.
    if (it == map_simHitInfoEE.end()) {
      map_simHitInfoEE[caloHit.id()] = SimHitInfo();
      map_simHitInfoEE[caloHit.id()].energy    = 0.f;
      map_simHitInfoEE[caloHit.id()].energyEM  = 0.f;
      map_simHitInfoEE[caloHit.id()].energyHAD = 0.f;
      map_simHitInfoEE[caloHit.id()].nCaloHits = 0;
    }
    map_simHitInfoEE[caloHit.id()].energy    += caloHit.energy();
    map_simHitInfoEE[caloHit.id()].energyEM  += caloHit.energyEM();
    map_simHitInfoEE[caloHit.id()].energyHAD += caloHit.energyHad();
    map_simHitInfoEE[caloHit.id()].nCaloHits += 1;
  }
  */
}
template <typename T, typename... Args>
T* ParticleFlowAnalysisHistogrammer::book(const Args &...args) const {
  T *t = fs_->make<T>(args...);
  return t;
}

void ParticleFlowAnalysisHistogrammer::bookTH1F(TString hName, TString hTitle, int nbinsx, float xmin, float xmax){
  m_Histos1D[hName] = book<TH1F>(hName, hTitle, nbinsx, xmin, xmax);
}
void ParticleFlowAnalysisHistogrammer::bookTH2F(TString hName, TString hTitle, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax){
  m_Histos2D[hName] = book<TH2F>(hName, hTitle, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
}
void ParticleFlowAnalysisHistogrammer::bookTH3F(TString hName, TString hTitle, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, int nbinsz, float zmin, float zmax){
  m_Histos3D[hName] = book<TH3F>(hName, hTitle, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax);
}
void ParticleFlowAnalysisHistogrammer::bookProfile(TString pName, TString pTitle, int nbinsx, float xmin, float xmax){
  m_Profiles[pName] = book<TProfile>(pName, pTitle, nbinsx, xmin, xmax);
}
void ParticleFlowAnalysisHistogrammer::bookProfile2D(TString pName, TString pTitle, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax){
  m_Profiles2D[pName] = book<TProfile2D>(pName, pTitle, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
}

// const SimTrack* ParticleFlowAnalysisHistogrammer::getSimTrack(unsigned fSimTrackId) {
//   for (unsigned i = 0; i < getSimTrackCollection()->size(); ++i) {
//     if ((*getSimTrackCollection())[i].trackId() == fSimTrackId)
//       return &(*getSimTrackCollection())[i];
//   }
//   return nullptr;
// }

// // https://cmssdt.cern.ch/lxr/source/RecoJets/JetProducers/src/JetMatchingTools.cc#0230
// // https://cmssdt.cern.ch/lxr/source/RecoJets/JetProducers/interface/JetMatchingTools.h
// const SimTrack* JetMatchingTools::getTrack(unsigned fSimTrackId) {
//   for (unsigned i = 0; i < getSimTrackCollection()->size(); ++i) {
//     if ((*getSimTrackCollection())[i].trackId() == fSimTrackId)
//       return &(*getSimTrackCollection())[i];
//   }
//   return nullptr;
// }

DEFINE_FWK_MODULE(ParticleFlowAnalysisHistogrammer);
