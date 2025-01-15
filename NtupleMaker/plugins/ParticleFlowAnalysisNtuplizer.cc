#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "ParticleFlowAnalysis/NtupleMaker/plugins/ParticleFlowAnalysisNtuplizer.h"

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

struct SimHitInfo {
  float energy;
  float energyEM;
  float energyHAD;
  int nCaloHits;
  float energy_25ns;
  float energyEM_25ns;
  float energyHAD_25ns;
  int nCaloHits_25ns;

  SimHitInfo(){
    energy=0; energyHAD=0; energyEM=0; nCaloHits=0;
    energy_25ns=0; energyHAD_25ns=0; energyEM_25ns=0; nCaloHits_25ns=0;
  }
};

ParticleFlowAnalysisNtuplizer::ParticleFlowAnalysisNtuplizer(const edm::ParameterSet& iConfig) {

  tokengenParticles_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<InputTag>("genParticles"));
  tokengenVertexXYZ_ = consumes<math::XYZPointF>(iConfig.getParameter<InputTag>("genVertex"));

  tokenPFCandidates_ = consumes<reco::PFCandidateCollection>(iConfig.getParameter<InputTag>("PFCandidates"));;

  pfClustersECALToken_ = consumes<reco::PFClusterCollection>(iConfig.getParameter<InputTag>("PFClustersECAL"));
  pfClustersPSToken_   = consumes<reco::PFClusterCollection>(iConfig.getParameter<InputTag>("PFClustersPS"));
  pfClustersHCALToken_ = consumes<reco::PFClusterCollection>(iConfig.getParameter<InputTag>("PFClustersHCAL"));

  pfRecHitsECALToken_ = consumes<std::vector<reco::PFRecHit>>(iConfig.getParameter<edm::InputTag>("PFRecHitsECAL"));
  pfRecHitsPSToken_   = consumes<std::vector<reco::PFRecHit>>(iConfig.getParameter<edm::InputTag>("PFRecHitsPS"));
  pfRecHitsHBHEToken_ = consumes<std::vector<reco::PFRecHit>>(iConfig.getParameter<edm::InputTag>("PFRecHitsHBHE"));
  pfRecHitsHFToken_   = consumes<std::vector<reco::PFRecHit>>(iConfig.getParameter<edm::InputTag>("PFRecHitsHF"));

  chargedHadronIsolationToken_ = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("chargedHadronIsolation"));

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

  // Smallest number of pixel hits
  nPixMin_ = iConfig.getParameter<int>("nPixMin");

  // Smallest number of track hits in different eta ranges
  nHitMin_ = iConfig.getParameter< std::vector<int> > ("nHitMin");
  etaMinForHitMin_ = iConfig.getParameter< std::vector<double> > ("etaMinForHitMin");
  etaMaxForHitMin_ = iConfig.getParameter< std::vector<double> > ("etaMaxForHitMin");

  dRMaxGenPartToPFCandChgHad_ = iConfig.getUntrackedParameter<double>("dRMaxGenPartToPFCandChgHad",0.4);
  dRMaxGenPartToPFCandNeuHad_ = iConfig.getUntrackedParameter<double>("dRMaxGenPartToPFCandNeuHad",0.4);
  dRMaxGenPartToPFCandPhoton_ = iConfig.getUntrackedParameter<double>("dRMaxGenPartToPFCandPhoton",0.4);

  savePFClustersECAL_ = iConfig.getUntrackedParameter<bool>("savePFClustersECAL",false);
  savePFClustersPS_ = iConfig.getUntrackedParameter<bool>("savePFClustersPS",false);
  savePFClustersHCAL_ = iConfig.getUntrackedParameter<bool>("savePFClustersHCAL",false);

  saveSimCaloHitHBHE_ = iConfig.getUntrackedParameter<bool>("saveSimCaloHitHBHE",false);
  saveSimCaloHitEB_ = iConfig.getUntrackedParameter<bool>("saveSimCaloHitEB",false);
  saveSimCaloHitEE_ = iConfig.getUntrackedParameter<bool>("saveSimCaloHitEE",false);
  saveSimHitHBHE_ = iConfig.getUntrackedParameter<bool>("saveSimHitHBHE",false);
  saveSimHitEB_ = iConfig.getUntrackedParameter<bool>("saveSimHitEB",false);
  saveSimHitEE_ = iConfig.getUntrackedParameter<bool>("saveSimHitEE",false);

  //
  //
  //
  theParameterMap = new HcalSimParameterMap(iConfig);

  // The root tuple
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("PFAnaTree", "PFAnaTree");

  tree->Branch("run",&orun,"orun/l");
  tree->Branch("evt",&oevt,"orun/l");
  tree->Branch("lumiBlock",&olumiBlock,"orun/l");
  tree->Branch("time",&otime,"orun/l");

  tree->Branch("genVertex_x",&genVertex_x);
  tree->Branch("genVertex_y",&genVertex_y);
  tree->Branch("genVertex_z",&genVertex_z);

  tree->Branch("genPart_pt", &genPart_pt);
  tree->Branch("genPart_eta", &genPart_eta);
  tree->Branch("genPart_phi", &genPart_phi);
  tree->Branch("genPart_mass", &genPart_mass);
  tree->Branch("genPart_energy", &genPart_energy);
  tree->Branch("genPart_p", &genPart_p);
  tree->Branch("genPart_pdgId", &genPart_pdgId);
  tree->Branch("genPart_charge", &genPart_charge);

  tree->Branch("genPart_extrapolated_okECAL",  &genPart_extrapolated_okECAL);
  tree->Branch("genPart_extrapolated_okHCAL",  &genPart_extrapolated_okHCAL);
  tree->Branch("genPart_extrapolated_EB_ieta", &genPart_extrapolated_EB_ieta);
  tree->Branch("genPart_extrapolated_EB_iphi", &genPart_extrapolated_EB_iphi);
  tree->Branch("genPart_extrapolated_EB_eta",  &genPart_extrapolated_EB_eta);
  tree->Branch("genPart_extrapolated_EB_phi",  &genPart_extrapolated_EB_phi);
  tree->Branch("genPart_extrapolated_PFRecHitEB_Idx",    &genPart_extrapolated_PFRecHitEB_Idx);
  tree->Branch("genPart_extrapolated_PFRecHitEB_detId",  &genPart_extrapolated_PFRecHitEB_detId);

  tree->Branch("genPart_extrapolated_EE_ix",     &genPart_extrapolated_EE_ix);
  tree->Branch("genPart_extrapolated_EE_iy",     &genPart_extrapolated_EE_iy);
  tree->Branch("genPart_extrapolated_EE_eta",    &genPart_extrapolated_EE_eta);
  tree->Branch("genPart_extrapolated_EE_phi",    &genPart_extrapolated_EE_phi);
  tree->Branch("genPart_extrapolated_PFRecHitEE_Idx",    &genPart_extrapolated_PFRecHitEE_Idx);
  tree->Branch("genPart_extrapolated_PFRecHitEE_detId",  &genPart_extrapolated_PFRecHitEE_detId);

  tree->Branch("genPart_extrapolated_EHCAL_ieta",              &genPart_extrapolated_EHCAL_ieta);
  tree->Branch("genPart_extrapolated_EHCAL_iphi",              &genPart_extrapolated_EHCAL_iphi);
  tree->Branch("genPart_extrapolated_EHCAL_depth",             &genPart_extrapolated_EHCAL_depth);
  tree->Branch("genPart_extrapolated_EHCAL_eta",               &genPart_extrapolated_EHCAL_eta);
  tree->Branch("genPart_extrapolated_EHCAL_phi",               &genPart_extrapolated_EHCAL_phi);
  tree->Branch("genPart_extrapolated_EHCAL_subdetId",          &genPart_extrapolated_EHCAL_subdetId);
  tree->Branch("genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx",  &genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx);
  tree->Branch("genPart_extrapolated_EHCAL_PFRecHitHBHE_detId",&genPart_extrapolated_EHCAL_PFRecHitHBHE_detId);
  tree->Branch("genPart_extrapolated_EHCAL_PFRecHitHF_Idx",    &genPart_extrapolated_EHCAL_PFRecHitHF_Idx);
  tree->Branch("genPart_extrapolated_EHCAL_PFRecHitHF_detId",  &genPart_extrapolated_EHCAL_PFRecHitHF_detId);

  tree->Branch("genPart_extrapolated_HCAL_ieta",              &genPart_extrapolated_HCAL_ieta);
  tree->Branch("genPart_extrapolated_HCAL_iphi",              &genPart_extrapolated_HCAL_iphi);
  tree->Branch("genPart_extrapolated_HCAL_depth",             &genPart_extrapolated_HCAL_depth);
  tree->Branch("genPart_extrapolated_HCAL_eta",               &genPart_extrapolated_HCAL_eta);
  tree->Branch("genPart_extrapolated_HCAL_phi",               &genPart_extrapolated_HCAL_phi);
  tree->Branch("genPart_extrapolated_HCAL_subdetId",          &genPart_extrapolated_HCAL_subdetId);
  tree->Branch("genPart_extrapolated_HCAL_PFRecHitHBHE_Idx",  &genPart_extrapolated_HCAL_PFRecHitHBHE_Idx);
  tree->Branch("genPart_extrapolated_HCAL_PFRecHitHBHE_detId",&genPart_extrapolated_HCAL_PFRecHitHBHE_detId);
  tree->Branch("genPart_extrapolated_HCAL_PFRecHitHF_Idx",    &genPart_extrapolated_HCAL_PFRecHitHF_Idx);
  tree->Branch("genPart_extrapolated_HCAL_PFRecHitHF_detId",  &genPart_extrapolated_HCAL_PFRecHitHF_detId);

  tree->Branch("genPart_extrapolated_pointECAL_x",&genPart_extrapolated_pointECAL_x);
  tree->Branch("genPart_extrapolated_pointECAL_y",&genPart_extrapolated_pointECAL_y);
  tree->Branch("genPart_extrapolated_pointECAL_z",&genPart_extrapolated_pointECAL_z);
  tree->Branch("genPart_extrapolated_pointECAL_eta",&genPart_extrapolated_pointECAL_eta);
  tree->Branch("genPart_extrapolated_pointECAL_phi",&genPart_extrapolated_pointECAL_phi);

  tree->Branch("genPart_extrapolated_pointHCAL_x",&genPart_extrapolated_pointHCAL_x);
  tree->Branch("genPart_extrapolated_pointHCAL_y",&genPart_extrapolated_pointHCAL_y);
  tree->Branch("genPart_extrapolated_pointHCAL_z",&genPart_extrapolated_pointHCAL_z);
  tree->Branch("genPart_extrapolated_pointHCAL_eta",&genPart_extrapolated_pointHCAL_eta);
  tree->Branch("genPart_extrapolated_pointHCAL_phi",&genPart_extrapolated_pointHCAL_phi);

  tree->Branch("nPFCand",&nPFCand);
  tree->Branch("PFCand_pt",&PFCand_pt);
  tree->Branch("PFCand_eta",&PFCand_eta);
  tree->Branch("PFCand_phi",&PFCand_phi);
  tree->Branch("PFCand_mass",&PFCand_mass);
  tree->Branch("PFCand_energy",&PFCand_energy);
  tree->Branch("PFCand_pdgId",&PFCand_pdgId);
  tree->Branch("PFCand_isChgHadIso",&PFCand_isChgHadIso);
  tree->Branch("PFCand_birthId",&PFCand_birthId);
  tree->Branch("PFCand_dRToGenPart",&PFCand_dRToGenPart);
  tree->Branch("PFCand_has_trk",&PFCand_has_trk);

  tree->Branch("PFCand_trk_pt",&PFCand_trk_pt);
  tree->Branch("PFCand_trk_p",&PFCand_trk_p);
  tree->Branch("PFCand_trk_eta",&PFCand_trk_eta);
  tree->Branch("PFCand_trk_phi",&PFCand_trk_phi);
  tree->Branch("PFCand_trk_algo",&PFCand_trk_algo);

  tree->Branch("PFCand_trk_charge",&PFCand_trk_charge);
  tree->Branch("PFCand_trk_ptError",&PFCand_trk_ptError);
  tree->Branch("PFCand_trk_dxy",&PFCand_trk_dxy);
  tree->Branch("PFCand_trk_dz",&PFCand_trk_dz);
  tree->Branch("PFCand_trk_dxyError",&PFCand_trk_dxyError);
  tree->Branch("PFCand_trk_dzError",&PFCand_trk_dzError);
  tree->Branch("PFCand_trk_chi2",&PFCand_trk_chi2);
  tree->Branch("PFCand_trk_ndof",&PFCand_trk_ndof);
  tree->Branch("PFCand_trk_vx",&PFCand_trk_vx);
  tree->Branch("PFCand_trk_vy",&PFCand_trk_vy);
  tree->Branch("PFCand_trk_vz",&PFCand_trk_vz);

  tree->Branch("PFCand_trk_numberOfValidPixelHits",&PFCand_trk_numberOfValidPixelHits);
  tree->Branch("PFCand_trk_numberOfValidTrackerHits",&PFCand_trk_numberOfValidTrackerHits);
  tree->Branch("PFCand_trk_pixelHitOK",&PFCand_trk_pixelHitOK);
  tree->Branch("PFCand_trk_trackerHitOK",&PFCand_trk_trackerHitOK);

  tree->Branch("PFCand_trk_positionAtECALEntrance_x",&PFCand_trk_positionAtECALEntrance_x);
  tree->Branch("PFCand_trk_positionAtECALEntrance_y",&PFCand_trk_positionAtECALEntrance_y);
  tree->Branch("PFCand_trk_positionAtECALEntrance_z",&PFCand_trk_positionAtECALEntrance_z);
  tree->Branch("PFCand_trk_positionAtECALEntrance_eta",&PFCand_trk_positionAtECALEntrance_eta);
  tree->Branch("PFCand_trk_positionAtECALEntrance_phi",&PFCand_trk_positionAtECALEntrance_phi);

  tree->Branch("PFCand_trk_closestECAL_eta",&PFCand_trk_closestECAL_eta);
  tree->Branch("PFCand_trk_closestECAL_phi",&PFCand_trk_closestECAL_phi);
  tree->Branch("PFCand_trk_closestECAL_detId",&PFCand_trk_closestECAL_detId);
  tree->Branch("PFCand_trk_closestECAL_EBieta",&PFCand_trk_closestECAL_EBieta);
  tree->Branch("PFCand_trk_closestECAL_EBiphi",&PFCand_trk_closestECAL_EBiphi);
  tree->Branch("PFCand_trk_closestECAL_EEix",&PFCand_trk_closestECAL_EEix);
  tree->Branch("PFCand_trk_closestECAL_EEiy",&PFCand_trk_closestECAL_EEiy);
  tree->Branch("PFCand_trk_closestHCAL_eta",&PFCand_trk_closestHCAL_eta);
  tree->Branch("PFCand_trk_closestHCAL_phi",&PFCand_trk_closestHCAL_phi);
  tree->Branch("PFCand_trk_closestHCAL_detId",&PFCand_trk_closestHCAL_detId);
  tree->Branch("PFCand_trk_closestHCAL_ieta",&PFCand_trk_closestHCAL_ieta);
  tree->Branch("PFCand_trk_closestHCAL_iphi",&PFCand_trk_closestHCAL_iphi);

  tree->Branch("PFCand_ecalEnergy",&PFCand_ecalEnergy);
  tree->Branch("PFCand_hcalEnergy",&PFCand_hcalEnergy);
  tree->Branch("PFCand_caloEnergy",&PFCand_caloEnergy);
  tree->Branch("PFCand_rawEcalEnergy",&PFCand_rawEcalEnergy);
  tree->Branch("PFCand_rawHcalEnergy",&PFCand_rawHcalEnergy);
  tree->Branch("PFCand_rawCaloEnergy",&PFCand_rawCaloEnergy);
  tree->Branch("PFCand_hoEnergy",&PFCand_hoEnergy);
  tree->Branch("PFCand_rawHoEnergy",&PFCand_rawHoEnergy);
  tree->Branch("PFCand_pS1Energy",&PFCand_pS1Energy);
  tree->Branch("PFCand_p21Energy",&PFCand_pS2Energy);
  tree->Branch("PFCand_hcalDepth1EFrac",&PFCand_hcalDepth1EFrac);
  tree->Branch("PFCand_hcalDepth2EFrac",&PFCand_hcalDepth2EFrac);
  tree->Branch("PFCand_hcalDepth3EFrac",&PFCand_hcalDepth3EFrac);
  tree->Branch("PFCand_hcalDepth4EFrac",&PFCand_hcalDepth4EFrac);
  tree->Branch("PFCand_hcalDepth5EFrac",&PFCand_hcalDepth5EFrac);
  tree->Branch("PFCand_hcalDepth6EFrac",&PFCand_hcalDepth6EFrac);
  tree->Branch("PFCand_hcalDepth7EFrac",&PFCand_hcalDepth7EFrac);
  tree->Branch("PFCand_PFBlock_Idx",&PFCand_PFBlock_Idx);

  tree->Branch("PFCand_nPFTrack",&PFCand_nPFTrack);
  tree->Branch("PFCand_nPFClusterHCAL",&PFCand_nPFClusterHCAL);
  tree->Branch("PFCand_nPFClusterECAL",&PFCand_nPFClusterECAL);
  tree->Branch("PFCand_nPFClusterPS",&PFCand_nPFClusterPS);
  tree->Branch("PFCand_nPFClusterHF",&PFCand_nPFClusterHF);
  tree->Branch("PFCand_nPFClusterSC",&PFCand_nPFClusterSC);
  tree->Branch("PFCand_nPFClusterGSF",&PFCand_nPFClusterGSF);
  tree->Branch("PFCand_nPFClusterBREM",&PFCand_nPFClusterBREM);
  tree->Branch("PFCand_PFClusterHCAL_Idx",&PFCand_PFClusterHCAL_Idx);
  tree->Branch("PFCand_PFClusterECAL_Idx",&PFCand_PFClusterECAL_Idx);
  tree->Branch("PFCand_PFClusterPS_Idx",&PFCand_PFClusterPS_Idx);
  tree->Branch("PFCand_PFClusterHF_Idx",&PFCand_PFClusterHF_Idx);
  tree->Branch("PFCand_nPFTrackInBlock",&PFCand_nPFTrackInBlock);
  tree->Branch("PFCand_nPFClusterHCALInBlock",&PFCand_nPFClusterHCALInBlock);
  tree->Branch("PFCand_nPFClusterECALInBlock",&PFCand_nPFClusterECALInBlock);
  tree->Branch("PFCand_nPFClusterPSInBlock",&PFCand_nPFClusterPSInBlock);
  tree->Branch("PFCand_nPFClusterHFInBlock",&PFCand_nPFClusterHFInBlock);
  tree->Branch("PFCand_nPFClusterSCInBlock",&PFCand_nPFClusterSCInBlock);
  tree->Branch("PFCand_nPFClusterGSFInBlock",&PFCand_nPFClusterGSFInBlock);
  tree->Branch("PFCand_nPFClusterBREMInBlock",&PFCand_nPFClusterBREMInBlock);
  tree->Branch("PFCand_PFClusterHCALInBlock_Idx",&PFCand_PFClusterHCALInBlock_Idx);
  tree->Branch("PFCand_PFClusterECALInBlock_Idx",&PFCand_PFClusterECALInBlock_Idx);
  tree->Branch("PFCand_PFClusterPSInBlock_Idx",&PFCand_PFClusterPSInBlock_Idx);
  tree->Branch("PFCand_PFClusterHFInBlock_Idx",&PFCand_PFClusterHFInBlock_Idx);

  tree->Branch("Idx_ClosestPFCandChgHad",&Idx_ClosestPFCandChgHad);
  tree->Branch("Idx_ClosestPFCandChgHadIso",&Idx_ClosestPFCandChgHadIso);
  tree->Branch("Idx_ClosestPFCandNeuHad",&Idx_ClosestPFCandNeuHad);
  tree->Branch("Idx_ClosestPFCandPhoton",&Idx_ClosestPFCandPhoton);

  tree->Branch("nPFRecHitHBHE",&nPFRecHitHBHE,"nPFRecHitHBHE/I");
  tree->Branch("PFRecHitHBHE_energy",&PFRecHitHBHE_energy);
  tree->Branch("PFRecHitHBHE_ieta",  &PFRecHitHBHE_ieta);
  tree->Branch("PFRecHitHBHE_iphi",  &PFRecHitHBHE_iphi);
  tree->Branch("PFRecHitHBHE_eta",   &PFRecHitHBHE_eta);
  tree->Branch("PFRecHitHBHE_phi",   &PFRecHitHBHE_phi);
  tree->Branch("PFRecHitHBHE_depth", &PFRecHitHBHE_depth);
  tree->Branch("PFRecHitHBHE_cutThreshold", &PFRecHitHBHE_cutThreshold);
  tree->Branch("PFRecHitHBHE_hcalRespCorr", &PFRecHitHBHE_hcalRespCorr);
  tree->Branch("PFRecHitHBHE_detId", &PFRecHitHBHE_detId);

  // tree->Branch("nPFRecHitHF",&nPFRecHitHF,"nPFRecHitHF/I");
  // tree->Branch("PFRecHitHF_energy",&PFRecHitHF_energy);
  // tree->Branch("PFRecHitHF_ieta",  &PFRecHitHF_ieta);
  // tree->Branch("PFRecHitHF_iphi",  &PFRecHitHF_iphi);
  // tree->Branch("PFRecHitHF_eta",   &PFRecHitHF_eta);
  // tree->Branch("PFRecHitHF_phi",   &PFRecHitHF_phi);
  // tree->Branch("PFRecHitHF_depth", &PFRecHitHF_depth);

  tree->Branch("nPFRecHitEB",&nPFRecHitEB,"nPFRecHitEB/I");
  tree->Branch("PFRecHitEB_energy",       &PFRecHitEB_energy);
  tree->Branch("PFRecHitEB_ieta",         &PFRecHitEB_ieta);
  tree->Branch("PFRecHitEB_iphi",         &PFRecHitEB_iphi);
  tree->Branch("PFRecHitEB_tower_ieta",   &PFRecHitEB_tower_ieta);
  tree->Branch("PFRecHitEB_tower_iphi",   &PFRecHitEB_tower_iphi);
  tree->Branch("PFRecHitEB_approxEta",    &PFRecHitEB_approxEta);
  tree->Branch("PFRecHitEB_eta",          &PFRecHitEB_eta);
  tree->Branch("PFRecHitEB_phi",          &PFRecHitEB_phi);
  tree->Branch("PFRecHitEB_cutThreshold", &PFRecHitEB_cutThreshold);
  tree->Branch("PFRecHitEB_detId", &PFRecHitEB_detId);

  tree->Branch("nPFRecHitEE",&nPFRecHitEE,"nPFRecHitEE/I");
  tree->Branch("PFRecHitEE_energy", &PFRecHitEE_energy);
  tree->Branch("PFRecHitEE_ix",     &PFRecHitEE_ix);
  tree->Branch("PFRecHitEE_iy",     &PFRecHitEE_iy);
  tree->Branch("PFRecHitEE_eta",    &PFRecHitEE_eta);
  tree->Branch("PFRecHitEE_phi",    &PFRecHitEE_phi);
  tree->Branch("PFRecHitEE_cutThreshold", &PFRecHitEE_cutThreshold);
  tree->Branch("PFRecHitEE_detId", &PFRecHitEE_detId);

  if (saveSimCaloHitHBHE_){
    tree->Branch("nG4SimHitPCaloHitHBHE",&nG4SimHitPCaloHitHBHE,"nG4SimHitPCaloHitHBHE/I");
    tree->Branch("G4SimHitPCaloHitHBHE_energy",&G4SimHitPCaloHitHBHE_energy);
    tree->Branch("G4SimHitPCaloHitHBHE_energyEM",&G4SimHitPCaloHitHBHE_energyEM);
    tree->Branch("G4SimHitPCaloHitHBHE_energyHAD",&G4SimHitPCaloHitHBHE_energyHAD);
    tree->Branch("G4SimHitPCaloHitHBHE_samplingFactor",&G4SimHitPCaloHitHBHE_samplingFactor);
    tree->Branch("G4SimHitPCaloHitHBHE_ieta",&G4SimHitPCaloHitHBHE_ieta);
    tree->Branch("G4SimHitPCaloHitHBHE_iphi",&G4SimHitPCaloHitHBHE_iphi);
    tree->Branch("G4SimHitPCaloHitHBHE_depth",&G4SimHitPCaloHitHBHE_depth);
    tree->Branch("G4SimHitPCaloHitHBHE_eta",&G4SimHitPCaloHitHBHE_eta);
    tree->Branch("G4SimHitPCaloHitHBHE_phi",&G4SimHitPCaloHitHBHE_phi);
    tree->Branch("G4SimHitPCaloHitHBHE_time",&G4SimHitPCaloHitHBHE_time);
    tree->Branch("G4SimHitPCaloHitHBHE_detId",&G4SimHitPCaloHitHBHE_detId);
    tree->Branch("G4SimHitPCaloHitHBHE_subdetId",&G4SimHitPCaloHitHBHE_subdetId);
  }

  if (saveSimCaloHitEB_){
    tree->Branch("nG4SimHitPCaloHitEB",&nG4SimHitPCaloHitEB,"nG4SimHitPCaloHitEB/I");
    tree->Branch("G4SimHitPCaloHitEB_energy",&G4SimHitPCaloHitEB_energy);
    tree->Branch("G4SimHitPCaloHitEB_energyEM",&G4SimHitPCaloHitEB_energyEM);
    tree->Branch("G4SimHitPCaloHitEB_energyHAD",&G4SimHitPCaloHitEB_energyHAD);
    tree->Branch("G4SimHitPCaloHitEB_ieta",&G4SimHitPCaloHitEB_ieta);
    tree->Branch("G4SimHitPCaloHitEB_iphi",&G4SimHitPCaloHitEB_iphi);
    tree->Branch("G4SimHitPCaloHitEB_depth",&G4SimHitPCaloHitEB_depth);
    tree->Branch("G4SimHitPCaloHitEB_eta",&G4SimHitPCaloHitEB_eta);
    tree->Branch("G4SimHitPCaloHitEB_phi",&G4SimHitPCaloHitEB_phi);
    tree->Branch("G4SimHitPCaloHitEB_time",&G4SimHitPCaloHitEB_time);
    tree->Branch("G4SimHitPCaloHitEB_detId",&G4SimHitPCaloHitEB_detId);
    tree->Branch("G4SimHitPCaloHitEB_subdetId",&G4SimHitPCaloHitEB_subdetId);
  }

  if (saveSimCaloHitEE_){
    tree->Branch("nG4SimHitPCaloHitEE",&nG4SimHitPCaloHitEE,"nG4SimHitPCaloHitEE/I");
    tree->Branch("G4SimHitPCaloHitEE_energy",&G4SimHitPCaloHitEE_energy);
    tree->Branch("G4SimHitPCaloHitEE_energyEM",&G4SimHitPCaloHitEE_energyEM);
    tree->Branch("G4SimHitPCaloHitEE_energyHAD",&G4SimHitPCaloHitEE_energyHAD);
    tree->Branch("G4SimHitPCaloHitEE_ix",&G4SimHitPCaloHitEE_ix);
    tree->Branch("G4SimHitPCaloHitEE_iy",&G4SimHitPCaloHitEE_iy);
    tree->Branch("G4SimHitPCaloHitEE_depth",&G4SimHitPCaloHitEE_depth);
    tree->Branch("G4SimHitPCaloHitEE_eta",&G4SimHitPCaloHitEE_eta);
    tree->Branch("G4SimHitPCaloHitEE_phi",&G4SimHitPCaloHitEE_phi);
    tree->Branch("G4SimHitPCaloHitEE_time",&G4SimHitPCaloHitEE_time);
    tree->Branch("G4SimHitPCaloHitEE_detId",&G4SimHitPCaloHitEE_detId);
    tree->Branch("G4SimHitPCaloHitEE_subdetId",&G4SimHitPCaloHitEE_subdetId);
  }

  if (saveSimHitHBHE_){
    tree->Branch("nSimHitHBHE",&nSimHitHBHE,"nSimHitHBHE/I");
    tree->Branch("SimHitHBHE_energy",&SimHitHBHE_energy);
    tree->Branch("SimHitHBHE_energyEM",&SimHitHBHE_energyEM);
    tree->Branch("SimHitHBHE_energyHAD",&SimHitHBHE_energyHAD);
    tree->Branch("SimHitHBHE_samplingFactor",&SimHitHBHE_samplingFactor);
    tree->Branch("SimHitHBHE_nCaloHits",&SimHitHBHE_nCaloHits);
    tree->Branch("SimHitHBHE_ieta",&SimHitHBHE_ieta);
    tree->Branch("SimHitHBHE_iphi",&SimHitHBHE_iphi);
    tree->Branch("SimHitHBHE_depth",&SimHitHBHE_depth);
    tree->Branch("SimHitHBHE_eta",&SimHitHBHE_eta);
    tree->Branch("SimHitHBHE_phi",&SimHitHBHE_phi);
    tree->Branch("SimHitHBHE_detId",&SimHitHBHE_detId);
    tree->Branch("SimHitHBHE_subdetId",&SimHitHBHE_subdetId);
    tree->Branch("SimHitHBHE_energy_25ns",&SimHitHBHE_energy_25ns);
    tree->Branch("SimHitHBHE_energyEM_25ns",&SimHitHBHE_energyEM_25ns);
    tree->Branch("SimHitHBHE_energyHAD_25ns",&SimHitHBHE_energyHAD_25ns);
  }
  if (saveSimHitEB_){
    tree->Branch("nSimHitEB",&nSimHitEB,"nSimHitEB/I");
    tree->Branch("SimHitEB_energy",&SimHitEB_energy);
    tree->Branch("SimHitEB_energyEM",&SimHitEB_energyEM);
    tree->Branch("SimHitEB_energyHAD",&SimHitEB_energyHAD);
    tree->Branch("SimHitEB_nCaloHits",&SimHitEB_nCaloHits);
    tree->Branch("SimHitEB_ieta",&SimHitEB_ieta);
    tree->Branch("SimHitEB_iphi",&SimHitEB_iphi);
    tree->Branch("SimHitEB_eta",&SimHitEB_eta);
    tree->Branch("SimHitEB_phi",&SimHitEB_phi);
    tree->Branch("SimHitEB_detId",&SimHitEB_detId);
    tree->Branch("SimHitEB_subdetId",&SimHitEB_subdetId);
  }

  if (saveSimHitEE_){
    tree->Branch("nSimHitEE",&nSimHitEE,"nSimHitEE/I");
    tree->Branch("SimHitEE_energy",&SimHitEE_energy);
    tree->Branch("SimHitEE_energyEM",&SimHitEE_energyEM);
    tree->Branch("SimHitEE_energyHAD",&SimHitEE_energyHAD);
    tree->Branch("SimHitEE_nCaloHits",&SimHitEE_nCaloHits);
    tree->Branch("SimHitEE_ix",&SimHitEE_ix);
    tree->Branch("SimHitEE_iy",&SimHitEE_iy);
    tree->Branch("SimHitEE_eta",&SimHitEE_eta);
    tree->Branch("SimHitEE_phi",&SimHitEE_phi);
    tree->Branch("SimHitEE_detId",&SimHitEE_detId);
    tree->Branch("SimHitEE_subdetId",&SimHitEE_subdetId);
  }

  if(savePFClustersECAL_){
    tree->Branch("nPFClusterECAL",&nPFClusterECAL);
    tree->Branch("PFClusterECAL_pt",&PFClusterECAL_pt);
    tree->Branch("PFClusterECAL_energy",&PFClusterECAL_energy);
    tree->Branch("PFClusterECAL_correctedEnergy",&PFClusterECAL_correctedEnergy);
    tree->Branch("PFClusterECAL_eta",&PFClusterECAL_eta);
    tree->Branch("PFClusterECAL_phi",&PFClusterECAL_phi);
    tree->Branch("PFClusterECAL_layer",&PFClusterECAL_layer);
    tree->Branch("PFClusterECAL_seedhit_detId",&PFClusterECAL_seedhit_detId);
    tree->Branch("PFClusterECAL_nhits",&PFClusterECAL_nhits);
    tree->Branch("PFClusterECAL_hits_detId",&PFClusterECAL_hits_detId);
    tree->Branch("PFClusterECAL_hits_fraction",&PFClusterECAL_hits_fraction);
    tree->Branch("PFClusterECAL_hits_PFRecHitEB_Idx",&PFClusterECAL_hits_PFRecHitEB_Idx);
    tree->Branch("PFClusterECAL_hits_PFRecHitEE_Idx",&PFClusterECAL_hits_PFRecHitEE_Idx);
    tree->Branch("PFClusterECAL_key",&PFClusterECAL_key);
  }

  if(savePFClustersPS_){
    tree->Branch("nPFClusterPS",&nPFClusterPS);
    tree->Branch("PFClusterPS_pt",&PFClusterPS_pt);
    tree->Branch("PFClusterPS_energy",&PFClusterPS_energy);
    tree->Branch("PFClusterPS_correctedEnergy",&PFClusterPS_correctedEnergy);
    tree->Branch("PFClusterPS_eta",&PFClusterPS_eta);
    tree->Branch("PFClusterPS_phi",&PFClusterPS_phi);
    tree->Branch("PFClusterPS_layer",&PFClusterPS_layer);
    tree->Branch("PFClusterPS_seedhit_detId",&PFClusterPS_seedhit_detId);
    tree->Branch("PFClusterPS_nhits",&PFClusterPS_nhits);
    tree->Branch("PFClusterPS_hits_detId",&PFClusterPS_hits_detId);
    tree->Branch("PFClusterPS_hits_fraction",&PFClusterPS_hits_fraction);
    tree->Branch("PFClusterPS_key",&PFClusterPS_key);
  }

  if(savePFClustersHCAL_){
    tree->Branch("nPFClusterHCAL",&nPFClusterHCAL);
    tree->Branch("PFClusterHCAL_pt",&PFClusterHCAL_pt);
    tree->Branch("PFClusterHCAL_energy",&PFClusterHCAL_energy);
    tree->Branch("PFClusterHCAL_correctedEnergy",&PFClusterHCAL_correctedEnergy);
    tree->Branch("PFClusterHCAL_eta",&PFClusterHCAL_eta);
    tree->Branch("PFClusterHCAL_phi",&PFClusterHCAL_phi);
    tree->Branch("PFClusterHCAL_layer",&PFClusterHCAL_layer);
    tree->Branch("PFClusterHCAL_seedhit_detId",&PFClusterHCAL_seedhit_detId);
    tree->Branch("PFClusterHCAL_nhits",&PFClusterHCAL_nhits);
    tree->Branch("PFClusterHCAL_hits_detId",&PFClusterHCAL_hits_detId);
    tree->Branch("PFClusterHCAL_hits_fraction",&PFClusterHCAL_hits_fraction);
    tree->Branch("PFClusterHCAL_hits_PFRecHitHBHE_Idx",&PFClusterHCAL_hits_PFRecHitHBHE_Idx);
    tree->Branch("PFClusterHCAL_key",&PFClusterHCAL_key);
  }
}

ParticleFlowAnalysisNtuplizer::~ParticleFlowAnalysisNtuplizer(){}

void ParticleFlowAnalysisNtuplizer::beginRun(const edm::Run& run, const edm::EventSetup & es) {}
void ParticleFlowAnalysisNtuplizer::endRun(const edm::Run& run, const edm::EventSetup & es) {}

void ParticleFlowAnalysisNtuplizer::analyze(const Event& iEvent, const EventSetup& iSetup) {
  LogDebug("ParticleFlowAnalysisNtuplizer") <<"START event: "<<iEvent.id().event() <<" in run "<<iEvent.id().run()<<endl;

  run  = iEvent.id().run();
  evt  = iEvent.id().event();
  lumiBlock = iEvent.id().luminosityBlock();
  time = iEvent.time();

  orun = (size_t)run;
  oevt = (size_t)evt;
  olumiBlock = (size_t)lumiBlock;
  otime = (size_t)((iEvent.time().value())>>32);

  //==========================================
  //
  // Reset all variables that we save in TTree
  // **VERY VERY** IMPORTANT to clear vectors.
  // ALWAYS CHECK THAT FOR EVERY VECTOR WE SAVE
  // IN TTree, that it is cleared here first.
  //
  //==========================================
  genPart_pt = -1.f;
  genPart_eta = -9.f;
  genPart_phi = -9.f;
  genPart_mass = -1.f;
  genPart_energy = -1.f;
  genPart_p = -1.f;
  genPart_pdgId = 0;
  genPart_charge = 0;

  genPart_extrapolated_okECAL = false;
  genPart_extrapolated_okHCAL = false;
  genPart_extrapolated_EB_ieta = 0;
  genPart_extrapolated_EB_iphi = 0;
  genPart_extrapolated_EB_eta = -9.f;
  genPart_extrapolated_EB_phi = -9.f;
  genPart_extrapolated_PFRecHitEB_Idx = -1;
  genPart_extrapolated_PFRecHitEB_detId = 0;

  genPart_extrapolated_EE_ix = 0;
  genPart_extrapolated_EE_iy = 0;
  genPart_extrapolated_EE_eta = -9.f;
  genPart_extrapolated_EE_phi = -9.f;
  genPart_extrapolated_PFRecHitEE_Idx = -1;
  genPart_extrapolated_PFRecHitEE_detId = 0;

  genPart_extrapolated_EHCAL_ieta = 0;
  genPart_extrapolated_EHCAL_iphi = 0;
  genPart_extrapolated_EHCAL_depth = 0;
  genPart_extrapolated_EHCAL_eta = -9.f;
  genPart_extrapolated_EHCAL_phi = -9.f;
  genPart_extrapolated_EHCAL_subdetId = -1;
  genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx = -1;
  genPart_extrapolated_EHCAL_PFRecHitHF_Idx = -1;
  genPart_extrapolated_EHCAL_PFRecHitHBHE_detId = 0;
  genPart_extrapolated_EHCAL_PFRecHitHF_detId = 0;

  genPart_extrapolated_HCAL_ieta = 0;
  genPart_extrapolated_HCAL_iphi = 0;
  genPart_extrapolated_HCAL_depth = 0;
  genPart_extrapolated_HCAL_eta = -9.f;
  genPart_extrapolated_HCAL_phi = -9.f;
  genPart_extrapolated_HCAL_subdetId = -1;
  genPart_extrapolated_HCAL_PFRecHitHBHE_Idx = -1;
  genPart_extrapolated_HCAL_PFRecHitHF_Idx = -1;
  genPart_extrapolated_HCAL_PFRecHitHBHE_detId = 0;
  genPart_extrapolated_HCAL_PFRecHitHF_detId = 0;

  genPart_extrapolated_pointECAL_x = 0.f;
  genPart_extrapolated_pointECAL_y = 0.f;
  genPart_extrapolated_pointECAL_z = 0.f;
  genPart_extrapolated_pointECAL_eta = 0.f;
  genPart_extrapolated_pointECAL_phi = 0.f;

  genPart_extrapolated_pointHCAL_x = 0.f;
  genPart_extrapolated_pointHCAL_y = 0.f;
  genPart_extrapolated_pointHCAL_z = 0.f;
  genPart_extrapolated_pointHCAL_eta = 0.f;
  genPart_extrapolated_pointHCAL_phi = 0.f;

  genVertex_x = 0.f;
  genVertex_y = 0.f;
  genVertex_z = 0.f;

  nPFCand = 0;
  PFCand_pt.clear();
  PFCand_eta.clear();
  PFCand_phi.clear();
  PFCand_mass.clear();
  PFCand_energy.clear();
  PFCand_pdgId.clear();
  PFCand_isChgHadIso.clear();
  PFCand_birthId.clear();
  PFCand_dRToGenPart.clear();

  PFCand_has_trk.clear();

  PFCand_trk_pt.clear();
  PFCand_trk_p.clear();
  PFCand_trk_eta.clear();
  PFCand_trk_phi.clear();
  PFCand_trk_algo.clear();

  PFCand_trk_charge.clear();
  PFCand_trk_ptError.clear();
  PFCand_trk_dxy.clear();
  PFCand_trk_dz.clear();
  PFCand_trk_dxyError.clear();
  PFCand_trk_dzError.clear();
  PFCand_trk_chi2.clear();
  PFCand_trk_ndof.clear();
  PFCand_trk_vx.clear();
  PFCand_trk_vy.clear();
  PFCand_trk_vz.clear();

  PFCand_trk_numberOfValidPixelHits.clear();
  PFCand_trk_numberOfValidTrackerHits.clear();
  PFCand_trk_pixelHitOK.clear();
  PFCand_trk_trackerHitOK.clear();

  PFCand_trk_positionAtECALEntrance_x.clear();
  PFCand_trk_positionAtECALEntrance_y.clear();
  PFCand_trk_positionAtECALEntrance_z.clear();
  PFCand_trk_positionAtECALEntrance_eta.clear();
  PFCand_trk_positionAtECALEntrance_phi.clear();

  PFCand_trk_closestECAL_eta.clear();
  PFCand_trk_closestECAL_phi.clear();
  PFCand_trk_closestECAL_detId.clear();
  PFCand_trk_closestECAL_EBieta.clear();
  PFCand_trk_closestECAL_EBiphi.clear();
  PFCand_trk_closestECAL_EEix.clear();
  PFCand_trk_closestECAL_EEiy.clear();
  PFCand_trk_closestHCAL_eta.clear();
  PFCand_trk_closestHCAL_phi.clear();
  PFCand_trk_closestHCAL_detId.clear();
  PFCand_trk_closestHCAL_ieta.clear();
  PFCand_trk_closestHCAL_iphi.clear();

  PFCand_ecalEnergy.clear();
  PFCand_hcalEnergy.clear();
  PFCand_caloEnergy.clear();
  PFCand_rawEcalEnergy.clear();
  PFCand_rawHcalEnergy.clear();
  PFCand_rawCaloEnergy.clear();
  PFCand_hoEnergy.clear();
  PFCand_rawHoEnergy.clear();
  PFCand_pS1Energy.clear();
  PFCand_pS2Energy.clear();
  PFCand_hcalDepth1EFrac.clear();
  PFCand_hcalDepth2EFrac.clear();
  PFCand_hcalDepth3EFrac.clear();
  PFCand_hcalDepth4EFrac.clear();
  PFCand_hcalDepth5EFrac.clear();
  PFCand_hcalDepth6EFrac.clear();
  PFCand_hcalDepth7EFrac.clear();
  PFCand_PFBlock_Idx.clear();
  PFCand_nPFTrack.clear();
  PFCand_nPFClusterHCAL.clear();
  PFCand_nPFClusterECAL.clear();
  PFCand_nPFClusterPS.clear();

  PFCand_nPFTrack.clear();
  PFCand_nPFClusterHCAL.clear();
  PFCand_nPFClusterECAL.clear();
  PFCand_nPFClusterPS.clear();
  PFCand_nPFClusterHF.clear();
  PFCand_nPFClusterSC.clear();
  PFCand_nPFClusterGSF.clear();
  PFCand_nPFClusterBREM.clear();

  PFCand_PFClusterHCAL_Idx.clear();
  PFCand_PFClusterECAL_Idx.clear();
  PFCand_PFClusterPS_Idx.clear();
  PFCand_PFClusterHF_Idx.clear();

  PFCand_nPFTrackInBlock.clear();
  PFCand_nPFClusterHCALInBlock.clear();
  PFCand_nPFClusterECALInBlock.clear();
  PFCand_nPFClusterPSInBlock.clear();
  PFCand_nPFClusterHFInBlock.clear();
  PFCand_nPFClusterSCInBlock.clear();
  PFCand_nPFClusterGSFInBlock.clear();
  PFCand_nPFClusterBREMInBlock.clear();

  PFCand_PFClusterHCALInBlock_Idx.clear();
  PFCand_PFClusterECALInBlock_Idx.clear();
  PFCand_PFClusterPSInBlock_Idx.clear();
  PFCand_PFClusterHFInBlock_Idx.clear();

  Idx_ClosestPFCandChgHad = -1;
  Idx_ClosestPFCandChgHadIso = -1;
  Idx_ClosestPFCandNeuHad = -1;
  Idx_ClosestPFCandPhoton = -1;

  nPFRecHitHBHE = 0;
  PFRecHitHBHE_energy.clear();
  PFRecHitHBHE_ieta.clear();
  PFRecHitHBHE_iphi.clear();
  PFRecHitHBHE_eta.clear();
  PFRecHitHBHE_phi.clear();
  PFRecHitHBHE_depth.clear();
  PFRecHitHBHE_cutThreshold.clear();
  PFRecHitHBHE_hcalRespCorr.clear();
  PFRecHitHBHE_detId.clear();

  nPFRecHitHF = 0;
  PFRecHitHF_energy.clear();
  PFRecHitHF_ieta.clear();
  PFRecHitHF_iphi.clear();
  PFRecHitHF_eta.clear();
  PFRecHitHF_phi.clear();
  PFRecHitHF_depth.clear();
  PFRecHitHF_cutThreshold.clear();
  PFRecHitHF_detId.clear();

  nPFRecHitEB = 0;
  PFRecHitEB_energy.clear();
  PFRecHitEB_ieta.clear();
  PFRecHitEB_iphi.clear();
  PFRecHitEB_tower_ieta.clear();
  PFRecHitEB_tower_iphi.clear();
  PFRecHitEB_approxEta.clear();
  PFRecHitEB_eta.clear();
  PFRecHitEB_phi.clear();
  PFRecHitEB_cutThreshold.clear();
  PFRecHitEB_detId.clear();

  nPFRecHitEE = 0;
  PFRecHitEE_energy.clear();
  PFRecHitEE_ix.clear();
  PFRecHitEE_iy.clear();
  PFRecHitEE_eta.clear();
  PFRecHitEE_phi.clear();
  PFRecHitEE_cutThreshold.clear();
  PFRecHitEE_detId.clear();

  nG4SimHitPCaloHitHBHE = 0;
  G4SimHitPCaloHitHBHE_energy.clear();
  G4SimHitPCaloHitHBHE_energyEM.clear();
  G4SimHitPCaloHitHBHE_energyHAD.clear();
  G4SimHitPCaloHitHBHE_samplingFactor.clear();
  G4SimHitPCaloHitHBHE_ieta.clear();
  G4SimHitPCaloHitHBHE_iphi.clear();
  G4SimHitPCaloHitHBHE_depth.clear();
  G4SimHitPCaloHitHBHE_eta.clear();
  G4SimHitPCaloHitHBHE_phi.clear();
  G4SimHitPCaloHitHBHE_time.clear();
  G4SimHitPCaloHitHBHE_detId.clear();
  G4SimHitPCaloHitHBHE_subdetId.clear();

  nG4SimHitPCaloHitEB = 0;
  G4SimHitPCaloHitEB_energy.clear();
  G4SimHitPCaloHitEB_energyEM.clear();
  G4SimHitPCaloHitEB_energyHAD.clear();
  G4SimHitPCaloHitEB_ieta.clear();
  G4SimHitPCaloHitEB_iphi.clear();
  G4SimHitPCaloHitEB_depth.clear();
  G4SimHitPCaloHitEB_eta.clear();
  G4SimHitPCaloHitEB_phi.clear();
  G4SimHitPCaloHitEB_time.clear();
  G4SimHitPCaloHitEB_detId.clear();
  G4SimHitPCaloHitEB_subdetId.clear();

  nG4SimHitPCaloHitEE = 0;
  G4SimHitPCaloHitEE_energy.clear();
  G4SimHitPCaloHitEE_energyEM.clear();
  G4SimHitPCaloHitEE_energyHAD.clear();
  G4SimHitPCaloHitEE_ix.clear();
  G4SimHitPCaloHitEE_iy.clear();
  G4SimHitPCaloHitEE_depth.clear();
  G4SimHitPCaloHitEE_eta.clear();
  G4SimHitPCaloHitEE_phi.clear();
  G4SimHitPCaloHitEE_time.clear();
  G4SimHitPCaloHitEE_detId.clear();
  G4SimHitPCaloHitEE_subdetId.clear();

  nSimHitHBHE = 0;
  SimHitHBHE_energy.clear();
  SimHitHBHE_energyEM.clear();
  SimHitHBHE_energyHAD.clear();
  SimHitHBHE_samplingFactor.clear();
  SimHitHBHE_nCaloHits.clear();
  SimHitHBHE_ieta.clear();
  SimHitHBHE_iphi.clear();
  SimHitHBHE_depth.clear();
  SimHitHBHE_eta.clear();
  SimHitHBHE_phi.clear();
  SimHitHBHE_detId.clear();
  SimHitHBHE_subdetId.clear();
  SimHitHBHE_energy_25ns.clear();
  SimHitHBHE_energyEM_25ns.clear();
  SimHitHBHE_energyHAD_25ns.clear();
  SimHitHBHE_nCaloHits_25ns.clear();

  nSimHitEB = 0;
  SimHitEB_energy.clear();
  SimHitEB_energyEM.clear();
  SimHitEB_energyHAD.clear();
  SimHitEB_nCaloHits.clear();
  SimHitEB_ieta.clear();
  SimHitEB_iphi.clear();
  SimHitEB_eta.clear();
  SimHitEB_phi.clear();
  SimHitEB_detId.clear();
  SimHitEB_subdetId.clear();

  nSimHitEE = 0;
  SimHitEE_energy.clear();
  SimHitEE_energyEM.clear();
  SimHitEE_energyHAD.clear();
  SimHitEE_nCaloHits.clear();
  SimHitEE_ix.clear();
  SimHitEE_iy.clear();
  SimHitEE_eta.clear();
  SimHitEE_phi.clear();
  SimHitEE_detId.clear();
  SimHitEE_subdetId.clear();

  nPFClusterECAL=0;
  PFClusterECAL_pt.clear();
  PFClusterECAL_energy.clear();
  PFClusterECAL_correctedEnergy.clear();
  PFClusterECAL_eta.clear();
  PFClusterECAL_phi.clear();
  PFClusterECAL_layer.clear();
  PFClusterECAL_nhits.clear();
  PFClusterECAL_seedhit_detId.clear();
  PFClusterECAL_hits_detId.clear();
  PFClusterECAL_hits_fraction.clear();
  PFClusterECAL_hits_PFRecHitEB_Idx.clear();
  PFClusterECAL_hits_PFRecHitEE_Idx.clear();
  PFClusterECAL_key.clear();

  nPFClusterPS=0;
  PFClusterPS_pt.clear();
  PFClusterPS_energy.clear();
  PFClusterPS_correctedEnergy.clear();
  PFClusterPS_eta.clear();
  PFClusterPS_phi.clear();
  PFClusterPS_layer.clear();
  PFClusterPS_nhits.clear();
  PFClusterPS_seedhit_detId.clear();
  PFClusterPS_hits_detId.clear();
  PFClusterPS_hits_fraction.clear();
  PFClusterPS_key.clear();

  nPFClusterHCAL=0;
  PFClusterHCAL_pt.clear();
  PFClusterHCAL_energy.clear();
  PFClusterHCAL_correctedEnergy.clear();
  PFClusterHCAL_eta.clear();
  PFClusterHCAL_phi.clear();
  PFClusterHCAL_layer.clear();
  PFClusterHCAL_nhits.clear();
  PFClusterHCAL_seedhit_detId.clear();
  PFClusterHCAL_hits_detId.clear();
  PFClusterHCAL_hits_fraction.clear();
  PFClusterHCAL_hits_PFRecHitHBHE_Idx.clear();
  PFClusterHCAL_key.clear();

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
  // Get GenParticle container and use spr::propagateCALO
  // (https://github.com/cms-sw/cmssw/blob/CMSSW_14_0_19_patch2/Calibration/IsolatedParticles/src/CaloPropagateTrack.cc#L429)
  // to estimate the point of entry of GenParticles with ECAL and HCAL.
  // Since we are analysis single particle guns, the container should *always* have EXACTLY one entry.
  //
  //===================================================================================================
  Handle<GenParticleCollection> genParticlesHandle;
  iEvent.getByToken(tokengenParticles_, genParticlesHandle);

  GlobalPoint posVec, posECAL;
  DetId genPart_detIdECAL, genPart_detIdEHCAL, genPart_detIdHCAL;

  std::vector<spr::propagatedGenParticleID> genPartIDs = spr::propagateCALO(genParticlesHandle, pdt, geometry, bField, /*etaMax=*/ 5.0, false);

  if (genPartIDs.size() == 1){
    reco::GenParticleCollection::const_iterator p = genPartIDs[0].trkItr;
    genPart_pt     = p->pt();
    genPart_eta    = p->eta();
    genPart_phi    = p->phi();
    genPart_mass   = p->mass();
    genPart_energy = p->energy();
    genPart_p      = p->p();
    genPart_pdgId  = p->pdgId();
    genPart_charge = p->charge();
    genPart_extrapolated_okECAL = genPartIDs[0].okECAL;
    genPart_extrapolated_okHCAL = genPartIDs[0].okHCAL;

    //
    // Propagate genpart trajectory to ECAL
    //
    if(genPartIDs[0].okECAL){
      genPart_detIdECAL = genPartIDs[0].detIdECAL;
      genPart_detIdEHCAL = genPartIDs[0].detIdEHCAL;

      //// if ECAL Barrel (EB)
      if (genPartIDs[0].detIdECAL.subdetId() == EcalSubdetector::EcalBarrel){
        EBDetId theECalEBDetId(genPart_detIdECAL);
        genPart_extrapolated_EB_ieta = theECalEBDetId.ieta();
        genPart_extrapolated_EB_iphi = theECalEBDetId.iphi();
        auto cellGeometry = geometry->getSubdetectorGeometry(theECalEBDetId)->getGeometry(theECalEBDetId);
        genPart_extrapolated_EB_eta = cellGeometry->getPosition().eta();
        genPart_extrapolated_EB_phi = cellGeometry->getPosition().phi();
        genPart_extrapolated_PFRecHitEB_detId = genPart_detIdECAL;
      }
      //// if ECAL Endcap (EE)
      else if(genPartIDs[0].detIdECAL.subdetId() == EcalSubdetector::EcalEndcap){
        EEDetId theECalEEDetId(genPart_detIdECAL);
        genPart_extrapolated_EE_ix = theECalEEDetId.ix();
        genPart_extrapolated_EE_iy = theECalEEDetId.iy();
        auto cellGeometry = geometry->getSubdetectorGeometry(theECalEEDetId)->getGeometry(theECalEEDetId);
        genPart_extrapolated_EE_eta = cellGeometry->getPosition().eta();
        genPart_extrapolated_EE_phi = cellGeometry->getPosition().phi();
        genPart_extrapolated_PFRecHitEE_detId = genPart_detIdECAL;
      }

      genPart_extrapolated_pointECAL_x = genPartIDs[0].pointECAL.x();
      genPart_extrapolated_pointECAL_y = genPartIDs[0].pointECAL.y();
      genPart_extrapolated_pointECAL_z = genPartIDs[0].pointECAL.z();
      genPart_extrapolated_pointECAL_eta = genPartIDs[0].pointECAL.eta();
      genPart_extrapolated_pointECAL_phi = genPartIDs[0].pointECAL.phi();

      // This would be the closest HCAL cell
      // with respect to the point where genparticle
      // crosses into ECAL
      HcalDetId theHcalDetId(genPart_detIdEHCAL);
      genPart_extrapolated_EHCAL_subdetId = genPartIDs[0].detIdEHCAL.subdetId();
      genPart_extrapolated_EHCAL_ieta = theHcalDetId.ieta();
      genPart_extrapolated_EHCAL_iphi = theHcalDetId.iphi();
      genPart_extrapolated_EHCAL_depth = theHcalDetId.depth();
       if (genPart_extrapolated_EHCAL_subdetId == HcalSubdetector::HcalBarrel){
        genPart_extrapolated_EHCAL_eta  = cellGeometryHB->getPosition(theHcalDetId).eta();
        genPart_extrapolated_EHCAL_phi  = cellGeometryHB->getPosition(theHcalDetId).phi();
      }
      else if (genPart_extrapolated_EHCAL_subdetId == HcalSubdetector::HcalEndcap){
        genPart_extrapolated_EHCAL_eta  = cellGeometryHE->getPosition(theHcalDetId).eta();
        genPart_extrapolated_EHCAL_phi  = cellGeometryHE->getPosition(theHcalDetId).phi();
      }
      else if (genPart_extrapolated_EHCAL_subdetId == HcalSubdetector::HcalForward){
        genPart_extrapolated_EHCAL_eta  = cellGeometryHF->getPosition(theHcalDetId).eta();
        genPart_extrapolated_EHCAL_phi  = cellGeometryHF->getPosition(theHcalDetId).phi();
      }

      if(genPart_extrapolated_EHCAL_subdetId == HcalSubdetector::HcalBarrel || genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalEndcap){
        genPart_extrapolated_EHCAL_PFRecHitHBHE_detId = genPart_detIdEHCAL;
      }
      else if(genPart_extrapolated_EHCAL_subdetId == HcalSubdetector::HcalForward){
        genPart_extrapolated_EHCAL_PFRecHitHF_detId = genPart_detIdEHCAL;
      }
    }

    //
    // Propagate genpart trajectory to HCAL HBHE
    //
    if(genPartIDs[0].okHCAL){
      genPart_detIdHCAL = genPartIDs[0].detIdHCAL;
      HcalDetId theHcalDetId(genPart_detIdHCAL);
      genPart_extrapolated_HCAL_subdetId = genPartIDs[0].detIdHCAL.subdetId();
      genPart_extrapolated_HCAL_ieta = theHcalDetId.ieta();
      genPart_extrapolated_HCAL_iphi = theHcalDetId.iphi();
      genPart_extrapolated_HCAL_depth = theHcalDetId.depth();
      if (genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalBarrel){
        genPart_extrapolated_HCAL_eta  = cellGeometryHB->getPosition(theHcalDetId).eta();
        genPart_extrapolated_HCAL_phi  = cellGeometryHB->getPosition(theHcalDetId).phi();
      }
      else if (genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalEndcap){
        genPart_extrapolated_HCAL_eta  = cellGeometryHE->getPosition(theHcalDetId).eta();
        genPart_extrapolated_HCAL_phi  = cellGeometryHE->getPosition(theHcalDetId).phi();
      }
      else if (genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalForward){
        genPart_extrapolated_EHCAL_eta  = cellGeometryHF->getPosition(theHcalDetId).eta();
        genPart_extrapolated_EHCAL_phi  = cellGeometryHF->getPosition(theHcalDetId).phi();
      }

      genPart_extrapolated_pointHCAL_x = genPartIDs[0].pointHCAL.x();
      genPart_extrapolated_pointHCAL_y = genPartIDs[0].pointHCAL.y();
      genPart_extrapolated_pointHCAL_z = genPartIDs[0].pointHCAL.z();
      genPart_extrapolated_pointHCAL_eta = genPartIDs[0].pointHCAL.eta();
      genPart_extrapolated_pointHCAL_phi = genPartIDs[0].pointHCAL.phi();

      if(genPart_extrapolated_EHCAL_subdetId == HcalSubdetector::HcalBarrel || genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalEndcap){
        genPart_extrapolated_HCAL_PFRecHitHBHE_detId = genPart_detIdHCAL;
      }
      else if(genPart_extrapolated_EHCAL_subdetId == HcalSubdetector::HcalForward){
        genPart_extrapolated_HCAL_PFRecHitHF_detId = genPart_detIdHCAL;
      }
    }
  }else{
    std::cout << "No genParts" << std::endl;
    return;
  }
  const reco::GenParticle& genPart = (*(genPartIDs[0].trkItr));

  //==============================================================
  //
  // FIKRI:
  // - Loop over PF candidates
  // - Consider PF ChgHads and PF NeutralHadrons and PFPhotons.
  // - Calculate DeltaR around the genpart and PF candidates,
  // apply cuts and save each of them in a vector of PF candidates
  // - Sort by dR2
  //
  //==============================================================
  // get PFCandidates
  Handle<PFCandidateCollection> pfCandidatesHandle;
  iEvent.getByToken(tokenPFCandidates_, pfCandidatesHandle);
  size_t n_pfCandidates = pfCandidatesHandle->size();
  std::vector<reco::PFCandidate> pfCandidates = *pfCandidatesHandle;
  //
  edm::Handle<edm::ValueMap<bool>> chargedHadronIsolationHandle;
  iEvent.getByToken(chargedHadronIsolationToken_, chargedHadronIsolationHandle);
  const edm::ValueMap<bool>& chargedHadronIsolation = *(chargedHadronIsolationHandle.product());

  //
  // Save in vector of pairs where the first value is the dR2 between PFCandidate and genPart
  //
  std::vector<std::pair<float, size_t>> pfCand_Selected_Idx;

  for (size_t icand = 0; icand < n_pfCandidates; ++icand){
    const reco::PFCandidate& pfCand = pfCandidates.at(icand);
    float dR2 = reco::deltaR2(genPart.eta(), genPart.phi(), pfCand.eta(), pfCand.phi());

    bool isChgHad = fabs(pfCand.pdgId()) == 211 && dR2 <= (dRMaxGenPartToPFCandChgHad_*dRMaxGenPartToPFCandChgHad_);
    bool isNeuHad = pfCand.pdgId() == 130 && dR2 <= (dRMaxGenPartToPFCandNeuHad_*dRMaxGenPartToPFCandNeuHad_);
    bool isPhoton = pfCand.pdgId() == 22 && dR2 <= (dRMaxGenPartToPFCandPhoton_*dRMaxGenPartToPFCandPhoton_);

    if (isChgHad || isNeuHad || isPhoton){
      pfCand_Selected_Idx.push_back(make_pair(dR2, icand));
    }
  }

  //
  // Sort by
  //
  auto dR2Comparator = [](const std::pair<float, size_t> lhs, const std::pair<float, size_t> rhs){
    return lhs.first < rhs.first;
  };
  std::stable_sort(pfCand_Selected_Idx.begin(), pfCand_Selected_Idx.end(), dR2Comparator);

  //==========================================
  //
  // Loop over PFCandidates
  //
  //==========================================
  nPFCand = pfCand_Selected_Idx.size();

  // std::vector<std::vector<unsigned int>> PFCand_PFElementTrack_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterHCAL_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterECAL_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterPS_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterHF_keys;
  // std::vector<std::vector<unsigned int>> PFCand_PFElementClusterSC_keys;
  // std::vector<std::vector<unsigned int>> PFCand_PFElementGSF_keys;
  // std::vector<std::vector<unsigned int>> PFCand_PFElementBREM_keys;

  // std::vector<std::vector<unsigned int>> PFCand_PFElementTrackInBlock_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterHCALInBlock_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterECALInBlock_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterPSInBlock_keys;
  std::vector<std::vector<unsigned int>> PFCand_PFElementClusterHFInBlock_keys;
  // std::vector<std::vector<unsigned int>> PFCand_PFElementClusterSCInBlock_keys;
  // std::vector<std::vector<unsigned int>> PFCand_PFElementGSFInBlock_keys;
  // std::vector<std::vector<unsigned int>> PFCand_PFElementBREMInBlock_keys;

  for(size_t icandsel=0 ; icandsel < pfCand_Selected_Idx.size(); icandsel++){
    const reco::PFCandidate& pfCand = pfCandidates.at(pfCand_Selected_Idx[icandsel].second);
    reco::PFCandidateRef pfCandRef(pfCandidatesHandle,pfCand_Selected_Idx[icandsel].second);


    PFCand_pt.push_back(pfCand.pt());
    PFCand_eta.push_back(pfCand.eta());
    PFCand_phi.push_back(pfCand.phi());
    PFCand_mass.push_back(pfCand.mass());
    PFCand_energy.push_back(pfCand.energy());
    PFCand_pdgId.push_back(pfCand.pdgId());
    unsigned birthId = 0;
    birthId = pfCand.getRecoLocationIdx();
    PFCand_birthId.push_back(birthId);
    PFCand_dRToGenPart.push_back(TMath::Sqrt(pfCand_Selected_Idx[icandsel].first));

    bool isChgHadIso = false;
    if (fabs(pfCand.pdgId()) == 211){
      isChgHadIso = chargedHadronIsolation[pfCandRef];
    }
    PFCand_isChgHadIso.push_back(int(isChgHadIso));

    if(Idx_ClosestPFCandChgHad == -1 && fabs(pfCand.pdgId()) == 211){
      Idx_ClosestPFCandChgHad = icandsel;
    }
    if(Idx_ClosestPFCandChgHadIso == -1 && fabs(pfCand.pdgId()) == 211 && isChgHadIso){
      Idx_ClosestPFCandChgHadIso = icandsel;
    }
    if(Idx_ClosestPFCandNeuHad == -1 && pfCand.pdgId() == 130){
      Idx_ClosestPFCandNeuHad = icandsel;
    }
    if(Idx_ClosestPFCandPhoton == -1 && pfCand.pdgId() == 22){
      Idx_ClosestPFCandPhoton = icandsel;
    }


    float track_pt = -1.f;
    float track_p = -1.f;
    float track_eta = -9.f;
    float track_phi = -9.f;
    int track_algo = -1;

    bool has_track = 0;
    int trk_charge = 0;
    float trk_ptError = -1.f;
    float trk_dxy = 0.f;
    float trk_dz = 0.f;
    float trk_dxyError = 0.f;
    float trk_dzError = 0.f;
    float trk_chi2 = 0.f;
    float trk_ndof = 0.f;
    float trk_vx = 0.f;
    float trk_vy = 0.f;
    float trk_vz = 0.f;

    // unsigned int tobN = 0;
    // unsigned int tecN = 0;
    // unsigned int tibN = 0;
    // unsigned int tidN = 0;
    // unsigned int pxbN = 0;
    // unsigned int pxdN = 0;
    int track_numberOfValidPixelHits = -1;
    int track_numberOfValidTrackerHits = -1;
    bool track_pixelHitOK = false;
    bool track_trackerHitOK = false;
 
    float trk_positionAtECALEntrance_x = 0.f;
    float trk_positionAtECALEntrance_y = 0.f;
    float trk_positionAtECALEntrance_z = 0.f;
    float trk_positionAtECALEntrance_eta = -9.f;
    float trk_positionAtECALEntrance_phi = -9.f;

    float trk_closestECAL_eta = -9.f;
    float trk_closestECAL_phi = -9.f;
    int trk_closestECAL_detId = 0;
    int trk_closestECAL_EBieta = 0;
    int trk_closestECAL_EBiphi = 0;
    int trk_closestECAL_EEix = 0;
    int trk_closestECAL_EEiy = 0;
    float trk_closestHCAL_eta = -9.f;
    float trk_closestHCAL_phi = -9.f;
    int trk_closestHCAL_detId = 0;
    int trk_closestHCAL_ieta = 0;
    int trk_closestHCAL_iphi = 0;

    const reco::Track* track = pfCand.bestTrack();

    if (track){
      has_track = true;
      track_pt = track->pt();
      track_p = track->p();
      track_eta = track->eta();
      track_phi = track->phi();
      track_algo = track->algo();

      trk_charge = track->charge();
      trk_ptError = track->ptError();
      trk_dxy = track->dxy();
      trk_dz = track->dz();
      trk_dxyError = track->dxyError();
      trk_dzError = track->dzError();
      trk_chi2 = track->chi2();
      trk_ndof = track->ndof();
      trk_vx = track->vx();
      trk_vy = track->vy();
      trk_vz = track->vz();

      //
      // TODO: Must use reco::PFRecTrack object
      //
      // reco::PFTrajectoryPoint::LayerType ecalEntrance = reco::PFTrajectoryPoint::ECALEntrance;
      // const reco::PFTrajectoryPoint& tpatecal = pfRecTrack.extrapolatedPoint( ecalEntrance );

      const reco::HitPattern& hp = track->hitPattern();
      track_numberOfValidPixelHits = hp.numberOfValidPixelHits();
      track_numberOfValidTrackerHits = hp.numberOfValidTrackerHits();
      // switch ( track_algo ) {
      //   case TrackBase::initialStep:
      //     tobN += hp.numberOfValidStripTOBHits();
      //     tecN += hp.numberOfValidStripTECHits();
      //     tibN += hp.numberOfValidStripTIBHits();
      //     tidN += hp.numberOfValidStripTIDHits();
      //     pxbN += hp.numberOfValidPixelBarrelHits();
      //     pxdN += hp.numberOfValidPixelEndcapHits();
      //     break;
      //   case TrackBase::lowPtQuadStep:
      //   case TrackBase::highPtTripletStep:
      //   case TrackBase::lowPtTripletStep:
      //   case TrackBase::detachedQuadStep:
      //   case TrackBase::detachedTripletStep:
      //   case TrackBase::pixelPairStep:
      //   case TrackBase::mixedTripletStep:
      //   case TrackBase::pixelLessStep:
      //   case TrackBase::tobTecStep:
      //   case TrackBase::jetCoreRegionalStep:
      //   default:
      //     break;
      // }

      // Number of pixel hits (eta-dependent cut)
      track_pixelHitOK = track_numberOfValidPixelHits >= nPixMin_;

      // Number of tracker hits (eta-dependent cut)
      for ( unsigned int ieta=0; ieta<etaMinForHitMin_.size(); ++ieta ) {
        track_trackerHitOK = fabs(track_eta) > etaMinForHitMin_[ieta] && fabs(track_eta) <= etaMaxForHitMin_[ieta] && track_numberOfValidPixelHits+track_numberOfValidTrackerHits > nHitMin_[ieta];
        if ( track_trackerHitOK ) break;
      }

      const math::XYZPointF& ecalPoint = pfCand.positionAtECALEntrance();
      GlobalPoint ecalGPoint(ecalPoint.X(), ecalPoint.Y(), ecalPoint.Z());

      trk_positionAtECALEntrance_x = ecalPoint.X();
      trk_positionAtECALEntrance_y = ecalPoint.Y();
      trk_positionAtECALEntrance_z = ecalPoint.Z();
      trk_positionAtECALEntrance_eta = ecalPoint.eta();
      trk_positionAtECALEntrance_phi = ecalPoint.phi();

      HcalDetId closestHCalDetId;
      if (fabs(track_eta) < 1.392f){
        closestHCalDetId = cellGeometryHB->getClosestCell(ecalGPoint);
        trk_closestHCAL_eta = cellGeometryHB->getPosition(closestHCalDetId).eta();
        trk_closestHCAL_phi = cellGeometryHB->getPosition(closestHCalDetId).phi();
      }
      else if (fabs(track_eta) < 3.0f){
        closestHCalDetId = cellGeometryHE->getClosestCell(ecalGPoint);
        trk_closestHCAL_eta = cellGeometryHE->getPosition(closestHCalDetId).eta();
        trk_closestHCAL_phi = cellGeometryHE->getPosition(closestHCalDetId).phi();
      }
      trk_closestHCAL_ieta = closestHCalDetId.ieta();
      trk_closestHCAL_iphi = closestHCalDetId.iphi();

      EBDetId closestEBDetId;
      EEDetId closestEEDetId;

      if (fabs(track_eta) < 1.479f){
        const CaloSubdetectorGeometry* barrelECALGeom = geometry->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
        closestEBDetId = barrelECALGeom->getClosestCell(ecalGPoint);
        auto cellGeometry = geometry->getSubdetectorGeometry(closestEBDetId)->getGeometry(closestEBDetId);
        trk_closestECAL_EBieta = closestEBDetId.ieta();
        trk_closestECAL_EBiphi = closestEBDetId.iphi();
        trk_closestECAL_detId = closestEBDetId;
      } else if (fabs(track_eta) < 3.0f){
        const CaloSubdetectorGeometry* endcapECALGeom = geometry->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
        closestEEDetId = endcapECALGeom->getClosestCell(ecalGPoint);
        auto cellGeometry = geometry->getSubdetectorGeometry(closestEEDetId)->getGeometry(closestEEDetId);
        trk_closestECAL_EEix = closestEEDetId.ix();
        trk_closestECAL_EEiy = closestEEDetId.iy();
        trk_closestECAL_detId = closestEEDetId;
      }
    }

    PFCand_has_trk.push_back(int(has_track));

    PFCand_trk_pt.push_back(track_pt);
    PFCand_trk_p.push_back(track_p);
    PFCand_trk_eta.push_back(track_eta);
    PFCand_trk_phi.push_back(track_phi);
    PFCand_trk_algo.push_back(track_algo);

    PFCand_trk_charge.push_back(trk_charge);
    PFCand_trk_ptError.push_back(trk_ptError);
    PFCand_trk_dxy.push_back(trk_dxy);
    PFCand_trk_dz.push_back(trk_dz);
    PFCand_trk_dxyError.push_back(trk_dxyError);
    PFCand_trk_dzError.push_back(trk_dzError);
    PFCand_trk_chi2.push_back(trk_chi2);
    PFCand_trk_ndof.push_back(trk_ndof);
    PFCand_trk_vx.push_back(trk_vx);
    PFCand_trk_vy.push_back(trk_vy);
    PFCand_trk_vz.push_back(trk_vz);

    PFCand_trk_numberOfValidPixelHits.push_back(track_numberOfValidPixelHits);
    PFCand_trk_numberOfValidTrackerHits.push_back(track_numberOfValidTrackerHits);
    PFCand_trk_pixelHitOK.push_back(int(track_pixelHitOK));
    PFCand_trk_trackerHitOK.push_back(int(track_trackerHitOK));
    PFCand_trk_positionAtECALEntrance_x.push_back(trk_positionAtECALEntrance_x);
    PFCand_trk_positionAtECALEntrance_y.push_back(trk_positionAtECALEntrance_y);
    PFCand_trk_positionAtECALEntrance_z.push_back(trk_positionAtECALEntrance_z);
    PFCand_trk_positionAtECALEntrance_eta.push_back(trk_positionAtECALEntrance_eta);
    PFCand_trk_positionAtECALEntrance_phi.push_back(trk_positionAtECALEntrance_phi);

    PFCand_trk_closestECAL_eta.push_back(trk_closestECAL_eta);
    PFCand_trk_closestECAL_phi.push_back(trk_closestECAL_phi);
    PFCand_trk_closestECAL_detId.push_back(trk_closestECAL_detId);
    PFCand_trk_closestECAL_EBieta.push_back(trk_closestECAL_EBieta);
    PFCand_trk_closestECAL_EBiphi.push_back(trk_closestECAL_EBiphi);
    PFCand_trk_closestECAL_EEix.push_back(trk_closestECAL_EEix);
    PFCand_trk_closestECAL_EEiy.push_back(trk_closestECAL_EEiy);
    PFCand_trk_closestHCAL_eta.push_back(trk_closestHCAL_eta);
    PFCand_trk_closestHCAL_phi.push_back(trk_closestHCAL_phi);
    PFCand_trk_closestHCAL_detId.push_back(trk_closestHCAL_detId);
    PFCand_trk_closestHCAL_ieta.push_back(trk_closestHCAL_ieta);
    PFCand_trk_closestHCAL_iphi.push_back(trk_closestHCAL_iphi);

    PFCand_ecalEnergy.push_back(pfCand.ecalEnergy());
    PFCand_hcalEnergy.push_back(pfCand.hcalEnergy());
    PFCand_rawEcalEnergy.push_back(pfCand.rawEcalEnergy());
    PFCand_rawHcalEnergy.push_back(pfCand.rawHcalEnergy());
    PFCand_hoEnergy.push_back(pfCand.hoEnergy());
    PFCand_rawHoEnergy.push_back(pfCand.rawHoEnergy());
    PFCand_caloEnergy.push_back(pfCand.ecalEnergy()+pfCand.hcalEnergy());
    PFCand_rawCaloEnergy.push_back(pfCand.rawEcalEnergy()+pfCand.rawHcalEnergy());
    PFCand_pS1Energy.push_back(pfCand.pS1Energy());
    PFCand_pS2Energy.push_back(pfCand.pS2Energy());
    float hcalDepth1EFrac = -1.f;
    float hcalDepth2EFrac = -1.f;
    float hcalDepth3EFrac = -1.f;
    float hcalDepth4EFrac = -1.f;
    float hcalDepth5EFrac = -1.f;
    float hcalDepth6EFrac = -1.f;
    float hcalDepth7EFrac = -1.f;
    if (pfCand.pdgId() == 130 || fabs(pfCand.pdgId()) == 211){
      hcalDepth1EFrac = pfCand.hcalDepthEnergyFraction(1);
      hcalDepth2EFrac = pfCand.hcalDepthEnergyFraction(2);
      hcalDepth3EFrac = pfCand.hcalDepthEnergyFraction(3);
      hcalDepth4EFrac = pfCand.hcalDepthEnergyFraction(4);
      hcalDepth5EFrac = pfCand.hcalDepthEnergyFraction(5);
      hcalDepth6EFrac = pfCand.hcalDepthEnergyFraction(6);
      hcalDepth7EFrac = pfCand.hcalDepthEnergyFraction(7);
    }
    PFCand_hcalDepth1EFrac.push_back(hcalDepth1EFrac);
    PFCand_hcalDepth2EFrac.push_back(hcalDepth2EFrac);
    PFCand_hcalDepth3EFrac.push_back(hcalDepth3EFrac);
    PFCand_hcalDepth4EFrac.push_back(hcalDepth4EFrac);
    PFCand_hcalDepth5EFrac.push_back(hcalDepth5EFrac);
    PFCand_hcalDepth6EFrac.push_back(hcalDepth6EFrac);
    PFCand_hcalDepth7EFrac.push_back(hcalDepth7EFrac);

    //=====================================================================
    //
    // Find the corresponding PF block elements for this pfCandidate
    //
    //=====================================================================
    int pfCand_blockIdx = -1;

    unsigned int nTrack = 0;
    unsigned int nClusterHCAL = 0;
    unsigned int nClusterECAL = 0;
    unsigned int nClusterPS = 0;
    unsigned int nClusterHF = 0;
    unsigned int nClusterSC = 0;
    unsigned int nGSF = 0;
    unsigned int nBREM = 0;

    unsigned int nTrackInBlock = 0;
    unsigned int nClusterHCALInBlock = 0;
    unsigned int nClusterECALInBlock = 0;
    unsigned int nClusterPSInBlock = 0;
    unsigned int nClusterHFInBlock = 0;
    unsigned int nClusterSCInBlock = 0;
    unsigned int nGSFInBlock = 0;
    unsigned int nBREMInBlock = 0;

    // std::vector<unsigned int> PFElementTrack_keys;
    std::vector<unsigned int> PFElementClusterHCAL_keys;
    std::vector<unsigned int> PFElementClusterECAL_keys;
    std::vector<unsigned int> PFElementClusterPS_keys;
    std::vector<unsigned int> PFElementClusterHF_keys;
    // std::vector<unsigned int> PFElementClusterSC_keys;
    // std::vector<unsigned int> PFElementGSF_keys;
    // std::vector<unsigned int> PFElementBREM_keys;

    // std::vector<unsigned int> PFElementTrackInBlock_keys;
    std::vector<unsigned int> PFElementClusterHCALInBlock_keys;
    std::vector<unsigned int> PFElementClusterECALInBlock_keys;
    std::vector<unsigned int> PFElementClusterPSInBlock_keys;
    std::vector<unsigned int> PFElementClusterHFInBlock_keys;
    // std::vector<unsigned int> PFElementClusterSCInBlock_keys;
    // std::vector<unsigned int> PFElementGSFInBlock_keys;
    // std::vector<unsigned int> PFElementBREMInBlock_keys;

    // Get the list of elements that this PF candidate consist of
    const PFCandidate::ElementsInBlocks& pfElements = pfCand.elementsInBlocks();

    //
    // Get the original block where this pfCandidate comes from and the elements in the block
    //
    if (pfElements.size() >= 1){
      pfCand_blockIdx = pfElements[0].first.index();
      const reco::PFBlockRef pfBlockRef = pfElements[0].first;
      const edm::OwnVector<reco::PFBlockElement>& elementsInOriBlock = pfBlockRef->elements();
      PFBlock::LinkData linkData =  pfBlockRef->linkData();

      //
      // Loop over elements in the block where it comes from
      //
      for(unsigned int eBlock=0; eBlock < elementsInOriBlock.size(); eBlock++) {
        PFBlockElement::Type typeInBlock = elementsInOriBlock[eBlock].type();
        switch( typeInBlock ) {
          case PFBlockElement::TRACK:
            nTrackInBlock++;
            break;
          case PFBlockElement::HCAL:
            nClusterHCALInBlock++;
            PFElementClusterHCAL_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
            break;
          case PFBlockElement::ECAL:
            nClusterECALInBlock++;
            PFElementClusterECALInBlock_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
            break;
          case PFBlockElement::PS1:
          case PFBlockElement::PS2:
            nClusterPSInBlock++;
            PFElementClusterPSInBlock_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
            break;
          case PFBlockElement::HFHAD:
          case PFBlockElement::HFEM:
            nClusterHFInBlock++;
            PFElementClusterHFInBlock_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
            break;
          case PFBlockElement::SC:
            nClusterSCInBlock++;
            break;
          case PFBlockElement::GSF:
            nGSFInBlock++;
            break;
          case PFBlockElement::BREM:
            nBREMInBlock++;
            break;
          default:
            break;
        }
        //
        // For a given element in the block, check if this corresponds to the element
        // that make up this PFcandidate
        //
        for(unsigned int e=0; e < pfElements.size(); e++) {
          if (elementsInOriBlock[eBlock].index() ==  pfElements[e].second){
            switch( typeInBlock ) {
              case PFBlockElement::TRACK:
                nTrack++;
                break;
              case PFBlockElement::HCAL:
                nClusterHCAL++;
                PFElementClusterHCAL_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
                break;
              case PFBlockElement::ECAL:
                nClusterECAL++;
                PFElementClusterECAL_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
                break;
              case PFBlockElement::PS1:
              case PFBlockElement::PS2:
                nClusterPS++;
                PFElementClusterPS_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
                break;
              case PFBlockElement::HFHAD:
              case PFBlockElement::HFEM:
                nClusterHF++;
                PFElementClusterHF_keys.push_back(elementsInOriBlock[eBlock].clusterRef().key());
                break;
              case PFBlockElement::SC:
                nClusterSC++;
                break;
              case PFBlockElement::GSF:
                nGSF++;
                break;
              case PFBlockElement::BREM:
                nBREM++;
                break;
              default:
                break;
            }
          }
        }
      }
    }
    PFCand_PFBlock_Idx.push_back(pfCand_blockIdx);

    PFCand_nPFTrack.push_back(nTrack);
    PFCand_nPFClusterHCAL.push_back(nClusterHCAL);
    PFCand_nPFClusterECAL.push_back(nClusterECAL);
    PFCand_nPFClusterPS.push_back(nClusterPS);
    PFCand_nPFClusterHF.push_back(nClusterHF);
    PFCand_nPFClusterSC.push_back(nClusterSC);
    PFCand_nPFClusterGSF.push_back(nGSF);
    PFCand_nPFClusterBREM.push_back(nBREM);

    PFCand_nPFTrackInBlock.push_back(nTrackInBlock);
    PFCand_nPFClusterHCALInBlock.push_back(nClusterHCALInBlock);
    PFCand_nPFClusterECALInBlock.push_back(nClusterECALInBlock);
    PFCand_nPFClusterPSInBlock.push_back(nClusterPSInBlock);
    PFCand_nPFClusterHFInBlock.push_back(nClusterHFInBlock);
    PFCand_nPFClusterSCInBlock.push_back(nClusterSCInBlock);
    PFCand_nPFClusterGSFInBlock.push_back(nGSFInBlock);
    PFCand_nPFClusterBREMInBlock.push_back(nBREMInBlock);

    // PFCand_PFElementTrack_keys.push_back(PFElementTrack_keys);
    PFCand_PFElementClusterHCAL_keys.push_back(PFElementClusterHCAL_keys);
    PFCand_PFElementClusterECAL_keys.push_back(PFElementClusterECAL_keys);
    PFCand_PFElementClusterPS_keys.push_back(PFElementClusterPS_keys);
    PFCand_PFElementClusterHF_keys.push_back(PFElementClusterHF_keys);
    // PFCand_PFElementClusterSC_keys.push_back(PFElementClusterSC_keys);
    // PFCand_PFElementGSF_keys.push_back(PFElementGSF_keys);
    // PFCand_PFElementBREM_keys.push_back(PFElementBREM_keys);

    // PFCand_PFElementTrackInBlock_keys.push_back(PFElementTrackInBlock_keys);
    PFCand_PFElementClusterHCALInBlock_keys.push_back(PFElementClusterHCALInBlock_keys);
    PFCand_PFElementClusterECALInBlock_keys.push_back(PFElementClusterECALInBlock_keys);
    PFCand_PFElementClusterPSInBlock_keys.push_back(PFElementClusterPSInBlock_keys);
    PFCand_PFElementClusterHFInBlock_keys.push_back(PFElementClusterHFInBlock_keys);
    // PFCand_PFElementClusterSCInBlock_keys.push_back(PFElementClusterSCInBlock_keys);
    // PFCand_PFElementGSFInBlock_keys.push_back(PFElementGSFInBlock_keys);
    // PFCand_PFElementBREMInBlock_keys.push_back(PFElementBREMInBlock_keys);

    PFCand_PFClusterHCAL_Idx.push_back(std::vector<int>());
    PFCand_PFClusterECAL_Idx.push_back(std::vector<int>());
    PFCand_PFClusterPS_Idx.push_back(std::vector<int>());
    PFCand_PFClusterHF_Idx.push_back(std::vector<int>());

    PFCand_PFClusterHCALInBlock_Idx.push_back(std::vector<int>());
    PFCand_PFClusterECALInBlock_Idx.push_back(std::vector<int>());
    PFCand_PFClusterPSInBlock_Idx.push_back(std::vector<int>());
    PFCand_PFClusterHFInBlock_Idx.push_back(std::vector<int>());
  }

  //==========================================
  //
  // PFRecHits HCAL: HBHE
  //
  //==========================================
  edm::Handle<std::vector<reco::PFRecHit>> pfRecHitsHBHEHandle;
  iEvent.getByToken(pfRecHitsHBHEToken_, pfRecHitsHBHEHandle);

  auto pfRecHitsHBHE = pfRecHitsHBHEHandle.product();
  size_t n_pfRecHitsHBHE = pfRecHitsHBHE->size();

  nPFRecHitHBHE=0;
  for (size_t idx = 0; idx < n_pfRecHitsHBHE; idx++){
    const reco::PFRecHit& pfrechit = pfRecHitsHBHE->at(idx);
    HcalDetId theHcalDetId(pfrechit.detId());
    double eta = cellGeometryHB->getPosition(theHcalDetId).eta();
    double phi = cellGeometryHB->getPosition(theHcalDetId).phi();
    PFRecHitHBHE_energy.push_back(pfrechit.energy());
    PFRecHitHBHE_ieta.push_back(theHcalDetId.ieta());
    PFRecHitHBHE_iphi.push_back(theHcalDetId.iphi());
    PFRecHitHBHE_eta.push_back(eta);
    PFRecHitHBHE_phi.push_back(phi);
    PFRecHitHBHE_depth.push_back(theHcalDetId.depth());
    const HcalPFCut* cutValue = hcalPFCuts->getValues(pfrechit.detId());
    float thresholdE = cutValue->noiseThreshold();
    PFRecHitHBHE_cutThreshold.push_back(thresholdE);
    float hcalCorr = 1.0f;
    if (hcalRespCorrs->exists(pfrechit.detId()))
      hcalCorr = hcalRespCorrs->getValues(pfrechit.detId())->getValue();
    PFRecHitHBHE_hcalRespCorr.push_back(hcalCorr);
    PFRecHitHBHE_detId.push_back(pfrechit.detId());
    nPFRecHitHBHE++;

    if (genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalBarrel || genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalEndcap){
      if (genPart_detIdEHCAL == pfrechit.detId()){
        genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx = PFRecHitHBHE_ieta.size()-1;
      }
      if (genPart_detIdHCAL == pfrechit.detId()){
        genPart_extrapolated_HCAL_PFRecHitHBHE_Idx = PFRecHitHBHE_ieta.size()-1;
      }
    }
  }

  //==========================================
  //
  // PFRecHits HCAL:HF
  //
  //==========================================
  edm::Handle<std::vector<reco::PFRecHit>> pfRecHitsHFHandle;
  iEvent.getByToken(pfRecHitsHFToken_, pfRecHitsHFHandle);

  auto pfRecHitsHF = pfRecHitsHFHandle.product();
  size_t n_pfRecHitsHF = pfRecHitsHF->size();

  // nPFRecHitHF=0;
  // for (size_t idx = 0; idx < n_pfRecHitsHF; idx++){
  //   const reco::PFRecHit& pfrechit = pfRecHitsHF->at(idx);
  //   HcalDetId theHcalDetId(pfrechit.detId());
  //   double eta = cellGeometryHF->getPosition(theHcalDetId).eta();
  //   double phi = cellGeometryHF->getPosition(theHcalDetId).phi();
  //   PFRecHitHF_energy.push_back(pfrechit.energy());
  //   PFRecHitHF_ieta.push_back(theHcalDetId.ieta());
  //   PFRecHitHF_iphi.push_back(theHcalDetId.iphi());
  //   PFRecHitHF_eta.push_back(eta);
  //   PFRecHitHF_phi.push_back(phi);
  //   PFRecHitHF_depth.push_back(theHcalDetId.depth());
  //   const HcalPFCut *cutValue = hcalPFCuts->getValues(pfrechit.detId());
  //   float thresholdE = cutValue->noiseThreshold();
  //   PFRecHitHF_cutThreshold.push_back(thresholdE);
  //   PFRecHitHF_detId.push_back(pfrechit.detId());
  //   nPFRecHitHF++;
  //   if (genPart_extrapolated_HCAL_subdetId == HcalSubdetector::HcalForward){
  //     if (genPart_detIdEHCAL == pfrechit.detId()){
  //       genPart_extrapolated_EHCAL_PFRecHitHF_Idx = PFRecHitHF_ieta.size()-1;
  //     }
  //     if (genPart_detIdHCAL == pfrechit.detId()){
  //       genPart_extrapolated_HCAL_PFRecHitHF_Idx = PFRecHitHF_ieta.size()-1;
  //     }
  //   }
  // }
  //==========================================
  //
  // PFRecHits ECAL
  //
  //==========================================
  edm::Handle<std::vector<reco::PFRecHit>> pfRecHitsECALHandle;
  iEvent.getByToken(pfRecHitsECALToken_, pfRecHitsECALHandle);

  auto pfRecHitsECAL = pfRecHitsECALHandle.product();
  size_t n_pfRecHitsECAL = pfRecHitsECAL->size();

  nPFRecHitEB=0;
  nPFRecHitEE=0;

  for (size_t idx = 0; idx < n_pfRecHitsECAL; idx++){
    const reco::PFRecHit& pfrechit = pfRecHitsECAL->at(idx);
    DetId detIdRaw(pfrechit.detId());
    //
    // ECAL Barrel (EB)
    //
    if (detIdRaw.det() == DetId::Ecal && detIdRaw.subdetId() == EcalSubdetector::EcalBarrel){
      EBDetId theECalEBDetId(pfrechit.detId());
      auto cellGeometry = geometry->getSubdetectorGeometry(theECalEBDetId)->getGeometry(theECalEBDetId);
      double eta = cellGeometry->getPosition().eta();
      double phi = cellGeometry->getPosition().phi();
      PFRecHitEB_energy.push_back(pfrechit.energy());
      PFRecHitEB_ieta.push_back(theECalEBDetId.ieta());
      PFRecHitEB_iphi.push_back(theECalEBDetId.iphi());
      PFRecHitEB_tower_ieta.push_back(theECalEBDetId.tower_ieta());
      PFRecHitEB_tower_iphi.push_back(theECalEBDetId.tower_iphi());
      PFRecHitEB_approxEta.push_back(theECalEBDetId.approxEta());
      PFRecHitEB_eta.push_back(eta);
      PFRecHitEB_phi.push_back(phi);
      PFRecHitEB_cutThreshold.push_back((*ecalPFCuts)[detIdRaw]);
      PFRecHitEB_detId.push_back(pfrechit.detId());
      nPFRecHitEB++;
      if (genPart_detIdECAL == pfrechit.detId()){
        genPart_extrapolated_PFRecHitEB_Idx = PFRecHitEB_ieta.size()-1;
      }
    }
    //
    // ECAL Endcap (EE)
    //
    else if(detIdRaw.det() == DetId::Ecal && detIdRaw.subdetId() == EcalSubdetector::EcalEndcap){
      EEDetId theECalEEDetId(pfrechit.detId());
      auto cellGeometry = geometry->getSubdetectorGeometry(theECalEEDetId)->getGeometry(theECalEEDetId);
      double eta = cellGeometry->getPosition().eta();
      double phi = cellGeometry->getPosition().phi();
      PFRecHitEE_energy.push_back(pfrechit.energy());
      PFRecHitEE_ix.push_back(theECalEEDetId.ix());
      PFRecHitEE_iy.push_back(theECalEEDetId.iy());
      PFRecHitEE_eta.push_back(eta);
      PFRecHitEE_phi.push_back(phi);
      PFRecHitEE_cutThreshold.push_back((*ecalPFCuts)[detIdRaw]);
      PFRecHitEE_detId.push_back(pfrechit.detId());
      nPFRecHitEE++;
      if (genPart_detIdECAL == pfrechit.detId()){
        genPart_extrapolated_PFRecHitEE_Idx = PFRecHitEE_ix.size()-1;
      }
    }
  }
  //==========================================
  //
  // PCaloHit: HCAL
  //
  //==========================================
  std::map<unsigned int, SimHitInfo> map_simHitInfoHBHE;

  if (saveSimCaloHitHBHE_ || saveSimHitHBHE_){
    auto hcalCond = &iSetup.getData(hcalDBToken_);
    theParameterMap->setDbService(hcalCond);

    edm::Handle<edm::PCaloHitContainer> g4SimHitsPCaloHitsHCALHandle;
    iEvent.getByToken(g4SimHitsPCaloHitsHCALToken_, g4SimHitsPCaloHitsHCALHandle);
    auto g4SimHitsPCaloHitsHCAL = g4SimHitsPCaloHitsHCALHandle.product();
    size_t n_g4SimHitsPCaloHitsHCAL = g4SimHitsPCaloHitsHCAL->size();

    nG4SimHitPCaloHitHBHE=0;
    // nG4SimHitPCaloHitHF=0;

    for (size_t idx = 0; idx < n_g4SimHitsPCaloHitsHCAL; idx++){
      const PCaloHit& caloHit = g4SimHitsPCaloHitsHCAL->at(idx);
      DetId newid = HcalHitRelabeller::relabel(caloHit.id(), hcaloConst);
      HcalDetId theHcalDetId(newid);
      double eta = -9.;
      double phi = -9.;
      if (newid.subdetId() == HcalSubdetector::HcalBarrel){
        eta = cellGeometryHB->getPosition(theHcalDetId).eta();
        phi = cellGeometryHB->getPosition(theHcalDetId).phi();
      }else if(newid.subdetId() == HcalSubdetector::HcalEndcap){
        eta = cellGeometryHE->getPosition(theHcalDetId).eta();
        phi = cellGeometryHE->getPosition(theHcalDetId).phi();
      }
      //
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

        if(saveSimCaloHitHBHE_){
          G4SimHitPCaloHitHBHE_energy.push_back(caloHit.energy());
          G4SimHitPCaloHitHBHE_energyEM.push_back(caloHit.energyEM());
          G4SimHitPCaloHitHBHE_energyHAD.push_back(caloHit.energyHad());
          G4SimHitPCaloHitHBHE_samplingFactor.push_back(samplingFactor);
          G4SimHitPCaloHitHBHE_ieta.push_back(theHcalDetId.ieta());
          G4SimHitPCaloHitHBHE_iphi.push_back(theHcalDetId.iphi());
          G4SimHitPCaloHitHBHE_depth.push_back(theHcalDetId.depth());
          G4SimHitPCaloHitHBHE_eta.push_back(eta);
          G4SimHitPCaloHitHBHE_phi.push_back(phi);
          G4SimHitPCaloHitHBHE_time.push_back(caloHit.time());
          G4SimHitPCaloHitHBHE_detId.push_back(newid.rawId());
          G4SimHitPCaloHitHBHE_subdetId.push_back(newid.subdetId());
          nG4SimHitPCaloHitHBHE++;
        }
      }
    }

    if (saveSimHitHBHE_){
      nSimHitHBHE=0;
      for (std::map<unsigned int, SimHitInfo>::const_iterator it = map_simHitInfoHBHE.cbegin(); it != map_simHitInfoHBHE.cend(); ++it) {
        DetId newid = HcalHitRelabeller::relabel(it->first, hcaloConst);
        HcalDetId theHcalDetId(newid);
        double eta = -9.;
        double phi = -9.;
        if (newid.subdetId() == HcalSubdetector::HcalBarrel){
          eta = cellGeometryHB->getPosition(theHcalDetId).eta();
          phi = cellGeometryHB->getPosition(theHcalDetId).phi();
        }else if(newid.subdetId() == HcalSubdetector::HcalEndcap){
          eta = cellGeometryHE->getPosition(theHcalDetId).eta();
          phi = cellGeometryHE->getPosition(theHcalDetId).phi();
        }
        const HcalSimParameters& pars = dynamic_cast<const HcalSimParameters&>(theParameterMap->simParameters(theHcalDetId));
        float samplingFactor = pars.samplingFactor(theHcalDetId);

        SimHitHBHE_energy.push_back((it->second).energy);
        SimHitHBHE_energyEM.push_back((it->second).energyEM);
        SimHitHBHE_energyHAD.push_back((it->second).energyHAD);
        SimHitHBHE_nCaloHits.push_back((it->second).nCaloHits);
        SimHitHBHE_samplingFactor.push_back(samplingFactor);
        SimHitHBHE_ieta.push_back(theHcalDetId.ieta());
        SimHitHBHE_iphi.push_back(theHcalDetId.iphi());
        SimHitHBHE_depth.push_back(theHcalDetId.depth());
        SimHitHBHE_eta.push_back(eta);
        SimHitHBHE_phi.push_back(phi);
        SimHitHBHE_detId.push_back(newid);
        SimHitHBHE_subdetId.push_back(newid.subdetId());
        SimHitHBHE_energy_25ns.push_back((it->second).energy_25ns);
        SimHitHBHE_energyEM_25ns.push_back((it->second).energyEM_25ns);
        SimHitHBHE_energyHAD_25ns.push_back((it->second).energyHAD_25ns);
        nSimHitHBHE++;
      }
    }
  }
  //==========================================
  //
  // PCaloHit: EB
  //
  //==========================================
  std::map<unsigned int, SimHitInfo> map_simHitInfoEB;

  if (saveSimCaloHitEB_ || saveSimHitEB_){
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

      if(saveSimCaloHitEB_){
        G4SimHitPCaloHitEB_energy.push_back(caloHit.energy());
        G4SimHitPCaloHitEB_energyEM.push_back(caloHit.energyEM());
        G4SimHitPCaloHitEB_energyHAD.push_back(caloHit.energyHad());
        G4SimHitPCaloHitEB_ieta.push_back(theECalEBDetId.ieta());
        G4SimHitPCaloHitEB_iphi.push_back(theECalEBDetId.iphi());
        G4SimHitPCaloHitEB_depth.push_back(caloHit.depth());
        G4SimHitPCaloHitEB_eta.push_back(eta);
        G4SimHitPCaloHitEB_phi.push_back(phi);
        G4SimHitPCaloHitEB_time.push_back(caloHit.time());
        G4SimHitPCaloHitEB_detId.push_back(detIdRaw);
        G4SimHitPCaloHitEB_subdetId.push_back(detIdRaw.subdetId());
        nG4SimHitPCaloHitEB++;
      }
    }

    if (saveSimHitEB_){
      nSimHitEB=0;
      for (std::map<unsigned int, SimHitInfo>::const_iterator it = map_simHitInfoEB.cbegin(); it != map_simHitInfoEB.cend(); ++it) {
        DetId detIdRaw(it->first);
        EBDetId theECalEBDetId(detIdRaw);
        auto cellGeometry = geometry->getSubdetectorGeometry(theECalEBDetId)->getGeometry(theECalEBDetId);
        double eta = cellGeometry->getPosition().eta();
        double phi = cellGeometry->getPosition().phi();
        SimHitEB_energy.push_back((it->second).energy);
        SimHitEB_energyEM.push_back((it->second).energyEM);
        SimHitEB_energyHAD.push_back((it->second).energyHAD);
        SimHitEB_nCaloHits.push_back((it->second).nCaloHits);
        SimHitEB_ieta.push_back(theECalEBDetId.ieta());
        SimHitEB_iphi.push_back(theECalEBDetId.iphi());
        SimHitEB_eta.push_back(eta);
        SimHitEB_phi.push_back(phi);
        SimHitEB_detId.push_back(detIdRaw);
        SimHitEB_subdetId.push_back(detIdRaw.subdetId());
        nSimHitEB++;
      }
    }
  }

  //==========================================
  //
  // PCaloHit: EE
  //
  //==========================================
  std::map<unsigned int, SimHitInfo> map_simHitInfoEE;

  if (saveSimCaloHitEE_ || saveSimHitEE_){
    edm::Handle<edm::PCaloHitContainer> g4SimHitsPCaloHitsEEHandle;
    iEvent.getByToken(g4SimHitsPCaloHitsEEToken_, g4SimHitsPCaloHitsEEHandle);
    auto g4SimHitsPCaloHitsEE = g4SimHitsPCaloHitsEEHandle.product();
    size_t n_g4SimHitsPCaloHitsEE = g4SimHitsPCaloHitsEE->size();

    nG4SimHitPCaloHitEE=0;
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

      if(saveSimCaloHitEE_){
        G4SimHitPCaloHitEE_energy.push_back(caloHit.energy());
        G4SimHitPCaloHitEE_energyEM.push_back(caloHit.energyEM());
        G4SimHitPCaloHitEE_energyHAD.push_back(caloHit.energyHad());
        G4SimHitPCaloHitEE_ix.push_back(theECalEEDetId.ix());
        G4SimHitPCaloHitEE_iy.push_back(theECalEEDetId.iy());
        G4SimHitPCaloHitEE_depth.push_back(caloHit.depth());
        G4SimHitPCaloHitEE_eta.push_back(eta);
        G4SimHitPCaloHitEE_phi.push_back(phi);
        G4SimHitPCaloHitEE_time.push_back(caloHit.time());
        G4SimHitPCaloHitEE_detId.push_back(detIdRaw);
        G4SimHitPCaloHitEE_subdetId.push_back(detIdRaw.subdetId());
        nG4SimHitPCaloHitEE++;
      }
    }

    if (saveSimHitEE_){
      nSimHitEE=0;
      for (std::map<unsigned int, SimHitInfo>::const_iterator it = map_simHitInfoEE.cbegin(); it != map_simHitInfoEE.cend(); ++it) {
        DetId detIdRaw(it->first);
        EEDetId theECalEEDetId(detIdRaw);
        auto cellGeometry = geometry->getSubdetectorGeometry(theECalEEDetId)->getGeometry(theECalEEDetId);
        double eta = cellGeometry->getPosition().eta();
        double phi = cellGeometry->getPosition().phi();
        SimHitEE_energy.push_back((it->second).energy);
        SimHitEE_energyEM.push_back((it->second).energyEM);
        SimHitEE_energyHAD.push_back((it->second).energyHAD);
        SimHitEE_nCaloHits.push_back((it->second).nCaloHits);
        SimHitEE_ix.push_back(theECalEEDetId.ix());
        SimHitEE_iy.push_back(theECalEEDetId.iy());
        SimHitEE_eta.push_back(eta);
        SimHitEE_phi.push_back(phi);
        SimHitEE_detId.push_back(detIdRaw);
        SimHitEE_subdetId.push_back(detIdRaw.subdetId());
        nSimHitEE++;
      }
    }
  }

  //===========================================================================
  //
  // TO DO: Setup to save PFTracks/PFRecTracks/PFBlockElementTrack
  //
  //===========================================================================
  // if (savePFTracks_){
  //  const auto& matched_pftrack = orig.trackRefPF();
  //   if (matched_pftrack.isNonnull()) {
  //     const auto& atECAL = matched_pftrack->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);
  //     const auto& atHCAL = matched_pftrack->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance);
  //     if (atECAL.isValid()) {
  //       eta_ecal = atECAL.positionREP().eta();
  //       phi_ecal = atECAL.positionREP().phi();
  //     }
  //     if (atHCAL.isValid()) {
  //       eta_hcal = atHCAL.positionREP().eta();
  //       phi_hcal = atHCAL.positionREP().phi();
  //     }
  //   }
  //   const auto& ref = ((const reco::PFBlockElementTrack*)&orig)->trackRef();
  //   pt = ref->pt();
  //   eta = ref->eta();
  //   phi = ref->phi();
  //   p = ref->p();
  //   px = ref->px();
  //   py = ref->py();
  //   pz = ref->pz();
  //   charge = ref->charge();
  //   nrechits = ref->recHitsSize();
  //   reco::MuonRef muonRef = orig.muonRef();
  //   if (muonRef.isNonnull()) {
  //     reco::TrackRef standAloneMu = muonRef->standAloneMuon();
  //     if (standAloneMu.isNonnull()) {
  //       trk_muon_dt_hits = standAloneMu->hitPattern().numberOfValidMuonDTHits();
  //       trk_muon_csc_hits = standAloneMu->hitPattern().numberOfValidMuonCSCHits();
  //     }
  //     trk_muon_type = muonRef->type();
  //   }
  // }

  //==========================================
  //
  // PFClustersECAL
  //
  //==========================================
  if (savePFClustersECAL_){
    Handle<reco::PFClusterCollection> pfClustersECALHandle;
    iEvent.getByToken(pfClustersECALToken_, pfClustersECALHandle);
    for( size_t icluster=0; icluster < pfClustersECALHandle->size(); ++icluster ){
      reco::PFClusterRef clusterRef( pfClustersECALHandle, icluster );
      PFClusterECAL_pt.push_back(clusterRef->pt());
      PFClusterECAL_energy.push_back(clusterRef->energy());
      PFClusterECAL_correctedEnergy.push_back(clusterRef->correctedEnergy());
      PFClusterECAL_eta.push_back(clusterRef->eta());
      PFClusterECAL_phi.push_back(clusterRef->phi());
      PFClusterECAL_layer.push_back(clusterRef->layer());
      PFClusterECAL_seedhit_detId.push_back(clusterRef->seed());

      const std::vector<std::pair<DetId, float> >& clusterHits = clusterRef->hitsAndFractions();
      PFClusterECAL_nhits.push_back(clusterHits.size());

      std::vector<unsigned int> hits_detId(clusterHits.size());
      std::vector<float> hits_fraction(clusterHits.size());
      std::vector<int> hits_indexEB(clusterHits.size());
      std::vector<int> hits_indexEE(clusterHits.size());

      for (size_t ihit=0; ihit < clusterHits.size(); ++ihit){
        hits_detId.push_back(clusterHits[ihit].first);
        hits_fraction.push_back(clusterHits[ihit].second);

        int idxEB = -1;
        auto itEB = std::find(PFRecHitEB_detId.begin(), PFRecHitEB_detId.end(), clusterHits[ihit].first);
        if (itEB != PFRecHitEB_detId.end()) {
          idxEB = static_cast<int>(std::distance(PFRecHitEB_detId.begin(), itEB));
        }
        hits_indexEB.push_back(idxEB);

        int idxEE = -1;
        if (idxEB == -1){
          auto itEE = std::find(PFRecHitEE_detId.begin(), PFRecHitEE_detId.end(), clusterHits[ihit].first);
          if (itEE != PFRecHitEE_detId.end()) {
            idxEE = static_cast<int>(std::distance(PFRecHitEE_detId.begin(), itEE));
          }
        }
        hits_indexEE.push_back(idxEE);
      }
      PFClusterECAL_hits_detId.push_back(hits_detId);
      PFClusterECAL_hits_fraction.push_back(hits_fraction);
      PFClusterECAL_hits_PFRecHitEB_Idx.push_back(hits_indexEB);
      PFClusterECAL_hits_PFRecHitEE_Idx.push_back(hits_indexEE);
      PFClusterECAL_key.push_back(clusterRef.key());
      nPFClusterECAL++;
    }
    //
    // Save indices mapping PFCluster to PFCandidates
    //
    for (size_t i=0; i < PFClusterECAL_key.size(); i++){
      for (size_t ii=0; ii < PFCand_PFElementClusterECAL_keys.size(); ii++){
      auto it = std::find(PFCand_PFElementClusterECAL_keys[ii].begin(), PFCand_PFElementClusterECAL_keys[ii].end(), PFClusterECAL_key[i]);
        if (it != PFCand_PFElementClusterECAL_keys[ii].end()){
          PFCand_PFClusterECAL_Idx[ii].push_back(i);
        }
      }
      for (size_t ii=0; ii < PFCand_PFElementClusterECALInBlock_keys.size(); ii++){
      auto it = std::find(PFCand_PFElementClusterECALInBlock_keys[ii].begin(), PFCand_PFElementClusterECALInBlock_keys[ii].end(), PFClusterECAL_key[i]);
        if (it != PFCand_PFElementClusterECALInBlock_keys[ii].end()){
          PFCand_PFClusterECALInBlock_Idx[ii].push_back(i);
        }
      }
    }
  }

  //==========================================
  //
  // PFClustersPS
  //
  //==========================================
  if (savePFClustersPS_){
    Handle<reco::PFClusterCollection> pfClustersPSHandle;
    iEvent.getByToken(pfClustersPSToken_, pfClustersPSHandle);
    for( size_t icluster=0; icluster < pfClustersPSHandle->size(); ++icluster ){
      reco::PFClusterRef clusterRef( pfClustersPSHandle, icluster );
      PFClusterPS_pt.push_back(clusterRef->pt());
      PFClusterPS_energy.push_back(clusterRef->energy());
      PFClusterPS_correctedEnergy.push_back(clusterRef->correctedEnergy());
      PFClusterPS_eta.push_back(clusterRef->eta());
      PFClusterPS_phi.push_back(clusterRef->phi());
      PFClusterPS_layer.push_back(clusterRef->layer());
      PFClusterPS_seedhit_detId.push_back(clusterRef->seed());

      const std::vector<std::pair<DetId, float> >& clusterHits = clusterRef->hitsAndFractions();
      PFClusterPS_nhits.push_back(clusterHits.size());

      std::vector<unsigned int> hits_detId(clusterHits.size());
      std::vector<float> hits_fraction(clusterHits.size());
      for (size_t ihit=0; ihit < clusterHits.size(); ++ihit){
        hits_detId.push_back(clusterHits[ihit].first);
        hits_fraction.push_back(clusterHits[ihit].second);
      }
      PFClusterPS_hits_detId.push_back(hits_detId);
      PFClusterPS_hits_fraction.push_back(hits_fraction);
      PFClusterPS_key.push_back(clusterRef.key());
      nPFClusterPS++;
    }

    //
    // Save indices mapping PFCluster to PFCandidates
    //
    for (size_t i=0; i < PFClusterPS_key.size(); i++){
      for (size_t ii=0; ii < PFCand_PFElementClusterPS_keys.size(); ii++){
      auto it = std::find(PFCand_PFElementClusterPS_keys[ii].begin(), PFCand_PFElementClusterPS_keys[ii].end(), PFClusterPS_key[i]);
        if (it != PFCand_PFElementClusterPS_keys[ii].end()){
          PFCand_PFClusterPS_Idx[ii].push_back(i);
        }
      }
      for (size_t ii=0; ii < PFCand_PFElementClusterPSInBlock_keys.size(); ii++){
      auto it = std::find(PFCand_PFElementClusterPSInBlock_keys[ii].begin(), PFCand_PFElementClusterPSInBlock_keys[ii].end(), PFClusterPS_key[i]);
        if (it != PFCand_PFElementClusterPSInBlock_keys[ii].end()){
          PFCand_PFClusterPSInBlock_Idx[ii].push_back(i);
        }
      }
    }
  }

  //==========================================
  //
  // PFClustersHCAL
  //
  //==========================================
  if (savePFClustersHCAL_){
    Handle<reco::PFClusterCollection> pfClustersHCALHandle;
    iEvent.getByToken(pfClustersHCALToken_, pfClustersHCALHandle);
    for( size_t icluster=0; icluster < pfClustersHCALHandle->size(); ++icluster ){
      reco::PFClusterRef clusterRef( pfClustersHCALHandle, icluster );
      PFClusterHCAL_pt.push_back(clusterRef->pt());
      PFClusterHCAL_energy.push_back(clusterRef->energy());
      PFClusterHCAL_correctedEnergy.push_back(clusterRef->correctedEnergy());
      PFClusterHCAL_eta.push_back(clusterRef->eta());
      PFClusterHCAL_phi.push_back(clusterRef->phi());
      PFClusterHCAL_layer.push_back(clusterRef->layer());
      PFClusterHCAL_seedhit_detId.push_back(clusterRef->seed());

      const std::vector<std::pair<DetId, float> >& clusterHits = clusterRef->hitsAndFractions();
      PFClusterHCAL_nhits.push_back(clusterHits.size());

      std::vector<unsigned int> hits_detId(clusterHits.size());
      std::vector<float> hits_fraction(clusterHits.size());
      std::vector<int> hits_index(clusterHits.size());
      for (size_t ihit=0; ihit < clusterHits.size(); ++ihit){
        hits_detId.push_back(clusterHits[ihit].first);
        hits_fraction.push_back(clusterHits[ihit].second);

        int idx =-1;
        auto it = std::find(PFRecHitHBHE_detId.begin(), PFRecHitHBHE_detId.end(), clusterHits[ihit].first);
        if (it != PFRecHitHBHE_detId.end()) {
          idx = static_cast<int>(std::distance(PFRecHitHBHE_detId.begin(), it));
        }
        hits_index.push_back(idx);
      }
      PFClusterHCAL_hits_detId.push_back(hits_detId);
      PFClusterHCAL_hits_fraction.push_back(hits_fraction);
      PFClusterHCAL_hits_PFRecHitHBHE_Idx.push_back(hits_index);
      PFClusterHCAL_key.push_back(clusterRef.key());
      nPFClusterHCAL++;
    }

    //
    // Save indices mapping PFCluster to PFCandidates
    //
    for (size_t i=0; i < PFClusterHCAL_key.size(); i++){
      for (size_t ii=0; ii < PFCand_PFElementClusterHCAL_keys.size(); ii++){
      auto it = std::find(PFCand_PFElementClusterHCAL_keys[ii].begin(), PFCand_PFElementClusterHCAL_keys[ii].end(), PFClusterHCAL_key[i]);
        if (it != PFCand_PFElementClusterHCAL_keys[ii].end()){
          PFCand_PFClusterHCAL_Idx[ii].push_back(i);
        }
      }
      for (size_t ii=0; ii < PFCand_PFElementClusterHCALInBlock_keys.size(); ii++){
      auto it = std::find(PFCand_PFElementClusterHCALInBlock_keys[ii].begin(), PFCand_PFElementClusterHCALInBlock_keys[ii].end(), PFClusterHCAL_key[i]);
        if (it != PFCand_PFElementClusterHCALInBlock_keys[ii].end()){
          PFCand_PFClusterHCALInBlock_Idx[ii].push_back(i);
        }
      }
    }
  }

  tree->Fill();
}

DEFINE_FWK_MODULE(ParticleFlowAnalysisNtuplizer);
