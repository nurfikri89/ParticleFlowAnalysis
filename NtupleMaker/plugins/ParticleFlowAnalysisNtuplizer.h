#ifndef ParticleFlowAnalysis_PFNtupleMaker_ParticleFlowAnalysisNtuplizer
#define ParticleFlowAnalysis_PFNtupleMaker_ParticleFlowAnalysisNtuplizer

#include <memory>
#include <string>
#include <iostream>
#include <numeric>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitDefs.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"

#include "SimCalorimetry/CaloSimAlgos/interface/CaloVSimParameterMap.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalSimParameterMap.h"

#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "CondFormats/DataRecord/interface/HcalPFCutsRcd.h"
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/HcalObjects/interface/HcalPFCuts.h"
#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TH1.h>
#include <TH2.h>
#include <math.h>

class ParticleFlowAnalysisNtuplizer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
 public:

  explicit ParticleFlowAnalysisNtuplizer(const edm::ParameterSet&);

  ~ParticleFlowAnalysisNtuplizer();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override;

 private:

  edm::EDGetTokenT<reco::GenParticleCollection> tokengenParticles_;
  edm::EDGetTokenT<math::XYZPointF> tokengenVertexXYZ_;

  edm::EDGetTokenT<reco::PFCandidateCollection> tokenPFCandidates_;

  edm::EDGetTokenT<edm::ValueMap<bool>> chargedHadronIsolationToken_;

  edm::EDGetTokenT<std::vector<reco::PFRecHit>>    pfRecHitsECALToken_;
  edm::EDGetTokenT<std::vector<reco::PFRecHit>>    pfRecHitsPSToken_;
  edm::EDGetTokenT<std::vector<reco::PFRecHit>>    pfRecHitsHBHEToken_;
  edm::EDGetTokenT<std::vector<reco::PFRecHit>>    pfRecHitsHFToken_;

  edm::EDGetTokenT<reco::PFClusterCollection>   pfClustersECALToken_;
  edm::EDGetTokenT<reco::PFClusterCollection>   pfClustersPSToken_;
  edm::EDGetTokenT<reco::PFClusterCollection>   pfClustersHCALToken_;

  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsEBToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsEEToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsPSToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsHCALToken_;
  // edm::EDGetTokenT<PCaloHitContainer>    g4SimPCaloHitsHFToken_;

  /// Min number of pixel hits for charged hadrons
  int nPixMin_;

  /// Min number of track hits for charged hadrons
  std::vector<int> nHitMin_;
  std::vector<double> etaMinForHitMin_;
  std::vector<double> etaMaxForHitMin_;

  TTree* tree;
  size_t orun,oevt,olumiBlock,otime;

  bool saveSimCaloHitHBHE_;
  bool saveSimCaloHitEB_;
  bool saveSimCaloHitEE_;
  bool saveSimHitHBHE_;
  bool saveSimHitEB_;
  bool saveSimHitEE_;

  bool savePFClustersECAL_;
  bool savePFClustersPS_;
  bool savePFClustersHCAL_;

  double dRMaxGenPartToPFCandChgHad_;
  double dRMaxGenPartToPFCandNeuHad_;
  double dRMaxGenPartToPFCandPhoton_;

  edm::RunNumber_t run;
  edm::EventNumber_t evt;
  edm::LuminosityBlockNumber_t lumiBlock;
  edm::Timestamp time;

  edm::ESGetToken<HepPDT::ParticleDataTable, PDTRecord> pdtToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopologyToken_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> hcalTopologyToken_;
  edm::ESGetToken<HcalDDDRecConstants, HcalRecNumberingRecord> hcalDDDrecToken_;
  edm::ESGetToken<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd> ecalPFRecHitThresholdsToken_;
  edm::ESGetToken<HcalPFCuts, HcalPFCutsRcd> hcalPFCutsToken_;
  edm::ESGetToken<HcalDbService, HcalDbRecord> hcalDBToken_;
  edm::ESGetToken<HcalRespCorrs, HcalRespCorrsRcd> hcalRespCorrsToken_;

  float genPart_pt;
  float genPart_eta;
  float genPart_phi;
  float genPart_mass;
  float genPart_energy;
  float genPart_p;
  int   genPart_pdgId;
  int   genPart_charge;

  bool  genPart_extrapolated_okECAL;
  bool  genPart_extrapolated_okHCAL;

  int   genPart_extrapolated_EB_ieta;
  int   genPart_extrapolated_EB_iphi;
  float genPart_extrapolated_EB_eta;
  float genPart_extrapolated_EB_phi;
  int   genPart_extrapolated_PFRecHitEB_Idx;
  unsigned int genPart_extrapolated_PFRecHitEB_detId;

  int   genPart_extrapolated_EE_ix;
  int   genPart_extrapolated_EE_iy;
  float genPart_extrapolated_EE_eta;
  float genPart_extrapolated_EE_phi;
  int   genPart_extrapolated_PFRecHitEE_Idx;
  unsigned int genPart_extrapolated_PFRecHitEE_detId;

  int   genPart_extrapolated_EHCAL_ieta;
  int   genPart_extrapolated_EHCAL_iphi;
  int   genPart_extrapolated_EHCAL_depth;
  float genPart_extrapolated_EHCAL_eta;
  float genPart_extrapolated_EHCAL_phi;
  int   genPart_extrapolated_EHCAL_subdetId;
  int   genPart_extrapolated_EHCAL_PFRecHitHBHE_Idx;
  int   genPart_extrapolated_EHCAL_PFRecHitHF_Idx;
  unsigned int genPart_extrapolated_EHCAL_PFRecHitHBHE_detId;
  unsigned int genPart_extrapolated_EHCAL_PFRecHitHF_detId;

  int   genPart_extrapolated_HCAL_ieta;
  int   genPart_extrapolated_HCAL_iphi;
  int   genPart_extrapolated_HCAL_depth;
  float genPart_extrapolated_HCAL_eta;
  float genPart_extrapolated_HCAL_phi;
  int   genPart_extrapolated_HCAL_subdetId;
  int   genPart_extrapolated_HCAL_PFRecHitHBHE_Idx;
  int   genPart_extrapolated_HCAL_PFRecHitHF_Idx;
  unsigned int genPart_extrapolated_HCAL_PFRecHitHBHE_detId;
  unsigned int genPart_extrapolated_HCAL_PFRecHitHF_detId;

  float genPart_extrapolated_pointECAL_x;
  float genPart_extrapolated_pointECAL_y;
  float genPart_extrapolated_pointECAL_z;
  float genPart_extrapolated_pointECAL_eta;
  float genPart_extrapolated_pointECAL_phi;

  float genPart_extrapolated_pointHCAL_x;
  float genPart_extrapolated_pointHCAL_y;
  float genPart_extrapolated_pointHCAL_z;
  float genPart_extrapolated_pointHCAL_eta;
  float genPart_extrapolated_pointHCAL_phi;

  float genVertex_x;
  float genVertex_y;
  float genVertex_z;

  int   nPFCand;
  std::vector<float> PFCand_pt;
  std::vector<float> PFCand_eta;
  std::vector<float> PFCand_phi;
  std::vector<float> PFCand_mass;
  std::vector<float> PFCand_energy;
  std::vector<int>   PFCand_pdgId;
  std::vector<int>   PFCand_isChgHadIso;
  std::vector<unsigned> PFCand_birthId;//Special
  std::vector<float> PFCand_dRToGenPart;

  std::vector<int>   PFCand_has_trk;
  std::vector<float> PFCand_trk_pt;
  std::vector<float> PFCand_trk_p;
  std::vector<float> PFCand_trk_eta;
  std::vector<float> PFCand_trk_phi;
  std::vector<int>   PFCand_trk_algo;

  std::vector<int>   PFCand_trk_charge;
  std::vector<float> PFCand_trk_ptError;
  std::vector<float> PFCand_trk_dxy;
  std::vector<float> PFCand_trk_dz;
  std::vector<float> PFCand_trk_dxyError;
  std::vector<float> PFCand_trk_dzError;
  std::vector<float> PFCand_trk_chi2;
  std::vector<float> PFCand_trk_ndof;
  std::vector<float> PFCand_trk_vx;
  std::vector<float> PFCand_trk_vy;
  std::vector<float> PFCand_trk_vz;

  std::vector<int> PFCand_trk_numberOfValidPixelHits;
  std::vector<int> PFCand_trk_numberOfValidTrackerHits;
  std::vector<int> PFCand_trk_pixelHitOK;
  std::vector<int> PFCand_trk_trackerHitOK;

  std::vector<float> PFCand_trk_positionAtECALEntrance_x;
  std::vector<float> PFCand_trk_positionAtECALEntrance_y;
  std::vector<float> PFCand_trk_positionAtECALEntrance_z;
  std::vector<float> PFCand_trk_positionAtECALEntrance_eta;
  std::vector<float> PFCand_trk_positionAtECALEntrance_phi;

  std::vector<float> PFCand_trk_closestECAL_eta;
  std::vector<float> PFCand_trk_closestECAL_phi;
  std::vector<int>   PFCand_trk_closestECAL_detId;
  std::vector<float> PFCand_trk_closestECAL_EBieta;
  std::vector<float> PFCand_trk_closestECAL_EBiphi;
  std::vector<float> PFCand_trk_closestECAL_EEix;
  std::vector<float> PFCand_trk_closestECAL_EEiy;
  std::vector<float> PFCand_trk_closestHCAL_eta;
  std::vector<float> PFCand_trk_closestHCAL_phi;
  std::vector<int>   PFCand_trk_closestHCAL_detId;
  std::vector<float> PFCand_trk_closestHCAL_ieta;
  std::vector<float> PFCand_trk_closestHCAL_iphi;

  std::vector<float> PFCand_ecalEnergy;
  std::vector<float> PFCand_hcalEnergy;
  std::vector<float> PFCand_caloEnergy;
  std::vector<float> PFCand_rawEcalEnergy;
  std::vector<float> PFCand_rawHcalEnergy;
  std::vector<float> PFCand_rawCaloEnergy;
  std::vector<float> PFCand_hoEnergy;
  std::vector<float> PFCand_rawHoEnergy;
  std::vector<float> PFCand_pS1Energy;
  std::vector<float> PFCand_pS2Energy;

  std::vector<float> PFCand_hcalDepth1EFrac;
  std::vector<float> PFCand_hcalDepth2EFrac;
  std::vector<float> PFCand_hcalDepth3EFrac;
  std::vector<float> PFCand_hcalDepth4EFrac;
  std::vector<float> PFCand_hcalDepth5EFrac;
  std::vector<float> PFCand_hcalDepth6EFrac;
  std::vector<float> PFCand_hcalDepth7EFrac;
  std::vector<int>   PFCand_PFBlock_Idx;
  std::vector<unsigned> PFCand_nPFTrack;
  std::vector<unsigned> PFCand_nPFClusterHCAL;
  std::vector<unsigned> PFCand_nPFClusterECAL;
  std::vector<unsigned> PFCand_nPFClusterPS;
  std::vector<unsigned> PFCand_nPFClusterHF;
  std::vector<unsigned> PFCand_nPFClusterSC;
  std::vector<unsigned> PFCand_nPFClusterGSF;
  std::vector<unsigned> PFCand_nPFClusterBREM;
  std::vector<std::vector<int>> PFCand_PFClusterHCAL_Idx;
  std::vector<std::vector<int>> PFCand_PFClusterECAL_Idx;
  std::vector<std::vector<int>> PFCand_PFClusterPS_Idx;
  std::vector<std::vector<int>> PFCand_PFClusterHF_Idx;
  std::vector<unsigned> PFCand_nPFTrackInBlock;
  std::vector<unsigned> PFCand_nPFClusterHCALInBlock;
  std::vector<unsigned> PFCand_nPFClusterECALInBlock;
  std::vector<unsigned> PFCand_nPFClusterPSInBlock;
  std::vector<unsigned> PFCand_nPFClusterHFInBlock;
  std::vector<unsigned> PFCand_nPFClusterSCInBlock;
  std::vector<unsigned> PFCand_nPFClusterGSFInBlock;
  std::vector<unsigned> PFCand_nPFClusterBREMInBlock;
  std::vector<std::vector<int>> PFCand_PFClusterHCALInBlock_Idx;
  std::vector<std::vector<int>> PFCand_PFClusterECALInBlock_Idx;
  std::vector<std::vector<int>> PFCand_PFClusterPSInBlock_Idx;
  std::vector<std::vector<int>> PFCand_PFClusterHFInBlock_Idx;

  int Idx_ClosestPFCandChgHad;
  int Idx_ClosestPFCandChgHadIso;
  int Idx_ClosestPFCandNeuHad;
  int Idx_ClosestPFCandPhoton;

  //
  //
  //
  int nPFRecHitHBHE;
  std::vector<float> PFRecHitHBHE_energy;
  std::vector<int>   PFRecHitHBHE_ieta;
  std::vector<int>   PFRecHitHBHE_iphi;
  std::vector<float> PFRecHitHBHE_eta;
  std::vector<float> PFRecHitHBHE_phi;
  std::vector<int>   PFRecHitHBHE_depth;
  std::vector<float> PFRecHitHBHE_cutThreshold;
  std::vector<float> PFRecHitHBHE_hcalRespCorr;
  std::vector<unsigned int>   PFRecHitHBHE_detId;
  //
  int nPFRecHitHF;
  std::vector<float> PFRecHitHF_energy;
  std::vector<int>   PFRecHitHF_ieta;
  std::vector<int>   PFRecHitHF_iphi;
  std::vector<float> PFRecHitHF_eta;
  std::vector<float> PFRecHitHF_phi;
  std::vector<int>   PFRecHitHF_depth;
  std::vector<float> PFRecHitHF_cutThreshold;
  std::vector<unsigned int>   PFRecHitHF_detId;
  //
  int nPFRecHitEB;
  std::vector<float>  PFRecHitEB_energy;
  std::vector<int>    PFRecHitEB_ieta;
  std::vector<int>    PFRecHitEB_iphi;
  std::vector<int>    PFRecHitEB_tower_ieta;
  std::vector<int>    PFRecHitEB_tower_iphi;
  std::vector<float>  PFRecHitEB_approxEta;
  std::vector<float>  PFRecHitEB_eta;
  std::vector<float>  PFRecHitEB_phi;
  std::vector<float>  PFRecHitEB_cutThreshold;
  std::vector<unsigned int>   PFRecHitEB_detId;
  //
  int nPFRecHitEE;
  std::vector<float>  PFRecHitEE_energy;
  std::vector<int>    PFRecHitEE_ix;
  std::vector<int>    PFRecHitEE_iy;
  std::vector<float>  PFRecHitEE_eta;
  std::vector<float>  PFRecHitEE_phi;
  std::vector<float>  PFRecHitEE_cutThreshold;
  std::vector<unsigned int>   PFRecHitEE_detId;
  //
  int nG4SimHitPCaloHitHBHE;
  std::vector<float> G4SimHitPCaloHitHBHE_energy;
  std::vector<float> G4SimHitPCaloHitHBHE_energyEM;
  std::vector<float> G4SimHitPCaloHitHBHE_energyHAD;
  std::vector<float> G4SimHitPCaloHitHBHE_samplingFactor;
  std::vector<int>   G4SimHitPCaloHitHBHE_ieta;
  std::vector<int>   G4SimHitPCaloHitHBHE_iphi;
  std::vector<int>   G4SimHitPCaloHitHBHE_depth;
  std::vector<float> G4SimHitPCaloHitHBHE_eta;
  std::vector<float> G4SimHitPCaloHitHBHE_phi;
  std::vector<float> G4SimHitPCaloHitHBHE_time;
  std::vector<float> G4SimHitPCaloHitHBHE_detId;
  std::vector<int>   G4SimHitPCaloHitHBHE_subdetId;
  //
  //
  //
  int nG4SimHitPCaloHitEB;
  std::vector<float> G4SimHitPCaloHitEB_energy;
  std::vector<float> G4SimHitPCaloHitEB_energyEM;
  std::vector<float> G4SimHitPCaloHitEB_energyHAD;
  std::vector<int>   G4SimHitPCaloHitEB_ieta;
  std::vector<int>   G4SimHitPCaloHitEB_iphi;
  std::vector<int>   G4SimHitPCaloHitEB_depth;
  std::vector<float> G4SimHitPCaloHitEB_eta;
  std::vector<float> G4SimHitPCaloHitEB_phi;
  std::vector<float> G4SimHitPCaloHitEB_time;
  std::vector<float> G4SimHitPCaloHitEB_detId;
  std::vector<int>   G4SimHitPCaloHitEB_subdetId;
  //
  //
  //
  int nG4SimHitPCaloHitEE;
  std::vector<float> G4SimHitPCaloHitEE_energy;
  std::vector<float> G4SimHitPCaloHitEE_energyEM;
  std::vector<float> G4SimHitPCaloHitEE_energyHAD;
  std::vector<int>   G4SimHitPCaloHitEE_ix;
  std::vector<int>   G4SimHitPCaloHitEE_iy;
  std::vector<int>   G4SimHitPCaloHitEE_depth;
  std::vector<float> G4SimHitPCaloHitEE_eta;
  std::vector<float> G4SimHitPCaloHitEE_phi;
  std::vector<float> G4SimHitPCaloHitEE_time;
  std::vector<float> G4SimHitPCaloHitEE_detId;
  std::vector<int>   G4SimHitPCaloHitEE_subdetId;
  //
  //
  //
  int nSimHitHBHE;
  std::vector<float> SimHitHBHE_energy;
  std::vector<float> SimHitHBHE_energyEM;
  std::vector<float> SimHitHBHE_energyHAD;
  std::vector<float> SimHitHBHE_samplingFactor;
  std::vector<int>   SimHitHBHE_nCaloHits;
  std::vector<int>   SimHitHBHE_ieta;
  std::vector<int>   SimHitHBHE_iphi;
  std::vector<int>   SimHitHBHE_depth;
  std::vector<float> SimHitHBHE_eta;
  std::vector<float> SimHitHBHE_phi;
  std::vector<float> SimHitHBHE_detId;
  std::vector<int>   SimHitHBHE_subdetId;
  std::vector<float> SimHitHBHE_energy_25ns;
  std::vector<float> SimHitHBHE_energyEM_25ns;
  std::vector<float> SimHitHBHE_energyHAD_25ns;
  std::vector<int>   SimHitHBHE_nCaloHits_25ns;

  //
  //
  //
  int nSimHitEB;
  std::vector<float> SimHitEB_energy;
  std::vector<float> SimHitEB_energyEM;
  std::vector<float> SimHitEB_energyHAD;
  std::vector<int>   SimHitEB_nCaloHits;
  std::vector<int>   SimHitEB_ieta;
  std::vector<int>   SimHitEB_iphi;
  std::vector<float> SimHitEB_eta;
  std::vector<float> SimHitEB_phi;
  std::vector<float> SimHitEB_detId;
  std::vector<int>   SimHitEB_subdetId;

  int nSimHitEE;
  std::vector<float> SimHitEE_energy;
  std::vector<float> SimHitEE_energyEM;
  std::vector<float> SimHitEE_energyHAD;
  std::vector<int>   SimHitEE_nCaloHits;
  std::vector<int>   SimHitEE_ix;
  std::vector<int>   SimHitEE_iy;
  std::vector<float> SimHitEE_eta;
  std::vector<float> SimHitEE_phi;
  std::vector<float> SimHitEE_detId;
  std::vector<int>   SimHitEE_subdetId;

  int nPFClusterECAL;
  std::vector<float> PFClusterECAL_pt;
  std::vector<float> PFClusterECAL_energy;
  std::vector<float> PFClusterECAL_correctedEnergy;
  std::vector<float> PFClusterECAL_eta;
  std::vector<float> PFClusterECAL_phi;
  std::vector<int>   PFClusterECAL_layer;
  std::vector<int>   PFClusterECAL_nhits;
  std::vector<unsigned int>   PFClusterECAL_seedhit_detId;
  std::vector<std::vector<unsigned int>> PFClusterECAL_hits_detId;
  std::vector<std::vector<float>>        PFClusterECAL_hits_fraction;
  std::vector<std::vector<int>> PFClusterECAL_hits_PFRecHitEB_Idx;
  std::vector<std::vector<int>> PFClusterECAL_hits_PFRecHitEE_Idx;
  std::vector<int>   PFClusterECAL_key;

  int nPFClusterPS;
  std::vector<float> PFClusterPS_pt;
  std::vector<float> PFClusterPS_energy;
  std::vector<float> PFClusterPS_correctedEnergy;
  std::vector<float> PFClusterPS_eta;
  std::vector<float> PFClusterPS_phi;
  std::vector<int>   PFClusterPS_layer;
  std::vector<int>   PFClusterPS_nhits;
  std::vector<unsigned int>   PFClusterPS_seedhit_detId;
  std::vector<std::vector<unsigned int>> PFClusterPS_hits_detId;
  std::vector<std::vector<float>>        PFClusterPS_hits_fraction;
  std::vector<int>   PFClusterPS_key;

  int nPFClusterHCAL;
  std::vector<float> PFClusterHCAL_pt;
  std::vector<float> PFClusterHCAL_energy;
  std::vector<float> PFClusterHCAL_correctedEnergy;
  std::vector<float> PFClusterHCAL_eta;
  std::vector<float> PFClusterHCAL_phi;
  std::vector<int>   PFClusterHCAL_layer;
  std::vector<int>   PFClusterHCAL_nhits;
  std::vector<unsigned int>   PFClusterHCAL_seedhit_detId;
  std::vector<std::vector<unsigned int>> PFClusterHCAL_hits_detId;
  std::vector<std::vector<float>>        PFClusterHCAL_hits_fraction;
  std::vector<std::vector<int>> PFClusterHCAL_hits_PFRecHitHBHE_Idx;
  std::vector<int>   PFClusterHCAL_key;

  //https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HcalSimAlgos/test/CaloSamplesAnalyzer.cc#0086
  HcalSimParameterMap* theParameterMap;

  bool   verbose_;
};

#endif