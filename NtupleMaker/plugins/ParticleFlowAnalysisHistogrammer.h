#ifndef ParticleFlowAnalysis_PFNtupleMaker_ParticleFlowAnalysisHistogrammer
#define ParticleFlowAnalysis_PFNtupleMaker_ParticleFlowAnalysisHistogrammer

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
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

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
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <math.h>

class ParticleFlowAnalysisHistogrammer : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
 public:

  explicit ParticleFlowAnalysisHistogrammer(const edm::ParameterSet&);

  ~ParticleFlowAnalysisHistogrammer();

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override;

 private:

  edm::EDGetTokenT<reco::GenParticleCollection> tokengenParticles_;
  edm::EDGetTokenT<math::XYZPointF> tokengenVertexXYZ_;

  edm::EDGetTokenT<edm::SimTrackContainer> input_simtrack_token_;
  edm::EDGetTokenT<edm::SimVertexContainer> input_simvertex_token_;

  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsEBToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsEEToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsPSToken_;
  edm::EDGetTokenT<edm::PCaloHitContainer>    g4SimHitsPCaloHitsHCALToken_;
  // edm::EDGetTokenT<PCaloHitContainer>    g4SimPCaloHitsHFToken_;

  template <typename T, typename... Args>
  T *book(const Args &...args) const;
  void bookTH1F(TString hName, TString hTitle, int nbinsx, float xmin, float xmax);
  void bookTH2F(TString hName, TString hTitle, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax);
  void bookTH3F(TString hName, TString hTitle, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, int nbinsz, float zmin, float zmax);
  void bookProfile(TString pName, TString pTitle, int nbinsx, float xmin, float xmax);
  void bookProfile2D(TString pName, TString pTitle, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax);

  float genPart_pt;
  float genPart_eta;
  float genPart_phi;
  float genPart_mass;
  float genPart_energy;
  float genPart_p;
  int   genPart_pdgId;
  int   genPart_charge;

  float genVertex_x;
  float genVertex_y;
  float genVertex_z;

  TTree* tree;
  size_t orun,oevt,olumiBlock,otime;

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

  //https://cmssdt.cern.ch/lxr/source/SimCalorimetry/HcalSimAlgos/test/CaloSamplesAnalyzer.cc#0086
  HcalSimParameterMap* theParameterMap;

  bool   verbose_;

  edm::Service<TFileService> fs_;

  std::map<TString, TH1F*>       m_Histos1D;
  std::map<TString, TH2F*>       m_Histos2D;
  std::map<TString, TH3F*>       m_Histos3D;
  std::map<TString, TProfile*>   m_Profiles;
  std::map<TString, TProfile2D*> m_Profiles2D;

  int eventCounter_;
};

#endif