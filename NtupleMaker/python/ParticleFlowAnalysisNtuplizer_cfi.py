import FWCore.ParameterSet.Config as cms


def AddParticleFlowAnalysisNtuplizer(process):
    process.load('SimCalorimetry.HcalSimProducers.hcalSimParameters_cfi')
    process.pfAnaNtuplizer = cms.EDAnalyzer( "ParticleFlowAnalysisNtuplizer",
        process.hcalSimParameters,
        genParticles = cms.InputTag("genParticles"),
        genVertex = cms.InputTag("genParticles","xyz0"),
        PFCandidates = cms.InputTag("particleFlow"),
        chargedHadronIsolation = cms.InputTag("chargedHadronPFTrackIsolation"),
        #
        dRMaxGenPartToPFCandChgHad = cms.untracked.double(0.4),
        dRMaxGenPartToPFCandNeuHad = cms.untracked.double(0.4),
        dRMaxGenPartToPFCandPhoton = cms.untracked.double(0.4),
        #
        PFClustersECAL = cms.InputTag("particleFlowClusterECAL"),
        PFClustersPS   = cms.InputTag("particleFlowClusterPS"),
        PFClustersHCAL = cms.InputTag("particleFlowClusterHCAL"),
        savePFClustersECAL = cms.untracked.bool(True),
        savePFClustersPS = cms.untracked.bool(True),
        savePFClustersHCAL = cms.untracked.bool(True),
        #
        PFRecHitsECAL = cms.InputTag("particleFlowRecHitECAL"),
        PFRecHitsPS  = cms.InputTag("particleFlowRecHitPS"),
        PFRecHitsHBHE = cms.InputTag("particleFlowRecHitHBHE"),
        PFRecHitsHF = cms.InputTag("particleFlowRecHitHF"),
        #
        g4SimHitsPCaloHitsEB   = cms.InputTag("g4SimHits","EcalHitsEB"),
        g4SimHitsPCaloHitsEE   = cms.InputTag("g4SimHits","EcalHitsEE"),
        g4SimHitsPCaloHitsPS   = cms.InputTag("g4SimHits","EcalHitsES"),
        g4SimHitsPCaloHitsHCAL = cms.InputTag("g4SimHits","HcalHits"),
        saveSimCaloHitHCAL = cms.untracked.bool(False),
        saveSimCaloHitEB   = cms.untracked.bool(False),
        saveSimCaloHitEE   = cms.untracked.bool(False),
        saveSimHitHBHE     = cms.untracked.bool(True),
        saveSimHitEB       = cms.untracked.bool(True),
        saveSimHitEE       = cms.untracked.bool(True),
        #
        #
        nPixMin = cms.int32(2),                     # Nb of pixel hits
        #
        nHitMin = cms.vint32(14,17,20,17,10),       # Nb of track hits
        etaMinForHitMin = cms.vdouble(0.0,1.4,1.6,2.0,2.6), # in these eta ranges
        etaMaxForHitMin = cms.vdouble(1.4,1.6,2.0,2.4,2.6), # in these eta ranges
        # hb = process.hcalSimParameters.hb.clone(),
        # he = process.hcalSimParameters.he.clone(),
        # hcalSimParameters = process.hcalSimParameters
    )
    return process