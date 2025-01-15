import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024
process = cms.Process('PFAna',Run3_2024)

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
# process.load('Configuration.StandardSequences.Generator_cff')
# process.load('GeneratorInterface.Core.genFilterSummary_cff')
# process.load('Configuration.StandardSequences.SimIdeal_cff')
# process.load('Configuration.StandardSequences.Digi_cff')
# process.load('Configuration.StandardSequences.SimL1Emulator_cff')
# process.load('Configuration.StandardSequences.DigiToRaw_cff')
# process.load('Configuration.StandardSequences.RawToDigi_cff')
# process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "140X_mcRun3_2024_realistic_v26")
# from FastSimulation.Event.ParticleFilter_cfi import *

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))
# process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
      'file:./SinglePionMinus_FEVTDEBUGHLT.root'
  ),
  duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

#
#
#
from ParticleFlowAnalysis.NtupleMaker.ParticleFlowAnalysisNtuplizer_cfi import AddParticleFlowAnalysisNtuplizer
process = AddParticleFlowAnalysisNtuplizer(process)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string("SinglePionMinus_PFAnalysis_NtuplesFromGENSIMRECO.root")
)
process.pAna = cms.Path(
  process.pfAnaNtuplizer
)

process.schedule = cms.Schedule(process.pAna)
# process.MessageLogger.cout.threshold = "INFO"
# process.MessageLogger.cout.enable    = True
# # process.MessageLogger.cout.default   = cms.untracked.PSet( limit = cms.untracked.int32(0) )
# # process.MessageLogger.cerr.threshold = "INFO"
# # process.MessageLogger.cerr.enable    = True
# # process.MessageLogger.destinations  = cms.untracked.vstring('cout','cerr')
# process.MessageLogger.debugModules = ['pfAnaNtuplizer']
