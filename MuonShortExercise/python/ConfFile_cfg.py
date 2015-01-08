import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/home/users/fgolf/cmsdas/CMSSW_7_3_0_pre3/src/cmsdas/MuonFilter/test/merged_qcd.root'
#        'file:/home/users/fgolf/cmsdas/CMSSW_7_3_0_pre3/src/cmsdas/MuonFilter/test/merged_dymm.root'
        "/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root"
    )
)

process.load("cmsdas.MuonShortExercise.CfiFile_cfi")
# process.demo = cms.EDAnalyzer('MuonShortExercise'
# )

process.TFileService = cms.Service("TFileService", fileName = cms.string("histos_ttbar.root"))

process.p = cms.Path(process.demo)
