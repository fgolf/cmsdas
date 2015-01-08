import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Phys14DR/QCD_Pt-30to50_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/00000/0A08C123-7274-E411-853B-001E6739753A.root'
    )
)

process.demo = cms.EDFilter('MuonFilter',
                            muonInputTag_ = cms.InputTag("slimmedMuons"),
                            nMus_ = cms.int32(1)
)


process.p = cms.Path(process.demo)

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        fileName     = cms.untracked.string('ntuple.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    ),
    dropMetaData = cms.untracked.string("NONE")
)
process.outpath      = cms.EndPath(process.out)
