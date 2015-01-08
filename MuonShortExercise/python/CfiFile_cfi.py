import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('MuonShortExercise',
                      muonInputTag_ = cms.InputTag("slimmedMuons"),
                      genInputTag_  = cms.InputTag("packedGenParticles"),
                      vertexInputTag_ = cms.InputTag("offlineSlimmedPrimaryVertices")
)
