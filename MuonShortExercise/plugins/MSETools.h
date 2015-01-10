#ifndef MSETOOLS_H
#define MSETOOLS_H

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef math::XYZTLorentzVector LorentzVector;

namespace mse
{
    enum MuonParentage {PROMPT, HF, LF, NOT_A_MUON, NMU_PAR_TYPES};
    extern const char* enumNames[4];
    
    bool isBHadron (int pdgId);
    bool decaysToB (const reco::GenParticle& gp);
    bool isCHadron (int pdgId);
    bool decaysToC (const reco::GenParticle& gp);

    const pat::PackedGenParticle getMatchedGenParticle(const pat::Muon&, const std::vector<pat::PackedGenParticle>&);
    const reco::GenParticle getMotherPacked(const pat::PackedGenParticle&);
    const reco::GenParticle getMother(const reco::GenParticle&);

    MuonParentage getParentType(const reco::GenParticle&);

    bool isGoodVertex(const reco::Vertex&);
};

#endif
