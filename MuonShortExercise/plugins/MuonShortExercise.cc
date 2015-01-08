// -*- C++ -*-
//
// Package:    cmsdas/MuonShortExercise
// Class:      MuonShortExercise
// 
/**\class MuonShortExercise MuonShortExercise.cc cmsdas/MuonShortExercise/plugins/MuonShortExercise.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 10 Dec 2014 02:49:55 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <iostream>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/VectorUtil.h"

#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"

typedef math::XYZTLorentzVector LorentzVector;

//
// some useful functions
//
//===============================================================================

bool isBHadron( int pdgId ) {
   bool retval = false ;
   int thirdDigit = (abs(pdgId) / 100) % 10 ;
   if ( thirdDigit == 5 ) retval = true ; // should catch all B mesons
   int fourthDigit = (abs(pdgId) / 1000) % 10 ;
   if ( fourthDigit == 5 ) retval = true ; // should catch all B baryons
   // int secondDigit = (abs(pdgId) / 10) % 10 ;
   // if ( secondDigit == 5 && thirdDigit == 5 ) retval = false ; // do not count bottomonium (b bbar mesons).
   return retval ;
} // isBHadron

//===============================================================================

bool decaysToB( const reco::GenParticle& gp ) {
   for ( size_t di=0; di<gp.numberOfDaughters(); di++ ) {
       const reco::Candidate* dau = gp.daughter( di ) ; 
      if ( isBHadron( dau->pdgId() ) ) return true ;
   } // di
   return false ;
} // decaysToB

//===============================================================================

bool isCHadron( int pdgId ) {
   bool retval = false ;
   int thirdDigit = (abs(pdgId) / 100) % 10 ;
   if ( thirdDigit == 4 ) retval = true ; // should catch all B mesons
   int fourthDigit = (abs(pdgId) / 1000) % 10 ;
   if ( fourthDigit == 4  ) retval = true ; // should catch all B baryons
   // int secondDigit = (abs(pdgId) / 10) % 10 ;
   // if ( secondDigit == 5 && thirdDigit == 5 ) retval = false ; // do not count bottomonium (b bbar mesons).
   return retval ;
} // isCHadron

//===============================================================================

bool decaysToC( const reco::GenParticle& gp ) {
   for ( size_t di=0; di<gp.numberOfDaughters(); di++ ) {
       const reco::Candidate* dau = gp.daughter( di ) ; 
      if ( isCHadron( dau->pdgId() ) ) return true ;
   } // di
   return false ;
} // decaysToC

//===============================================================================


//
// class declaration
//

class MuonShortExercise : public edm::EDAnalyzer {
public:
    explicit MuonShortExercise(const edm::ParameterSet&);
    ~MuonShortExercise();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    
    
    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
    enum MuonParentage {PROMPT, HF, LF, NOT_A_MUON, NMU_PAR_TYPES};
    static const char* enumNames[];
    
    const pat::PackedGenParticle getMatchedGenParticle(const pat::Muon&, const std::vector<pat::PackedGenParticle>&);
    const reco::GenParticle getMotherPacked(const pat::PackedGenParticle&);
    const reco::GenParticle getMother(const reco::GenParticle&);
    MuonParentage getParentType(const reco::GenParticle&);

    edm::InputTag muonInputTag;
    edm::InputTag genInputTag;
    edm::InputTag vertexInputTag;

    TH1F* h_pt[4];
    TH1F* h_eta[4];
    TH1F* h_nchi2[4];
    TH1F* h_nhits[4];
    TH1F* h_nstations[4];
    TH1F* h_npixels[4];
    TH1F* h_nlayers[4];
    TH1F* h_d0[4];
    TH1F* h_chiso[4];
    TH1F* h_nhiso[4];
    TH1F* h_emiso[4];
    TH1F* h_puiso[4];
    TH1F* h_coriso[4];
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
const char* MuonShortExercise::enumNames[] = {"prompt", "hf", "lf", "other"};


//
// constructors and destructor
//
MuonShortExercise::MuonShortExercise(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    muonInputTag   = iConfig.getParameter<edm::InputTag>("muonInputTag_");
    genInputTag    = iConfig.getParameter<edm::InputTag>("genInputTag_");
    vertexInputTag = iConfig.getParameter<edm::InputTag>("vertexInputTag_");

    edm::Service<TFileService> fs;
    
    for (unsigned int idx = 0; idx < MuonParentage::NMU_PAR_TYPES; idx++)
    {
        h_pt[idx]        = fs->make<TH1F>(Form("h_pt_%s", enumNames[idx]), Form("h_pt_%s", enumNames[idx]), 20, 0, 100);
        h_eta[idx]       = fs->make<TH1F>(Form("h_eta_%s", enumNames[idx]), Form("h_eta_%s", enumNames[idx]), 25, -2.5, 2.5);
        h_nchi2[idx]     = fs->make<TH1F>(Form("h_nchi2_%s", enumNames[idx]), Form("h_nchi2_%s", enumNames[idx]), 55, 0, 11);
        h_nhits[idx]     = fs->make<TH1F>(Form("h_nhits_%s", enumNames[idx]), Form("h_nhits_%s", enumNames[idx]), 40, -0.5, 39.5);
        h_nstations[idx] = fs->make<TH1F>(Form("h_nstations_%s", enumNames[idx]), Form("h_nstations_%s", enumNames[idx]), 6, -0.5, 5.5);
        h_npixels[idx]   = fs->make<TH1F>(Form("h_npixels_%s", enumNames[idx]), Form("h_npixels_%s", enumNames[idx]), 6, -0.5, 5.5);
        h_nlayers[idx]   = fs->make<TH1F>(Form("h_nlayers_%s", enumNames[idx]), Form("h_nlayers_%s", enumNames[idx]), 25, -0.5, 24.5);
        h_d0[idx]        = fs->make<TH1F>(Form("h_d0_%s", enumNames[idx]), Form("h_d0_%s", enumNames[idx]), 40, -0.10, 0.11);
        h_chiso[idx]     = fs->make<TH1F>(Form("h_chiso_%s", enumNames[idx]), Form("h_chiso_%s", enumNames[idx]), 20, 0, 10);
        h_nhiso[idx]     = fs->make<TH1F>(Form("h_nhiso_%s", enumNames[idx]), Form("h_nhiso_%s", enumNames[idx]), 20, 0, 10);
        h_emiso[idx]     = fs->make<TH1F>(Form("h_emiso_%s", enumNames[idx]), Form("h_emiso_%s", enumNames[idx]), 20, 0, 10);
        h_puiso[idx]     = fs->make<TH1F>(Form("h_puiso_%s", enumNames[idx]), Form("h_puiso_%s", enumNames[idx]), 20, 0, 10);
        h_coriso[idx]    = fs->make<TH1F>(Form("h_coriso_%s", enumNames[idx]), Form("h_coriso_%s", enumNames[idx]), 30, 0, 0.3);
    }
}


MuonShortExercise::~MuonShortExercise()
{
 
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    // for (unsigned int idx = 0; idx < MuonParentage::NMU_PAR_TYPES; idx++)
    // {
    //     delete h_pt[idx];
    //     delete h_eta[idx];
    //     delete h_nchi2[idx];
    //     delete h_nhits[idx];
    //     delete h_nstations[idx];
    //     delete h_npixels[idx];
    //     delete h_nlayers[idx];
    //     delete h_d0[idx];
    //     delete h_chiso[idx];
    //     delete h_nhiso[idx];
    //     delete h_emiso[idx];
    //     delete h_puiso[idx];
    //     delete h_coriso[idx];
    // }
}


//
// member functions
//

const pat::PackedGenParticle MuonShortExercise::getMatchedGenParticle(const pat::Muon& muon, const std::vector<pat::PackedGenParticle>& gpcol)
{
    LorentzVector mup4(muon.p4());    
    double minDR = 999.;
    std::vector<pat::PackedGenParticle>::const_iterator gen_end = gpcol.end();
    pat::PackedGenParticle matched_gen;
    for (std::vector<pat::PackedGenParticle>::const_iterator it = gpcol.begin(); it != gen_end; it++)
    {
        LorentzVector gp4(it->p4());
        double dr = ROOT::Math::VectorUtil::DeltaR(mup4, gp4);
        if (dr < 0.1 && abs(it->pdgId()) == 13) return (*it);
        if (dr < minDR)
        {
            minDR = dr;
            matched_gen = *it;
        }
    }

    if (minDR < 0.1) return matched_gen;
    return pat::PackedGenParticle();
}

const reco::GenParticle MuonShortExercise::getMotherPacked(const pat::PackedGenParticle& pgp)
{
    if (pgp.numberOfMothers() > 0)
    {        
        const reco::GenParticle* firstMother = (const reco::GenParticle*)pgp.mother(0);
        if (firstMother != 0)
        {
            if (firstMother->pdgId() != pgp.pdgId()) return (*firstMother);
            const reco::GenParticle newMother = MuonShortExercise::getMother(*firstMother);
            return newMother;
        }
        else return reco::GenParticle();
    }

    return reco::GenParticle();
}

const reco::GenParticle MuonShortExercise::getMother(const reco::GenParticle& gp)
{
    const reco::GenParticle* mom = &gp;
    while (mom->numberOfMothers() > 0)
    {
        for (unsigned int idx = 0; idx < mom->numberOfMothers(); idx++)
        {
            mom = dynamic_cast<const reco::GenParticle*>(mom->mother(idx));
            if (mom->pdgId() != gp.pdgId())
                return (*mom);
        }
    }

    return (*mom);
}

MuonShortExercise::MuonParentage MuonShortExercise::getParentType(const reco::GenParticle& gp)
{
    unsigned int pdgId = abs(gp.pdgId());
    if (pdgId == 15 || pdgId == 23 || pdgId == 24 || pdgId == 25) return MuonParentage::PROMPT;
    if (isBHadron(pdgId) || isCHadron(pdgId)) return MuonParentage::HF;    
    return MuonShortExercise::MuonParentage::LF;
}

// ------------ method called for each event  ------------
void
MuonShortExercise::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<pat::Muon> > muon_h;
    iEvent.getByLabel(muonInputTag, muon_h);

    edm::Handle<pat::PackedGenParticleCollection> genp_h;
    iEvent.getByLabel(genInputTag, genp_h);
    std::vector<pat::PackedGenParticle> gpcol = (*genp_h);

    edm::Handle<reco::VertexCollection> vtx_h;
    iEvent.getByLabel(vertexInputTag, vtx_h);
    reco::VertexCollection::const_iterator firstGoodVertex = vtx_h->end();
    for (reco::VertexCollection::const_iterator it = vtx_h->begin(); it != firstGoodVertex; it++)
    {
        if (it->isFake()) continue;
        if (it->ndof() < 4) continue;
        if (it->position().Rho() > 2.) continue;
        if (fabs(it->position().Z()) > 24.) continue;
        firstGoodVertex = it;
        break;
    }

    // require a good vertex
    if (firstGoodVertex == vtx_h->end()) return;
    
    //
    // loop over muons
    //
    if (!muon_h.isValid() || !genp_h.isValid()) return;
    edm::View<pat::Muon>::const_iterator muon_end = muon_h->end();
    for (edm::View<pat::Muon>::const_iterator it = muon_h->begin(); it != muon_end; it++)
    {
        // require that muon is a global muon
        if (!it->globalTrack().isNonnull()) continue;

        // require muon to have a silicon track and be matched to the PV (dz < 0.2)
        if (!it->innerTrack().isNonnull()) continue;
        if (fabs(it->innerTrack()->dz(firstGoodVertex->position())) > 0.2) continue;
        
        // require that muon is within the tracker volume and has pt > 20 GeV
        if (fabs(it->eta()) > 2.5) continue;
        if (it->pt() < 20) continue;

        MuonParentage parentage = MuonParentage::NOT_A_MUON;
        int momid = 0;
        const pat::PackedGenParticle matchedPackedGenParticle = getMatchedGenParticle(*it, gpcol);        
        if (abs(matchedPackedGenParticle.pdgId()) == 13)
        {
            const reco::GenParticle momgp = getMotherPacked(matchedPackedGenParticle);
            if (momgp.pdgId() != 0)
            {
                momid = momgp.pdgId();
                parentage = getParentType(momgp);                
            }
        }

        h_pt[parentage]->Fill(std::min(it->pt(),(float)99.9));
        h_eta[parentage]->Fill(std::max(std::min(it->eta(), (float)2.49), (float)-2.49));
        h_nchi2[parentage]->Fill(std::min(it->globalTrack()->normalizedChi2(), 10.99));
        h_nhits[parentage]->Fill(std::min(it->globalTrack()->hitPattern().numberOfValidMuonHits(), (int)39));
        h_nstations[parentage]->Fill(std::min(it->numberOfMatchedStations(), (int)5));
        h_npixels[parentage]->Fill(std::min(it->innerTrack()->hitPattern().numberOfValidPixelHits(), (int)5));
        h_nlayers[parentage]->Fill(std::min(it->innerTrack()->hitPattern().trackerLayersWithMeasurement(), (int)14));
        h_d0[parentage]->Fill(std::min(std::max(it->muonBestTrack()->dxy(firstGoodVertex->position()), -0.299), 0.299));

        if (it->isIsolationValid())
        {
            reco::MuonPFIsolation pfR03 = it->pfIsolationR03();
            h_chiso[parentage]->Fill(std::min(pfR03.sumChargedHadronPt, (float)9.9));
            h_nhiso[parentage]->Fill(std::min(pfR03.sumNeutralHadronEt, (float)9.9));
            h_emiso[parentage]->Fill(std::min(pfR03.sumPhotonEt, (float)9.9));
            h_puiso[parentage]->Fill(std::min(pfR03.sumPUPt, (float)9.9));

            double coriso = pfR03.sumChargedHadronPt + std::max(0., pfR03.sumNeutralHadronEt+pfR03.sumPhotonEt-0.5*pfR03.sumPUPt);
            h_coriso[parentage]->Fill(std::min(coriso/it->pt(), 0.299));
        }
        
        // with parentage of muon in hand, let's fill some histograms
        if (parentage < MuonParentage::LF && parentage > MuonParentage::LF)
        // if (parentage < MuonParentage::LF)
            printf("pt, eta, pdgId, momId, parentage: %4.2f, %4.2f, %d, %d, %d\n", it->pt(), it->eta(), it->pdgId(), momid, parentage);
    }
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonShortExercise::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonShortExercise::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  MuonShortExercise::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  MuonShortExercise::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  MuonShortExercise::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  MuonShortExercise::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonShortExercise::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonShortExercise);
