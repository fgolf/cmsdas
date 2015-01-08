// -*- C++ -*-
//
// Package:    cmsdas/MuonFilter
// Class:      MuonFilter
// 
/**\class MuonFilter MuonFilter.cc cmsdas/MuonFilter/plugins/MuonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Mon, 05 Jan 2015 16:56:34 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

//
// class declaration
//

class MuonFilter : public edm::EDFilter {
   public:
      explicit MuonFilter(const edm::ParameterSet&);
      ~MuonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
    edm::InputTag muonInputTag;
    unsigned int nMus;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonFilter::MuonFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
    muonInputTag = iConfig.getParameter<edm::InputTag>("muonInputTag_");
    nMus         = iConfig.getParameter<int>("nMus_");
}


MuonFilter::~MuonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
MuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<pat::Muon> > mu_h;
    iEvent.getByLabel(muonInputTag, mu_h);

    unsigned int numMus = 0;
    edm::View<pat::Muon>::const_iterator mu_end = mu_h->end();
    for (edm::View<pat::Muon>::const_iterator it = mu_h->begin(); it != mu_end; it++)
    {
        if (!it->globalTrack().isNonnull()) continue;
        if (it->pt() < 20.) continue;
        if (fabs(it->eta()) > 2.5) continue;
        ++numMus;
    }
    
    return (numMus >= nMus);
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
MuonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MuonFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonFilter);
