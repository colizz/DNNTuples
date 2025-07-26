/*
 * DeepNtuplizer.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "Compression.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DeepNTuples/NtupleCommons/interface/TreeWriter.h"

#include "DeepNTuples/Ntupler/interface/JetInfoFiller.h"
#include "DeepNTuples/Ntupler/interface/FatJetInfoFiller.h"
#include "DeepNTuples/Ntupler/interface/SVFiller.h"
#include "DeepNTuples/Ntupler/interface/PFCompleteFiller.h"


using namespace deepntuples;

class DeepNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DeepNtuplizer(const edm::ParameterSet&);
  ~DeepNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  double jetR = -1;

  edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> candToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetWithNuMatchToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetWithNuSoftDropMatchToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetNoNuMatchToken_;
  edm::EDGetTokenT<edm::Association<reco::GenJetCollection>> genJetNoNuSoftDropMatchToken_;

  edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>> electronToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticleToken_;

  bool addLowLevel_;

  edm::Service<TFileService> fs;
  TreeWriter *treeWriter = nullptr;

  NtupleBase* addModule(NtupleBase *m){
    modules_.push_back(m);
    return m;
  }
  std::vector<NtupleBase*> modules_;
};

DeepNtuplizer::DeepNtuplizer(const edm::ParameterSet& iConfig):
    jetR(iConfig.getParameter<double>("jetR")),
    jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
    candToken_(consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pfcands"))),
    genJetWithNuMatchToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetsWithNuMatch"))),
    genJetWithNuSoftDropMatchToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetsWithNuSoftDropMatch"))),
    genJetNoNuMatchToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetsNoNuMatch"))),
    genJetNoNuSoftDropMatchToken_(consumes<edm::Association<reco::GenJetCollection>>(iConfig.getParameter<edm::InputTag>("genJetsNoNuSoftDropMatch"))),
    muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("METs"))),
    vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    genParticleToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    addLowLevel_(iConfig.getUntrackedParameter<bool>("addLowLevel", true))
{

  // register modules
  JetInfoFiller *jetinfo = new JetInfoFiller("", jetR);
  addModule(jetinfo);

  FatJetInfoFiller *fjinfo = new FatJetInfoFiller("", jetR);
  addModule(fjinfo);

  if (addLowLevel_) {
    SVFiller *sv = new SVFiller("", jetR);
    addModule(sv);

    PFCompleteFiller *parts = new PFCompleteFiller("", jetR);
    addModule(parts);
  }

  // read config and init modules
  for(auto& m: modules_)
    m->readConfig(iConfig, consumesCollector());

}

DeepNtuplizer::~DeepNtuplizer()
{
  for(auto *m : modules_)
    delete m;
}


// ------------ method called for each event  ------------
void DeepNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  for(auto *m : modules_){
    m->readEvent(iEvent, iSetup);
  }

  edm::Handle<edm::View<pat::Jet>> jets;
  iEvent.getByToken(jetToken_, jets);

  edm::Handle<edm::View<reco::Candidate>> candHandle;
  iEvent.getByToken(candToken_, candHandle);

  edm::Handle<edm::Association<reco::GenJetCollection>> genJetWithNuMatchHandle;
  iEvent.getByToken(genJetWithNuMatchToken_, genJetWithNuMatchHandle);

  edm::Handle<edm::Association<reco::GenJetCollection>> genJetWithNuSoftDropMatchHandle;
  iEvent.getByToken(genJetWithNuSoftDropMatchToken_, genJetWithNuSoftDropMatchHandle);

  edm::Handle<edm::Association<reco::GenJetCollection>> genJetNoNuMatchHandle;
  iEvent.getByToken(genJetNoNuMatchToken_, genJetNoNuMatchHandle);

  edm::Handle<edm::Association<reco::GenJetCollection>> genJetNoNuSoftDropMatchHandle;
  iEvent.getByToken(genJetNoNuSoftDropMatchToken_, genJetNoNuSoftDropMatchHandle);

  // Get muons, electrons, MET, and vertices
  edm::Handle<edm::View<pat::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<edm::View<pat::Electron>> electrons;
  iEvent.getByToken(electronToken_, electrons);

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);

  // Get gen particles
  edm::Handle<edm::View<reco::GenParticle>> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);

  // // Get primary vertex
  // const reco::Vertex& primaryVertex = vertices->at(0);

  // Build V boson from reco-level leptons + MET (**not used**)
  /*
  reco::Candidate::LorentzVector maxLepP4(0,0,0,0);
  float maxLepPt = -1;
  bool hasLep = false;

  // Select muons
  for (const auto& mu : *muons) {
    // muon selection requirement from VHcc 1L channel:
    // lep.pt > 25 and abs(lep.dxy) < 0.05 and abs(lep.dz) < 0.2 and lep.tightId and lep.pfRelIso04_all < 0.06
    bool isTightId = mu.isLooseMuon(primaryVertex);
    double pfRelIso04_all = (mu.pfIsolationR04().sumChargedHadronPt + 
                             std::max(mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - mu.pfIsolationR04().sumPUPt/2, 0.0f)) / mu.pt();
    if (mu.pt() > 25 && std::abs(mu.dB(pat::Muon::PV2D)) < 0.05 && std::abs(mu.dB(pat::Muon::PVDZ)) < 0.2
        && isTightId && pfRelIso04_all < 0.06) {
      if (mu.pt() > maxLepPt) {
        maxLepPt = mu.pt();
        maxLepP4 = mu.p4();
        hasLep = true;
      }
    }
  }

  // Select electrons
  for (const auto& ele : *electrons) {
    // electron selection requirement from VHcc 1L channel:
    // lep.pt > 30 and lep.mvaFall17V2Iso_WP80
    if (ele.pt() > 30 && ele.electronID("mvaEleID-Fall17-iso-V2-wp80")) {
      if (ele.pt() > maxLepPt) {
        maxLepPt = ele.pt();
        maxLepP4 = ele.p4();
        hasLep = true;
      }
    }
  }

  // Build V boson four-momentum
  reco::Candidate::LorentzVector vBosonP4(0,0,0,0);
  if (hasLep && mets->size() > 0) {
    vBosonP4 = maxLepP4 + mets->at(0).p4();
  }
  */

  // Build V boson from GEN level particles
  std::vector<const reco::GenParticle*> leptonicWbosons;
  
  for (const auto& genParticle : *genParticles) {
    // Look for W bosons (PDG ID = ±24)
    if (std::abs(genParticle.pdgId()) == 24 && genParticle.status() == 62) {
      bool isLeptonicDecay = false;
      reco::Candidate::LorentzVector lepP4(0,0,0,0);
      
      // Check daughters for leptonic decay
      for (unsigned int i = 0; i < genParticle.numberOfDaughters(); ++i) {
        const reco::GenParticle* daughter = dynamic_cast<const reco::GenParticle*>(genParticle.daughter(i));
        if (!daughter) continue;
        
        int daughterPdgId = std::abs(daughter->pdgId());
        // Check for electron (11), muon (13), or tau (15) and their neutrinos (12, 14, 16)
        if ((daughterPdgId >= 11 && daughterPdgId <= 16) && daughterPdgId % 2 == 1) {
          // Found a charged lepton
          isLeptonicDecay = true;
          break;
        }
      }
      
      if (isLeptonicDecay) {
        leptonicWbosons.push_back(&genParticle);
      }
    }
  }

  // Warning if multiple leptonic W bosons found
  if (leptonicWbosons.size() > 1) {
    std::cout << "WARNING: Found " << leptonicWbosons.size() << " leptonic W bosons in event!" << std::endl;
    for (size_t i = 0; i < leptonicWbosons.size(); ++i) {
      std::cout << "  W boson " << i << ": pt = " << leptonicWbosons[i]->pt() 
                << ", eta = " << leptonicWbosons[i]->eta() << std::endl;
    }
  }
  if (leptonicWbosons.empty()) {
      std::cout << "WARNING: No leptonic W boson found in this event" << std::endl;
      return;
  }

  for (unsigned idx=0; idx<jets->size(); ++idx){
    bool write_ = true;

    // Check deltaR cut between jet and V boson
    double dR = reco::deltaR(jets->at(idx).p4(), leptonicWbosons[0]->p4());
    if (dR <= 2.5) continue; // Skip jets too close to V boson

    const auto& jet = jets->at(idx); // need to keep the JEC for puppi sdmass corr
    JetHelper jet_helper(&jet, candHandle);
    jet_helper.setGenjetWithNu((*genJetWithNuMatchHandle)[jets->refAt(idx)]);
    jet_helper.setGenjetWithNuSoftDrop((*genJetWithNuSoftDropMatchHandle)[jets->refAt(idx)]);
    jet_helper.setGenjetNoNu((*genJetNoNuMatchHandle)[jets->refAt(idx)]);
    jet_helper.setGenjetNoNuSoftDrop((*genJetNoNuSoftDropMatchHandle)[jets->refAt(idx)]);

    for (auto *m : modules_){
      if (!m->fillBranches(jet.correctedJet("Uncorrected"), idx, jet_helper)){
        write_ = false;
        break;
      }
    }

    if (write_) {
      treeWriter->fill();
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void DeepNtuplizer::beginJob() {
  if( !fs ){
    throw edm::Exception( edm::errors::Configuration,
        "TFile Service is not registered in cfg file" );
  }
  fs->file().SetCompressionAlgorithm(ROOT::kLZ4);
  fs->file().SetCompressionLevel(4);
  treeWriter = new TreeWriter(fs->make<TTree>("tree" ,"tree"));

  for(auto *m : modules_)
    m->initBranches(treeWriter);

}

// ------------ method called once each job just after ending the event loop  ------------
void DeepNtuplizer::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DeepNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DeepNtuplizer);
