/*
 * JetInfoAK8.h
 *
 *  Created on: May 23, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_JETINFOFILLER_H_
#define NTUPLER_INTERFACE_JETINFOFILLER_H_

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"
#include "DeepNTuples/BTagHelpers/interface/FlavorDefinition.h"

namespace deepntuples {

class JetInfoFiller: public NtupleBase {
public:
  JetInfoFiller() : JetInfoFiller("") {}
  JetInfoFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR), flavorDef(jetR),
      jetIdTight(PFJetIDSelectionFunctor::SUMMER18PUPPI, PFJetIDSelectionFunctor::TIGHT) {}
  virtual ~JetInfoFiller() {}

  // get input parameters from the cfg file
  virtual void readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) override;

  // read event content or event setup for each event
  virtual void readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) override;

protected:
  // declare the data branches (name, type, default values)
  virtual void book() override;
  // fill the branches
  virtual bool fill(const pat::Jet &jet, size_t jetidx, const JetHelper &jet_helper) override;

private:
  FlavorDefinition flavorDef;

  double minPt_ = 0;
  double maxPt_ = 0;
  double maxAbsEta_ = 0;
  bool isQCDSample_ = false;
  bool isTTBarSample_ = false;
  bool isTrainSample_ = false;
  std::vector<std::string> btag_discriminators_;
  std::vector<std::string> bDiscriminatorsCompactSave1_;
  std::vector<std::string> bDiscriminatorsCompactSave2_;
  std::vector<std::string> bDiscriminatorsCompactSave3_;

  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
  edm::EDGetTokenT<double> rhoToken_;

  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
  edm::Handle<double> rhoInfo;

  unsigned event_ = 0;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;

  PFJetIDSelectionFunctor jetIdTight;

  bool addMET_ = false;
  edm::EDGetTokenT<std::vector<pat::MET>> metToken_;
  edm::Handle<std::vector<pat::MET>> mets;

  bool addHLT_ = false;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  std::vector<std::string> triggerList_;
  std::map<std::string, bool> triggerPassedMap_;
};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_JetInfoFiller_H_ */
