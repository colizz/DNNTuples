/*
 * JetInfoAK8.cc
 *
 *  Created on: May 23, 2017
 *      Author: hqu
 */

#include "DeepNTuples/Ntupler/interface/JetInfoFiller.h"

#include <algorithm>


namespace deepntuples {

void JetInfoFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector && cc) {
  minPt_ = iConfig.getUntrackedParameter<double>("jetPtMin", 150);
  maxPt_ = iConfig.getUntrackedParameter<double>("jetPtMax", -1);
  maxAbsEta_ = iConfig.getUntrackedParameter<double>("jetAbsEtaMax", 2.4);
  isQCDSample_ = iConfig.getUntrackedParameter<bool>("isQCDSample", false);
  isTTBarSample_ = iConfig.getUntrackedParameter<bool>("isTTBarSample", false);
  isTrainSample_ = iConfig.getUntrackedParameter<bool>("isTrainSample", false);
  btag_discriminators_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminators");
  bDiscriminatorsCompactSave1_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminatorsCompactSave1");
  bDiscriminatorsCompactSave2_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminatorsCompactSave2");
  bDiscriminatorsCompactSave3_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminatorsCompactSave3");
  bDiscriminatorsCompactSave4_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminatorsCompactSave4");
  bDiscriminatorsCompactSave5_ = iConfig.getParameter<std::vector<std::string>>("bDiscriminatorsCompactSave5");

  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  puToken_ = cc.consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puInfo"));
  rhoToken_ = cc.consumes<double>(iConfig.getParameter<edm::InputTag>("rhoInfo"));
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  metToken_ = cc.consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("METs"));
  addMET_ = iConfig.getUntrackedParameter<bool>("addMET", false);

  triggerBits_ = cc.consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTBits")),
  triggerObjects_ = cc.consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("HLTObjects")),
  triggerPrescales_ = cc.consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("HLTPrescales")),
  triggerList_ = iConfig.getUntrackedParameter<std::vector<std::string>>("HLTList"),
  addHLT_ = iConfig.getUntrackedParameter<bool>("addHLT", false);
}

void JetInfoFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(puToken_, puInfo);
  iEvent.getByToken(rhoToken_, rhoInfo);
  event_ = iEvent.id().event();
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
  flavorDef.setGenParticles(*genParticlesHandle);
  if (addMET_) {
    iEvent.getByToken(metToken_, mets);
  }
  if (addHLT_) {
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);

    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    std::string triggersPassed = "";
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      std::string triggerName = names.triggerName(i);
      if (strstr(triggerName.c_str(), "HLT_AK8PFHT") || strstr(triggerName.c_str(), "HLT_AK8PFJet") || strstr(triggerName.c_str(), "HLT_CaloJet")
          || strstr(triggerName.c_str(), "HLT_DiPFJetAve") || strstr(triggerName.c_str(), "HLT_PFHT") || strstr(triggerName.c_str(), "HLT_PFJet")
          || strstr(triggerName.c_str(), "HLT_QuadPFJet")
      ) {
        if (triggerBits->accept(i)) {
          triggersPassed += triggerName;
        }
      }
    }
    // std::cout << triggersPassed << std::endl;
    for (unsigned int i = 0; i < triggerList_.size(); ++i) {
      triggerPassedMap_[triggerList_.at(i)] = strstr(triggersPassed.c_str(), triggerList_.at(i).c_str());
    }
  }
}

bool JetInfoFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {
  // pv selection
  if (vertices->empty()) return false;

  // jet selection
  if (jet.pt() < minPt_) return false;
  if (maxPt_ > 0 && jet.pt() > maxPt_) return false;
  if (std::abs(jet.eta()) > maxAbsEta_) return false;

  // QCD and ttbar samples for inference: keep only 1/7 of the events
  if (!isTrainSample_ && (isQCDSample_ || isTTBarSample_)) {
    if (event_ %7 != 0) return false;
  }

  // event information
  data.fill<unsigned>("event_no", event_);
  data.fill<unsigned>("jet_no", jetidx);
  data.fill<float>("npv", vertices->size());
  data.fill<float>("rho", *rhoInfo);
  for (const auto &v : *puInfo) {
    int bx = v.getBunchCrossing();
    if (bx == 0) {
      data.fill<float>("ntrueInt", v.getTrueNumInteractions());
    }
  }

  // MET information
  if (addMET_) {
    data.fill<float>("met_pt", mets->front().pt());
    data.fill<float>("met_phi", mets->front().phi());
    data.fill<float>("met_sumEt", mets->front().sumEt());
    data.fill<float>("met_significance", mets->front().significance());
  }

  // HLT information
  if (addHLT_) {
    for (unsigned int i = 0; i < triggerList_.size(); ++i) {
      data.fill<bool>(triggerList_.at(i).c_str(), triggerPassedMap_[triggerList_.at(i)]);
    }
  }

  // truth labels
  float gen_pt = jet.genJet() ? jet.genJet()->pt() : 0;
  data.fill<float>("gen_pt", gen_pt);
  data.fill<float>("Delta_gen_pt", gen_pt - jet.correctedJet("Uncorrected").pt());

  auto flavor = flavorDef.jet_flavour(jet);
  data.fill<int>("isB", flavor==JetFlavor::B);
  data.fill<int>("isBB", flavor==JetFlavor::BB);
  data.fill<int>("isLeptonicB", flavor==JetFlavor::LeptonicB);
  data.fill<int>("isLeptonicB_C", flavor==JetFlavor::LeptonicB_C);
  data.fill<int>("isC", flavor==JetFlavor::C);
  data.fill<int>("isUD", flavor==JetFlavor::UD);
  data.fill<int>("isS", flavor==JetFlavor::S);
  data.fill<int>("isG", flavor==JetFlavor::G);
  data.fill<int>("isUndefined", flavor==JetFlavor::UNDEFINED);

  // jet variables
  data.fill<float>("jet_pt", jet.correctedJet("Uncorrected").pt());
  data.fill<float>("jet_corr_pt", jet.pt());
  data.fill<float>("jet_eta", jet.eta());
  data.fill<float>("jet_phi", jet.phi());

  // jet id
  data.fill<float>("jet_tightId", jetIdTight(jet));

  // discriminators
  if (!isTrainSample_) {
    for(const auto& disc : btag_discriminators_) {
      std::string name(disc);
      std::replace(name.begin(), name.end(), ':', '_');
      data.fill<float>(name, catchInfs(jet.bDiscriminator(disc), -99));
    }
  }

  // special jet scores
  if (!bDiscriminatorsCompactSave1_.empty()) {
    for (const auto& disc : bDiscriminatorsCompactSave1_) {
      data.fillMulti<float>("jet_custom_discs", catchInfs(jet.bDiscriminator(disc), -99));
    }
  }
  if (!bDiscriminatorsCompactSave2_.empty()) {
    for (const auto& disc : bDiscriminatorsCompactSave2_) {
      data.fillMulti<float>("jet_custom_discs_2", catchInfs(jet.bDiscriminator(disc), -99));
    }
  }
  if (!bDiscriminatorsCompactSave3_.empty()) {
    for (const auto& disc : bDiscriminatorsCompactSave3_) {
      data.fillMulti<float>("jet_custom_discs_3", catchInfs(jet.bDiscriminator(disc), -99));
    }
  }
  if (!bDiscriminatorsCompactSave4_.empty()) {
    for (const auto& disc : bDiscriminatorsCompactSave4_) {
      data.fillMulti<float>("jet_custom_discs_4", catchInfs(jet.bDiscriminator(disc), -99));
    }
  }
  if (!bDiscriminatorsCompactSave5_.empty()) {
    for (const auto& disc : bDiscriminatorsCompactSave5_) {
      data.fillMulti<float>("jet_custom_discs_5", catchInfs(jet.bDiscriminator(disc), -99));
    }
  }

  return true;
}

void JetInfoFiller::book() {
  // event information
  data.add<float>("npv", 0);
  data.add<float>("rho", 0);
  data.add<float>("ntrueInt", 0);
  data.add<unsigned>("event_no", 0);
  data.add<unsigned>("jet_no", 0);

  // MET information
  if (addMET_) {
    data.add<float>("met_pt", 0);
    data.add<float>("met_phi", 0);
    data.add<float>("met_sumEt", 0);
    data.add<float>("met_significance", 0);
  }

  // HLT information
  if (addHLT_) {
    for (unsigned int i = 0; i < triggerList_.size(); ++i) {
      data.add<bool>(triggerList_.at(i).c_str(), false);
    }
  }

  // truth labels
  data.add<float>("gen_pt", 0);
  data.add<float>("Delta_gen_pt", 0);

  data.add<int>("isB", 0);
  data.add<int>("isBB", 0);
  data.add<int>("isLeptonicB", 0);
  data.add<int>("isLeptonicB_C", 0);
  data.add<int>("isC", 0);
  data.add<int>("isUD", 0);
  data.add<int>("isS", 0);
  data.add<int>("isG", 0);
  data.add<int>("isUndefined", 0);

  // jet variables
  data.add<float>("jet_pt", 0);
  data.add<float>("jet_corr_pt", 0);
  data.add<float>("jet_eta", 0);
  data.add<float>("jet_phi", 0);

  // jet id
  data.add<float>("jet_tightId", 0);

  // jet discriminators
  if (!isTrainSample_) {
    for(auto name : btag_discriminators_) {
      std::replace(name.begin(), name.end(), ':', '_');
      data.add<float>(name, 0);
    }
  }

  if (!bDiscriminatorsCompactSave1_.empty()) {
    data.addMulti<float>("jet_custom_discs");
  }
  if (!bDiscriminatorsCompactSave2_.empty()) {
    data.addMulti<float>("jet_custom_discs_2");
  }
  if (!bDiscriminatorsCompactSave3_.empty()) {
    data.addMulti<float>("jet_custom_discs_3");
  }
  if (!bDiscriminatorsCompactSave4_.empty()) {
    data.addMulti<float>("jet_custom_discs_4");
  }
  if (!bDiscriminatorsCompactSave5_.empty()) {
    data.addMulti<float>("jet_custom_discs_5");
  }
}


} /* namespace deepntuples */

