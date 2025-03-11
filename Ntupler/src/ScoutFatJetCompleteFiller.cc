/*
 * ScoutFatJetCompleteFiller.cc
 *
 *  Created on: Mar 2, 2025
 *      Author: colizz
 */

#include <TVector3.h>

#include "DeepNTuples/Ntupler/interface/ScoutFatJetCompleteFiller.h"
#include <string>
#include <algorithm>

namespace deepntuples {

void ScoutFatJetCompleteFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  isHVV2DVarMassSample_ = iConfig.getUntrackedParameter<bool>("isHVV2DVarMassSample", false);
  isMDTagger_ = iConfig.getUntrackedParameter<bool>("isMDTagger", true);

  // value maps (for aux features of scouting PF cands)
  for (auto& name: value_map_float_names_) {
    value_map_float_tokens_[name] = cc.consumes<edm::ValueMap<float>>(edm::InputTag(name));
  }
  for (auto& name: value_map_int_names_) {
    value_map_int_tokens_[name] = cc.consumes<edm::ValueMap<int>>(edm::InputTag(name));
  }
}


void ScoutFatJetCompleteFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
  for (auto& name: value_map_float_names_) {
    iEvent.getByToken(value_map_float_tokens_[name], value_map_float_handles_[name]);
  }
  for (auto& name: value_map_int_names_) {
    iEvent.getByToken(value_map_int_tokens_[name], value_map_int_handles_[name]);
  }
}

void ScoutFatJetCompleteFiller::book() {
  // assign labels individually to scouting jet:

  // truth labels
  data.add<int>("scoutfj_isTop", 0);
  data.add<int>("scoutfj_isW", 0);
  data.add<int>("scoutfj_isZ", 0);
  data.add<int>("scoutfj_isH2p", 0);
  data.add<int>("scoutfj_isHWW", 0);
  data.add<int>("scoutfj_isHZZ", 0);
  data.add<int>("scoutfj_isQCD", 0);

  data.add<int>("scoutfj_label", 0);
  if (labels_.empty()) {
    // initialize labels
    if (isMDTagger_) {
      for (auto& l: labelTop_)  labels_.push_back("Top_" + l);
      for (auto& l: labelH2p_)  labels_.push_back("H_" + l);
      for (auto& l: labelHWW_)  labels_.push_back("H_WW_" + l);
      for (auto& l: labelHWW_)  labels_.push_back("H_WxWx_" + l);
      for (auto& l: labelHWW_)  labels_.push_back("H_WxWxStar_" + l);
      for (auto& l: labelHZZ_)  labels_.push_back("H_ZZ_" + l);
      for (auto& l: labelHZZ_)  labels_.push_back("H_ZxZx_" + l);
      for (auto& l: labelHZZ_)  labels_.push_back("H_ZxZxStar_" + l);
      for (auto& l: labelQCD_)  labels_.push_back("QCD_" + l);
      for (auto& l: labelH2pExt_) labels_.push_back("Hext_" + l);
      for (auto& l: labelHHVExt_) labels_.push_back("H_HVext_" + l);
    } else {
      for (auto& l: labelTop_)  labels_.push_back("Top_" + l);
      for (auto& l: labelH2p_)  labels_.push_back("H_" + l);
      for (auto& l: labelW_)    labels_.push_back("W_" + l);
      for (auto& l: labelZ_)    labels_.push_back("Z_" + l);
      for (auto& l: labelQCD_)  labels_.push_back("QCD_" + l);
    }

    if (debug_) {
      std::cout << "All labels (isMDTagger = " << isMDTagger_ << "):" << std::endl;
      for (auto& l: labels_)  std::cout << "'" << l << "', ";
      std::cout << std::endl;
    }
  }

  // gen-matched particle (top/W/etc.)
  data.add<float>("scoutfj_gen_pt", 0);
  data.add<float>("scoutfj_gen_eta", 0);
  data.add<float>("scoutfj_gen_phi", 0);
  data.add<float>("scoutfj_gen_mass", 0);
  data.add<float>("scoutfj_gen_pid", 0);
  data.add<float>("scoutfj_gen_deltaR", 999);
  data.add<float>("scoutfj_gendau1_pt", 0);
  data.add<float>("scoutfj_gendau1_eta", 0);
  data.add<float>("scoutfj_gendau1_phi", 0);
  data.add<float>("scoutfj_gendau1_mass", 0);
  data.add<float>("scoutfj_gendau1_pid", 0);
  data.add<float>("scoutfj_gendau1_deltaR", 999);
  data.add<float>("scoutfj_gendau2_pt", 0);
  data.add<float>("scoutfj_gendau2_eta", 0);
  data.add<float>("scoutfj_gendau2_phi", 0);
  data.add<float>("scoutfj_gendau2_mass", 0);
  data.add<float>("scoutfj_gendau2_pid", 0);
  data.add<float>("scoutfj_gendau2_deltaR", 999);

  // sum of all hard particles inside a jet 
  data.add<float>("scoutfj_genparts_pt", 0);
  data.add<float>("scoutfj_genparts_eta", 0);
  data.add<float>("scoutfj_genparts_phi", 0);
  data.add<float>("scoutfj_genparts_mass", 0);
  data.add<float>("scoutfj_genpart1_pid", 0);
  data.add<float>("scoutfj_genpart2_pid", 0);
  data.add<float>("scoutfj_genpart3_pid", 0);
  data.add<float>("scoutfj_genpart4_pid", 0);

  // fatjet kinematics
  data.add<float>("scoutfj_pt", 0);
  data.add<float>("scoutfj_eta", 0);
  data.add<float>("scoutfj_phi", 0);
  data.add<float>("scoutfj_mass", 0);
  data.add<float>("scoutfj_energy", 0);

  // substructure
  data.add<float>("scoutfj_tau1", 0);
  data.add<float>("scoutfj_tau2", 0);
  data.add<float>("scoutfj_tau3", 0);
  data.add<float>("scoutfj_tau21", 0);
  data.add<float>("scoutfj_tau32", 0);

  // soft drop
  data.add<float>("scoutfj_sdmass", 0);

  // jet tagging probabilities
  data.add<float>("scoutfj_probQCD", 0);
  data.add<float>("scoutfj_probHbb", 0);
  data.add<float>("scoutfj_probHcc", 0);
  data.add<float>("scoutfj_probHqq", 0);

  // regressed mass
  data.add<float>("scoutfj_massreg", 0);

  // ----------------------------------------------------------------
  // scouting pf features
  data.add<int>("n_scoutpfcands", 0);

  data.addMulti<float>("scoutpfcand_px");
  data.addMulti<float>("scoutpfcand_py");
  data.addMulti<float>("scoutpfcand_pz");
  data.addMulti<float>("scoutpfcand_energy");
  data.addMulti<float>("scoutpfcand_pt");

  data.addMulti<float>("scoutpfcand_charge");
  data.addMulti<float>("scoutpfcand_isEl");
  data.addMulti<float>("scoutpfcand_isMu");
  data.addMulti<float>("scoutpfcand_isChargedHad");
  data.addMulti<float>("scoutpfcand_isGamma");
  data.addMulti<float>("scoutpfcand_isNeutralHad");
  data.addMulti<float>("scoutpfcand_phirel");
  data.addMulti<float>("scoutpfcand_etarel");
  data.addMulti<float>("scoutpfcand_deltaR");
  data.addMulti<float>("scoutpfcand_abseta");
  data.addMulti<float>("scoutpfcand_ptrel_log");
  data.addMulti<float>("scoutpfcand_ptrel");
  data.addMulti<float>("scoutpfcand_erel_log");
  data.addMulti<float>("scoutpfcand_erel");
  data.addMulti<float>("scoutpfcand_pt_log");
  data.addMulti<float>("scoutpfcand_e_log");          // Added

  data.addMulti<float>("scoutpfcand_normchi2");
  data.addMulti<float>("scoutpfcand_lostInnerHits");
  data.addMulti<float>("scoutpfcand_quality");
  data.addMulti<float>("scoutpfcand_dz");
  data.addMulti<float>("scoutpfcand_dzsig");
  data.addMulti<float>("scoutpfcand_dxy");
  data.addMulti<float>("scoutpfcand_dxysig");
  data.addMulti<float>("scoutpfcand_btagEtaRel");
  data.addMulti<float>("scoutpfcand_btagPtRatio");
  data.addMulti<float>("scoutpfcand_btagPParRatio");

}

bool ScoutFatJetCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  if (jet_helper.getScoutJet() == nullptr) { // this is the case when no scouting jet exists in an event, so ScoutJet is not set
    return false;
  }

  const edm::RefToBase<reco::Jet>& scoutjetRef = *jet_helper.getScoutJet();

  if (scoutjetRef.isNonnull()) {

    // ----------------------------------------------------------------
    const pat::Jet ajet(scoutjetRef);
    if (debug_) {
      std::cout << "matched: " << "jet pt = " << jet.pt()
      << "; scouting jet pt = " << ajet.pt() << std::endl;
    }

    // gen-matching
    fjmatch_.flavorLabel(&ajet, *genParticlesHandle, jetR_, isMDTagger_);
    std::string fjlabel = fjmatch_.getResult().label;
    auto& resparts = fjmatch_.getResult().resParticles;
    auto& parts = fjmatch_.getResult().particles;

    // update the label for HWW and HZZ
    if (fjlabel.rfind("H_WW", 0) == 0 && isHVV2DVarMassSample_) {
      float mass_asymm = std::abs(resparts[1]->mass() - resparts[2]->mass()) / (resparts[1]->mass() + resparts[2]->mass());
      if (mass_asymm < 0.1) {
        fjlabel.replace(fjlabel.find("H_WW"), 4, "H_WxWx");
      } else {
        fjlabel.replace(fjlabel.find("H_WW"), 4, "H_WxWxStar");
      }
      if (debug_) std::cout << "V masses: " << resparts[1]->mass() << " " << resparts[2]->mass() << " Asymm: " << mass_asymm << std::endl;
    }
    if (fjlabel.rfind("H_ZZ", 0) == 0 && isHVV2DVarMassSample_) {
      float mass_asymm = std::abs(resparts[1]->mass() - resparts[2]->mass()) / (resparts[1]->mass() + resparts[2]->mass());
      if (mass_asymm < 0.1) {
        fjlabel.replace(fjlabel.find("H_ZZ"), 4, "H_ZxZx");
      } else {
        fjlabel.replace(fjlabel.find("H_ZZ"), 4, "H_ZxZxStar");
      }
      if (debug_) std::cout << "V masses: " << resparts[1]->mass() << " " << resparts[2]->mass() << " Asymm: " << mass_asymm << std::endl;
    }
    if (debug_) {
      std::cout << ">> debug fjlabel: " << fjlabel << "  " << std::endl;
      std::cout << "   resonance parts: "; for (auto& p: resparts) {std::cout << p->pdgId() << " ";} std::cout << std::endl;
      std::cout << "   parts: "; for (auto& p: parts) {std::cout << p->pdgId() << " ";} std::cout << std::endl;
    }
  
    data.fill<int>("scoutfj_isTop", fjlabel.rfind("Top_", 0) == 0);
    data.fill<int>("scoutfj_isW",   fjlabel.rfind("W_", 0) == 0);
    data.fill<int>("scoutfj_isZ",   fjlabel.rfind("Z_", 0) == 0);
    data.fill<int>("scoutfj_isH2p", fjlabel.rfind("H_", 0) == 0 && !fjlabel.rfind("H_WW_", 0) == 0 && !fjlabel.rfind("H_WxWx_", 0) == 0 && !fjlabel.rfind("H_WxWxStar_", 0) == 0 && !fjlabel.rfind("H_ZZ_", 0) == 0 && !fjlabel.rfind("H_ZxZx_", 0) == 0 && !fjlabel.rfind("H_ZxZxStar_", 0) == 0);
    data.fill<int>("scoutfj_isHWW", fjlabel.rfind("H_WW_", 0) == 0 || fjlabel.rfind("H_WxWx_", 0) == 0 || fjlabel.rfind("H_WxWxStar_", 0) == 0);
    data.fill<int>("scoutfj_isHZZ", fjlabel.rfind("H_ZZ_", 0) == 0 || fjlabel.rfind("H_ZxZx_", 0) == 0 || fjlabel.rfind("H_ZxZxStar_", 0) == 0);
    data.fill<int>("scoutfj_isQCD", fjlabel.rfind("QCD_", 0) == 0);
  
    // find the label index
    int label_index = -1;
    auto it = std::find(labels_.begin(), labels_.end(), fjlabel);
    if (it != labels_.end()) {
      label_index = std::distance(labels_.begin(), it);
    }else {
      throw std::logic_error("[ScoutFatJetCompleteFiller::fill]: unexpected label " + fjlabel);
    }
    // if (debug_) std::cout << "   label_index: " << label_index << std::endl;
  
    // fill the label index
    data.fill<int>("scoutfj_label", label_index);
  

    // gen-matched particle (top/W/etc.)
    int resparts_size = resparts.size();
    data.fill<float>("scoutfj_gen_pt", resparts_size > 0 ? resparts[0]->pt() : -999);
    data.fill<float>("scoutfj_gen_eta", resparts_size > 0 ? resparts[0]->eta() : -999);
    data.fill<float>("scoutfj_gen_phi", resparts_size > 0 ? resparts[0]->phi() : -999);
    data.fill<float>("scoutfj_gen_mass", resparts_size > 0 ? resparts[0]->mass() : 0);
    data.fill<float>("scoutfj_gen_pid", resparts_size > 0 ? resparts[0]->pdgId() : 0);
    data.fill<float>("scoutfj_gen_deltaR", resparts_size > 0 ? reco::deltaR(jet, resparts[0]->p4()) : 999);
    data.fill<float>("scoutfj_gendau1_pt", resparts_size > 1 ? resparts[1]->pt() : -999);
    data.fill<float>("scoutfj_gendau1_eta", resparts_size > 1 ? resparts[1]->eta() : -999);
    data.fill<float>("scoutfj_gendau1_phi", resparts_size > 1 ? resparts[1]->phi() : -999);
    data.fill<float>("scoutfj_gendau1_mass", resparts_size > 1 ? resparts[1]->mass() : 0);
    data.fill<float>("scoutfj_gendau1_pid", resparts_size > 1 ? resparts[1]->pdgId() : 0);
    data.fill<float>("scoutfj_gendau1_deltaR", resparts_size > 1 ? reco::deltaR(jet, resparts[1]->p4()) : 999);
    data.fill<float>("scoutfj_gendau2_pt", resparts_size > 2 ? resparts[2]->pt() : -999);
    data.fill<float>("scoutfj_gendau2_eta", resparts_size > 2 ? resparts[2]->eta() : -999);
    data.fill<float>("scoutfj_gendau2_phi", resparts_size > 2 ? resparts[2]->phi() : -999);
    data.fill<float>("scoutfj_gendau2_mass", resparts_size > 2 ? resparts[2]->mass() : 0);
    data.fill<float>("scoutfj_gendau2_pid", resparts_size > 2 ? resparts[2]->pdgId() : 0);
    data.fill<float>("scoutfj_gendau2_deltaR", resparts_size > 2 ? reco::deltaR(jet, resparts[2]->p4()) : 999);

    // sum of all hard particles inside a jet
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > sumP4(0,0,0,0);
    for (auto& p: parts) sumP4 += p->p4();
    int parts_size = parts.size();
    data.fill<float>("scoutfj_genparts_pt", sumP4.pt());
    data.fill<float>("scoutfj_genparts_eta", sumP4.eta());
    data.fill<float>("scoutfj_genparts_phi", sumP4.phi());
    data.fill<float>("scoutfj_genparts_mass", sumP4.mass());
    data.fill<float>("scoutfj_genpart1_pid", parts_size > 0 ? parts[0]->pdgId() : 0);
    data.fill<float>("scoutfj_genpart2_pid", parts_size > 1 ? parts[1]->pdgId() : 0);
    data.fill<float>("scoutfj_genpart3_pid", parts_size > 2 ? parts[2]->pdgId() : 0);
    data.fill<float>("scoutfj_genpart4_pid", parts_size > 3 ? parts[3]->pdgId() : 0);

    if (debug_) {
      std::cout << "   gen resonance mass " << (resparts_size > 0 ? resparts[0]->mass() : 0) << std::endl;
      std::cout << "   gen particle mass  " << sumP4.mass() << std::endl;
    }

    // ----------------------------------

    // fatjet kinematics
    data.fill<float>("scoutfj_pt", ajet.pt());
    data.fill<float>("scoutfj_eta", ajet.eta());
    data.fill<float>("scoutfj_phi", ajet.phi());
    data.fill<float>("scoutfj_mass", ajet.mass());
    data.fill<float>("scoutfj_energy", ajet.energy());

    // substructure
    float tau1 = (*value_map_float_handles_["scoutingFatPFJetReclusterNjettiness:tau1"])[scoutjetRef];
    float tau2 = (*value_map_float_handles_["scoutingFatPFJetReclusterNjettiness:tau2"])[scoutjetRef];
    float tau3 = (*value_map_float_handles_["scoutingFatPFJetReclusterNjettiness:tau3"])[scoutjetRef];
    data.fill<float>("scoutfj_tau1", tau1);
    data.fill<float>("scoutfj_tau2", tau2);
    data.fill<float>("scoutfj_tau3", tau3);
    data.fill<float>("scoutfj_tau21", tau1 > 0 ? tau2/tau1 : 1.01);
    data.fill<float>("scoutfj_tau32", tau2 > 0 ? tau3/tau2 : 1.01);

    // soft drop
    auto msd_uncorr = (*value_map_float_handles_["scoutingFatPFJetReclusterSoftDropMass"])[scoutjetRef];
    data.fill<float>("scoutfj_sdmass", msd_uncorr);

    // jet tagging probs
    data.fill<float>("scoutfj_probQCD", (*value_map_float_handles_["scoutingFatPFJetReclusterParticleNetJetTags:probQCDall"])[scoutjetRef]);
    data.fill<float>("scoutfj_probHbb", (*value_map_float_handles_["scoutingFatPFJetReclusterParticleNetJetTags:probHbb"])[scoutjetRef]);
    data.fill<float>("scoutfj_probHcc", (*value_map_float_handles_["scoutingFatPFJetReclusterParticleNetJetTags:probHcc"])[scoutjetRef]);
    data.fill<float>("scoutfj_probHqq", (*value_map_float_handles_["scoutingFatPFJetReclusterParticleNetJetTags:probHqq"])[scoutjetRef]);

    // mass regression
    data.fill<float>("scoutfj_massreg", (*value_map_float_handles_["scoutingFatPFJetReclusterParticleNetMassRegressionJetTags:mass"])[scoutjetRef]);
  
    // ----------------------------------------------------------------

    // scouting pf candidates
    std::vector<reco::CandidatePtr> daughters;
    for (const auto &dau : ajet.daughterPtrVector()) {
      // filling daughters
      daughters.push_back(dau);
    }
    // sort daughters by pt
    std::sort(daughters.begin(), daughters.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });

    data.fill<int>("n_scoutpfcands", daughters.size());

    // fill daughters
    math::XYZVector jet_dir = ajet.momentum().Unit();
    TVector3 jet_direction(ajet.momentum().Unit().x(), ajet.momentum().Unit().y(), ajet.momentum().Unit().z());
    const float etasign = jet.eta() > 0 ? 1 : -1;
  
    for (const auto &cand : daughters) {
      const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));
      auto candP4 = cand->p4();

      data.fillMulti<float>("scoutpfcand_px", candP4.px());
      data.fillMulti<float>("scoutpfcand_py", candP4.py());
      data.fillMulti<float>("scoutpfcand_pz", candP4.pz());
      data.fillMulti<float>("scoutpfcand_energy", candP4.energy());
      data.fillMulti<float>("scoutpfcand_pt", candP4.pt());
      
      data.fillMulti<float>("scoutpfcand_charge",  reco_cand->charge());
      data.fillMulti<float>("scoutpfcand_isEl",    std::abs(reco_cand->pdgId()) == 11);
      data.fillMulti<float>("scoutpfcand_isMu",    std::abs(reco_cand->pdgId()) == 13);
      data.fillMulti<float>("scoutpfcand_isChargedHad", std::abs(reco_cand->pdgId()) == 211);
      data.fillMulti<float>("scoutpfcand_isGamma", std::abs(reco_cand->pdgId()) == 22);
      data.fillMulti<float>("scoutpfcand_isNeutralHad", std::abs(reco_cand->pdgId()) == 130);
      data.fillMulti<float>("scoutpfcand_phirel", reco::deltaPhi(candP4, ajet));
      data.fillMulti<float>("scoutpfcand_etarel", etasign * (candP4.eta() - ajet.eta()));   // Edited
      data.fillMulti<float>("scoutpfcand_deltaR", reco::deltaR(candP4, ajet));
      data.fillMulti<float>("scoutpfcand_abseta", std::abs(candP4.eta()));
      data.fillMulti<float>("scoutpfcand_ptrel_log", std::log(candP4.pt() / ajet.pt()));
      data.fillMulti<float>("scoutpfcand_ptrel", candP4.pt() / ajet.pt());
      data.fillMulti<float>("scoutpfcand_erel_log", std::log(candP4.energy() / ajet.energy()));
      data.fillMulti<float>("scoutpfcand_erel", candP4.energy() / ajet.energy());
      data.fillMulti<float>("scoutpfcand_pt_log", std::log(candP4.pt()));
      data.fillMulti<float>("scoutpfcand_e_log", std::log(candP4.energy()));          // Added

      if ((*value_map_float_handles_["scoutingPFCandidate:normchi2"])[cand] > 900) {
        data.fillMulti<float>("scoutpfcand_normchi2", 0);
        data.fillMulti<float>("scoutpfcand_lostInnerHits", 0);
        data.fillMulti<float>("scoutpfcand_quality", 0);
        data.fillMulti<float>("scoutpfcand_dz", 0);
        data.fillMulti<float>("scoutpfcand_dzsig", 0);
        data.fillMulti<float>("scoutpfcand_dxy", 0);
        data.fillMulti<float>("scoutpfcand_dxysig", 0);
        data.fillMulti<float>("scoutpfcand_btagEtaRel", 0);
        data.fillMulti<float>("scoutpfcand_btagPtRatio", 0);
        data.fillMulti<float>("scoutpfcand_btagPParRatio", 0);
      } else {
        data.fillMulti<float>("scoutpfcand_normchi2", (*value_map_float_handles_["scoutingPFCandidate:normchi2"])[cand]);
        data.fillMulti<float>("scoutpfcand_lostInnerHits", (*value_map_int_handles_["scoutingPFCandidate:lostInnerHits"])[cand]);
        data.fillMulti<float>("scoutpfcand_quality", (*value_map_int_handles_["scoutingPFCandidate:quality"])[cand]);
        data.fillMulti<float>("scoutpfcand_dz", (*value_map_float_handles_["scoutingPFCandidate:dz"])[cand]);
        data.fillMulti<float>("scoutpfcand_dzsig", (*value_map_float_handles_["scoutingPFCandidate:dzsig"])[cand]);
        data.fillMulti<float>("scoutpfcand_dxy", (*value_map_float_handles_["scoutingPFCandidate:dxy"])[cand]);
        data.fillMulti<float>("scoutpfcand_dxysig", (*value_map_float_handles_["scoutingPFCandidate:dxysig"])[cand]);
        float trk_px = (*value_map_float_handles_["scoutingPFCandidate:trkPt"])[cand] * std::cos((*value_map_float_handles_["scoutingPFCandidate:trkPhi"])[cand]);
        float trk_py = (*value_map_float_handles_["scoutingPFCandidate:trkPt"])[cand] * std::sin((*value_map_float_handles_["scoutingPFCandidate:trkPhi"])[cand]);
        float trk_pz = (*value_map_float_handles_["scoutingPFCandidate:trkPt"])[cand] * std::sinh((*value_map_float_handles_["scoutingPFCandidate:trkEta"])[cand]);
        math::XYZVector track_mom(trk_px, trk_py, trk_pz);
        TVector3 track_direction(trk_px, trk_py, trk_pz);
        double track_mag = sqrt(trk_px * trk_px + trk_py * trk_py + trk_pz * trk_pz);
        data.fillMulti<float>("scoutpfcand_btagEtaRel", reco::btau::etaRel(jet_dir, track_mom));
        data.fillMulti<float>("scoutpfcand_btagPtRatio", track_direction.Perp(jet_direction) / track_mag);
        data.fillMulti<float>("scoutpfcand_btagPParRatio", jet_dir.Dot(track_mom) / track_mag);
      }

    }

  } else {
    if (debug_)
      std::cout << "!!!not matched: " << "jet.pt() = " << jet.pt() << std::endl;
  }


  return true;
}

}