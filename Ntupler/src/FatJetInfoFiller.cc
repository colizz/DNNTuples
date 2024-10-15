/*
 * FatJetInfoFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "DeepNTuples/Ntupler/interface/FatJetInfoFiller.h"
#include <string>
#include <algorithm>

namespace deepntuples {

void FatJetInfoFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  genParticlesToken_ = cc.consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  fjTagInfoName = iConfig.getParameter<std::string>("fjTagInfoName");
  useReclusteredJets_ = iConfig.getParameter<bool>("useReclusteredJets");
  isQCDSample_ = iConfig.getUntrackedParameter<bool>("isQCDSample", false);
  isTTBarSample_ = iConfig.getUntrackedParameter<bool>("isTTBarSample", false);
  isHVV2DVarMassSample_ = iConfig.getUntrackedParameter<bool>("isHVV2DVarMassSample", false);
  sample_use_pythia_ = iConfig.getParameter<bool>("isPythia");
  sample_use_herwig_ = iConfig.getParameter<bool>("isHerwig");
  sample_use_madgraph_ = iConfig.getParameter<bool>("isMadGraph");
  isTrainSample_ = iConfig.getUntrackedParameter<bool>("isTrainSample", false);
  isMDTagger_ = iConfig.getUntrackedParameter<bool>("isMDTagger", true);
  fillSeparateLabels_ = iConfig.getUntrackedParameter<bool>("fillSeparateLabels", false);
  adhocFixMode_ = iConfig.getUntrackedParameter<int>("adhocFixMode", 0);
  fjName = iConfig.getParameter<std::string>("jetType") + std::to_string(int(10*jetR_));
}

void FatJetInfoFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
}

void FatJetInfoFiller::book() {
  // truth labels
  data.add<int>("fj_isTop", 0);
  data.add<int>("fj_isW", 0);
  data.add<int>("fj_isZ", 0);
  data.add<int>("fj_isH2p", 0);
  data.add<int>("fj_isHWW", 0);
  data.add<int>("fj_isHZZ", 0);
  data.add<int>("fj_isQCD", 0);

  data.add<int>("fj_label", 0);
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

  if (fillSeparateLabels_) {
    for (auto& label: labels_) {
      data.add<int>("fj_label_" + label, 0);
    }
  }

  data.add<int>("sample_isQCD", 0);
  data.add<int>("sample_isTTBar", 0);
  data.add<int>("sample_isHVV2DVarMass", 0);
  data.add<int>("sample_use_pythia", 0);
  data.add<int>("sample_use_herwig", 0);
  data.add<int>("sample_use_madgraph", 0);

  data.add<int>("sample_useReclusteredJets", useReclusteredJets_);

//  // legacy labels
//  data.add<int>("fj_labelLegacy", 0);

  // gen-matched particle (top/W/etc.)
  data.add<float>("fj_gen_pt", 0);
  data.add<float>("fj_gen_eta", 0);
  data.add<float>("fj_gen_phi", 0);
  data.add<float>("fj_gen_mass", 0);
  data.add<float>("fj_gen_pid", 0);
  data.add<float>("fj_gen_deltaR", 999);
  data.add<float>("fj_gendau1_pt", 0);
  data.add<float>("fj_gendau1_eta", 0);
  data.add<float>("fj_gendau1_phi", 0);
  data.add<float>("fj_gendau1_mass", 0);
  data.add<float>("fj_gendau1_pid", 0);
  data.add<float>("fj_gendau1_deltaR", 999);
  data.add<float>("fj_gendau2_pt", 0);
  data.add<float>("fj_gendau2_eta", 0);
  data.add<float>("fj_gendau2_phi", 0);
  data.add<float>("fj_gendau2_mass", 0);
  data.add<float>("fj_gendau2_pid", 0);
  data.add<float>("fj_gendau2_deltaR", 999);

  // sum of all hard particles inside a jet 
  data.add<float>("fj_genparts_pt", 0);
  data.add<float>("fj_genparts_eta", 0);
  data.add<float>("fj_genparts_phi", 0);
  data.add<float>("fj_genparts_mass", 0);
  data.add<float>("fj_genpart1_pid", 0);
  data.add<float>("fj_genpart2_pid", 0);
  data.add<float>("fj_genpart3_pid", 0);
  data.add<float>("fj_genpart4_pid", 0);

  // --- jet energy/mass regression ---
  data.add<float>("fj_genjet_pt", 0);
  data.add<float>("fj_genOverReco_pt", 1); // default to 1 if not gen-matched
  data.add<float>("fj_genOverReco_pt_null", 0); // default to 0 if not gen-matched
  data.add<float>("fj_genjet_mass", 0);
  data.add<float>("fj_genOverReco_mass", 1); // default to 1 if not gen-matched
  data.add<float>("fj_genOverReco_mass_null", 0); // default to 0 if not gen-matched
  data.add<float>("fj_genjet_sdmass", 0);
  data.add<float>("fj_genjet_sdmass_sqrt", 0);
  data.add<float>("fj_genjet_targetmass", 0);
  data.add<float>("fj_genOverReco_sdmass", 1); // default to 1 if not gen-matched
  data.add<float>("fj_genOverReco_sdmass_null", 0); // default to 0 if not gen-matched
  data.add<float>("fj_genjet_nomu_pt", 0);
  data.add<float>("fj_genOverReco_nomu_pt_null", 0); // default to 0 if not gen-matched
  data.add<float>("fj_genjet_nomu_mass", 0);
  data.add<float>("fj_genOverReco_nomu_mass_null", 0); // default to 0 if not gen-matched
  data.add<float>("fj_genjet_nomu_sdmass", 0);
  data.add<float>("fj_genOverReco_nomu_sdmass_null", 0); // default to 0 if not gen-matched
  // ----------------------------------

  // fatjet kinematics
  data.add<float>("fj_pt", 0);
  data.add<float>("fj_eta", 0);
  data.add<float>("fj_phi", 0);
  data.add<float>("fj_mass", 0);
  data.add<float>("fj_energy", 0);

  // substructure
  data.add<float>("fj_tau1", 0);
  data.add<float>("fj_tau2", 0);
  data.add<float>("fj_tau3", 0);
  data.add<float>("fj_tau21", 0);
  data.add<float>("fj_tau32", 0);

  // soft drop
  data.add<float>("fj_sdmass", 0);
  data.add<float>("fj_sdmass_fromsubjets", 0);
  data.add<float>("fj_rho", 0);
  data.add<float>("fj_uncorrsdmass", 0);
  data.add<float>("fj_corrsdmass", 0);

  // subjets: soft drop gives up to 2 subjets
  data.add<float>("fj_n_sdsubjets", 0);

  data.add<float>("fj_sdsj1_pt", 0);
  data.add<float>("fj_sdsj1_eta", 0);
  data.add<float>("fj_sdsj1_phi", 0);
  data.add<float>("fj_sdsj1_mass", 0);
  data.add<float>("fj_sdsj1_csv", 0);
  data.add<float>("fj_sdsj1_ptD", 0);
  data.add<float>("fj_sdsj1_axis1", 0);
  data.add<float>("fj_sdsj1_axis2", 0);
  data.add<float>("fj_sdsj1_mult", 0);

  data.add<float>("fj_sdsj2_pt", 0);
  data.add<float>("fj_sdsj2_eta", 0);
  data.add<float>("fj_sdsj2_phi", 0);
  data.add<float>("fj_sdsj2_mass", 0);
  data.add<float>("fj_sdsj2_csv", 0);
  data.add<float>("fj_sdsj2_ptD", 0);
  data.add<float>("fj_sdsj2_axis1", 0);
  data.add<float>("fj_sdsj2_axis2", 0);
  data.add<float>("fj_sdsj2_mult", 0);

  // some variables used in a baseline tagger
  data.add<float>("fj_ptDR", 0);
  data.add<float>("fj_relptdiff", 0);
  data.add<float>("fj_sdn2", 0);


  //double-b
  data.add<float>("fj_doubleb", 0);

  //flavor info
  data.add<int>("fj_isBB", 0);
  data.add<int>("fj_isNonBB", 0);
  data.add<int>("fj_nbHadrons", 0);
  data.add<int>("fj_ncHadrons", 0);

  //double-b inputs
  // data.add<float>("fj_z_ratio", 0);
  // data.add<float>("fj_trackSipdSig_3", 0);
  // data.add<float>("fj_trackSipdSig_2", 0);
  // data.add<float>("fj_trackSipdSig_1", 0);
  // data.add<float>("fj_trackSipdSig_0", 0);
  // data.add<float>("fj_trackSipdSig_1_0", 0);
  // data.add<float>("fj_trackSipdSig_0_0", 0);
  // data.add<float>("fj_trackSipdSig_1_1", 0);
  // data.add<float>("fj_trackSipdSig_0_1", 0);
  // data.add<float>("fj_trackSip2dSigAboveCharm_0", 0);
  // data.add<float>("fj_trackSip2dSigAboveBottom_0", 0);
  // data.add<float>("fj_trackSip2dSigAboveBottom_1", 0);
  // data.add<float>("fj_tau1_trackEtaRel_0", 0);
  // data.add<float>("fj_tau1_trackEtaRel_1", 0);
  // data.add<float>("fj_tau1_trackEtaRel_2", 0);
  // data.add<float>("fj_tau0_trackEtaRel_0", 0);
  // data.add<float>("fj_tau0_trackEtaRel_1", 0);
  // data.add<float>("fj_tau0_trackEtaRel_2", 0);
  // data.add<float>("fj_tau_vertexMass_0", 0);
  // data.add<float>("fj_tau_vertexEnergyRatio_0", 0);
  // data.add<float>("fj_tau_vertexDeltaR_0", 0);
  // data.add<float>("fj_tau_flightDistance2dSig_0", 0);
  // data.add<float>("fj_tau_vertexMass_1", 0);
  // data.add<float>("fj_tau_vertexEnergyRatio_1", 0);
  // data.add<float>("fj_tau_flightDistance2dSig_1", 0);
  // data.add<float>("fj_jetNTracks", 0);
  // data.add<float>("fj_nSV", 0);

}

bool FatJetInfoFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  // ----------------------------------------------------------------
  fjmatch_.flavorLabel(&jet, *genParticlesHandle, jetR_, isMDTagger_);
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

  // veto unmatched jets in signal samples for training
  if (isTrainSample_ && !isQCDSample_ && fjlabel.rfind("QCD_", 0) == 0)
    return false;

  // ad-hoc fixes
  if (adhocFixMode_ == 1) {
    if (resparts.size() > 0 && resparts[0]->mass() > 255) {
      // std::cout << ">> ad-hoc fix: high-mass resonance " << resparts[0]->mass() << std::endl;
      return false;
    }
  }

  data.fill<int>("fj_isTop", fjlabel.rfind("Top_", 0) == 0);
  data.fill<int>("fj_isW",   fjlabel.rfind("W_", 0) == 0);
  data.fill<int>("fj_isZ",   fjlabel.rfind("Z_", 0) == 0);
  data.fill<int>("fj_isH2p", fjlabel.rfind("H_", 0) == 0 && !fjlabel.rfind("H_WW_", 0) == 0 && !fjlabel.rfind("H_WxWx_", 0) == 0 && !fjlabel.rfind("H_WxWxStar_", 0) == 0 && !fjlabel.rfind("H_ZZ_", 0) == 0 && !fjlabel.rfind("H_ZxZx_", 0) == 0 && !fjlabel.rfind("H_ZxZxStar_", 0) == 0);
  data.fill<int>("fj_isHWW", fjlabel.rfind("H_WW_", 0) == 0 || fjlabel.rfind("H_WxWx_", 0) == 0 || fjlabel.rfind("H_WxWxStar_", 0) == 0);
  data.fill<int>("fj_isHZZ", fjlabel.rfind("H_ZZ_", 0) == 0 || fjlabel.rfind("H_ZxZx_", 0) == 0 || fjlabel.rfind("H_ZxZxStar_", 0) == 0);
  data.fill<int>("fj_isQCD", fjlabel.rfind("QCD_", 0) == 0);

  // find the label index
  int label_index = -1;
  auto it = std::find(labels_.begin(), labels_.end(), fjlabel);
  if (it != labels_.end()) {
    label_index = std::distance(labels_.begin(), it);
  }else {
    throw std::logic_error("[FatJetInfoFiller::fill]: unexpected label " + fjlabel);
  }
  // if (debug_) std::cout << "   label_index: " << label_index << std::endl;

  // fill the label index
  data.fill<int>("fj_label", label_index);

  // fill separate labels
  if (fillSeparateLabels_) {
    for (int i=0; i < (int)labels_.size(); ++i) {
      data.fill<int>("label_" + labels_[i], i == label_index);
    }
  }

  data.fill<int>("sample_isQCD",  isQCDSample_);
  data.fill<int>("sample_isTTBar",  isTTBarSample_);
  data.fill<int>("sample_isHVV2DVarMass",  isHVV2DVarMassSample_);
  data.fill<int>("sample_use_pythia", sample_use_pythia_);
  data.fill<int>("sample_use_herwig", sample_use_herwig_);
  data.fill<int>("sample_use_madgraph", sample_use_madgraph_); // MG can be interfaced w/ either pythia or herwig

  // gen-matched particle (top/W/etc.)
  int resparts_size = resparts.size();
  data.fill<float>("fj_gen_pt", resparts_size > 0 ? resparts[0]->pt() : -999);
  data.fill<float>("fj_gen_eta", resparts_size > 0 ? resparts[0]->eta() : -999);
  data.fill<float>("fj_gen_phi", resparts_size > 0 ? resparts[0]->phi() : -999);
  data.fill<float>("fj_gen_mass", resparts_size > 0 ? resparts[0]->mass() : 0);
  data.fill<float>("fj_gen_pid", resparts_size > 0 ? resparts[0]->pdgId() : 0);
  data.fill<float>("fj_gen_deltaR", resparts_size > 0 ? reco::deltaR(jet, resparts[0]->p4()) : 999);
  data.fill<float>("fj_gendau1_pt", resparts_size > 1 ? resparts[1]->pt() : -999);
  data.fill<float>("fj_gendau1_eta", resparts_size > 1 ? resparts[1]->eta() : -999);
  data.fill<float>("fj_gendau1_phi", resparts_size > 1 ? resparts[1]->phi() : -999);
  data.fill<float>("fj_gendau1_mass", resparts_size > 1 ? resparts[1]->mass() : 0);
  data.fill<float>("fj_gendau1_pid", resparts_size > 1 ? resparts[1]->pdgId() : 0);
  data.fill<float>("fj_gendau1_deltaR", resparts_size > 1 ? reco::deltaR(jet, resparts[1]->p4()) : 999);
  data.fill<float>("fj_gendau2_pt", resparts_size > 2 ? resparts[2]->pt() : -999);
  data.fill<float>("fj_gendau2_eta", resparts_size > 2 ? resparts[2]->eta() : -999);
  data.fill<float>("fj_gendau2_phi", resparts_size > 2 ? resparts[2]->phi() : -999);
  data.fill<float>("fj_gendau2_mass", resparts_size > 2 ? resparts[2]->mass() : 0);
  data.fill<float>("fj_gendau2_pid", resparts_size > 2 ? resparts[2]->pdgId() : 0);
  data.fill<float>("fj_gendau2_deltaR", resparts_size > 2 ? reco::deltaR(jet, resparts[2]->p4()) : 999);

  // sum of all hard particles inside a jet
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > sumP4(0,0,0,0);
  int parts_size = parts.size();
  for (auto& p: parts) sumP4 += p->p4();
  data.fill<float>("fj_genparts_pt", sumP4.pt());
  data.fill<float>("fj_genparts_eta", sumP4.eta());
  data.fill<float>("fj_genparts_phi", sumP4.phi());
  data.fill<float>("fj_genparts_mass", sumP4.mass());
  data.fill<float>("fj_genpart1_pid", parts_size > 0 ? parts[0]->pdgId() : 0);
  data.fill<float>("fj_genpart2_pid", parts_size > 1 ? parts[1]->pdgId() : 0);
  data.fill<float>("fj_genpart3_pid", parts_size > 2 ? parts[2]->pdgId() : 0);
  data.fill<float>("fj_genpart4_pid", parts_size > 3 ? parts[3]->pdgId() : 0);

  if (debug_) {
    std::cout << "   gen resonance mass " << (resparts_size > 0 ? resparts[0]->mass() : 0) << std::endl;
    std::cout << "   gen particle mass  " << sumP4.mass() << std::endl;
  }

  // ----------------------------------

  // fatjet kinematics
  data.fill<float>("fj_pt", jet.pt());
  data.fill<float>("fj_eta", jet.eta());
  data.fill<float>("fj_phi", jet.phi());
  data.fill<float>("fj_mass", jet.mass());
  data.fill<float>("fj_energy", jet.energy());

  // substructure
  float tau1 = jet.userFloat("Njettiness" + fjName +"Puppi:tau1");
  float tau2 = jet.userFloat("Njettiness" + fjName +"Puppi:tau2");
  float tau3 = jet.userFloat("Njettiness" + fjName +"Puppi:tau3");
  data.fill<float>("fj_tau1", tau1);
  data.fill<float>("fj_tau2", tau2);
  data.fill<float>("fj_tau3", tau3);
  data.fill<float>("fj_tau21", tau1 > 0 ? tau2/tau1 : 1.01);
  data.fill<float>("fj_tau32", tau2 > 0 ? tau3/tau2 : 1.01);

  // soft drop
  std::string name(fjName);
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  auto msd_uncorr = jet.userFloat(name +"PFJetsPuppiSoftDropMass");
  data.fill<float>("fj_sdmass", msd_uncorr);
  data.fill<float>("fj_sdmass_fromsubjets", jet.groomedMass());
  data.fill<float>("fj_rho", 2 * std::log(std::max(jet.groomedMass(), 0.01) / jet_helper.jet().pt())); // use corrected pt

  // subjets: soft drop gives up to 2 subjets
  const auto& subjets = jet_helper.getSubJets();

  data.fill<float>("fj_n_sdsubjets", subjets.size());
  auto sdcorr = jet_helper.getCorrectedPuppiSoftDropMass(subjets);
  data.fill<float>("fj_uncorrsdmass", sdcorr.first);
  data.fill<float>("fj_corrsdmass", sdcorr.second);

  if (subjets.size() > 0){
    const auto &sj1 = subjets.at(0);
    data.fill<float>("fj_sdsj1_pt", sj1->pt());
    data.fill<float>("fj_sdsj1_eta", sj1->eta());
    data.fill<float>("fj_sdsj1_phi", sj1->phi());
    data.fill<float>("fj_sdsj1_mass", sj1->mass());
    data.fill<float>("fj_sdsj1_csv", sj1->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    if (subjets.size() > 1){
      const auto &sj2 = subjets.at(1);
      data.fill<float>("fj_sdsj2_pt", sj2->pt());
      data.fill<float>("fj_sdsj2_eta", sj2->eta());
      data.fill<float>("fj_sdsj2_phi", sj2->phi());
      data.fill<float>("fj_sdsj2_mass", sj2->mass());
      data.fill<float>("fj_sdsj2_csv", sj2->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

      // some variables used in a baseline tagger
      float deltaR = reco::deltaR(*sj1, *sj2);
      float var_sd_0 = sj2->pt()/(sj1->pt()+sj2->pt());
      data.fill<float>("fj_ptDR", jet.pt() * deltaR);
      data.fill<float>("fj_relptdiff", std::abs(sj1->pt()-sj2->pt()) / jet.pt());
      data.fill<float>("fj_sdn2", var_sd_0/std::pow(deltaR,-2));
    }
  }

  // ----------------------------------------------------------------

  // --- jet energy/mass regression ---
  const auto *genjet = jet_helper.genjetWithNu();
  if (genjet){
    // jet here points to the uncorrected jet
    data.fill<float>("fj_genjet_pt", genjet->pt());
    data.fill<float>("fj_genOverReco_pt", catchInfs(genjet->pt() / jet.pt(), 1));
    data.fill<float>("fj_genOverReco_pt_null", catchInfs(genjet->pt() / jet.pt(), 0));
    data.fill<float>("fj_genjet_mass", genjet->mass());
    data.fill<float>("fj_genOverReco_mass", catchInfs(genjet->mass() / jet.mass(), 1));
    data.fill<float>("fj_genOverReco_mass_null", catchInfs(genjet->mass() / jet.mass(), 0));
    if (debug_) std::cout << "   gen jet mass        " << genjet->mass() << std::endl;
  }
  const auto *sdgenjet = jet_helper.genjetWithNuSoftDrop();
  if (sdgenjet){
    // jet here points to the uncorrected jet
    auto pos = [](double x){ return x<0 ? 0 : x; };
    data.fill<float>("fj_genjet_sdmass", pos(sdgenjet->mass()));
    data.fill<float>("fj_genjet_sdmass_sqrt", std::sqrt(pos(sdgenjet->mass())));
    data.fill<float>("fj_genjet_targetmass", (!fjlabel.rfind("QCD_", 0) == 0 && !resparts.empty()) ? resparts[0]->mass() : pos(sdgenjet->mass()) );
    data.fill<float>("fj_genOverReco_sdmass", catchInfs(pos(sdgenjet->mass()) / pos(msd_uncorr), 1));
    data.fill<float>("fj_genOverReco_sdmass_null", catchInfs(pos(sdgenjet->mass()) / pos(msd_uncorr), 0));
    if (debug_) std::cout << "   gen jet sdmass      " << pos(sdgenjet->mass()) << std::endl;
  }
  const auto *genjetnomu = jet_helper.genjetNoNu();
  if (genjetnomu){
    // jet here points to the uncorrected jet
    data.fill<float>("fj_genjet_nomu_pt", genjetnomu->pt());
    data.fill<float>("fj_genOverReco_nomu_pt_null", catchInfs(genjetnomu->pt() / jet.pt(), 0));
    data.fill<float>("fj_genjet_nomu_mass", genjetnomu->mass());
    data.fill<float>("fj_genOverReco_nomu_mass_null", catchInfs(genjetnomu->mass() / jet.mass(), 0));
    if (debug_) std::cout << "   gen jet nomu mass   " << genjetnomu->mass() << std::endl;
  }
  const auto *sdgenjetnomu = jet_helper.genjetNoNuSoftDrop();
  if (sdgenjetnomu){
    // jet here points to the uncorrected jet
    auto pos = [](double x){ return x<0 ? 0 : x; };
    data.fill<float>("fj_genjet_nomu_sdmass", pos(sdgenjetnomu->mass()));
    data.fill<float>("fj_genOverReco_nomu_sdmass_null", catchInfs(pos(sdgenjetnomu->mass()) / pos(msd_uncorr), 0));
    if (debug_) std::cout << "   gen jet nomu sdmass " << pos(sdgenjetnomu->mass()) << std::endl;
  }


  // --------
  // double-b

  // const auto *bdsvTagInfo = jet.tagInfoBoostedDoubleSV(fjTagInfoName);
  // assert(bdsvTagInfo);
  // const auto &vars = bdsvTagInfo->taggingVariables();

  data.fill<float>("fj_doubleb", jet.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));

  //flavor info
  data.fill<int>("fj_isBB", jet.jetFlavourInfo().getbHadrons().size() >= 2);
  data.fill<int>("fj_isNonBB", jet.jetFlavourInfo().getbHadrons().size() < 2);
  data.fill<int>("fj_nbHadrons", jet.jetFlavourInfo().getbHadrons().size());
  data.fill<int>("fj_ncHadrons", jet.jetFlavourInfo().getcHadrons().size());

  //double-b inputs
  // data.fill<float>("fj_z_ratio", vars.get(reco::btau::z_ratio));
  // data.fill<float>("fj_trackSipdSig_3", vars.get(reco::btau::trackSip3dSig_3));
  // data.fill<float>("fj_trackSipdSig_2", vars.get(reco::btau::trackSip3dSig_2));
  // data.fill<float>("fj_trackSipdSig_1", vars.get(reco::btau::trackSip3dSig_1));
  // data.fill<float>("fj_trackSipdSig_0", vars.get(reco::btau::trackSip3dSig_0));
  // data.fill<float>("fj_trackSipdSig_1_0", vars.get(reco::btau::tau2_trackSip3dSig_0));
  // data.fill<float>("fj_trackSipdSig_0_0", vars.get(reco::btau::tau1_trackSip3dSig_0));
  // data.fill<float>("fj_trackSipdSig_1_1", vars.get(reco::btau::tau2_trackSip3dSig_1));
  // data.fill<float>("fj_trackSipdSig_0_1", vars.get(reco::btau::tau1_trackSip3dSig_1));
  // data.fill<float>("fj_trackSip2dSigAboveCharm_0", vars.get(reco::btau::trackSip2dSigAboveCharm));
  // data.fill<float>("fj_trackSip2dSigAboveBottom_0", vars.get(reco::btau::trackSip2dSigAboveBottom_0));
  // data.fill<float>("fj_trackSip2dSigAboveBottom_1", vars.get(reco::btau::trackSip2dSigAboveBottom_1));
  // data.fill<float>("fj_tau1_trackEtaRel_0", vars.get(reco::btau::tau2_trackEtaRel_0));
  // data.fill<float>("fj_tau1_trackEtaRel_1", vars.get(reco::btau::tau2_trackEtaRel_1));
  // data.fill<float>("fj_tau1_trackEtaRel_2", vars.get(reco::btau::tau2_trackEtaRel_2));
  // data.fill<float>("fj_tau0_trackEtaRel_0", vars.get(reco::btau::tau1_trackEtaRel_0));
  // data.fill<float>("fj_tau0_trackEtaRel_1", vars.get(reco::btau::tau1_trackEtaRel_1));
  // data.fill<float>("fj_tau0_trackEtaRel_2", vars.get(reco::btau::tau1_trackEtaRel_2));
  // data.fill<float>("fj_tau_vertexMass_0", vars.get(reco::btau::tau1_vertexMass));
  // data.fill<float>("fj_tau_vertexEnergyRatio_0", vars.get(reco::btau::tau1_vertexEnergyRatio));
  // data.fill<float>("fj_tau_vertexDeltaR_0", vars.get(reco::btau::tau1_vertexDeltaR));
  // data.fill<float>("fj_tau_flightDistance2dSig_0", vars.get(reco::btau::tau1_flightDistance2dSig));
  // data.fill<float>("fj_tau_vertexMass_1", vars.get(reco::btau::tau2_vertexMass));
  // data.fill<float>("fj_tau_vertexEnergyRatio_1", vars.get(reco::btau::tau2_vertexEnergyRatio));
  // data.fill<float>("fj_tau_flightDistance2dSig_1", vars.get(reco::btau::tau2_flightDistance2dSig));
  // data.fill<float>("fj_jetNTracks", vars.get(reco::btau::jetNTracks));
  // data.fill<float>("fj_nSV", vars.get(reco::btau::jetNSecondaryVertices));

  return true;
}

}