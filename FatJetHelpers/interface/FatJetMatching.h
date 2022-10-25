/*
 * FatJetMatching.hh
 *
 *  Created on: Feb 1, 2017
 *      Author: hqu
 */

#ifndef FATJETHELPERS_INTERFACE_FATJETMATCHING_H_
#define FATJETHELPERS_INTERFACE_FATJETMATCHING_H_

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include <unordered_set>
#include <utility>

namespace deepntuples {

namespace ParticleID{
enum PdgId { p_unknown, p_d, p_u, p_s, p_c, p_b, p_t, p_bprime, p_tprime,
  p_eminus = 11, p_nu_e, p_muminus, p_nu_mu, p_tauminus, p_nu_tau,
  p_tauprimeminus, p_nu_tauprime, p_g = 21, p_gamma, p_Z0,
  p_Wplus, p_h0, p_Zprime0 = 32, p_Zpprime0, p_Wprimeplus, p_H0,
  p_A0, p_Hplus, p_G = 39, p_R0 = 41, p_H30 = 45, p_A20 = 46,
  p_LQ, p_cluster = 91, p_string,
  p_pi0 = 111, p_rho0 = 113, p_klong = 130, p_piplus = 211, p_rhoplus = 213, p_eta = 221, p_omega = 223,
  p_kshort = 310, p_k0, p_kstar0 = 313, p_kplus = 321, p_kstarplus = 323, p_phi = 333,
  p_dplus = 411, p_d0 = 421, p_dsplus = 431, p_b0 =511, p_bplus = 521,
  p_bs0 = 531, p_bcplus = 541,
  p_neutron = 2112, p_proton = 2212,
  p_sigmaminus = 3112, p_lambda0 = 3122,
  p_sigma0 = 3212, p_sigmaplus = 3222, p_ximinus = 3312, p_xi0 = 3322, p_omegaminus = 3334,
  p_sigmac0 = 4112, p_lambdacplus = 4122, p_xic0 = 4132,
  p_sigmacplus = 4212, p_sigmacpp = 4222, p_xicplus = 4232, p_omegac0 = 4332,
  p_sigmabminus = 5112, p_lambdab0 = 5122, p_xibminus = 5132, p_sigmab0 = 5212, p_sigmabplus = 5222,
  p_xib0 = 5232, p_omegabminus = 5332,
  p_hc = 5000003,
};
}

class FatJetMatching {
public:
  enum FatJetFlavor {
    Default = 0,
    Top = 1,
    W = 2,
    Z = 3,
    H = 4,
  };

  enum FatJetLabel {
    Invalid=0,
    Top_all=100, Top_bcq, Top_bqq, Top_bc, Top_bq, Top_bev, Top_bmv, Top_bleptauev, Top_bleptaumv, Top_bhadtauv,
    W_all=200, W_cq_b, W_qq_b, W_ev_b, W_mv_b, W_leptauev_b, W_leptaumv_b, W_hadtauv_b,
      W_cq_c, W_qq_c, W_ev_c, W_mv_c, W_leptauev_c, W_leptaumv_c, W_hadtauv_c,
      W_cq, W_qq, W_ev, W_mv, W_leptauev, W_leptaumv, W_hadtauv,
    Z_all=300, Z_bb, Z_cc, Z_qq,
    H_all=400, H_bb, H_bc, H_cc, H_ss, H_qq, H_qqqq, H_leptauehadtau, H_leptaumhadtau, H_hadtauhadtau,
    H_wwall=500, H_ww4q_2c, H_ww4q_1c, H_ww4q_0c, H_ww3q_2c, H_ww3q_1c, H_ww3q_0c, H_ww2qsame, H_ww2qsep, 
      H_wwevqq_1c, H_wwevqq_0c, H_wwmvqq_1c, H_wwmvqq_0c, H_wwleptauevqq_1c, H_wwleptauevqq_0c, H_wwleptaumvqq_1c, H_wwleptaumvqq_0c, H_wwhadtauvqq_1c, H_wwhadtauvqq_0c,
    QCD_all=900, QCD_bb, QCD_cc, QCD_b, QCD_c, QCD_others
  };

public:
  FatJetMatching() {}
  FatJetMatching(double jet_R, bool matchQuarks) : jetR_(jet_R), requiresQuarksContained_(matchQuarks) {}

  virtual ~FatJetMatching() {}

  std::pair<FatJetFlavor, const reco::GenParticle*> flavorJMAR(const pat::Jet *jet, const reco::GenParticleCollection& genParticles, double genRadius = 0.6);

  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > flavorLabel(const pat::Jet *jet, const reco::GenParticleCollection& genParticles, double distR);

private:
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > top_label(const pat::Jet *jet, const reco::GenParticle *parton, const reco::GenParticleCollection& genParticles, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > w_label(const pat::Jet *jet, const reco::GenParticle *parton, const reco::GenParticleCollection& genParticles, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > z_label(const pat::Jet *jet, const reco::GenParticle *parton, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > higgs_label(const pat::Jet *jet, const reco::GenParticle *parton, double distR);
  std::pair<FatJetLabel, std::vector<const reco::GenParticle*> > qcd_label(const pat::Jet *jet, const reco::GenParticleCollection& genParticles, double distR);


private:
  void printGenInfoHeader() const;
  void printGenParticleInfo(const reco::GenParticle* genParticle, const int idx) const;
  const reco::GenParticle* getFinal(const reco::GenParticle* particle);
  bool isHadronic(const reco::GenParticle* particle) const;
  std::vector<const reco::GenParticle*> getDaughterQuarks(const reco::GenParticle* particle);
  std::vector<const reco::GenParticle*> getDaughters(const reco::GenParticle* particle);
  const reco::GenParticle* getHeavyHadronAncestor(const reco::GenParticle* particle);
  std::pair<std::vector<const reco::GenParticle*>, int> getTauDaughters(const reco::GenParticle* particle);
  template <typename T>
  double maxDeltaRToDaughterQuarks(const T *center, const reco::GenParticle* mother) const {
    // mother particle needs to be the final version before decay
    double maxDeltaR = -1;
    for (const auto &q : mother->daughterRefVector()){
      if (std::abs(q->pdgId()) > ParticleID::p_b) continue;
      double deltaR = reco::deltaR(q->p4(), center->p4());
      if (deltaR > maxDeltaR) maxDeltaR = deltaR;
    }
    return maxDeltaR > 0 ? maxDeltaR : 1e9;
  }

private:
  double jetR_ = 0.8;
  bool   requiresQuarksContained_ = true;

  bool debug_ = false;
  std::unordered_set<const reco::GenParticle*> processed_;


};

}

#endif /* FATJETHELPERS_INTERFACE_FATJETMATCHING_H_ */
