/*
 * FatJetMatching.cc
 *
 *  Created on: Feb 1, 2017
 *      Author: hqu
 */

#include "DeepNTuples/FatJetHelpers/interface/FatJetMatching.h"

#include <unordered_set>
#include "TString.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace deepntuples;

void FatJetMatching::flavorLabel(const pat::Jet* jet,
    const reco::GenParticleCollection& genParticles, double distR, bool isMDTagger) {

  processed_.clear();
  clearResult();

  if (debug_print_genparts_table_) {
    std::cout << "\n=======\nJet (energy, pT, eta, phi) = "
        << jet->energy() << ", " << jet->pt() << ", " << jet->eta() << ", " << jet->phi()
        << std::endl << std::endl;
    printGenInfoHeader();
    for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
      printGenParticleInfo(&genParticles[ipart], ipart);
    }
  }

  bool found_higgs = false;
  bool found_reslike_higgs = false;
  for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
    const auto *gp = &genParticles[ipart];

    if (processed_.count(gp)) continue;
    processed_.insert(gp);

    auto pdgid = std::abs(gp->pdgId());
    if (pdgid == ParticleID::p_t){
      clearResult();
      top_label(jet, gp, genParticles, distR);
      if (getResult().label != "Invalid"){
        return;
      }
    }else if (pdgid == ParticleID::p_H0 || pdgid == ParticleID::p_Hplus || pdgid == ParticleID::p_Hbsm || pdgid == ParticleID::p_LQbsm || (!found_reslike_higgs && pdgid == ParticleID::p_h0)){
      found_higgs = true;
      if (pdgid != ParticleID::p_h0){
        found_reslike_higgs = true;
      }
      // Higgs found in the record so we'll stop recongnizing
      // following W or Zs as Higgs (in case of HWW/ZZ)
      clearResult();
      higgs_label(jet, gp, distR);
      if (getResult().label != "Invalid"){
        return;
      }
    }else if (!found_higgs && isMDTagger && (pdgid == ParticleID::p_Wplus || pdgid == ParticleID::p_Z0)){
      // for MD tagger, treat W/Z boson also as a Higgs boson
      clearResult();
      higgs_label(jet, gp, distR);
      if (getResult().label != "Invalid"){
        return;
      }
    }
    // special W/Z labels for non-MD tagger
    else if (!isMDTagger && pdgid == ParticleID::p_Wplus){
      clearResult();
      w_label(jet, gp, distR, /*is_from_top=*/ false);
      if (getResult().label != "Invalid"){
        return;
      }
    }else if (!isMDTagger && pdgid == ParticleID::p_Z0){
      clearResult();
      z_label(jet, gp, distR);
      if (getResult().label != "Invalid"){
        return;
      }
    }
  }

  if (genParticles.size() != processed_.size())
    throw std::logic_error("[FatJetMatching::flavor] Not all genParticles are processed!");

  clearResult();
  qcd_label(jet, genParticles, distR);

}


void FatJetMatching::printGenInfoHeader() const {
  using namespace std;
  cout    << right << setw(6) << "#" << " " << setw(10) << "pdgId"
      << "  " << "Chg" << "  " << setw(10) << "Mass" << "  " << setw(48) << " Momentum"
      << left << "  " << setw(10) << "Mothers" << " " << setw(30) << "Daughters" << endl;
}

void FatJetMatching::printGenParticleInfo(const reco::GenParticle* genParticle, const int idx) const {
  using namespace std;
  cout  << right << setw(3) << genParticle->status();
  cout  << right << setw(3) << idx << " " << setw(10) << genParticle->pdgId() << "  ";
  cout  << right << "  " << setw(3) << genParticle->charge() << "  " << TString::Format("%10.3g", genParticle->mass() < 1e-5 ? 0 : genParticle->mass());
  cout  << left << setw(50) << TString::Format("  (E=%6.4g pT=%6.4g eta=%7.3g phi=%7.3g)", genParticle->energy(), genParticle->pt(), genParticle->eta(), genParticle->phi());

  TString                     mothers;
  for (unsigned int iMom = 0; iMom < genParticle->numberOfMothers(); ++iMom) {
    if (mothers.Length())     mothers        += ",";
    mothers   += genParticle->motherRef(iMom).key();
  }
  cout << "  " << setw(10) << mothers;
  TString                     daughters;
  for (unsigned int iDau = 0; iDau < genParticle->numberOfDaughters(); ++iDau) {
    if (daughters.Length())   daughters      += ",";
    daughters += genParticle->daughterRef(iDau).key();
  }
  cout << " " << setw(30) << daughters << endl;
}

const reco::GenParticle* FatJetMatching::getFinal(const reco::GenParticle* particle) {
  // will mark intermediate particles as processed
  if (!particle) return nullptr;
  processed_.insert(particle);
  const reco::GenParticle *final = particle;

  while (final->numberOfDaughters()) {
    const reco::GenParticle *chain = nullptr;
    for (unsigned idau = 0; idau < final->numberOfDaughters(); ++idau){
      if (final->daughter(idau)->pdgId() == particle->pdgId()) {
        chain = dynamic_cast<const reco::GenParticle*>(final->daughter(idau));
        processed_.insert(chain);
        break;
      }
    }
    if (!chain) break;
    final = chain;
  }
  return final;
}

bool FatJetMatching::isHadronic(const reco::GenParticle* particle) const {
  // particle needs to be the final version before decay
  if (!particle) throw std::invalid_argument("[FatJetMatching::isHadronic()] Null particle!");
  for(const auto &dau : particle->daughterRefVector()){
    auto pdgid = std::abs(dau->pdgId());
    if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b) return true;
  }
  return false;
}

bool FatJetMatching::isNeutrino(const reco::GenParticle* particle) const {
  // particle needs to be the final version before decay
  if (!particle) throw std::invalid_argument("[FatJetMatching::isHadronic()] Null particle!");
  auto pdgid = std::abs(particle->pdgId());
  return pdgid == ParticleID::p_nu_e || pdgid == ParticleID::p_nu_mu || pdgid == ParticleID::p_nu_tau;
}

std::vector<const reco::GenParticle*> FatJetMatching::getDaughterQuarks(const reco::GenParticle* particle) {
  std::vector<const reco::GenParticle*> daughters;

  for (unsigned i=0; i<particle->numberOfDaughters(); ++i){
    const auto *dau = dynamic_cast<const reco::GenParticle*>(particle->daughter(i));
    auto pdgid = std::abs(dau->pdgId());
    if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b){
      daughters.push_back(dau);
    }
  }

  return daughters;
}

std::vector<const reco::GenParticle*> FatJetMatching::getDaughters(const reco::GenParticle* particle) {
  std::vector<const reco::GenParticle*> daughters;

  for (unsigned i=0; i<particle->numberOfDaughters(); ++i){
    const auto *dau = dynamic_cast<const reco::GenParticle*>(particle->daughter(i));
    auto pdgid = std::abs(dau->pdgId());
    if (pdgid != ParticleID::p_nu_e && pdgid != ParticleID::p_nu_mu && pdgid != ParticleID::p_nu_tau){
      daughters.push_back(dau);
    }
  }
  return daughters;
}

const reco::GenParticle* FatJetMatching::getHeavyHadronAncestor(const reco::GenParticle* particle) {
  auto pdgid = std::abs(particle->pdgId());
  int type;
  if ((pdgid >= 400 && pdgid < 500) || (pdgid >= 4000 && pdgid < 5000))  type = 4;
  else if ((pdgid >= 500 && pdgid < 600) || (pdgid >= 5000 && pdgid < 6000))  type = 5;
  else  return nullptr;

  auto final = particle;
  const reco::GenParticle *ancestor = nullptr;
  while (final->numberOfMothers()){
    // std::cout << "New iteration" << std::endl;
    for (unsigned imom = 0; imom < final->numberOfMothers(); ++imom){
      // std::cout << "checking mom: "; printGenParticleInfo(dynamic_cast<const reco::GenParticle*>(final->mother(imom)), -1);
      auto pdgid_m = std::abs(final->mother(imom)->pdgId());
      if (pdgid_m == type || (pdgid_m >= type * 100 && pdgid_m < (type+1) * 100) || (pdgid_m >= type * 1000 && pdgid_m < (type+1) * 1000)){
        final = dynamic_cast<const reco::GenParticle*>(final->mother(imom));
        break;
      }
      if (pdgid_m == ParticleID::p_proton || pdgid_m == ParticleID::p_t || pdgid_m == ParticleID::p_Z0 || pdgid_m == ParticleID::p_Wplus || pdgid_m == ParticleID::p_h0){
        ancestor = dynamic_cast<const reco::GenParticle*>(final->mother(imom));
        break;
      }
      if (imom == final->numberOfMothers() - 1) return nullptr;
    }
    if (ancestor) break;
  }
  return ancestor;
}

std::pair<std::vector<const reco::GenParticle*>, int> FatJetMatching::getTauDaughters(const reco::GenParticle* particle) {
  auto tau = getFinal(particle);
  auto daughters = getDaughters(tau);
  for (const auto & dau: daughters){
    auto pdgid = std::abs(dau->pdgId());
    if (pdgid == ParticleID::p_eminus)  return std::make_pair(daughters, 0);
    if (pdgid == ParticleID::p_muminus)  return std::make_pair(daughters, 1);
  }
  return std::make_pair(daughters, 2); // hadronic mode
}

void FatJetMatching::top_label(const pat::Jet* jet, const reco::GenParticle *parton, const reco::GenParticleCollection& genParticles, double distR)
{
  // top
  auto top = getFinal(parton);
  getResult().resParticles.push_back(top); // push a genparticle top quark

  // find the W and test if it's hadronic
  const reco::GenParticle *w_from_top = nullptr, *b_from_top = nullptr;
  for (const auto &dau : top->daughterRefVector()){
    if (std::abs(dau->pdgId()) == ParticleID::p_Wplus){
      w_from_top = getFinal(&(*dau));
      getResult().resParticles.push_back(w_from_top); // push the W from top
    }else if (std::abs(dau->pdgId()) <= ParticleID::p_b){
      // ! use <= p_b ! -- can also have charms etc.
      b_from_top = dynamic_cast<const reco::GenParticle*>(&(*dau));
    }
  }
  if (!w_from_top || !b_from_top) throw std::logic_error("[FatJetMatching::top_label] Cannot find b or W from top decay!");

  // detect W decay mode
  enum WDecay {W_cs, W_qq, W_ev, W_mv, W_tauev, W_taumv, W_tauhv, W_null};
  WDecay wdecay = W_null;
  std::vector<const reco::GenParticle*> tau_daus;
  auto wdaus = getDaughters(w_from_top);
  auto pdgid_0 = std::abs(wdaus.at(0)->pdgId());
  if (pdgid_0 == ParticleID::p_eminus)  wdecay = W_ev;
  else if (pdgid_0 == ParticleID::p_muminus)  wdecay = W_mv;
  else if (pdgid_0 == ParticleID::p_tauminus){
    auto tau_daus_info = getTauDaughters(wdaus.at(0));
    tau_daus = tau_daus_info.first;
    if (tau_daus_info.second == 0)  wdecay = W_tauev;
    else if (tau_daus_info.second == 1)  wdecay = W_taumv;
    else if (tau_daus_info.second == 2)  wdecay = W_tauhv;
  }
  else { // hadronic decay
    if (std::abs(wdaus.at(0)->pdgId()) == ParticleID::p_c || std::abs(wdaus.at(1)->pdgId()) == ParticleID::p_c){
      wdecay = W_cs;
    }
    else wdecay = W_qq;
  }

  if (wdecay == W_cs || wdecay == W_qq) {
    if (debug_){
      using namespace std;
      cout << "jet: " << jet->polarP4() << endl;
      cout << "top: "; printGenParticleInfo(top, -1);
      cout << "b:   "; printGenParticleInfo(b_from_top, -1);
      cout << "W:   "; printGenParticleInfo(w_from_top, -1);
    }

    if (wdaus.size() < 2) throw std::logic_error("[FatJetMatching::top_label] W decay has less than 2 quarks!");
//    if (wdaus.size() >= 2)
    {
      double dr_b     = reco::deltaR(jet->p4(), b_from_top->p4());
      double dr_q1    = reco::deltaR(jet->p4(), wdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), wdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(wdaus.at(0), wdaus.at(1));
      }

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, b)     : " << dr_b << endl;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
      }

      if (dr_b < distR){
        getResult().particles.push_back(b_from_top); // push to particle list
        auto pdgid_q1 = std::abs(wdaus.at(0)->pdgId());
        auto pdgid_q2 = std::abs(wdaus.at(1)->pdgId());
        if (debug_){
          using namespace std;
          cout << "pdgid(q1)        : " << pdgid_q1 << endl;
          cout << "pdgid(q2)        : " << pdgid_q2 << endl;
        }

        if (dr_q1<distR && dr_q2<distR){
          getResult().particles.push_back(wdaus.at(0));
          getResult().particles.push_back(wdaus.at(1));
          if (pdgid_q1 >= ParticleID::p_c || pdgid_q2 >= ParticleID::p_c) {
            getResult().label = "Top_bWcs";
          }
          else {
            getResult().label = "Top_bWqq";
          }
          return;
        }else if (dr_q1<distR && dr_q2>=distR){
          getResult().particles.push_back(wdaus.at(0)); // only push the particle within the jet
          if (pdgid_q1 >= ParticleID::p_c){
            getResult().label = "Top_bWc";
          }else if (pdgid_q1 == ParticleID::p_s){
            getResult().label = "Top_bWs";
          }else{
            getResult().label = "Top_bWq";
          }
          return;
        }
      }else{
        // test for W if dr(b, jet) > distR
        w_label(jet, w_from_top, distR, /*is_from_top=*/ true);
        return;
      }
    }
  } else {
    // leptonic W
    if (debug_){
      using namespace std;
      cout << "jet: " << jet->polarP4() << endl;
      cout << "top: "; printGenParticleInfo(top, -1);
      cout << "b:   "; printGenParticleInfo(b_from_top, -1);
      cout << "W:   "; printGenParticleInfo(w_from_top, -1);
    }

    const reco::GenParticle* lep = nullptr;
    for (unsigned i=0; i<w_from_top->numberOfDaughters(); ++i){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(w_from_top->daughter(i));
      auto pdgid = std::abs(dau->pdgId());
      if (pdgid == ParticleID::p_eminus || pdgid == ParticleID::p_muminus || pdgid == ParticleID::p_tauminus){
        // use final version here!
        lep = getFinal(dau); break;
      }
    }

    if (!lep) throw std::logic_error("[FatJetMatching::top_label] Cannot find charged lepton from leptonic W decay!");

    double dr_b     = reco::deltaR(jet->p4(), b_from_top->p4());
    double dr_l     = reco::deltaR(jet->p4(), lep->p4());
    if (debug_){
      using namespace std;
      cout << "deltaR(jet, b)     : " << dr_b << endl;
      cout << "deltaR(jet, l)     : " << dr_l << endl;
      cout << "pdgid(l)           : " << lep->pdgId() << endl;
      if (!tau_daus.empty()) {
        cout << "delta(jet, tau-daus) : " << reco::deltaR(jet->p4(), tau_daus.at(0)->p4()) << endl;
        cout << "pdgid(tau-daus)      : " << tau_daus.at(0)->pdgId() << endl;
      }
    }

    if (dr_b < distR){
      getResult().particles.push_back(b_from_top);
      if ((wdecay == W_ev || wdecay == W_mv || wdecay == W_tauhv) && dr_l < distR) {
        if (wdecay == W_ev) {
          getResult().label = "Top_bWev";
          getResult().particles.push_back(lep);
        }
        else if (wdecay == W_mv) {
          getResult().label = "Top_bWmv";
          getResult().particles.push_back(lep);
        }
        else if (wdecay == W_tauhv) {
          getResult().label = "Top_bWtauhv";
          // push only visable tau daughters (exclude tau neutrino anyway)
          getResult().particles.insert(getResult().particles.end(), tau_daus.begin(), tau_daus.end());
        }
        return;
      }
      if ((wdecay == W_tauev || wdecay == W_taumv) && reco::deltaR(jet->p4(), tau_daus.at(0)->p4()) < distR) {
        if (wdecay == W_tauev)  getResult().label = "Top_bWtauev";
        else if (wdecay == W_taumv)  getResult().label = "Top_bWtaumv";
        getResult().particles.push_back(tau_daus.at(0));
        return;
      }
    } else {
      // test for W if dr(b, jet) > distR
      w_label(jet, w_from_top, distR, /*is_from_top=*/ true);
      return;
    }
  }

}

void FatJetMatching::w_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR, bool is_from_top)
{

  auto w = getFinal(parton);
  if (!is_from_top) {
    // if the function is not called from the top matching, push the W to genparts
    getResult().resParticles.push_back(w);
  }
  std::vector<const reco::GenParticle*> w_v = {w}, empty_v = {};
  if (debug_){
    using namespace std;
    cout << "jet: " << jet->polarP4() << endl;
    cout << "W:   "; printGenParticleInfo(w, -1);
  }

  // detect W decay mode
  enum WDecay {W_cs, W_qq, W_ev, W_mv, W_tauev, W_taumv, W_tauhv, W_null};
  WDecay wdecay = W_null;
  std::vector<const reco::GenParticle*> tau_daus;
  auto wdaus = getDaughters(w);
  auto pdgid_0 = std::abs(wdaus.at(0)->pdgId());
  if (pdgid_0 == ParticleID::p_eminus)  wdecay = W_ev;
  else if (pdgid_0 == ParticleID::p_muminus)  wdecay = W_mv;
  else if (pdgid_0 == ParticleID::p_tauminus){
    auto tau_daus_info = getTauDaughters(wdaus.at(0));
    tau_daus = tau_daus_info.first;
    if (tau_daus_info.second == 0)  wdecay = W_tauev;
    else if (tau_daus_info.second == 1)  wdecay = W_taumv;
    else if (tau_daus_info.second == 2)  wdecay = W_tauhv;
  }
  else { // hadronic decay
    if (std::abs(wdaus.at(0)->pdgId()) == ParticleID::p_c || std::abs(wdaus.at(1)->pdgId()) == ParticleID::p_c){
      wdecay = W_cs;
    }
    else wdecay = W_qq;
  }

  if (debug_){
    using namespace std;
    cout << "deltaR(jet, W-dau1)    : " << reco::deltaR(jet->p4(), wdaus.at(0)->p4()) << endl;
    cout << "pdgid(W-dau1)          : " << wdaus.at(0)->pdgId() << endl;
    if (wdaus.size() > 1) {
      cout << "deltaR(jet, W-dau2)    : " << reco::deltaR(jet->p4(), wdaus.at(1)->p4()) << endl;
      cout << "pdgid(W-dau2)          : " << wdaus.at(1)->pdgId() << endl;
    }
    if (!tau_daus.empty()) {
      cout << "deltaR(jet, tau-daus)  : " << reco::deltaR(jet->p4(), tau_daus.at(0)->p4()) << endl;
      cout << "pdgid(tau-daus)        : " << tau_daus.at(0)->pdgId() << endl;
    }
  }

  // invalid cases when decayed product to not match to jets
  if ((wdecay == W_ev || wdecay == W_mv || wdecay == W_tauhv) && reco::deltaR(jet->p4(), wdaus.at(0)->p4()) >= distR){
    getResult().label = "Invalid";
    return;
  }
  if ((wdecay == W_tauev || wdecay == W_taumv) && reco::deltaR(jet->p4(), tau_daus.at(0)->p4()) >= distR){
    getResult().label = "Invalid";
    return; 
  }
  if ((wdecay == W_cs || wdecay == W_qq) && (reco::deltaR(jet->p4(), wdaus.at(0)->p4()) >= distR || reco::deltaR(jet->p4(), wdaus.at(1)->p4()) >= distR)){
    getResult().label = "Invalid";
    return;
  }

  // designate labels depending on whether it's called from the top label function
  if (!is_from_top) {
    if (wdecay == W_cs)  getResult().label = "W_cs";
    else if (wdecay == W_qq)  getResult().label = "W_qq";
    else if (wdecay == W_ev)  getResult().label = "W_ev";
    else if (wdecay == W_mv)  getResult().label = "W_mv";
    else if (wdecay == W_tauev)  getResult().label = "W_tauev";
    else if (wdecay == W_taumv)  getResult().label = "W_taumv";
    else if (wdecay == W_tauhv)  getResult().label = "W_tauhv";
  }
  else {
    if (wdecay == W_cs)  getResult().label = "Top_Wcs";
    else if (wdecay == W_qq)  getResult().label = "Top_Wqq";
    else if (wdecay == W_ev)  getResult().label = "Top_Wev";
    else if (wdecay == W_mv)  getResult().label = "Top_Wmv";
    else if (wdecay == W_tauev)  getResult().label = "Top_Wtauev";
    else if (wdecay == W_taumv)  getResult().label = "Top_Wtaumv";
    else if (wdecay == W_tauhv)  getResult().label = "Top_Wtauhv";
  }
  // push gen particles of W decay products
  if (wdecay == W_cs || wdecay == W_qq){
    getResult().particles.push_back(wdaus.at(0));
    getResult().particles.push_back(wdaus.at(1));
  }
  else if (wdecay == W_ev || wdecay == W_mv){
    getResult().particles.push_back(wdaus.at(0));
  }
  else if (wdecay == W_tauev || wdecay == W_taumv){
    getResult().particles.push_back(tau_daus.at(0));
  }
  else if (wdecay == W_tauhv){
    getResult().particles.insert(getResult().particles.end(), tau_daus.begin(), tau_daus.end());
  }
}

void FatJetMatching::z_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{

  auto z = getFinal(parton);
  getResult().resParticles.push_back(z);
  if (isHadronic(z)) {
    if (debug_){
      using namespace std;
      cout << "jet: " << jet->polarP4() << endl;
      cout << "Z:   "; printGenParticleInfo(z, -1);
    }

    auto zdaus = getDaughterQuarks(z);
    if (zdaus.size() < 2) throw std::logic_error("[FatJetMatching::z_label] Z decay has less than 2 quarks!");
//    if (zdaus.size() >= 2)
    {
      double dr_q1    = reco::deltaR(jet->p4(), zdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), zdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(zdaus.at(0), zdaus.at(1));
      }
      auto pdgid_q1 = std::abs(zdaus.at(0)->pdgId());
      auto pdgid_q2 = std::abs(zdaus.at(1)->pdgId());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        cout << "pdgid(q1)        : " << pdgid_q1 << endl;
        cout << "pdgid(q2)        : " << pdgid_q2 << endl;
      }

      if (dr_q1<distR && dr_q2<distR){
        // both quarks are matched
        getResult().particles.push_back(zdaus.at(0));
        getResult().particles.push_back(zdaus.at(1));

        if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
          getResult().label = "Z_bb";
        }else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
          getResult().label = "Z_cc";
        }else if (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_s) {
          getResult().label = "Z_ss";
        }else {
          getResult().label = "Z_qq";
        }
      }
    }
  }

}

void FatJetMatching::higgs_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{

  auto higgs = getFinal(parton);
  getResult().resParticles.push_back(higgs);

  if (debug_){
    using namespace std;
    cout << "jet: " << jet->polarP4() << endl;
    cout << "H:   "; printGenParticleInfo(higgs, -1);
  }

  enum HDecay {h_2p, h_tautau, h_qtau, h_WW, h_ZZ, h_WHorZH, h_null};
  HDecay hdecay = h_null;
  auto hdaus = getDaughters(higgs);
  if (higgs->numberOfDaughters() >= 3) {
    // e.g., h->Vqq or h->qqqq
    throw std::runtime_error("[FatJetMatching::higgs_label] H decays to 3/4 objects: not implemented");
  }else {
    auto pdgid1 = std::abs(hdaus.at(0)->pdgId()), pdgid2 = std::abs(hdaus.at(1)->pdgId());
    if (pdgid1 == ParticleID::p_Wplus && pdgid2 == ParticleID::p_Wplus) {
      hdecay = h_WW;
    }else if (pdgid1 == ParticleID::p_Z0 && pdgid2 == ParticleID::p_Z0) {
      hdecay = h_ZZ;
    }else if (((pdgid1 == ParticleID::p_Wplus || pdgid1 == ParticleID::p_Z0) && pdgid2 == ParticleID::p_h0) || 
              ((pdgid2 == ParticleID::p_Wplus || pdgid2 == ParticleID::p_Z0) && pdgid1 == ParticleID::p_h0)) {
      hdecay = h_WHorZH;
    }else if (pdgid1 == ParticleID::p_tauminus && pdgid2 == ParticleID::p_tauminus) {
      hdecay = h_tautau;
    }else if (((pdgid1 >= ParticleID::p_u && pdgid1 <= ParticleID::p_b) && pdgid2 == ParticleID::p_tauminus) || 
              ((pdgid2 >= ParticleID::p_u && pdgid2 <= ParticleID::p_b) && pdgid1 == ParticleID::p_tauminus)) {
      hdecay = h_qtau;
    }else {
      hdecay = h_2p;
    }
  }

  if (hdecay == h_WW || hdecay == h_ZZ || hdecay == h_WHorZH){
    // h->WW or h->ZZ
    enum VMode {v_had, v_lep, v_aa, v_null};
    VMode hVV_modes[2] = {v_null, v_null};
    std::vector<const reco::GenParticle*> hVV_daus;
    // found daughters of WW or ZZ, and determine the W/Z decay (had or lep),
    // then switch order to make sure had daughters go first
    for (unsigned idau=0; idau<higgs->numberOfDaughters(); ++idau){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(higgs->daughter(idau));
      auto daufinal = getFinal(dau);
      getResult().resParticles.push_back(daufinal); // push the Z or W to the list

      for (unsigned j=0; j<daufinal->numberOfDaughters(); ++j){
        const auto *ddau = dynamic_cast<const reco::GenParticle*>(daufinal->daughter(j));
        // determine the V decay mode
        if (j == 0) {
          auto dpdgid = std::abs(ddau->pdgId());
          if (dpdgid >= ParticleID::p_d && dpdgid <= ParticleID::p_b)  hVV_modes[idau] = v_had;
          else if (dpdgid == ParticleID::p_gamma)  hVV_modes[idau] = v_aa;
          else  hVV_modes[idau] = v_lep;
        }
        hVV_daus.push_back(ddau);
      }
    }
    if (hVV_modes[0] == v_lep && hVV_modes[1] == v_lep) {
      // require not both V are leptonic
      return;
    }
    else if (hVV_modes[0] == v_lep && hVV_modes[1] == v_had) {
      // hadronic V goes first
      std::swap(hVV_daus.at(0), hVV_daus.at(2));
      std::swap(hVV_daus.at(1), hVV_daus.at(3));
      hVV_modes[0] = v_had;
      hVV_modes[1] = v_lep;
    }
    else if (hVV_modes[1] == v_aa) {
      // H->gamma gamma goes first
      std::swap(hVV_daus.at(0), hVV_daus.at(2));
      std::swap(hVV_daus.at(1), hVV_daus.at(3));
      std::swap(hVV_modes[0], hVV_modes[1]);
    }

    // let e/mu/tau appears before neutrinos
    if (hVV_modes[1] == v_lep) {
      auto pdgid = std::abs(hVV_daus.at(2)->pdgId());
      if (pdgid == ParticleID::p_nu_e || pdgid == ParticleID::p_nu_mu || pdgid == ParticleID::p_nu_tau) {
        std::swap(hVV_daus.at(2), hVV_daus.at(3));
      }
    }

    // go to dedicated HWW and HZZ label functions
    // result.resParticles is completed; result.label and result.particles to be assigned inside those functions
    if (hdecay == h_WW) {
      higgs_WW_label(jet, hVV_daus, distR);
    }else if (hdecay == h_ZZ) {
      higgs_ZZ_label(jet, hVV_daus, distR);
    }else if (hdecay == h_WHorZH) {
      higgs_WHorZH_label(jet, hVV_daus, distR);
    }
    return;

  }else if (hdecay == h_2p) {
    // direct h->qq

    if (hdaus.size() < 2) throw std::logic_error("[FatJetMatching::higgs_label] Higgs decay has less than 2 quarks!");
//    if (zdaus.size() >= 2)
    {
      double dr_q1    = reco::deltaR(jet->p4(), hdaus.at(0)->p4());
      double dr_q2    = reco::deltaR(jet->p4(), hdaus.at(1)->p4());
      if (dr_q1 > dr_q2){
        // swap q1 and q2 so that dr_q1<=dr_q2
        std::swap(dr_q1, dr_q2);
        std::swap(hdaus.at(0), hdaus.at(1));
      }
      auto pdgid_q1 = std::abs(hdaus.at(0)->pdgId());
      auto pdgid_q2 = std::abs(hdaus.at(1)->pdgId());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, q1)    : " << dr_q1 << endl;
        cout << "deltaR(jet, q2)    : " << dr_q2 << endl;
        cout << "pdgid(q1)        : " << pdgid_q1 << endl;
        cout << "pdgid(q2)        : " << pdgid_q2 << endl;
      }

      if (dr_q1<distR && dr_q2<distR){
        getResult().particles.push_back(hdaus.at(0));
        getResult().particles.push_back(hdaus.at(1));

        if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
          getResult().label = "H_bb";
        }else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
          getResult().label = "H_cc";
        }else if (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_s) {
          getResult().label = "H_ss";
        }else if ((pdgid_q1 == ParticleID::p_u && pdgid_q2 == ParticleID::p_u) || (pdgid_q1 == ParticleID::p_d && pdgid_q2 == ParticleID::p_d) || 
                  (pdgid_q1 == ParticleID::p_u && pdgid_q2 == ParticleID::p_d) || (pdgid_q1 == ParticleID::p_d && pdgid_q2 == ParticleID::p_u)) {
          getResult().label = "H_qq";
        }else if ((pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_c) || (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_b)) {
          getResult().label = "H_bc";
        }else if ((pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_s) || (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_b)) {
          getResult().label = "H_bs";
        }else if ((pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_s) || (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_c)) {
          getResult().label = "H_cs";
        }else if (pdgid_q1 == ParticleID::p_g && pdgid_q2 == ParticleID::p_g) {
          getResult().label = "H_gg";
        }else if (pdgid_q1 == ParticleID::p_eminus && pdgid_q2 == ParticleID::p_eminus) {
          getResult().label = "H_ee";
        }else if (pdgid_q1 == ParticleID::p_muminus && pdgid_q2 == ParticleID::p_muminus) {
          getResult().label = "H_mm";
        }else if (pdgid_q1 == ParticleID::p_gamma && pdgid_q2 == ParticleID::p_gamma) {
          getResult().label = "Hext_aa";
        }
      }
      return;
    }
  }else if (hdecay == h_tautau) {
    if (hdaus.size() == 2){
      // higgs -> tautau
      // use first version or last version of the tau in dr?
      double dr_tau1    = reco::deltaR(jet->p4(), hdaus.at(0)->p4());
      double dr_tau2    = reco::deltaR(jet->p4(), hdaus.at(1)->p4());

      auto tau1_daus_info = getTauDaughters(hdaus.at(0));
      auto tau2_daus_info = getTauDaughters(hdaus.at(1));

      if (debug_){
        using namespace std;
        cout << "tau1 decay ID: " << tau1_daus_info.second << endl;
        cout << "tau2 decay ID: " << tau2_daus_info.second << endl;
        cout << "deltaR(jet, tau1)    : " << dr_tau1 << endl;
        cout << "deltaR(jet, tau1-dau): " << reco::deltaR(tau1_daus_info.first.at(0)->p4(), jet->p4()) << endl;
        cout << "deltaR(jet, tau2)    : " << dr_tau2 << endl;
        cout << "deltaR(jet, tau2-dau): " << reco::deltaR(tau2_daus_info.first.at(0)->p4(), jet->p4()) << endl;
      }

      // let hadronic tau be the first one
      if (tau1_daus_info.second < 2 && tau2_daus_info.second == 2){
        std::swap(dr_tau1, dr_tau2);
        std::swap(hdaus.at(0), hdaus.at(1));
        std::swap(tau1_daus_info, tau2_daus_info);
      }
      // the first tau must be hadronic
      if (tau1_daus_info.second == 2 && dr_tau1 < distR){
        // push the first tauh
        getResult().particles.insert(getResult().particles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());

        // inspect the second tau
        if ((tau2_daus_info.second == 0 || tau2_daus_info.second == 1) && reco::deltaR(tau2_daus_info.first.at(0)->p4(), jet->p4()) < distR){
          getResult().particles.push_back(tau2_daus_info.first.at(0));
          if (tau2_daus_info.second == 0) {
            getResult().label = "H_tauhtaue";
          }else {
            getResult().label = "H_tauhtaum";
          }
        }else if (tau2_daus_info.second == 2 && dr_tau2 < distR){
          getResult().particles.insert(getResult().particles.end(), tau2_daus_info.first.begin(), tau2_daus_info.first.end());
          getResult().label = "H_tauhtauh";
        }
      }
      return;
    }
  }else if (hdecay == h_qtau) { // customized categories
    if (hdaus.size() == 2){
      // resonance -> qtau
      if (std::abs(hdaus.at(0)->pdgId()) == ParticleID::p_tauminus){
        std::swap(hdaus.at(0), hdaus.at(1));
      }
      double dr_q      = reco::deltaR(jet->p4(), hdaus.at(0)->p4());
      double dr_tau    = reco::deltaR(jet->p4(), hdaus.at(1)->p4());

      auto tau_daus_info = getTauDaughters(hdaus.at(1));

      if (debug_){
        using namespace std;
        cout << "tau decay ID: " << tau_daus_info.second << endl;
        cout << "deltaR(jet, q)      : " << dr_q << endl;
        cout << "deltaR(jet, tau)    : " << dr_tau << endl;
        cout << "deltaR(jet, tau-dau): " << reco::deltaR(tau_daus_info.first.at(0)->p4(), jet->p4()) << endl;
      }

      if (dr_q < distR){
        // check the quark type
        std::string qtype = "";
        if (std::abs(hdaus.at(0)->pdgId()) == ParticleID::p_b) {
          qtype = "b";
        }else if (std::abs(hdaus.at(0)->pdgId()) == ParticleID::p_c) {
          qtype = "c";
        }else if (std::abs(hdaus.at(0)->pdgId()) == ParticleID::p_s) {
          qtype = "s";
        }else if (std::abs(hdaus.at(0)->pdgId()) == ParticleID::p_u || std::abs(hdaus.at(0)->pdgId()) == ParticleID::p_d) {
          qtype = "q";
        }else if (std::abs(hdaus.at(0)->pdgId()) == ParticleID::p_g) {
          qtype = "g";
        }
        // push the first quark
        getResult().particles.push_back(hdaus.at(0));

        // inspect the tau
        if ((tau_daus_info.second == 0 || tau_daus_info.second == 1) && reco::deltaR(tau_daus_info.first.at(0)->p4(), jet->p4()) < distR){
          getResult().particles.push_back(tau_daus_info.first.at(0));
          if (tau_daus_info.second == 0) {
            getResult().label = "Cust_" + qtype + "taue";
          }else {
            getResult().label = "Cust_" + qtype + "taum";
          }
        }else if (tau_daus_info.second == 2 && dr_tau < distR){
          getResult().particles.insert(getResult().particles.end(), tau_daus_info.first.begin(), tau_daus_info.first.end());
          getResult().label = "Cust_" + qtype + "tauh";
        }
      }
      return;
    }
  }
}

void FatJetMatching::higgs_WW_label(const pat::Jet* jet, std::vector<const reco::GenParticle*>& hVV_daughters, double distR)
{
  enum WDecay {W_cs, W_qq, W_ev, W_mv, W_tauev, W_taumv, W_tauhv, W_null};

  // determine the decay mode of two Ws
  WDecay wsdecay[2] = {W_null, W_null};
  int daus_matched[4] = {0, 0, 0, 0};
  for (int i=0; i<2; ++i) {
    const reco::GenParticle* daus[2] = {hVV_daughters.at(i*2), hVV_daughters.at(i*2+1)};

    if (std::abs(daus[0]->pdgId()) >= ParticleID::p_d && std::abs(daus[0]->pdgId()) <= ParticleID::p_b) {
      // W->qq
      if (std::abs(daus[0]->pdgId()) == ParticleID::p_c || std::abs(daus[1]->pdgId()) == ParticleID::p_c) {
        wsdecay[i] = W_cs;
        // the c quark goes first
        if (std::abs(daus[0]->pdgId()) != ParticleID::p_c) {
          std::swap(daus[0], daus[1]);
        }
      }
      else {
        wsdecay[i] = W_qq;
      }

      // check if the quarks are matched to the jet
      for (int j=0; j<2; ++j) {
        if (reco::deltaR(jet->p4(), daus[j]->p4()) < distR) {
          daus_matched[i*2+j] = 1;
          getResult().particles.push_back(daus[j]);
        }
      }
    }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_eminus) {
      wsdecay[i] = W_ev;
      if (reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
        daus_matched[i*2] = 1;
        getResult().particles.push_back(daus[0]);
      }
    }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_muminus) {
      wsdecay[i] = W_mv;
      if (reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
        daus_matched[i*2] = 1;
        getResult().particles.push_back(daus[0]);
      }
    }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_tauminus) {
      auto tau_daus_info = getTauDaughters(daus[0]);
      auto tau_daus = tau_daus_info.first;
      int tau_decay = tau_daus_info.second;
      if (tau_decay == 0)
        wsdecay[i] = W_tauev;
      else if (tau_decay == 1)
        wsdecay[i] = W_taumv;
      else
        wsdecay[i] = W_tauhv;

      // check if the tau (or daughter e/mu, if leptonic decay) is matched to the jet
      if ((tau_decay == 0 || tau_decay == 1) && reco::deltaR(jet->p4(), tau_daus.at(0)->p4()) < distR) {
        daus_matched[i*2] = 1;
        getResult().particles.push_back(tau_daus.at(0));
      }else if (tau_decay == 2 && reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
        daus_matched[i*2] = 1;
        getResult().particles.insert(getResult().particles.end(), tau_daus.begin(), tau_daus.end());
      }
    }

    if (debug_){
      using namespace std;
      cout << "W_" << i+1 << " decay mode  :" << wsdecay[i] << endl;
      cout <<  "  pdgid(dau1)   : " << daus[0]->pdgId() << endl;
      cout <<  "  pdgid(dau2)   : " << daus[1]->pdgId() << endl;
      cout <<  "  deltaR(jet, dau1)  : " << reco::deltaR(jet->p4(), daus[0]->p4()) << "  dau/tau-dau is matched: " << daus_matched[i*2] << endl;
      cout <<  "  deltaR(jet, dau2)  : " << reco::deltaR(jet->p4(), daus[1]->p4()) << "  dau/tau-dau is matched: " << daus_matched[i*2+1] << endl;
      if (wsdecay[i] == W_tauev || wsdecay[i] == W_taumv || wsdecay[i] == W_tauhv) {
        auto tau_daus_info = getTauDaughters(daus[0]);
        auto tau_daus = tau_daus_info.first;
        int tau_decay = tau_daus_info.second;
        cout << "  tau decay mode: " << tau_decay << endl;
        cout << "  tau daughters: " << endl;
        for (auto dau : tau_daus) {
          cout << "    pdgid: " << dau->pdgId() << "  deltaR(jet, dau): " << reco::deltaR(jet->p4(), dau->p4()) << endl;
        }
      }
    }
  }

  // make labels
  // be sure not using hVV_daughters following on (order not switched)
  if (wsdecay[0] == W_cs && wsdecay[1] == W_cs) {
    if (daus_matched[0] && daus_matched[1] && daus_matched[2] && daus_matched[3])
      getResult().label = "H_WW_cscs";
    else if (daus_matched[0] && daus_matched[2] && (daus_matched[1] || daus_matched[3]))
      getResult().label = "H_WW_csc";
    else if (daus_matched[1] && daus_matched[3] && (daus_matched[0] || daus_matched[2]))
      getResult().label = "H_WW_css";
    else if (daus_matched[0] + daus_matched[1] + daus_matched[2] + daus_matched[3] >= 3)
      throw std::logic_error("[FatJetMatching::higgs_WW_label]: All 3-q matched cases should be handled");
  }
  else if (wsdecay[0] == W_cs && wsdecay[1] == W_qq) {
    if (daus_matched[0] && daus_matched[1] && daus_matched[2] && daus_matched[3])
      getResult().label = "H_WW_csqq";
    else if (daus_matched[0] && daus_matched[1] && (daus_matched[2] || daus_matched[3]))
      getResult().label = "H_WW_csq";
    else if (daus_matched[0] && !daus_matched[1] && daus_matched[2] && daus_matched[3])
      getResult().label = "H_WW_qqc";
    else if (!daus_matched[0] && daus_matched[1] && daus_matched[2] && daus_matched[3])
      getResult().label = "H_WW_qqs";
    else if (daus_matched[0] + daus_matched[1] + daus_matched[2] + daus_matched[3] >= 3)
      throw std::logic_error("[FatJetMatching::higgs_WW_label]: All 3-q matched cases should be handled");
  }
  else if (wsdecay[0] == W_qq && wsdecay[1] == W_cs) {
    if (daus_matched[2] && daus_matched[3] && daus_matched[0] && daus_matched[1])
      getResult().label = "H_WW_csqq";
    else if (daus_matched[2] && daus_matched[3] && (daus_matched[0] || daus_matched[1]))
      getResult().label = "H_WW_csq";
    else if (daus_matched[2] && !daus_matched[3] && daus_matched[0] && daus_matched[1])
      getResult().label = "H_WW_qqc";
    else if (!daus_matched[2] && daus_matched[3] && daus_matched[0] && daus_matched[1])
      getResult().label = "H_WW_qqs";
    else if (daus_matched[0] + daus_matched[1] + daus_matched[2] + daus_matched[3] >= 3)
      throw std::logic_error("[FatJetMatching::higgs_WW_label]: All 3-q matched cases should be handled");
  }
  else if (wsdecay[0] == W_qq && wsdecay[1] == W_qq) {
    if (daus_matched[0] && daus_matched[1] && daus_matched[2] && daus_matched[3])
      getResult().label = "H_WW_qqqq";
    else if (daus_matched[0] + daus_matched[1] + daus_matched[2] + daus_matched[3] == 3)
      getResult().label = "H_WW_qqq";
  }
  // then will be qqlv case
  else {
    if (daus_matched[0] && daus_matched[1] && daus_matched[2]) {
      if (wsdecay[1] == W_ev)  getResult().label = (wsdecay[0] == W_cs) ? "H_WW_csev" : "H_WW_qqev";
      else if (wsdecay[1] == W_mv)  getResult().label = (wsdecay[0] == W_cs) ? "H_WW_csmv" : "H_WW_qqmv";
      else if (wsdecay[1] == W_tauev)  getResult().label = (wsdecay[0] == W_cs) ? "H_WW_cstauev" : "H_WW_qqtauev";
      else if (wsdecay[1] == W_taumv)  getResult().label = (wsdecay[0] == W_cs) ? "H_WW_cstaumv" : "H_WW_qqtaumv";
      else if (wsdecay[1] == W_tauhv)  getResult().label = (wsdecay[0] == W_cs) ? "H_WW_cstauhv" : "H_WW_qqtauhv";
    }
  }
}

void FatJetMatching::higgs_ZZ_label(const pat::Jet* jet, std::vector<const reco::GenParticle*>& hVV_daughters, double distR)
{
  enum ZDecay {Z_bb, Z_cc, Z_ss, Z_qq, Z_ee, Z_mm, Z_tauhtaue, Z_tauhtaum, Z_tauhtauh, Z_null};

  // determine the decay mode of two Ws
  ZDecay zsdecay[2] = {Z_null, Z_null};
  int daus_matched[4] = {0, 0, 0, 0};
  for (int i=0; i<2; ++i) {
    const reco::GenParticle* daus[2] = {hVV_daughters.at(i*2), hVV_daughters.at(i*2+1)};

    if (std::abs(daus[0]->pdgId()) >= ParticleID::p_d && std::abs(daus[0]->pdgId()) <= ParticleID::p_nu_mu) {
      // Z->qq/ee/mumu
      if (std::abs(daus[0]->pdgId()) == ParticleID::p_b && std::abs(daus[1]->pdgId()) == ParticleID::p_b) {
        zsdecay[i] = Z_bb;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_c && std::abs(daus[1]->pdgId()) == ParticleID::p_c) {
        zsdecay[i] = Z_cc;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_s && std::abs(daus[1]->pdgId()) == ParticleID::p_s) {
        zsdecay[i] = Z_ss;
      }else if ((std::abs(daus[0]->pdgId()) == ParticleID::p_u && std::abs(daus[1]->pdgId()) == ParticleID::p_u) || 
                (std::abs(daus[0]->pdgId()) == ParticleID::p_d && std::abs(daus[1]->pdgId()) == ParticleID::p_d)) {
        zsdecay[i] = Z_qq;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_eminus && std::abs(daus[1]->pdgId()) == ParticleID::p_eminus) {
        zsdecay[i] = Z_ee;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_muminus && std::abs(daus[1]->pdgId()) == ParticleID::p_muminus) {
        zsdecay[i] = Z_mm;
      }
      // check if the quarks are matched to the jet
      for (int j=0; j<2; ++j) {
        if (reco::deltaR(jet->p4(), daus[j]->p4()) < distR) {
          daus_matched[i*2+j] = 1;
          getResult().particles.push_back(daus[j]);
        }
      }
    }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_tauminus && std::abs(daus[1]->pdgId()) == ParticleID::p_tauminus) {
      // Z->tautau

      auto tau1_daus_info = getTauDaughters(daus[0]);
      auto tau2_daus_info = getTauDaughters(daus[1]);

      // let hadronic tau be the first one
      if (tau1_daus_info.second < 2 && tau2_daus_info.second == 2){
        std::swap(daus[0], daus[1]);
        std::swap(tau1_daus_info, tau2_daus_info);
      }
      if (tau1_daus_info.second == 2){ // first tau must be tauh
        if (tau2_daus_info.second == 0 || tau2_daus_info.second == 1){
          if (tau2_daus_info.second == 0) {
            zsdecay[i] = Z_tauhtaue;
          }else {
            zsdecay[i] = Z_tauhtaum;
          }
          // matching of the first tauh and the second taue/taum
          if (reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
            daus_matched[i*2] = 1;
            getResult().particles.insert(getResult().particles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());
          }
          if (reco::deltaR(jet->p4(), tau2_daus_info.first.at(0)->p4()) < distR) {
            daus_matched[i*2+1] = 1;
            getResult().particles.push_back(tau2_daus_info.first.at(0));
          }
        }else if (tau2_daus_info.second == 2){
          zsdecay[i] = Z_tauhtauh;
          // matching of both tauh
          if (reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
            daus_matched[i*2] = 1;
            getResult().particles.insert(getResult().particles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());
          }
          if (reco::deltaR(jet->p4(), daus[1]->p4()) < distR) {
            daus_matched[i*2+1] = 1;
            getResult().particles.insert(getResult().particles.end(), tau2_daus_info.first.begin(), tau2_daus_info.first.end());
          }
        }
      }
    }
    if (debug_){
      using namespace std;
      cout << "Z_" << i+1 << " decay mode  :" << zsdecay[i] << endl;
      cout <<  "  pdgid(dau1)   : " << daus[0]->pdgId() << endl;
      cout <<  "  pdgid(dau2)   : " << daus[1]->pdgId() << endl;
      cout <<  "  deltaR(jet, dau1)  : " << reco::deltaR(jet->p4(), daus[0]->p4()) << "  dau/tau-dau is matched: " << daus_matched[i*2] << endl;
      cout <<  "  deltaR(jet, dau2)  : " << reco::deltaR(jet->p4(), daus[1]->p4()) << "  dau/tau-dau is matched: " << daus_matched[i*2+1] << endl;
      if (zsdecay[i] == Z_tauhtaue || zsdecay[i] == Z_tauhtaum || zsdecay[i] == Z_tauhtauh) {
        for (int j=0; j<2; j++) {
          auto tau_daus_info = getTauDaughters(daus[j]);
          auto tau_daus = tau_daus_info.first;
          int tau_decay = tau_daus_info.second;
          cout << "  tau_" << j+1 << " :" << endl;
          cout << "    tau decay mode: " << tau_decay << endl;
          cout << "    tau daughters: " << endl;
          for (auto dau : tau_daus) {
            cout << "      pdgid: " << dau->pdgId() << "  deltaR(jet, dau): " << reco::deltaR(jet->p4(), dau->p4()) << endl;
          }
        }
      }
    }

  }
  if (zsdecay[0] == Z_null || zsdecay[1] == Z_null) {
    // other Z modes will not be handled (e.g Z->tauetaue)
    return;
  }
  if (zsdecay[1] == Z_tauhtaue || zsdecay[1] == Z_tauhtaum || zsdecay[1] == Z_tauhtauh) {
    // requires both taus should be matched
    if (!daus_matched[2] || !daus_matched[3]) {
      return;
    }
  }
  if (daus_matched[0] + daus_matched[1] + daus_matched[2] + daus_matched[3] < 3) {
    // should have at least has three matched daughters to desinate the ZZ label
    return;
  }

  // make labels
  // be sure not using hVV_daughters following on (order not switched)

  std::string matched_parts_str = "";
  for (int i=0; i<2; ++i) {
    if (zsdecay[i] == Z_bb)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "bb" : "b";
    else if (zsdecay[i] == Z_cc)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "cc" : "c";
    else if (zsdecay[i] == Z_ss)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "ss" : "s";
    else if (zsdecay[i] == Z_qq)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "qq" : "q";
    else if (zsdecay[i] == Z_ee)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "ee" : "e";
    else if (zsdecay[i] == Z_mm)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "mm" : "m";
    else if (zsdecay[i] == Z_tauhtaue)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "tauhtaue" : "NULL";
    else if (zsdecay[i] == Z_tauhtaum)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "tauhtaum" : "NULL";
    else if (zsdecay[i] == Z_tauhtauh)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "tauhtauh" : "NULL";
  }
  // desinate the label
  std::map<std::vector<std::string>, std::string> acceptable_strs_map = {
    {{"bbbb"}, "H_ZZ_bbbb"}, {{"bbcc", "ccbb"}, "H_ZZ_bbcc"}, {{"bbss", "ssbb"}, "H_ZZ_bbss"}, {{"bbqq", "qqbb"}, "H_ZZ_bbqq"}, {{"bbee", "eebb"}, "H_ZZ_bbee"}, {{"bbmm", "mmbb"}, "H_ZZ_bbmm"},
    {{"bbtauhtaue", "tauhtauebb"}, "H_ZZ_bbtauhtaue"}, {{"bbtauhtaum", "tauhtaumbb"}, "H_ZZ_bbtauhtaum"}, {{"bbtauhtauh", "tauhtauhbb"}, "H_ZZ_bbtauhtauh"},
    {{"bbb"}, "H_ZZ_bbb"}, {{"bbc", "cbb"}, "H_ZZ_bbc"}, {{"bbs", "sbb"}, "H_ZZ_bbs"}, {{"bbq", "qbb"}, "H_ZZ_bbq"}, {{"bbe", "ebb"}, "H_ZZ_bbe"}, {{"bbm", "mbb"}, "H_ZZ_bbm"},
    {{"cccc"}, "H_ZZ_cccc"}, {{"ccss", "sscc"}, "H_ZZ_ccss"}, {{"ccqq", "qqcc"}, "H_ZZ_ccqq"}, {{"ccee", "eecc"}, "H_ZZ_ccee"}, {{"ccmm", "mmcc"}, "H_ZZ_ccmm"},
    {{"cctauhtaue", "tauhtauecc"}, "H_ZZ_cctauhtaue"}, {{"cctauhtaum", "tauhtaumcc"}, "H_ZZ_cctauhtaum"}, {{"cctauhtauh", "tauhtauhcc"}, "H_ZZ_cctauhtauh"},
    {{"ccb", "bcc"}, "H_ZZ_ccb"}, {{"ccc"}, "H_ZZ_ccc"}, {{"ccs", "scc"}, "H_ZZ_ccs"}, {{"ccq", "qcc"}, "H_ZZ_ccq"}, {{"cce", "ecc"}, "H_ZZ_cce"}, {{"ccm", "mcc"}, "H_ZZ_ccm"},
    {{"ssss"}, "H_ZZ_ssss"}, {{"ssqq", "qqss"}, "H_ZZ_ssqq"}, {{"ssee", "eess"}, "H_ZZ_ssee"}, {{"ssmm", "mmss"}, "H_ZZ_ssmm"},
    {{"sstauhtaue", "tauhtauess"}, "H_ZZ_sstauhtaue"}, {{"sstauhtaum", "tauhtaumss"}, "H_ZZ_sstauhtaum"}, {{"sstauhtauh", "tauhtauhss"}, "H_ZZ_sstauhtauh"},
    {{"ssb", "bss"}, "H_ZZ_ssb"}, {{"ssc", "css"}, "H_ZZ_ssc"}, {{"sss"}, "H_ZZ_sss"}, {{"ssq", "qss"}, "H_ZZ_ssq"}, {{"sse", "ess"}, "H_ZZ_sse"}, {{"ssm", "mss"}, "H_ZZ_ssm"},
    {{"qqqq"}, "H_ZZ_qqqq"}, {{"qqee", "eeqq"}, "H_ZZ_qqee"}, {{"qqmm", "mmqq"}, "H_ZZ_qqmm"},
    {{"qqtauhtaue", "tauhtaueqq"}, "H_ZZ_qqtauhtaue"}, {{"qqtauhtaum", "tauhtaumqq"}, "H_ZZ_qqtauhtaum"}, {{"qqtauhtauh", "tauhtauhqq"}, "H_ZZ_qqtauhtauh"},
    {{"qqb", "bqq"}, "H_ZZ_qqb"}, {{"qqc", "cqq"}, "H_ZZ_qqc"}, {{"qqs", "sqq"}, "H_ZZ_qqs"}, {{"qqq"}, "H_ZZ_qqq"}, {{"qqe", "eqq"}, "H_ZZ_qqe"}, {{"qqm", "mqq"}, "H_ZZ_qqm"},
    {{"eeb", "bee"}, "H_ZZ_bee"}, {{"eec", "cee"}, "H_ZZ_cee"}, {{"ees", "see"}, "H_ZZ_see"}, {{"eeq", "qee"}, "H_ZZ_qee"},
    {{"mmb", "bmm"}, "H_ZZ_bmm"}, {{"mmc", "cmm"}, "H_ZZ_cmm"}, {{"mms", "smm"}, "H_ZZ_smm"}, {{"mmq", "qmm"}, "H_ZZ_qmm"},
    {{"tauhtaueb", "btauhtaue"}, "H_ZZ_btauhtaue"}, {{"tauhtauec", "ctauhtaue"}, "H_ZZ_ctauhtaue"}, {{"tauhtaues", "stauhtaue"}, "H_ZZ_stauhtaue"}, {{"tauhtaueq", "qtauhtaue"}, "H_ZZ_qtauhtaue"},
    {{"tauhtaumb", "btauhtaum"}, "H_ZZ_btauhtaum"}, {{"tauhtaumc", "ctauhtaum"}, "H_ZZ_ctauhtaum"}, {{"tauhtaums", "stauhtaum"}, "H_ZZ_stauhtaum"}, {{"tauhtaumq", "qtauhtaum"}, "H_ZZ_qtauhtaum"},
    {{"tauhtauhb", "btauhtauh"}, "H_ZZ_btauhtauh"}, {{"tauhtauhc", "ctauhtauh"}, "H_ZZ_ctauhtauh"}, {{"tauhtauhs", "stauhtauh"}, "H_ZZ_stauhtauh"}, {{"tauhtauhq", "qtauhtauh"}, "H_ZZ_qtauhtauh"},
  };
  for (auto it = acceptable_strs_map.begin(); it != acceptable_strs_map.end(); ++it) {
    auto& acceptable_strs = it->first;
    std::string& label = it->second;

    // iterate over all string options
    for (auto it_str = acceptable_strs.begin(); it_str != acceptable_strs.end(); ++it_str) {
      if (matched_parts_str == (*it_str)) {
        getResult().label = label;
        return;
      }
    }
  }
  throw std::logic_error("[FatJetMatching::higgs_ZZ_label] Unmatched HZZ label: " + matched_parts_str);

}

void FatJetMatching::higgs_WHorZH_label(const pat::Jet* jet, std::vector<const reco::GenParticle*>& hVV_daughters, double distR)
{
  // for W/Z->anything, H->gamma gamma
  enum RDecay {H_aa, Z_bb, Z_cc, Z_ss, Z_qq, W_bc, W_cs, W_bq, W_cq, W_sq, Z_gg, Z_ee, Z_mm, W_ev, W_mv, Z_tauhtaue, Z_tauhtaum, Z_tauhtauh, W_tauev, W_taumv, W_tauhv, R_null}; // here "Z" represents both W/Z

  // determine the decay mode of WH/ZH
  RDecay rsdecay[2] = {R_null, R_null};
  int daus_matched[4] = {0, 0, 0, 0};
  for (int i=0; i<2; ++i) {
    const reco::GenParticle* daus[2] = {hVV_daughters.at(i*2), hVV_daughters.at(i*2+1)};

    if ((std::abs(daus[0]->pdgId()) >= ParticleID::p_d && std::abs(daus[0]->pdgId()) <= ParticleID::p_nu_mu) || 
        std::abs(daus[0]->pdgId()) == ParticleID::p_g || 
        std::abs(daus[0]->pdgId()) == ParticleID::p_gamma) {
      // non-tau cases
      if (std::abs(daus[0]->pdgId()) == ParticleID::p_gamma && std::abs(daus[1]->pdgId()) == ParticleID::p_gamma) {
        rsdecay[i] = H_aa;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_b && std::abs(daus[1]->pdgId()) == ParticleID::p_b) {
        rsdecay[i] = Z_bb;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_c && std::abs(daus[1]->pdgId()) == ParticleID::p_c) {
        rsdecay[i] = Z_cc;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_s && std::abs(daus[1]->pdgId()) == ParticleID::p_s) {
        rsdecay[i] = Z_ss;
      }else if ((std::abs(daus[0]->pdgId()) == ParticleID::p_u && std::abs(daus[1]->pdgId()) == ParticleID::p_u) || 
                (std::abs(daus[0]->pdgId()) == ParticleID::p_d && std::abs(daus[1]->pdgId()) == ParticleID::p_d)) {
        rsdecay[i] = Z_qq;
      }else if ((std::abs(daus[0]->pdgId()) == ParticleID::p_b && std::abs(daus[1]->pdgId()) == ParticleID::p_c)) {
        rsdecay[i] = W_bc;
      }else if ((std::abs(daus[0]->pdgId()) == ParticleID::p_c && std::abs(daus[1]->pdgId()) == ParticleID::p_b)) {
        std::swap(daus[0], daus[1]);
        rsdecay[i] = W_bc;
      }else if ((std::abs(daus[0]->pdgId()) == ParticleID::p_c && std::abs(daus[1]->pdgId()) == ParticleID::p_s)) {
        rsdecay[i] = W_cs;
      }else if ((std::abs(daus[0]->pdgId()) == ParticleID::p_s && std::abs(daus[1]->pdgId()) == ParticleID::p_c)) {
        std::swap(daus[0], daus[1]);
        rsdecay[i] = W_cs;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_b && (std::abs(daus[1]->pdgId()) == ParticleID::p_u || std::abs(daus[1]->pdgId()) == ParticleID::p_d)) {
        rsdecay[i] = W_bq;
      }else if (std::abs(daus[1]->pdgId()) == ParticleID::p_b && (std::abs(daus[0]->pdgId()) == ParticleID::p_u || std::abs(daus[0]->pdgId()) == ParticleID::p_d)) {
        std::swap(daus[0], daus[1]);
        rsdecay[i] = W_bq;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_c && (std::abs(daus[1]->pdgId()) == ParticleID::p_u || std::abs(daus[1]->pdgId()) == ParticleID::p_d)) {
        rsdecay[i] = W_cq;
      }else if (std::abs(daus[1]->pdgId()) == ParticleID::p_c && (std::abs(daus[0]->pdgId()) == ParticleID::p_u || std::abs(daus[0]->pdgId()) == ParticleID::p_d)) {
        std::swap(daus[0], daus[1]);
        rsdecay[i] = W_cq;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_s && (std::abs(daus[1]->pdgId()) == ParticleID::p_u || std::abs(daus[1]->pdgId()) == ParticleID::p_d)) {
        rsdecay[i] = W_sq;
      }else if (std::abs(daus[1]->pdgId()) == ParticleID::p_s && (std::abs(daus[0]->pdgId()) == ParticleID::p_u || std::abs(daus[0]->pdgId()) == ParticleID::p_d)) {
        std::swap(daus[0], daus[1]);
        rsdecay[i] = W_sq;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_g && std::abs(daus[1]->pdgId()) == ParticleID::p_g) {
        rsdecay[i] = Z_gg;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_eminus && std::abs(daus[1]->pdgId()) == ParticleID::p_eminus) {
        rsdecay[i] = Z_ee;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_muminus && std::abs(daus[1]->pdgId()) == ParticleID::p_muminus) {
        rsdecay[i] = Z_mm;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_eminus && std::abs(daus[1]->pdgId()) == ParticleID::p_nu_e) {
        rsdecay[i] = W_ev;
      }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_muminus && std::abs(daus[1]->pdgId()) == ParticleID::p_nu_mu) {
        rsdecay[i] = W_mv;
      }
      // check if the quarks are matched to the jet
      for (int j=0; j<2; ++j) {
        if (!isNeutrino(daus[j]) && (reco::deltaR(jet->p4(), daus[j]->p4()) < distR)) {
          daus_matched[i*2+j] = 1;
          getResult().particles.push_back(daus[j]);
        }
      }
    }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_tauminus && std::abs(daus[1]->pdgId()) == ParticleID::p_tauminus) {
      // Z->tautau

      auto tau1_daus_info = getTauDaughters(daus[0]);
      auto tau2_daus_info = getTauDaughters(daus[1]);

      // let hadronic tau be the first one
      if (tau1_daus_info.second < 2 && tau2_daus_info.second == 2){
        std::swap(daus[0], daus[1]);
        std::swap(tau1_daus_info, tau2_daus_info);
      }
      if (tau1_daus_info.second == 2){ // first tau must be tauh
        if (tau2_daus_info.second == 0 || tau2_daus_info.second == 1){
          if (tau2_daus_info.second == 0) {
            rsdecay[i] = Z_tauhtaue;
          }else {
            rsdecay[i] = Z_tauhtaum;
          }
          // matching of the first tauh and the second taue/taum
          if (reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
            daus_matched[i*2] = 1;
            getResult().particles.insert(getResult().particles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());
          }
          if (reco::deltaR(jet->p4(), tau2_daus_info.first.at(0)->p4()) < distR) {
            daus_matched[i*2+1] = 1;
            getResult().particles.push_back(tau2_daus_info.first.at(0));
          }
        }else if (tau2_daus_info.second == 2){
          rsdecay[i] = Z_tauhtauh;
          // matching of both tauh
          if (reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
            daus_matched[i*2] = 1;
            getResult().particles.insert(getResult().particles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());
          }
          if (reco::deltaR(jet->p4(), daus[1]->p4()) < distR) {
            daus_matched[i*2+1] = 1;
            getResult().particles.insert(getResult().particles.end(), tau2_daus_info.first.begin(), tau2_daus_info.first.end());
          }
        }
      }
    }else if (std::abs(daus[0]->pdgId()) == ParticleID::p_tauminus && std::abs(daus[1]->pdgId()) == ParticleID::p_nu_tau) {
      // W->tauv
      auto tau1_daus_info = getTauDaughters(daus[0]);

      if (tau1_daus_info.second == 0 || tau1_daus_info.second == 1){
        if (tau1_daus_info.second == 0) {
          rsdecay[i] = W_tauev;
        }else {
          rsdecay[i] = W_taumv;
        }
        // matching of taue/taum
        if (reco::deltaR(jet->p4(), tau1_daus_info.first.at(0)->p4()) < distR) {
          daus_matched[i*2] = 1;
          getResult().particles.push_back(tau1_daus_info.first.at(0));
        }
      }else if (tau1_daus_info.second == 2){
        rsdecay[i] = W_tauhv;
        // matching of tauh
        if (reco::deltaR(jet->p4(), daus[0]->p4()) < distR) {
          daus_matched[i*2] = 1;
          getResult().particles.insert(getResult().particles.end(), tau1_daus_info.first.begin(), tau1_daus_info.first.end());
        }
      }
    }
    if (debug_){
      using namespace std;
      cout << "Z/W/H " << i+1 << " decay mode  :" << rsdecay[i] << endl;
      cout <<  "  pdgid(dau1)   : " << daus[0]->pdgId() << endl;
      cout <<  "  pdgid(dau2)   : " << daus[1]->pdgId() << endl;
      cout <<  "  deltaR(jet, dau1)  : " << reco::deltaR(jet->p4(), daus[0]->p4()) << "  dau/tau-dau is matched: " << daus_matched[i*2] << endl;
      cout <<  "  deltaR(jet, dau2)  : " << reco::deltaR(jet->p4(), daus[1]->p4()) << "  dau/tau-dau is matched: " << daus_matched[i*2+1] << endl;
      if (rsdecay[i] == Z_tauhtaue || rsdecay[i] == Z_tauhtaum || rsdecay[i] == Z_tauhtauh || rsdecay[i] == W_tauev || rsdecay[i] == W_taumv || rsdecay[i] == W_tauhv) {
        int num_tau = 0;
        if (rsdecay[i] == Z_tauhtaue || rsdecay[i] == Z_tauhtaum || rsdecay[i] == Z_tauhtauh)  num_tau = 2;
        else if (rsdecay[i] == W_tauev || rsdecay[i] == W_taumv || rsdecay[i] == W_tauhv)  num_tau = 1;
        for (int j=0; j<num_tau; j++) {
          auto tau_daus_info = getTauDaughters(daus[j]);
          auto tau_daus = tau_daus_info.first;
          int tau_decay = tau_daus_info.second;
          cout << "  tau_" << j+1 << " :" << endl;
          cout << "    tau decay mode: " << tau_decay << endl;
          cout << "    tau daughters: " << endl;
          for (auto dau : tau_daus) {
            cout << "      pdgid: " << dau->pdgId() << "  deltaR(jet, dau): " << reco::deltaR(jet->p4(), dau->p4()) << endl;
          }
        }
      }
    }

  }
  if (rsdecay[0] == R_null || rsdecay[1] == R_null) {
    // other Z modes will not be handled (e.g Z->tauetaue)
    return;
  }
  if (rsdecay[1] == Z_tauhtaue || rsdecay[1] == Z_tauhtaum || rsdecay[1] == Z_tauhtauh) {
    // requires both taus should be matched
    if (!daus_matched[2] || !daus_matched[3]) {
      return;
    }
  }
  if (rsdecay[1] == W_tauev || rsdecay[1] == W_taumv || rsdecay[1] == W_tauhv) {
    // requires tau should be matched
    if (!daus_matched[2]) {
      return;
    }
  }
  if (daus_matched[0] + daus_matched[1] + daus_matched[2] + daus_matched[3] < 3) {
    // should have at least has three matched daughters to desinate the ZZ label
    return;
  }

  // make labels
  // be sure not using hVV_daughters following on (order not switched)

  std::string matched_parts_str = "";
  for (int i=0; i<2; ++i) {
    if (rsdecay[i] == Z_bb)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "bb" : "b";
    else if (rsdecay[i] == Z_cc)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "cc" : "c";
    else if (rsdecay[i] == Z_ss)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "ss" : "s";
    else if (rsdecay[i] == Z_qq)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "qq" : "q";
    else if (rsdecay[i] == W_bc)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "bc" : daus_matched[i*2] ? "b" : "c";
    else if (rsdecay[i] == W_cs)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "cs" : daus_matched[i*2] ? "c" : "s";
    else if (rsdecay[i] == W_bq)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "bq" : daus_matched[i*2] ? "b" : "q";
    else if (rsdecay[i] == W_cq)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "cq" : daus_matched[i*2] ? "c" : "q";
    else if (rsdecay[i] == W_sq)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "sq" : daus_matched[i*2] ? "s" : "q";
    else if (rsdecay[i] == Z_gg)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "gg" : "g";
    else if (rsdecay[i] == H_aa)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "aa" : "a";
    else if (rsdecay[i] == Z_ee)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "ee" : "e";
    else if (rsdecay[i] == Z_mm)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "mm" : "m";
    else if (rsdecay[i] == W_ev)  matched_parts_str += (daus_matched[i*2]) ? "e" : "NULL";
    else if (rsdecay[i] == W_mv)  matched_parts_str += (daus_matched[i*2]) ? "m" : "NULL";
    else if (rsdecay[i] == Z_tauhtaue)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "tauhtaue" : "NULL"; // NULL is not possible
    else if (rsdecay[i] == Z_tauhtaum)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "tauhtaum" : "NULL";
    else if (rsdecay[i] == Z_tauhtauh)  matched_parts_str += (daus_matched[i*2] && daus_matched[i*2+1]) ? "tauhtauh" : "NULL";
    else if (rsdecay[i] == W_tauev)     matched_parts_str += (daus_matched[i*2]) ? "taue" : "NULL";
    else if (rsdecay[i] == W_taumv)     matched_parts_str += (daus_matched[i*2]) ? "taum" : "NULL";
    else if (rsdecay[i] == W_tauhv)     matched_parts_str += (daus_matched[i*2]) ? "tauh" : "NULL";
  }
  // desinate the label
  std::map<std::vector<std::string>, std::string> acceptable_strs_map = {
    {{"aabb"}, "H_HVext_aabb"}, {{"aacc"}, "H_HVext_aacc"}, {{"aass"}, "H_HVext_aass"}, {{"aaqq"}, "H_HVext_aaqq"}, {{"aabc"}, "H_HVext_aabc"}, {{"aacs"}, "H_HVext_aacs"}, {{"aabq"}, "H_HVext_aabq"}, {{"aacq"}, "H_HVext_aacq"}, {{"aasq"}, "H_HVext_aasq"}, {{"aagg"}, "H_HVext_aagg"}, {{"aaee"}, "H_HVext_aaee"}, {{"aamm"}, "H_HVext_aamm"}, {{"aatauhtaue"}, "H_HVext_aatauhtaue"}, {{"aatauhtaum"}, "H_HVext_aatauhtaum"}, {{"aatauhtauh"}, "H_HVext_aatauhtauh"},
    {{"aab"}, "H_HVext_aab"}, {{"aac"}, "H_HVext_aac"}, {{"aas"}, "H_HVext_aas"}, {{"aaq"}, "H_HVext_aaq"}, {{"aag"}, "H_HVext_aag"}, {{"aae"}, "H_HVext_aae"}, {{"aam"}, "H_HVext_aam"}, {{"aataue"}, "H_HVext_aataue"}, {{"aataum"}, "H_HVext_aataum"}, {{"aatauh"}, "H_HVext_aatauh"},
    {{"abb"}, "H_HVext_abb"}, {{"acc"}, "H_HVext_acc"}, {{"ass"}, "H_HVext_ass"}, {{"aqq"}, "H_HVext_aqq"}, {{"abc"}, "H_HVext_abc"}, {{"acs"}, "H_HVext_acs"}, {{"abq"}, "H_HVext_abq"}, {{"acq"}, "H_HVext_acq"}, {{"asq"}, "H_HVext_asq"}, {{"agg"}, "H_HVext_agg"}, {{"aee"}, "H_HVext_aee"}, {{"amm"}, "H_HVext_amm"}, {{"atauhtaue"}, "H_HVext_atauhtaue"}, {{"atauhtaum"}, "H_HVext_atauhtaum"}, {{"atauhtauh"}, "H_HVext_atauhtauh"},
  };
  for (auto it = acceptable_strs_map.begin(); it != acceptable_strs_map.end(); ++it) {
    auto& acceptable_strs = it->first;
    std::string& label = it->second;

    // iterate over all string options
    for (auto it_str = acceptable_strs.begin(); it_str != acceptable_strs.end(); ++it_str) {
      if (matched_parts_str == (*it_str)) {
        getResult().label = label;
        return;
      }
    }
  }
  throw std::logic_error("[FatJetMatching::higgs_WHorZH_label] Unmatched label: " + matched_parts_str);

}

void FatJetMatching::qcd_label(const pat::Jet* jet, const reco::GenParticleCollection& genParticles, double distR)
{
  const reco::GenParticle *parton = nullptr;
  double minDR = 999;
  for (const auto &gp : genParticles){
    if (gp.status() != 23) continue;
    auto pdgid = std::abs(gp.pdgId());
    if (!(pdgid<ParticleID::p_t || pdgid==ParticleID::p_g)) continue;
    auto dr = reco::deltaR(gp, *jet);
    if (dr<distR && dr<minDR){
      minDR = dr;
      parton = &gp;
    }
  }
  if (debug_){
    using namespace std;
    if (parton){
      cout << "parton"; printGenParticleInfo(parton, -1);
      cout << "dr(jet, parton): " << minDR << endl;
    }
  }

  auto n_bHadrons = jet->jetFlavourInfo().getbHadrons().size();
  auto n_cHadrons = jet->jetFlavourInfo().getcHadrons().size();

  if (n_bHadrons>=2) {
    getResult().label = "QCD_bb";
  }else if (n_bHadrons==1){
    getResult().label = "QCD_b";
  }else if (n_cHadrons>=2){
    getResult().label = "QCD_cc";
  }else if (n_cHadrons==1){
    getResult().label = "QCD_c";
  }else{
    getResult().label = "QCD_others";
  }

}