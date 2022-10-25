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

std::pair<FatJetMatching::FatJetFlavor, const reco::GenParticle*> FatJetMatching::flavorJMAR(const pat::Jet* jet,
    const reco::GenParticleCollection& genParticles, double genRadius) {

  processed_.clear();

  if (debug_) {
    std::cout << "\n=======\nJet (energy, pT, eta, phi) = "
        << jet->energy() << ", " << jet->pt() << ", " << jet->eta() << ", " << jet->phi()
        << std::endl << std::endl;
    printGenInfoHeader();
    for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
      printGenParticleInfo(&genParticles[ipart], ipart);
    }
  }

  for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
    const auto *gp = &genParticles[ipart];

    if (processed_.count(gp)) continue;
    processed_.insert(gp);

    auto pdgid = std::abs(gp->pdgId());
    if (pdgid == ParticleID::p_t){
      // top
      auto top = getFinal(gp);
      // find the W and test if it's hadronic
      const reco::GenParticle *w_from_top = nullptr, *b_from_top = nullptr;
      for (const auto &dau : top->daughterRefVector()){
        if (std::abs(dau->pdgId()) == ParticleID::p_Wplus){
          w_from_top = getFinal(&(*dau));
        }else if (std::abs(dau->pdgId()) <= ParticleID::p_b){
          // ! use <= p_b ! -- can also have charms etc.
          // for quarks use the first one in the decay chain
          b_from_top = dynamic_cast<const reco::GenParticle*>(&(*dau));
        }
      }
      if (!w_from_top || !b_from_top) throw std::logic_error("[FatJetMatching::flavor] Cannot find b or W from top decay: "+std::to_string(ipart));
      if (isHadronic(w_from_top)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "top: "; printGenParticleInfo(top, -1);
          cout << "b:   "; printGenParticleInfo(b_from_top, -1);
          cout << "W:   "; printGenParticleInfo(w_from_top, -1);
        }

        double dr_jet_top = reco::deltaR(jet->p4(), top->p4());
        double dr_top_wdaus = maxDeltaRToDaughterQuarks(top, w_from_top);
        double dr_top_b     = reco::deltaR(top->p4(), b_from_top->p4());
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, top)   : " << dr_jet_top << endl;
          cout << "deltaR(top, b)     : " << dr_top_b << endl;
          cout << "deltaR(top, w daus): " << dr_top_wdaus << endl;
        }
        // top
        if (dr_top_wdaus < genRadius && dr_top_b < genRadius && dr_jet_top < genRadius) return std::make_pair(FatJetFlavor::Top, top);

        double dr_jet_w = reco::deltaR(jet->p4(), w_from_top->p4());
        double dr_w = maxDeltaRToDaughterQuarks(w_from_top, w_from_top);
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, w)     : " << dr_jet_w << endl;
          cout << "deltaR(w, w_daus)  : " << dr_w << endl;
        }
        if (dr_w < genRadius && dr_jet_w < genRadius) return std::make_pair(FatJetFlavor::W, w_from_top);
      }
    }else if (pdgid == ParticleID::p_h0) {
      // Higgs
      auto h = getFinal(gp);
      if (isHadronic(h)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "H:   "; printGenParticleInfo(h, -1);
        }
        double dr_jet_h = reco::deltaR(jet->p4(), h->p4());
        double dr_hdaus = maxDeltaRToDaughterQuarks(h, h); // only works for h->bb??
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, H)   : " << dr_jet_h << endl;
          cout << "deltaR(h, h daus): " << dr_hdaus << endl;
        }
        if (dr_hdaus < genRadius && dr_jet_h < genRadius) return std::make_pair(FatJetFlavor::H, h);
      }
    }else if (pdgid == ParticleID::p_Wplus){
      // W: not from top, or top not in jet cone
      auto w = getFinal(gp);
      if (isHadronic(w)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "W:   "; printGenParticleInfo(w, -1);
        }
        double dr_jet_w = reco::deltaR(jet->p4(), w->p4());
        double dr_wdaus = maxDeltaRToDaughterQuarks(w, w);
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, w)   : " << dr_jet_w << endl;
          cout << "deltaR(w, w daus): " << dr_wdaus << endl;
        }
        if (dr_wdaus < genRadius && dr_jet_w < genRadius) return std::make_pair(FatJetFlavor::W, w);

      }
    }else if (pdgid == ParticleID::p_Z0) {
      // Z
      auto z = getFinal(gp);
      if (isHadronic(z)) {
        if (debug_){
          using namespace std;
          cout << "jet: " << jet->polarP4() << endl;
          cout << "Z:   "; printGenParticleInfo(z, -1);
        }
        double dr_jet_z = reco::deltaR(jet->p4(), z->p4());
        double dr_zdaus = maxDeltaRToDaughterQuarks(z, z);
        if (debug_){
          using namespace std;
          cout << "deltaR(jet, Z)   : " << dr_jet_z << endl;
          cout << "deltaR(Z, Z daus): " << dr_zdaus << endl;
        }
        if (dr_zdaus < genRadius && dr_jet_z < genRadius) return std::make_pair(FatJetFlavor::Z, z);
      }
    }else {
      // ?
    }
  }

  if (genParticles.size() != processed_.size())
    throw std::logic_error("[FatJetMatching::flavor] Not all genParticles are processed!");

  const reco::GenParticle *parton = nullptr;
  double minDR = 999;
  for (const auto &gp : genParticles){
    if (gp.status() != 23) continue;
    auto pdgid = std::abs(gp.pdgId());
    if (!(pdgid<ParticleID::p_t || pdgid==ParticleID::p_g)) continue;
    auto dr = reco::deltaR(gp, *jet);
    if (dr<genRadius && dr<minDR){
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

  return std::make_pair(FatJetFlavor::Default, parton);

}



std::pair<FatJetMatching::FatJetLabel, std::vector<const reco::GenParticle*> > FatJetMatching::flavorLabel(const pat::Jet* jet,
    const reco::GenParticleCollection& genParticles, double distR) {

  processed_.clear();

  if (debug_) {
    std::cout << "\n=======\nJet (energy, pT, eta, phi) = "
        << jet->energy() << ", " << jet->pt() << ", " << jet->eta() << ", " << jet->phi()
        << std::endl << std::endl;
    printGenInfoHeader();
    for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
      printGenParticleInfo(&genParticles[ipart], ipart);
    }
  }

  for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
    const auto *gp = &genParticles[ipart];

    if (processed_.count(gp)) continue;
    processed_.insert(gp);

    auto pdgid = std::abs(gp->pdgId());
    if (pdgid == ParticleID::p_t){
      auto result = top_label(jet, gp, genParticles, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      }
    }else if (pdgid == ParticleID::p_h0 || pdgid == ParticleID::p_hc){
      auto result = higgs_label(jet, gp, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      }
    }else if (pdgid == ParticleID::p_Wplus){
      auto result = w_label(jet, gp, genParticles, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      }
    }else if (pdgid == ParticleID::p_Z0){
      auto result = z_label(jet, gp, distR);
      if (result.first != FatJetLabel::Invalid){
        return result;
      }
    }
  }

  if (genParticles.size() != processed_.size())
    throw std::logic_error("[FatJetMatching::flavor] Not all genParticles are processed!");

  return qcd_label(jet, genParticles, distR);

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

std::pair<FatJetMatching::FatJetLabel, std::vector<const reco::GenParticle*> > FatJetMatching::top_label(const pat::Jet* jet, const reco::GenParticle *parton, const reco::GenParticleCollection& genParticles, double distR)
{

  // top
  auto top = getFinal(parton);
  std::vector<const reco::GenParticle*> top_v = {top}, empty_v = {};
  // find the W and test if it's hadronic
  const reco::GenParticle *w_from_top = nullptr, *b_from_top = nullptr;
  for (const auto &dau : top->daughterRefVector()){
    if (std::abs(dau->pdgId()) == ParticleID::p_Wplus){
      w_from_top = getFinal(&(*dau));
      top_v.push_back(w_from_top);
    }else if (std::abs(dau->pdgId()) <= ParticleID::p_b){
      // ! use <= p_b ! -- can also have charms etc.
      b_from_top = dynamic_cast<const reco::GenParticle*>(&(*dau));
    }
  }
  if (!w_from_top || !b_from_top) throw std::logic_error("[FatJetMatching::top_label] Cannot find b or W from top decay!");

  // detect W decay mode
  enum WDecay {W_cq, W_qq, W_ev, W_mv, W_leptauev, W_leptaumv, W_hadtauv, W_null};
  WDecay wdecay = W_null;
  std::vector<const reco::GenParticle*> tau_daus;
  auto wdaus = getDaughters(w_from_top);
  auto pdgid_0 = std::abs(wdaus.at(0)->pdgId());
  if (pdgid_0 == ParticleID::p_eminus)  wdecay = W_ev;
  else if (pdgid_0 == ParticleID::p_muminus)  wdecay = W_mv;
  else if (pdgid_0 == ParticleID::p_tauminus){
    auto tau_daus_info = getTauDaughters(wdaus.at(0));
    tau_daus = tau_daus_info.first;
    if (tau_daus_info.second == 0)  wdecay = W_leptauev;
    else if (tau_daus_info.second == 1)  wdecay = W_leptaumv;
    else if (tau_daus_info.second == 2)  wdecay = W_hadtauv;
  }
  else { // hadronic decay
    if (std::abs(wdaus.at(0)->pdgId()) == ParticleID::p_c || std::abs(wdaus.at(1)->pdgId()) == ParticleID::p_c){
      wdecay = W_cq;
    }
    else wdecay = W_qq;
  }

  if (wdecay == W_cq || wdecay == W_qq) {
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
        auto pdgid_q1 = std::abs(wdaus.at(0)->pdgId());
        auto pdgid_q2 = std::abs(wdaus.at(1)->pdgId());
        if (debug_){
          using namespace std;
          cout << "pdgid(q1)        : " << pdgid_q1 << endl;
          cout << "pdgid(q2)        : " << pdgid_q2 << endl;
        }

        if (dr_q1<distR && dr_q2<distR){
          if (pdgid_q1 >= ParticleID::p_c || pdgid_q2 >= ParticleID::p_c) {
            return std::make_pair(FatJetLabel::Top_bcq, top_v);
          }
          else {
            return std::make_pair(FatJetLabel::Top_bqq, top_v);
          }
        }else if (dr_q1<distR && dr_q2>=distR){
          if (pdgid_q1 >= ParticleID::p_c){
            return std::make_pair(FatJetLabel::Top_bc, top_v);
          }else{
            return std::make_pair(FatJetLabel::Top_bq, top_v);
          }
        }
      }else{
        // test for W if dr(b, jet) > distR
        return w_label(jet, w_from_top, genParticles, distR);
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
    }

    if (dr_b < distR){
      if ((wdecay == W_ev || wdecay == W_mv || wdecay == W_hadtauv) && dr_l < distR) {
        if (wdecay == W_ev) return std::make_pair(FatJetLabel::Top_bev, top_v);
        else if (wdecay == W_mv) return std::make_pair(FatJetLabel::Top_bmv, top_v);
        else if (wdecay == W_hadtauv) return std::make_pair(FatJetLabel::Top_bhadtauv, top_v);
      }
      if ((wdecay == W_leptauev || wdecay == W_leptaumv) && reco::deltaR(jet->p4(), tau_daus.at(0)->p4()) < distR) {
        if (wdecay == W_leptauev) return std::make_pair(FatJetLabel::Top_bleptauev, top_v);
        else if (wdecay == W_leptaumv) return std::make_pair(FatJetLabel::Top_bleptaumv, top_v);
      }
    } else {
      // test for W if dr(b, jet) > distR
      return w_label(jet, w_from_top, genParticles, distR);
    }
  }

  return std::make_pair(FatJetLabel::Invalid, empty_v);

}

std::pair<FatJetMatching::FatJetLabel, std::vector<const reco::GenParticle*> > FatJetMatching::w_label(const pat::Jet* jet, const reco::GenParticle *parton, const reco::GenParticleCollection& genParticles, double distR)
{

  auto w = getFinal(parton);
  std::vector<const reco::GenParticle*> w_v = {w}, empty_v = {};
  if (debug_){
    using namespace std;
    cout << "jet: " << jet->polarP4() << endl;
    cout << "W:   "; printGenParticleInfo(w, -1);
  }

  // detect W decay mode
  enum WDecay {W_cq, W_qq, W_ev, W_mv, W_leptauev, W_leptaumv, W_hadtauv, W_null};
  WDecay wdecay = W_null;
  std::vector<const reco::GenParticle*> tau_daus;
  auto wdaus = getDaughters(w);
  auto pdgid_0 = std::abs(wdaus.at(0)->pdgId());
  if (pdgid_0 == ParticleID::p_eminus)  wdecay = W_ev;
  else if (pdgid_0 == ParticleID::p_muminus)  wdecay = W_mv;
  else if (pdgid_0 == ParticleID::p_tauminus){
    auto tau_daus_info = getTauDaughters(wdaus.at(0));
    tau_daus = tau_daus_info.first;
    if (tau_daus_info.second == 0)  wdecay = W_leptauev;
    else if (tau_daus_info.second == 1)  wdecay = W_leptaumv;
    else if (tau_daus_info.second == 2)  wdecay = W_hadtauv;
  }
  else { // hadronic decay
    if (std::abs(wdaus.at(0)->pdgId()) == ParticleID::p_c || std::abs(wdaus.at(1)->pdgId()) == ParticleID::p_c){
      wdecay = W_cq;
    }
    else wdecay = W_qq;
  }

  // invalid cases when decayed product to not match to jets
  if ((wdecay == W_ev || wdecay == W_mv || wdecay == W_hadtauv) && reco::deltaR(jet->p4(), wdaus.at(0)->p4()) >= distR){
    return std::make_pair(FatJetLabel::Invalid, empty_v); 
  }
  if ((wdecay == W_leptauev || wdecay == W_leptaumv) && reco::deltaR(jet->p4(), tau_daus.at(0)->p4()) >= distR){
    return std::make_pair(FatJetLabel::Invalid, empty_v); 
  }
  if ((wdecay == W_cq || wdecay == W_qq) && (reco::deltaR(jet->p4(), wdaus.at(0)->p4()) >= distR || reco::deltaR(jet->p4(), wdaus.at(1)->p4()) >= distR)){
    return std::make_pair(FatJetLabel::Invalid, empty_v); 
  }

  // comment out code for finding additional b/c-hadrons
  /*
  // now we find a valid W jets. Next step is to see if additional b/c-hadrons fall into the jet
  enum WMatchQuark {match_b, match_c, match_none};
  WMatchQuark wmatchq = match_none;
  reco::GenParticleRefVector hadrons;
  for (int icheck = 0; icheck <= 1; ++icheck) {
    if (icheck == 0) hadrons = jet->jetFlavourInfo().getbHadrons();
    else             hadrons = jet->jetFlavourInfo().getcHadrons();

    for(reco::GenParticleRefVector::const_iterator it = hadrons.begin(); it != hadrons.end(); ++it){
      if (reco::deltaR(jet->p4(), (*it)->p4()) < distR){
        // find the corresponding hadrons in the genParticle list
        const reco::GenParticle* had_from_gp = nullptr;
        for (const auto & gp: genParticles){
          if ((*it)->pdgId() == gp.pdgId() && std::abs((*it)->pt() - gp.pt())/gp.pt() <= 1e-5 && std::abs((*it)->eta() - gp.eta())/gp.eta() <= 1e-5){
            had_from_gp = &gp;
            if (debug_){
              std::cout << "find hadron: "; printGenParticleInfo(&gp, -1);
            }
          }
        }
        if (!had_from_gp) continue;

        auto ances = getHeavyHadronAncestor(had_from_gp);
        if (ances){
          if (debug_){
            std::cout << "match ances: "; printGenParticleInfo(ances, -1);
          }
          if (std::abs(ances->pdgId()) == ParticleID::p_proton){
            // std::cout << "=================" << std::endl;
            // std::cout << "match ances: "; printGenParticleInfo(ances, -1);
            // std::cout << "hadron: "; printGenParticleInfo(had_from_gp, -1);
            // std::cout << "*****************" << std::endl;
            // for (unsigned ipart = 0; ipart<genParticles.size(); ++ipart){
            //   printGenParticleInfo(&genParticles[ipart], ipart);
            // }
            if (icheck == 0)  wmatchq = match_b;
            else              wmatchq = match_c;
          }
        }
      }
    }
    if (wmatchq == match_b) break; // if find an additional b, stop finding c
  }

  if (wmatchq == match_b) {
    if (wdecay == W_cq)  return std::make_pair(FatJetLabel::W_cq_b, w_v);
    else if (wdecay == W_qq)  return std::make_pair(FatJetLabel::W_qq_b, w_v);
    else if (wdecay == W_ev)  return std::make_pair(FatJetLabel::W_ev_b, w_v);
    else if (wdecay == W_mv)  return std::make_pair(FatJetLabel::W_mv_b, w_v);
    else if (wdecay == W_leptauev)  return std::make_pair(FatJetLabel::W_leptauev_b, w_v);
    else if (wdecay == W_leptaumv)  return std::make_pair(FatJetLabel::W_leptaumv_b, w_v);
    else if (wdecay == W_hadtauv)  return std::make_pair(FatJetLabel::W_hadtauv_b, w_v);
  }
  if (wmatchq == match_c) {
    if (wdecay == W_cq)  return std::make_pair(FatJetLabel::W_cq_c, w_v);
    else if (wdecay == W_qq)  return std::make_pair(FatJetLabel::W_qq_c, w_v);
    else if (wdecay == W_ev)  return std::make_pair(FatJetLabel::W_ev_c, w_v);
    else if (wdecay == W_mv)  return std::make_pair(FatJetLabel::W_mv_c, w_v);
    else if (wdecay == W_leptauev)  return std::make_pair(FatJetLabel::W_leptauev_c, w_v);
    else if (wdecay == W_leptaumv)  return std::make_pair(FatJetLabel::W_leptaumv_c, w_v);
    else if (wdecay == W_hadtauv)  return std::make_pair(FatJetLabel::W_hadtauv_c, w_v);
  }
  else */
  {
    if (wdecay == W_cq)  return std::make_pair(FatJetLabel::W_cq, w_v);
    else if (wdecay == W_qq)  return std::make_pair(FatJetLabel::W_qq, w_v);
    else if (wdecay == W_ev)  return std::make_pair(FatJetLabel::W_ev, w_v);
    else if (wdecay == W_mv)  return std::make_pair(FatJetLabel::W_mv, w_v);
    else if (wdecay == W_leptauev)  return std::make_pair(FatJetLabel::W_leptauev, w_v);
    else if (wdecay == W_leptaumv)  return std::make_pair(FatJetLabel::W_leptaumv, w_v);
    else if (wdecay == W_hadtauv)  return std::make_pair(FatJetLabel::W_hadtauv, w_v);
  }
  return std::make_pair(FatJetLabel::Invalid, empty_v);

}

std::pair<FatJetMatching::FatJetLabel, std::vector<const reco::GenParticle*> > FatJetMatching::z_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{

  auto z = getFinal(parton);
  std::vector<const reco::GenParticle*> z_v = {z}, empty_v = {};
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
        if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
          return std::make_pair(FatJetLabel::Z_bb, z_v);
        }else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
          return std::make_pair(FatJetLabel::Z_cc, z_v);
        }else {
          return std::make_pair(FatJetLabel::Z_qq, z_v);
        }
      }
    }
  }

  return std::make_pair(FatJetLabel::Invalid, empty_v);

}

std::pair<FatJetMatching::FatJetLabel, std::vector<const reco::GenParticle*> > FatJetMatching::higgs_label(const pat::Jet* jet, const reco::GenParticle *parton, double distR)
{

  auto higgs = getFinal(parton);
  std::vector<const reco::GenParticle*> higgs_v = {higgs}, empty_v = {};

  if (debug_){
    using namespace std;
    cout << "jet: " << jet->polarP4() << endl;
    cout << "H:   "; printGenParticleInfo(higgs, -1);
  }

  bool is_hVV = false;
  if (higgs->numberOfDaughters() >= 3) {
    // e.g., h->Vqq or h->qqqq
    is_hVV = true;
  }else {
    // e.g., h->VV*
    for (const auto &p : higgs->daughterRefVector()){
      auto pdgid = std::abs(p->pdgId());
      if (pdgid == ParticleID::p_Wplus || pdgid == ParticleID::p_Z0){
        is_hVV = true;
        break;
      }
    }
  }

  if (is_hVV){
    // h->WW or h->ZZ
    std::vector<const reco::GenParticle*> hVV_daus;
    int n_el = 0, n_mu = 0, n_tau = 0, n_quarks = 0;
    bool is_lhe_hvv = false;
    // determine the decay channel according to the higgs daughters
    for (unsigned idau=0; idau<higgs->numberOfDaughters(); ++idau){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(higgs->daughter(idau));
      auto pdgid = std::abs(higgs->daughter(idau)->pdgId());
      if (pdgid >= ParticleID::p_d && pdgid <= ParticleID::p_b) {
        ++n_quarks;
        hVV_daus.push_back(dau);
      }else if (pdgid >= ParticleID::p_eminus && pdgid <= ParticleID::p_nu_tau){
        if (pdgid == ParticleID::p_eminus) ++n_el;
        else if (pdgid == ParticleID::p_muminus) ++n_mu;
        else if (pdgid == ParticleID::p_tauminus) ++n_tau;
        hVV_daus.push_back(dau);
      }else{
        // const auto d = getDaughterQuarks(getFinal(dau));
        // hVV_daus.insert(hVV_daus.end(), d.begin(), d.end());
        higgs_v.push_back(dau); // also store the gen W from Higgs
        is_lhe_hvv = true;
        auto daufinal = getFinal(dau);
        for (unsigned j=0; j<daufinal->numberOfDaughters(); ++j){
          const auto *ddau = dynamic_cast<const reco::GenParticle*>(daufinal->daughter(j));
          auto dpdgid = std::abs(ddau->pdgId());
          if (dpdgid >= ParticleID::p_d && dpdgid <= ParticleID::p_b){
            ++n_quarks;
            hVV_daus.push_back(ddau);
          }else if (dpdgid >= ParticleID::p_eminus && dpdgid <= ParticleID::p_nu_tau){
            if (dpdgid == ParticleID::p_eminus) ++n_el;
            else if (dpdgid == ParticleID::p_muminus) ++n_mu;
            else if (dpdgid == ParticleID::p_tauminus) ++n_tau;
            hVV_daus.push_back(ddau);
          }
        }
      }
    }
    auto is_neutrino = [&](int pdgid){
      pdgid = std::abs(pdgid);
      return pdgid == ParticleID::p_nu_e || pdgid == ParticleID::p_nu_mu || pdgid == ParticleID::p_nu_tau;
    };
    // auto printDaus = [&](std::vector<const reco::GenParticle*>& parts, const pat::Jet* jet){
    //   using namespace std;
    //   cout << "Found " << parts.size() << " quarks from Higgs decay" << endl;
    //   for (const auto * gp : parts){
    //     using namespace std;
    //     printGenParticleInfo(gp, -1);
    //     cout << " ... dR(q, jet) = " << reco::deltaR(*gp, *jet) << endl;
    //   }
    // };
    // if four daughters are all quarks
    if (n_quarks == 4){
      int n_daus_in_jet = 0, n_charms = 0, idx = 0;
      std::vector<int> idx_daus_in_jet;
      for (const auto *gp : hVV_daus){
        auto dr = reco::deltaR(*gp, *jet);
        if (dr < distR){
          ++n_daus_in_jet;
          idx_daus_in_jet.push_back(idx);
          if (std::abs(gp->pdgId()) == ParticleID::p_c){
            ++n_charms;
          }
        }
        ++idx;
      }
      if (n_daus_in_jet == 4){
        if (n_charms == 2){
          return std::make_pair(FatJetLabel::H_ww4q_2c, higgs_v);
        }else if (n_charms == 1){
          return std::make_pair(FatJetLabel::H_ww4q_1c, higgs_v);
        }else if (n_charms == 0){
          return std::make_pair(FatJetLabel::H_ww4q_0c, higgs_v);
        }
      }else if (n_daus_in_jet == 3){
        if (n_charms == 2){
          return std::make_pair(FatJetLabel::H_ww3q_2c, higgs_v);
        }else if (n_charms == 1){
          return std::make_pair(FatJetLabel::H_ww3q_1c, higgs_v);
        }else if (n_charms == 0){
          return std::make_pair(FatJetLabel::H_ww3q_0c, higgs_v);
        }
      }else if (n_daus_in_jet >= 2 && is_lhe_hvv){
        if (idx_daus_in_jet == std::vector<int>({0,1}) || idx_daus_in_jet == std::vector<int>({2,3})){
          // std::cout << "2q from same W!" << std::endl;
          return std::make_pair(FatJetLabel::H_ww2qsame, higgs_v);
        }else {
          // std::cout << "2q from separate Ws!" << std::endl;
          return std::make_pair(FatJetLabel::H_ww2qsep, higgs_v);
        }
      }
    }else if (n_quarks == 2 && n_tau == 0) { // lvqq
      int n_daus_in_jet = 0, n_charms = 0;
      for (const auto *gp : hVV_daus){
        if (is_neutrino(gp->pdgId())){
          continue;
        }
        auto dr = reco::deltaR(*gp, *jet);
        if (dr < distR){
          ++n_daus_in_jet;
          if (std::abs(gp->pdgId()) == ParticleID::p_c){
            ++n_charms;
          }
        }
      }
      if (n_daus_in_jet >= 3){
        // std::cout << "lvqq!" << std::endl;
        if (n_el > 0){
          if (n_charms == 1){
            return std::make_pair(FatJetLabel::H_wwevqq_1c, higgs_v);
          }else if (n_charms == 0){
            return std::make_pair(FatJetLabel::H_wwevqq_0c, higgs_v);
          }
        }else if (n_mu > 0){
          if (n_charms == 1){
            return std::make_pair(FatJetLabel::H_wwmvqq_1c, higgs_v);
          }else if (n_charms == 0){
            return std::make_pair(FatJetLabel::H_wwmvqq_0c, higgs_v);
          }
        }
      }
    }else if (n_quarks == 2 && n_tau == 1) { // tauvqq
      int tau_pos = 0;
      int n_daus_in_jet = 0, n_charms = 0;
      for(std::size_t igp = 0; igp < hVV_daus.size(); ++igp) {
        auto gp = hVV_daus[igp];
        if (is_neutrino(gp->pdgId())){
          continue;
        }
        if (std::abs(gp->pdgId()) == ParticleID::p_tauminus){
          tau_pos = igp;
          continue;
        }
        auto dr = reco::deltaR(*gp, *jet);
        if (dr < distR){
          ++n_daus_in_jet;
          if (std::abs(gp->pdgId()) == ParticleID::p_c){
            ++n_charms;
          }
        }
      }
      auto tau_daus_info = getTauDaughters(hVV_daus[tau_pos]);
      auto tau_daus = tau_daus_info.first;
      int tau_decay = tau_daus_info.second;
      auto tau = getFinal(hVV_daus[tau_pos]);
      if (n_daus_in_jet >= 2){ // both two quarks are in
        if ((tau_decay == 0 || tau_decay == 1) && reco::deltaR(tau_daus.at(0)->p4(), jet->p4()) < distR){ // tau to e/mu
          if (tau_decay == 0){
            if (n_charms == 1) return std::make_pair(FatJetLabel::H_wwleptauevqq_1c, higgs_v);
            else if (n_charms == 0) return std::make_pair(FatJetLabel::H_wwleptauevqq_0c, higgs_v);
          }else {
            if (n_charms == 1) return std::make_pair(FatJetLabel::H_wwleptaumvqq_1c, higgs_v);
            else if (n_charms == 0) return std::make_pair(FatJetLabel::H_wwleptaumvqq_0c, higgs_v);
          }
        }
        if (tau_decay == 2 && reco::deltaR(tau->p4(), jet->p4()) < distR){ // hadronic taus
          if (n_charms == 1) return std::make_pair(FatJetLabel::H_wwhadtauvqq_1c, higgs_v);
          else if (n_charms == 0) return std::make_pair(FatJetLabel::H_wwhadtauvqq_0c, higgs_v);
        }
      }
    }else {
      // std::cout << "**undefined!" << std::endl;
      return std::make_pair(FatJetLabel::Invalid, empty_v);
    }
    // std::cout << "**unmatch!" << std::endl;
    return std::make_pair(FatJetLabel::Invalid, empty_v);

    // if (debug_){
    //   using namespace std;
    //   cout << "Found " << hVV_daus.size() << " quarks from Higgs decay" << endl;
    //   for (const auto * gp : hVV_daus){
    //     using namespace std;
    //     printGenParticleInfo(gp, -1);
    //     cout << " ... dR(q, jet) = " << reco::deltaR(*gp, *jet) << endl;
    //   }
    // }

    // unsigned n_quarks_in_jet = 0;
    // for (const auto *gp : hVV_daus){
    //   auto dr = reco::deltaR(*gp, *jet);
    //   if (dr < distR){
    //     ++n_quarks_in_jet;
    //   }
    // }
    // if (n_quarks_in_jet >= 4){
    //   return std::make_pair(FatJetLabel::H_qqqq, higgs);
    // }

  }else if (isHadronic(higgs)) {
    // direct h->qq

    auto hdaus = getDaughterQuarks(higgs);
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
        if (pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_b) {
          return std::make_pair(FatJetLabel::H_bb, higgs_v);
        }else if ((pdgid_q1 == ParticleID::p_b && pdgid_q2 == ParticleID::p_c) || (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_b)) {
          return std::make_pair(FatJetLabel::H_bc, higgs_v);
        }else if (pdgid_q1 == ParticleID::p_c && pdgid_q2 == ParticleID::p_c) {
          return std::make_pair(FatJetLabel::H_cc, higgs_v);
        }else if (pdgid_q1 == ParticleID::p_s && pdgid_q2 == ParticleID::p_s) {
          return std::make_pair(FatJetLabel::H_ss, higgs_v);
        }else {
          return std::make_pair(FatJetLabel::H_qq, higgs_v);
        }
      }
    }
  }else {
    // test h->tautau
    std::vector<const reco::GenParticle*> taus;
    for (unsigned i=0; i<higgs->numberOfDaughters(); ++i){
      const auto *dau = dynamic_cast<const reco::GenParticle*>(higgs->daughter(i));
      if (std::abs(dau->pdgId()) == ParticleID::p_tauminus){
        taus.push_back(dau);
      }
    }
    if (taus.size() == 2){
      // higgs -> tautau
      // use first version or last version of the tau in dr?
      double dr_tau1    = reco::deltaR(jet->p4(), taus.at(0)->p4());
      double dr_tau2    = reco::deltaR(jet->p4(), taus.at(1)->p4());

      if (debug_){
        using namespace std;
        cout << "deltaR(jet, tau1)    : " << dr_tau1 << endl;
        cout << "deltaR(jet, tau2)    : " << dr_tau2 << endl;
      }

      auto tau1_daus_info = getTauDaughters(taus.at(0));
      auto tau2_daus_info = getTauDaughters(taus.at(1));

      // let hadronic tau be the last
      if (tau1_daus_info.second == 2 && tau2_daus_info.second < 2){
        std::swap(dr_tau1, dr_tau2);
        std::swap(tau1_daus_info, tau2_daus_info);
      }
      if (tau2_daus_info.second == 2 && dr_tau2 < distR){
        if (tau1_daus_info.second == 0 && reco::deltaR(tau1_daus_info.first.at(0)->p4(), jet->p4()) < distR){
          return std::make_pair(FatJetLabel::H_leptauehadtau, higgs_v);
        }else if (tau1_daus_info.second == 1 && reco::deltaR(tau1_daus_info.first.at(0)->p4(), jet->p4()) < distR){
          return std::make_pair(FatJetLabel::H_leptaumhadtau, higgs_v);
        }else if (tau1_daus_info.second == 2 && dr_tau1 < distR){
          return std::make_pair(FatJetLabel::H_hadtauhadtau, higgs_v);
        }
      }
    }
  }

  return std::make_pair(FatJetLabel::Invalid, empty_v);

}

std::pair<FatJetMatching::FatJetLabel, std::vector<const reco::GenParticle*> > FatJetMatching::qcd_label(const pat::Jet* jet, const reco::GenParticleCollection& genParticles, double distR)
{

  const reco::GenParticle *parton = nullptr;
  std::vector<const reco::GenParticle*> parton_v;
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

  if (parton)  parton_v = {parton};
  else  parton_v = {};

  if (n_bHadrons>=2) {
    return std::make_pair(FatJetLabel::QCD_bb, parton_v);
  }else if (n_bHadrons==1){
    return std::make_pair(FatJetLabel::QCD_b, parton_v);
  }else if (n_cHadrons>=2){
    return std::make_pair(FatJetLabel::QCD_cc, parton_v);
  }else if (n_cHadrons==1){
    return std::make_pair(FatJetLabel::QCD_c, parton_v);
  }

  return std::make_pair(FatJetLabel::QCD_others, parton_v);
}
