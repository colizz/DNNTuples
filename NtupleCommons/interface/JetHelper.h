/*
 * JetHelper.hh
 *
 *  Created on: Jan 27, 2017
 *      Author: hqu
 */

#ifndef NTUPLECOMMONS_INTERFACE_JETHELPER_H_
#define NTUPLECOMMONS_INTERFACE_JETHELPER_H_

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include <map>

namespace deepntuples {

class JetHelper {
public:
  JetHelper() {}
  JetHelper(const pat::Jet *jet, const edm::Handle<reco::CandidateView> &pfcands);

  virtual ~JetHelper() {}

  // set genjet (clustered w/ and w/o neutrino)
  void setGenjetWithNu(const reco::GenJetRef &genjetRef) { genjetWithNu_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  void setGenjetWithNuSoftDrop(const reco::GenJetRef &genjetRef) { genjetWithNuSoftDrop_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  void setGenjetNoNu(const reco::GenJetRef &genjetRef) { genjetNoNu_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  void setGenjetNoNuSoftDrop(const reco::GenJetRef &genjetRef) { genjetNoNuSoftDrop_ = (genjetRef.isNull() ? nullptr : &(*genjetRef)); }
  void setScoutJet(const edm::RefToBase<reco::Jet> &scoutJetRef) { scoutJetRef_ = &scoutJetRef; }
  // ------

  // return jet constituents (PF candidates)
  const std::vector<reco::CandidatePtr>& getJetConstituents() const { return daughters_; }
  unsigned int numberOfDaughters() const { return daughters_.size(); }
  float getPuppiWeight(const reco::CandidatePtr &cand) const {
    auto iter = puppi_wgt_cache_.find(cand.key());
    if (iter == puppi_wgt_cache_.end()){
      throw cms::Exception("[JetHelper::getPuppiWeight] Cannot get puppi wgt!");
    }
    return iter->second;
  }

  const pat::Jet& jet() const { return *jet_; }
  const std::vector<const pat::Jet*>& getSubJets() const { return subjets_; }
  const std::vector<const pat::Jet*>& getUncorrSubJets() const { return uncorr_subjets_; }

  const reco::GenJet* genjetWithNu() const { return genjetWithNu_; }
  const reco::GenJet* genjetWithNuSoftDrop() const { return genjetWithNuSoftDrop_; }
  const reco::GenJet* genjetNoNu() const { return genjetNoNu_; }
  const reco::GenJet* genjetNoNuSoftDrop() const { return genjetNoNuSoftDrop_; }
  const edm::RefToBase<reco::Jet>* getScoutJet() const { return scoutJetRef_; }

  std::pair<double, double> getCorrectedPuppiSoftDropMass(const std::vector<const pat::Jet*> &puppisubjets) const; // tmp


private:
  void initializeConstituents(const edm::Handle<reco::CandidateView> &pfcands);


private:
  // data members
  const pat::Jet *jet_ = nullptr;
  const reco::GenJet *genjetWithNu_ = nullptr;
  const reco::GenJet *genjetWithNuSoftDrop_ = nullptr;
  const reco::GenJet *genjetNoNu_ = nullptr;
  const reco::GenJet *genjetNoNuSoftDrop_ = nullptr;
  const edm::RefToBase<reco::Jet>* scoutJetRef_ = nullptr;
  std::vector<const pat::Jet*> subjets_;
  std::vector<const pat::Jet*> uncorr_subjets_;
  std::vector<reco::CandidatePtr> daughters_;
  std::map<reco::CandidatePtr::key_type, float> puppi_wgt_cache_;

};

} /* namespace deepntuples */

#endif /* NTUPLECOMMONS_INTERFACE_JETHELPER_H_ */
