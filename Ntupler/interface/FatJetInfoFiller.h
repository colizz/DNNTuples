/*
 * FatJetInfoFiller.h
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#ifndef NTUPLER_INTERFACE_FATJETINFOFILLER_H_
#define NTUPLER_INTERFACE_FATJETINFOFILLER_H_

#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"
#include "DataFormats/BTauReco/interface/BoostedDoubleSVTagInfo.h"

#include "DeepNTuples/NtupleCommons/interface/NtupleBase.h"
#include "DeepNTuples/FatJetHelpers/interface/FatJetMatching.h"


namespace deepntuples {

class FatJetInfoFiller: public NtupleBase {
public:
  FatJetInfoFiller() : FatJetInfoFiller("") {}
  FatJetInfoFiller(std::string branchName, double jetR=0.8) : NtupleBase(branchName, jetR), fjmatch_(jetR, true) {}
  virtual ~FatJetInfoFiller() {}

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
  FatJetMatching fjmatch_;
  bool useReclusteredJets_ = false;
  bool isQCDSample_ = false;
  bool isTTBarSample_ = false;
  bool isHVV2DVarMassSample_ = false;
  bool isTrainSample_ = false;
  bool isMDTagger_ = true;
  bool fillSeparateLabels_ = false;

  bool sample_use_pythia_ = false;
  bool sample_use_herwig_ = false;
  bool sample_use_madgraph_ = false;

  std::string fjTagInfoName;
  std::string fjName;

  bool debug_ = false;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;

  std::vector<std::string> labelTop_ = {
    "bWcs", "bWqq", "bWc", "bWs", "bWq", "bWev", "bWmv", "bWtauev", "bWtaumv", "bWtauhv",
    "Wcs", "Wqq", "Wev", "Wmv", "Wtauev", "Wtaumv", "Wtauhv"};
  std::vector<std::string> labelW_ = {"cs", "qq", "ev", "mv", "tauev", "taumv", "tauhv"};
  std::vector<std::string> labelZ_ = {"bb", "cc", "ss", "qq"};
  std::vector<std::string> labelH2p_ = {
    "bb", "cc", "ss", "qq", "bc", "bs", "cs", "gg",
    "ee", "mm", "tauhtaue", "tauhtaum", "tauhtauh"};
  std::vector<std::string> labelHWW_ = {
    "cscs", "csqq", "qqqq",
    "csc", "css", "csq", "qqc", "qqs", "qqq",
    "csev", "qqev", "csmv", "qqmv", "cstauev", "qqtauev", "cstaumv", "qqtaumv", "cstauhv", "qqtauhv"};
  std::vector<std::string> labelHZZ_ = {
    "bbbb", "bbcc", "bbss", "bbqq", "cccc", "ccss", "ccqq", "ssss", "ssqq", "qqqq",
    "bbb", "bbc", "bbs", "bbq", "ccb", "ccc", "ccs", "ccq", "ssb", "ssc", "sss", "ssq", "qqb", "qqc", "qqs", "qqq",
    "bbee", "bbmm", "bbe", "bbm", "bee", "bmm", "bbtauhtaue", "bbtauhtaum", "bbtauhtauh", "btauhtaue", "btauhtaum", "btauhtauh", "ccee", "ccmm", "cce", "ccm", "cee", "cmm", "cctauhtaue", "cctauhtaum", "cctauhtauh", "ctauhtaue", "ctauhtaum", "ctauhtauh", "ssee", "ssmm", "sse", "ssm", "see", "smm", "sstauhtaue", "sstauhtaum", "sstauhtauh", "stauhtaue", "stauhtaum", "stauhtauh", "qqee", "qqmm", "qqe", "qqm", "qee", "qmm", "qqtauhtaue", "qqtauhtaum", "qqtauhtauh", "qtauhtaue", "qtauhtaum", "qtauhtauh"};
  std::vector<std::string> labelQCD_ = {"bb", "cc", "b", "c", "others"};
  std::vector<std::string> labelH2pExt_ = {
    "btaue", "btaum", "btauh", "aa", "x5", "x6", "x7", "x8", "x9", "x10",
    "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19", "x20",
    "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30"}; // reserve 30 places
  std::vector<std::string> labelHHVExt_ = {
    "aabb", "aacc", "aass", "aaqq", "aabc", "aacs", "aabq", "aacq", "aasq", "aagg", "aaee", "aamm", "aatauhtaue", "aatauhtaum", "aatauhtauh",
    "aab", "aac", "aas", "aaq", "aag", "aae", "aam", "aataue", "aataum", "aatauh",
    "abb", "acc", "ass", "aqq", "abc", "acs", "abq", "acq", "asq", "agg", "aee", "amm", "atauhtaue", "atauhtaum", "atauhtauh"};
  std::vector<std::string> labels_ = {};

};

} /* namespace deepntuples */

#endif /* NTUPLER_INTERFACE_FATJETINFOFILLER_H_ */
