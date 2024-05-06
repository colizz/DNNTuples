#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/BTauReco/interface/DeepBoostedJetTagInfo.h"

#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"

using namespace btagbtvdeep;

class ParticleTransformerAK8TagInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit ParticleTransformerAK8TagInfoProducer(const edm::ParameterSet &);
  ~ParticleTransformerAK8TagInfoProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  typedef std::vector<reco::DeepBoostedJetTagInfo> DeepBoostedJetTagInfoCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;
  typedef edm::View<reco::Candidate> CandidateView;

  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override {}

  void fillParticleFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet, const double &k_scale, const reco::Jet &jetScaled);
  void fillSVFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet, const double &k_scale, const reco::Jet &jetScaled);

  float findClosestEdge(const std::vector<double>& edges, float target, int i);

  const double jet_radius_;
  const double min_jet_pt_;
  const double max_jet_eta_;
  const double min_pt_for_track_properties_;
  const bool use_puppiP4_;
  const bool include_neutrals_;
  const bool sort_by_sip2dsig_;
  const double min_puppi_wgt_;
  const bool flip_ip_sign_;
  const double max_sip3dsig_;
  const std::vector<double> mass_edges_;
  const int mass_edge_index_;
  const bool remove_lead_lepton_;

  edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;
  edm::EDGetTokenT<CandidateView> pfcand_token_;
  edm::EDGetTokenT<CandidateView> lt_token_;

  bool use_puppi_value_map_;
  bool use_pvasq_value_map_;

  edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  edm::EDGetTokenT<edm::ValueMap<int>> pvasq_value_map_token_;
  edm::EDGetTokenT<edm::Association<VertexCollection>> pvas_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<SVCollection> svs_;
  edm::Handle<CandidateView> pfcands_;
  edm::Handle<CandidateView> lts_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;
  edm::Handle<edm::ValueMap<float>> puppi_value_map_;
  edm::Handle<edm::ValueMap<int>> pvasq_value_map_;
  edm::Handle<edm::Association<VertexCollection>> pvas_;

  const static std::vector<std::string> charged_particle_features_;
  const static std::vector<std::string> neutral_particle_features_;
  const static std::vector<std::string> sv_features_;
  const reco::Vertex *pv_ = nullptr;
};

const std::vector<std::string> ParticleTransformerAK8TagInfoProducer::charged_particle_features_{
    "cpfcandlt_puppiw",        "cpfcandlt_hcalFrac",       "cpfcandlt_VTX_ass",      "cpfcandlt_lostInnerHits",
    "cpfcandlt_quality",       "cpfcandlt_charge",         "cpfcandlt_isEl",         "cpfcandlt_isMu",
    "cpfcandlt_isChargedHad",  "cpfcandlt_phirel",
    "cpfcandlt_etarel",        "cpfcandlt_deltaR",         "cpfcandlt_abseta",       "cpfcandlt_ptrel_log",
    "cpfcandlt_erel_log",      "cpfcandlt_pt_log",         "cpfcandlt_drminsv",      "cpfcandlt_drsubjet1",
    "cpfcandlt_drsubjet2",     "cpfcandlt_normchi2",       "cpfcandlt_dz",           "cpfcandlt_dzsig",
    "cpfcandlt_dxy",           "cpfcandlt_dxysig",         "cpfcandlt_dptdpt",       "cpfcandlt_detadeta",
    "cpfcandlt_dphidphi",      "cpfcandlt_dxydxy",         "cpfcandlt_dzdz",         "cpfcandlt_dxydz",
    "cpfcandlt_dphidxy",       "cpfcandlt_dlambdadz",      "cpfcandlt_btagEtaRel",   "cpfcandlt_btagPtRatio",
    "cpfcandlt_btagPParRatio", "cpfcandlt_btagSip2dVal",   "cpfcandlt_btagSip2dSig", "cpfcandlt_btagSip3dVal",
    "cpfcandlt_btagSip3dSig",  "cpfcandlt_btagJetDistVal", "cpfcandlt_mask",         "cpfcandlt_pt_log_nopuppi",
    "cpfcandlt_e_log_nopuppi", "cpfcandlt_ptrel",          "cpfcandlt_erel",
    "cpfcandlt_isLostTrack",
    "cpfcandlt_pixelBarrelLayersWithMeasurement", "cpfcandlt_pixelEndcapLayersWithMeasurement",
    "cpfcandlt_stripTECLayersWithMeasurement",    "cpfcandlt_stripTIBLayersWithMeasurement",
    "cpfcandlt_stripTIDLayersWithMeasurement",    "cpfcandlt_stripTOBLayersWithMeasurement",
    "cpfcandlt_px",            "cpfcandlt_py",             "cpfcandlt_pz",           "cpfcandlt_energy"};

const std::vector<std::string> ParticleTransformerAK8TagInfoProducer::neutral_particle_features_{
    "npfcand_puppiw",        "npfcand_hcalFrac",       "npfcand_isGamma",      "npfcand_isNeutralHad",
    "npfcand_phirel",
    "npfcand_etarel",        "npfcand_deltaR",         "npfcand_abseta",       "npfcand_ptrel_log",
    "npfcand_erel_log",      "npfcand_pt_log",         "npfcand_mask",         "npfcand_pt_log_nopuppi",
    "npfcand_e_log_nopuppi", "npfcand_ptrel",          "npfcand_erel",
    "npfcand_px",            "npfcand_py",             "npfcand_pz",           "npfcand_energy"};

const std::vector<std::string> ParticleTransformerAK8TagInfoProducer::sv_features_{
    "sv_mask", "sv_ptrel",     "sv_erel",     "sv_phirel", "sv_etarel",       "sv_deltaR",  "sv_abseta",
    "sv_mass", "sv_ptrel_log", "sv_erel_log", "sv_pt_log", "sv_pt",           "sv_ntracks", "sv_normchi2",
    "sv_dxy",  "sv_dxysig",    "sv_d3d",      "sv_d3dsig", "sv_costhetasvpv",
    "sv_px",   "sv_py",        "sv_pz",       "sv_energy"
};

ParticleTransformerAK8TagInfoProducer::ParticleTransformerAK8TagInfoProducer(const edm::ParameterSet &iConfig)
    : jet_radius_(iConfig.getParameter<double>("jet_radius")),
      min_jet_pt_(iConfig.getParameter<double>("min_jet_pt")),
      max_jet_eta_(iConfig.getParameter<double>("max_jet_eta")),
      min_pt_for_track_properties_(iConfig.getParameter<double>("min_pt_for_track_properties")),
      use_puppiP4_(iConfig.getParameter<bool>("use_puppiP4")),
      include_neutrals_(iConfig.getParameter<bool>("include_neutrals")),
      sort_by_sip2dsig_(iConfig.getParameter<bool>("sort_by_sip2dsig")),
      min_puppi_wgt_(iConfig.getParameter<double>("min_puppi_wgt")),
      flip_ip_sign_(iConfig.getParameter<bool>("flip_ip_sign")),
      max_sip3dsig_(iConfig.getParameter<double>("sip3dSigMax")),
      mass_edges_(iConfig.getParameter<std::vector<double> >("mass_edges")),
      mass_edge_index_(iConfig.getParameter<int>("mass_edge_index")),
      remove_lead_lepton_(iConfig.getParameter<bool>("remove_lead_lepton")),
      jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
      pfcand_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pf_candidates"))),
      lt_token_(consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("lost_tracks"))),
      use_puppi_value_map_(false),
      use_pvasq_value_map_(false) {
  const auto &puppi_value_map_tag = iConfig.getParameter<edm::InputTag>("puppi_value_map");
  if (!puppi_value_map_tag.label().empty()) {
    puppi_value_map_token_ = consumes<edm::ValueMap<float>>(puppi_value_map_tag);
    use_puppi_value_map_ = true;
  }

  const auto &pvas_tag = iConfig.getParameter<edm::InputTag>("vertex_associator");
  if (!pvas_tag.label().empty()) {
    pvasq_value_map_token_ = consumes<edm::ValueMap<int>>(pvas_tag);
    pvas_token_ = consumes<edm::Association<VertexCollection>>(pvas_tag);
    use_pvasq_value_map_ = true;
  }

  produces<DeepBoostedJetTagInfoCollection>();
}

ParticleTransformerAK8TagInfoProducer::~ParticleTransformerAK8TagInfoProducer() {}


float ParticleTransformerAK8TagInfoProducer::findClosestEdge(const std::vector<double>& edges, float target, int i) {

  // Find the closest edge to the left of the target
  auto lower = std::lower_bound(edges.begin(), edges.end(), target, std::less_equal<float>());
  int closestIndex = std::distance(edges.begin(), lower) - 1;

  // Adjust the index based on 'i'
  int finalIndex = closestIndex + i;

  // Check if finalIndex is within the range
  if (finalIndex >= 0 && finalIndex < int(edges.size())) {
    // std::cout << "target mass = " << target << "; i = " << i << "; returned edge = " << edges[finalIndex] << std::endl;
    return edges[finalIndex];
  } else {
    // throw std::out_of_range("Index is out of range of the edge vector.");
    return (finalIndex < 0) ? edges[0] : edges[edges.size() - 1];
  }
}


void ParticleTransformerAK8TagInfoProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // pfParticleTransformerAK8TagInfos
  edm::ParameterSetDescription desc;
  desc.add<double>("jet_radius", 0.8);
  desc.add<double>("min_jet_pt", 150);
  desc.add<double>("max_jet_eta", 99);
  desc.add<double>("min_pt_for_track_properties", -1);
  desc.add<bool>("use_puppiP4", true);
  desc.add<bool>("include_neutrals", true);
  desc.add<bool>("sort_by_sip2dsig", false);
  desc.add<double>("min_puppi_wgt", 0.01);
  desc.add<bool>("flip_ip_sign", false);
  desc.add<double>("sip3dSigMax", -1);
  desc.add<std::vector<double> >("mass_edges", {});
  desc.add<int>("mass_edge_index", 0);
  desc.add<bool>("remove_lead_lepton", false);
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
  desc.add<edm::InputTag>("pf_candidates", edm::InputTag("particleFlow"));
  desc.add<edm::InputTag>("lost_tracks", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("jets", edm::InputTag("ak8PFJetsPuppi"));
  desc.add<edm::InputTag>("puppi_value_map", edm::InputTag("puppi"));
  desc.add<edm::InputTag>("vertex_associator", edm::InputTag("primaryVertexAssociation", "original"));
  descriptions.add("pfParticleTransformerAK8TagInfos", desc);
}

void ParticleTransformerAK8TagInfoProducer::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  auto output_tag_infos = std::make_unique<DeepBoostedJetTagInfoCollection>();

  auto jets = iEvent.getHandle(jet_token_);

  iEvent.getByToken(vtx_token_, vtxs_);
  if (vtxs_->empty()) {
    // produce empty TagInfos in case no primary vertex
    iEvent.put(std::move(output_tag_infos));
    return;  // exit event
  }
  // primary vertex
  pv_ = &vtxs_->at(0);

  iEvent.getByToken(sv_token_, svs_);
  iEvent.getByToken(pfcand_token_, pfcands_);
  iEvent.getByToken(lt_token_, lts_);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);

  if (use_puppi_value_map_) {
    iEvent.getByToken(puppi_value_map_token_, puppi_value_map_);
  }

  if (use_pvasq_value_map_) {
    iEvent.getByToken(pvasq_value_map_token_, pvasq_value_map_);
    iEvent.getByToken(pvas_token_, pvas_);
  }

  for (std::size_t jet_n = 0; jet_n < jets->size(); jet_n++) {
    const auto &jet = (*jets)[jet_n];
    edm::RefToBase<reco::Jet> jet_ref(jets, jet_n);

    // create jet features
    DeepBoostedJetFeatures features;
    // declare all the feature variables (init as empty vector)
    for (const auto &name : charged_particle_features_) {
      features.add(name);
    }
    for (const auto &name : neutral_particle_features_) {
      features.add(name);
    }
    for (const auto &name : sv_features_) {
      features.add(name);
    }
    if (remove_lead_lepton_) {
      for (const auto &name : {"removed_lep_pt", "removed_lep_eta", "removed_lep_phi", "removed_lep_pid"}) {
        features.add(name);
      }
    }

    // fill values only if above pt threshold and has daughters, otherwise left
    // empty
    bool fill_vars = true;
    if (jet.pt() < min_jet_pt_ || std::abs(jet.eta()) > max_jet_eta_)
      fill_vars = false;
    if (jet.numberOfDaughters() == 0)
      fill_vars = false;

    // calculate the scaling factor for the jet
    const auto *patJet = dynamic_cast<const pat::Jet *>(&jet);
    if (!patJet) {
      throw edm::Exception(edm::errors::InvalidReference) << "Input is not a pat::Jet.";
    }
    double msd = patJet->userFloat("ak8PFJetsPuppiSoftDropMass");
    const float k_scale = (mass_edges_.empty() || msd <= 0) ? 1. : findClosestEdge(mass_edges_, msd, mass_edge_index_) / msd;

    // scaling jet p4 by iterating its daughters
    math::XYZTLorentzVector jetP4;
    for (const auto &dau : jet.daughterPtrVector()) {
      jetP4 += math::XYZTLorentzVector(dau->px() * k_scale, dau->py() * k_scale, dau->pz() * k_scale, std::sqrt(dau->massSqr() + std::pow(dau->p() * k_scale, 2)));
    }
    reco::Jet jetScaled(jet);
    jetScaled.setP4(jetP4);

    if (fill_vars) {
      fillParticleFeatures(features, jet, k_scale, jetScaled);
      fillSVFeatures(features, jet, k_scale, jetScaled);

      features.check_consistency(charged_particle_features_);
      features.check_consistency(neutral_particle_features_);
      features.check_consistency(sv_features_);
    }

    // this should always be done even if features are not filled
    output_tag_infos->emplace_back(features, jet_ref);
  }

  iEvent.put(std::move(output_tag_infos));
}

void ParticleTransformerAK8TagInfoProducer::fillParticleFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet, const double &k_scale, const reco::Jet &jetScaled) {
  // require the input to be a pat::Jet
  const auto *patJet = dynamic_cast<const pat::Jet *>(&jet);
  if (!patJet) {
    throw edm::Exception(edm::errors::InvalidReference) << "Input is not a pat::Jet.";
  }

  // do nothing if jet does not have constituents
  if (jet.numberOfDaughters() == 0)
    return;

  // some jet properties
  math::XYZVector jet_dir = jet.momentum().Unit();
  GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
  const float etasign = jet.eta() > 0 ? 1 : -1;

  TrackInfoBuilder trkinfo(track_builder_);

  std::map<reco::CandidatePtr::key_type, float> puppi_wgt_cache;
  auto puppiWgt = [&](const reco::CandidatePtr &cand) {
    const auto *pack_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));
    float wgt = 1.;
    if (pack_cand) {
      wgt = pack_cand->puppiWeight();
    } else if (reco_cand) {
      if (use_puppi_value_map_) {
        wgt = (*puppi_value_map_)[cand];
      } else {
        throw edm::Exception(edm::errors::InvalidReference) << "Puppi value map is missing";
      }
    } else {
      throw edm::Exception(edm::errors::InvalidReference) << "Cannot convert to either pat::PackedCandidate or "
                                                             "reco::PFCandidate";
    }
    puppi_wgt_cache[cand.key()] = wgt;
    return wgt;
  };

  std::vector<reco::CandidatePtr> cpfPtrs, npfPtrs;
  std::map<reco::CandidatePtr::key_type, bool> isLostTrackMap;

  for (const auto &dau : jet.daughterPtrVector()) {
    // remove particles w/ extremely low puppi weights
    // [Note] use jet daughters here to get the puppiWgt correctly
    if ((puppiWgt(dau)) < min_puppi_wgt_)
      continue;
    // from here: get the original reco/packed candidate not scaled by the puppi weight
    auto cand = pfcands_->ptrAt(dau.key());
    // charged candidate selection (for Higgs Interaction Net)
    if (!include_neutrals_ && (cand->charge() == 0 || cand->pt() < min_pt_for_track_properties_))
      continue;
    // only when computing the nagative tagger: remove charged candidates with high sip3d
    if (flip_ip_sign_ && cand->charge()) {
      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      if (trkinfo.getTrackSip3dSig() > max_sip3dsig_)
        continue;
    }
    if (cand->charge() != 0) {
      cpfPtrs.push_back(cand);
      isLostTrackMap[cand.key()] = false;
    }else {
      npfPtrs.push_back(cand);
    }
  }
  // lost tracks: fill to cpfcands
  for (size_t i = 0; i < lts_->size(); ++i) {
    auto cand = lts_->ptrAt(i);
    if (reco::deltaR(*cand, jet) < jet_radius_) {
      cpfPtrs.push_back(cand);
      isLostTrackMap[cand.key()] = true;
      puppi_wgt_cache[cand.key()] = 1.; // set puppi weight to 1 for lost tracks
    }
  }

  std::vector<btagbtvdeep::SortingClass<reco::CandidatePtr>> c_sorted;
  if (sort_by_sip2dsig_) {
    // sort charged pf candidates by 2d impact parameter significance
    for (const auto &cand : cpfPtrs) {
      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      c_sorted.emplace_back(cand,
                            trkinfo.getTrackSip2dSig(),
                            -btagbtvdeep::mindrsvpfcand(*svs_, &(*cand), jet_radius_),
                            cand->pt() / jet.pt());
      std::sort(c_sorted.begin(), c_sorted.end(), btagbtvdeep::SortingClass<reco::CandidatePtr>::compareByABCInv);
    }
    for (unsigned int i = 0; i < c_sorted.size(); i++) {
      const auto &c = c_sorted.at(i);
      const auto &cand = c.get();
      cpfPtrs.at(i) = cand;
    }
  } else {
    if (use_puppiP4_) {
      // sort by Puppi-weighted pt
      std::sort(cpfPtrs.begin(),
                cpfPtrs.end(),
                [&puppi_wgt_cache](const reco::CandidatePtr &a, const reco::CandidatePtr &b) {
                  return puppi_wgt_cache.at(a.key()) * a->pt() > puppi_wgt_cache.at(b.key()) * b->pt();
                });
      std::sort(npfPtrs.begin(),
                npfPtrs.end(),
                [&puppi_wgt_cache](const reco::CandidatePtr &a, const reco::CandidatePtr &b) {
                  return puppi_wgt_cache.at(a.key()) * a->pt() > puppi_wgt_cache.at(b.key()) * b->pt();
                });
    } else {
      // sort by original pt (not Puppi-weighted)
      std::sort(cpfPtrs.begin(), cpfPtrs.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });
      std::sort(npfPtrs.begin(), npfPtrs.end(), [](const auto &a, const auto &b) { return a->pt() > b->pt(); });
    }
  }

  // remove leptons
  if (remove_lead_lepton_) {
    if (sort_by_sip2dsig_) {
      throw edm::Exception(edm::errors::InvalidReference) << "Cannot remove leptons when sorting by sip2dsig";
    }
    // check cpfPtrs from the begining, remove the leading lepton (pid = +-11 or +-13)
    for (auto it = cpfPtrs.begin(); it != cpfPtrs.end(); ++it) {
      if (std::abs((*it)->pdgId()) == 11 || std::abs((*it)->pdgId()) == 13) {
        cpfPtrs.erase(it);
        // fill the removed lepton features
        fts.fill("removed_lep_pt", (*it)->pt());
        fts.fill("removed_lep_eta", (*it)->eta());
        fts.fill("removed_lep_phi", (*it)->phi());
        fts.fill("removed_lep_pid", (*it)->pdgId());
        break;
      }
    }
  }

  // reserve space
  for (const auto &name : charged_particle_features_) {
    fts.reserve(name, cpfPtrs.size());
  }

  auto useTrackProperties = [&](const reco::PFCandidate *reco_cand) {
    const auto *trk = reco_cand->bestTrack();
    return trk != nullptr && trk->pt() > min_pt_for_track_properties_;
  };

  for (const auto &cand : cpfPtrs) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));

    if (!include_neutrals_ &&
        ((packed_cand && !packed_cand->hasTrackDetails()) || (reco_cand && !useTrackProperties(reco_cand))))
      continue;

    const float ip_sign = flip_ip_sign_ ? -1 : 1;

    auto candP4 = use_puppiP4_ ? puppi_wgt_cache.at(cand.key()) * cand->p4() : cand->p4();
    math::XYZTLorentzVector candP4Scaled(candP4.px() * k_scale, candP4.py() * k_scale, candP4.pz() * k_scale, std::sqrt(candP4.M2() + candP4.P2() * k_scale * k_scale));
    if (packed_cand) {
      float hcal_fraction = 0.;
      if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
        hcal_fraction = packed_cand->hcalFraction();
      } else if (packed_cand->isIsolatedChargedHadron()) {
        hcal_fraction = packed_cand->rawHcalFraction();
      }

      fts.fill("cpfcandlt_hcalFrac", hcal_fraction);
      fts.fill("cpfcandlt_VTX_ass", packed_cand->pvAssociationQuality());
      fts.fill("cpfcandlt_lostInnerHits", packed_cand->lostInnerHits());
      fts.fill("cpfcandlt_quality", packed_cand->bestTrack() ? packed_cand->bestTrack()->qualityMask() : 0);

      fts.fill("cpfcandlt_charge", packed_cand->charge());
      fts.fill("cpfcandlt_isEl", std::abs(packed_cand->pdgId()) == 11);
      fts.fill("cpfcandlt_isMu", std::abs(packed_cand->pdgId()) == 13);
      fts.fill("cpfcandlt_isChargedHad", std::abs(packed_cand->pdgId()) == 211);
      fts.fill("cpfcandlt_isLostTrack", isLostTrackMap[cand.key()]);

      // impact parameters
      fts.fill("cpfcandlt_dz", ip_sign * packed_cand->dz());
      fts.fill("cpfcandlt_dxy", ip_sign * packed_cand->dxy());
      fts.fill("cpfcandlt_dzsig", packed_cand->bestTrack() ? ip_sign * packed_cand->dz() / packed_cand->dzError() : 0);
      fts.fill("cpfcandlt_dxysig", packed_cand->bestTrack() ? ip_sign * packed_cand->dxy() / packed_cand->dxyError() : 0);

    } else if (reco_cand) {
      // get vertex association quality
      int pv_ass_quality = 0;  // fallback value
      float vtx_ass = 0;
      if (use_pvasq_value_map_) {
        pv_ass_quality = (*pvasq_value_map_)[cand];
        const reco::VertexRef &PV_orig = (*pvas_)[cand];
        vtx_ass = vtx_ass_from_pfcand(*reco_cand, pv_ass_quality, PV_orig);
      } else {
        throw edm::Exception(edm::errors::InvalidReference) << "Vertex association missing";
      }

      fts.fill("cpfcandlt_hcalFrac", reco_cand->hcalEnergy() / (reco_cand->ecalEnergy() + reco_cand->hcalEnergy()));
      fts.fill("cpfcandlt_VTX_ass", vtx_ass);
      fts.fill("cpfcandlt_lostInnerHits", useTrackProperties(reco_cand) ? lost_inner_hits_from_pfcand(*reco_cand) : 0);
      fts.fill("cpfcandlt_quality", useTrackProperties(reco_cand) ? quality_from_pfcand(*reco_cand) : 0);

      fts.fill("cpfcandlt_charge", reco_cand->charge());
      fts.fill("cpfcandlt_isEl", std::abs(reco_cand->pdgId()) == 11);
      fts.fill("cpfcandlt_isMu", std::abs(reco_cand->pdgId()) == 13);
      fts.fill("cpfcandlt_isChargedHad", std::abs(reco_cand->pdgId()) == 211);
      fts.fill("cpfcandlt_isLostTrack", isLostTrackMap[cand.key()]);

      // impact parameters
      const auto *trk = reco_cand->bestTrack();
      float dz = trk ? ip_sign * trk->dz(pv_->position()) : 0;
      float dxy = trk ? ip_sign * trk->dxy(pv_->position()) : 0;
      fts.fill("cpfcandlt_dz", dz);
      fts.fill("cpfcandlt_dzsig", trk ? dz / trk->dzError() : 0);
      fts.fill("cpfcandlt_dxy", dxy);
      fts.fill("cpfcandlt_dxysig", trk ? dxy / trk->dxyError() : 0);
    }

    // basic kinematics
    fts.fill("cpfcandlt_px", candP4Scaled.px());
    fts.fill("cpfcandlt_py", candP4Scaled.py());
    fts.fill("cpfcandlt_pz", candP4Scaled.pz());
    fts.fill("cpfcandlt_energy", candP4Scaled.energy());

    fts.fill("cpfcandlt_puppiw", puppi_wgt_cache.at(cand.key()));
    fts.fill("cpfcandlt_phirel", reco::deltaPhi(candP4Scaled, jetScaled));
    fts.fill("cpfcandlt_etarel", etasign * (candP4Scaled.eta() - jetScaled.eta()));
    fts.fill("cpfcandlt_deltaR", reco::deltaR(candP4Scaled, jetScaled));
    fts.fill("cpfcandlt_abseta", std::abs(candP4Scaled.eta()));

    fts.fill("cpfcandlt_ptrel_log", std::log(candP4Scaled.pt() / jetScaled.pt()));
    fts.fill("cpfcandlt_ptrel", candP4Scaled.pt() / jetScaled.pt());
    fts.fill("cpfcandlt_erel_log", std::log(candP4Scaled.energy() / jetScaled.energy()));
    fts.fill("cpfcandlt_erel", candP4Scaled.energy() / jetScaled.energy());
    fts.fill("cpfcandlt_pt_log", std::log(candP4Scaled.pt()));

    fts.fill("cpfcandlt_mask", 1);
    fts.fill("cpfcandlt_pt_log_nopuppi", std::log(cand->pt() * k_scale));
    fts.fill("cpfcandlt_e_log_nopuppi", std::log(std::sqrt(cand->massSqr() + std::pow(cand->p() * k_scale, 2))));

    float drminpfcandsv = btagbtvdeep::mindrsvpfcand(*svs_, &(*cand), std::numeric_limits<float>::infinity());
    fts.fill("cpfcandlt_drminsv", drminpfcandsv);

    // subjets
    if (patJet->nSubjetCollections() > 0) {
      auto subjets = patJet->subjets();
      std::sort(subjets.begin(), subjets.end(), [](const edm::Ptr<pat::Jet> &p1, const edm::Ptr<pat::Jet> &p2) {
        return p1->pt() > p2->pt();
      });  // sort by pt
      fts.fill("cpfcandlt_drsubjet1", !subjets.empty() ? reco::deltaR(*cand, *subjets.at(0)) : -1);
      fts.fill("cpfcandlt_drsubjet2", subjets.size() > 1 ? reco::deltaR(*cand, *subjets.at(1)) : -1);
    } else {
      fts.fill("cpfcandlt_drsubjet1", -1);
      fts.fill("cpfcandlt_drsubjet2", -1);
    }

    const reco::Track *trk = nullptr;
    if (packed_cand) {
      trk = packed_cand->bestTrack();
    } else if (reco_cand && useTrackProperties(reco_cand)) {
      trk = reco_cand->bestTrack();
    }
    if (trk) {
      fts.fill("cpfcandlt_normchi2", std::floor(trk->normalizedChi2()));

      // track covariance
      auto cov = [&](unsigned i, unsigned j) { return trk->covariance(i, j); };
      fts.fill("cpfcandlt_dptdpt", cov(0, 0));
      fts.fill("cpfcandlt_detadeta", cov(1, 1));
      fts.fill("cpfcandlt_dphidphi", cov(2, 2));
      fts.fill("cpfcandlt_dxydxy", cov(3, 3));
      fts.fill("cpfcandlt_dzdz", cov(4, 4));
      fts.fill("cpfcandlt_dxydz", cov(3, 4));
      fts.fill("cpfcandlt_dphidxy", cov(2, 3));
      fts.fill("cpfcandlt_dlambdadz", cov(1, 4));

      trkinfo.buildTrackInfo(&(*cand), jet_dir, jet_ref_track_dir, *pv_);
      fts.fill("cpfcandlt_btagEtaRel", trkinfo.getTrackEtaRel());
      fts.fill("cpfcandlt_btagPtRatio", trkinfo.getTrackPtRatio());
      fts.fill("cpfcandlt_btagPParRatio", trkinfo.getTrackPParRatio());
      fts.fill("cpfcandlt_btagSip2dVal", ip_sign * trkinfo.getTrackSip2dVal());
      fts.fill("cpfcandlt_btagSip2dSig", ip_sign * trkinfo.getTrackSip2dSig());
      fts.fill("cpfcandlt_btagSip3dVal", ip_sign * trkinfo.getTrackSip3dVal());
      fts.fill("cpfcandlt_btagSip3dSig", ip_sign * trkinfo.getTrackSip3dSig());
      fts.fill("cpfcandlt_btagJetDistVal", trkinfo.getTrackJetDistVal());

      fts.fill("cpfcandlt_pixelBarrelLayersWithMeasurement", trk->hitPattern().pixelBarrelLayersWithMeasurement());
      fts.fill("cpfcandlt_pixelEndcapLayersWithMeasurement", trk->hitPattern().pixelEndcapLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTIBLayersWithMeasurement", trk->hitPattern().stripTIBLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTIDLayersWithMeasurement", trk->hitPattern().stripTIDLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTOBLayersWithMeasurement", trk->hitPattern().stripTOBLayersWithMeasurement());
      fts.fill("cpfcandlt_stripTECLayersWithMeasurement", trk->hitPattern().stripTECLayersWithMeasurement());

    } else {
      fts.fill("cpfcandlt_normchi2", 999);

      fts.fill("cpfcandlt_dptdpt", 0);
      fts.fill("cpfcandlt_detadeta", 0);
      fts.fill("cpfcandlt_dphidphi", 0);
      fts.fill("cpfcandlt_dxydxy", 0);
      fts.fill("cpfcandlt_dzdz", 0);
      fts.fill("cpfcandlt_dxydz", 0);
      fts.fill("cpfcandlt_dphidxy", 0);
      fts.fill("cpfcandlt_dlambdadz", 0);

      fts.fill("cpfcandlt_btagEtaRel", 0);
      fts.fill("cpfcandlt_btagPtRatio", 0);
      fts.fill("cpfcandlt_btagPParRatio", 0);
      fts.fill("cpfcandlt_btagSip2dVal", 0);
      fts.fill("cpfcandlt_btagSip2dSig", 0);
      fts.fill("cpfcandlt_btagSip3dVal", 0);
      fts.fill("cpfcandlt_btagSip3dSig", 0);
      fts.fill("cpfcandlt_btagJetDistVal", 0);

      fts.fill("cpfcandlt_pixelBarrelLayersWithMeasurement", 0); 
      fts.fill("cpfcandlt_pixelEndcapLayersWithMeasurement", 0); 
      fts.fill("cpfcandlt_stripTIBLayersWithMeasurement", 0); 
      fts.fill("cpfcandlt_stripTIDLayersWithMeasurement", 0); 
      fts.fill("cpfcandlt_stripTOBLayersWithMeasurement", 0); 
      fts.fill("cpfcandlt_stripTECLayersWithMeasurement", 0); 
    }

    // pixel hits pattern variables

  }

  // fill neutral candidate features
  for (const auto &cand : npfPtrs) {
    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));
    const auto *reco_cand = dynamic_cast<const reco::PFCandidate *>(&(*cand));

    if (!include_neutrals_ &&
        ((packed_cand && !packed_cand->hasTrackDetails()) || (reco_cand && !useTrackProperties(reco_cand))))
      continue;

    auto candP4 = use_puppiP4_ ? puppi_wgt_cache.at(cand.key()) * cand->p4() : cand->p4();
    math::XYZTLorentzVector candP4Scaled(candP4.px() * k_scale, candP4.py() * k_scale, candP4.pz() * k_scale, std::sqrt(candP4.M2() + candP4.P2() * k_scale * k_scale));
    if (packed_cand) {
      float hcal_fraction = 0.;
      if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
        hcal_fraction = packed_cand->hcalFraction();
      } else if (packed_cand->isIsolatedChargedHadron()) {
        hcal_fraction = packed_cand->rawHcalFraction();
      }

      fts.fill("npfcand_hcalFrac", hcal_fraction);

      fts.fill("npfcand_isGamma", std::abs(packed_cand->pdgId()) == 22);
      fts.fill("npfcand_isNeutralHad", std::abs(packed_cand->pdgId()) == 130);

    } else if (reco_cand) {

      fts.fill("npfcand_hcalFrac", reco_cand->hcalEnergy() / (reco_cand->ecalEnergy() + reco_cand->hcalEnergy()));

      fts.fill("npfcand_isGamma", std::abs(reco_cand->pdgId()) == 22);
      fts.fill("npfcand_isNeutralHad", std::abs(reco_cand->pdgId()) == 130);
    }

    // basic kinematics
    fts.fill("npfcand_px", candP4Scaled.px());
    fts.fill("npfcand_py", candP4Scaled.py());
    fts.fill("npfcand_pz", candP4Scaled.pz());
    fts.fill("npfcand_energy", candP4Scaled.energy());

    fts.fill("npfcand_puppiw", puppi_wgt_cache.at(cand.key()));
    fts.fill("npfcand_phirel", reco::deltaPhi(candP4Scaled, jetScaled));
    fts.fill("npfcand_etarel", etasign * (candP4Scaled.eta() - jetScaled.eta()));
    fts.fill("npfcand_deltaR", reco::deltaR(candP4Scaled, jetScaled));
    fts.fill("npfcand_abseta", std::abs(candP4Scaled.eta()));

    fts.fill("npfcand_ptrel_log", std::log(candP4Scaled.pt() / jetScaled.pt()));
    fts.fill("npfcand_ptrel", candP4Scaled.pt() / jetScaled.pt());
    fts.fill("npfcand_erel_log", std::log(candP4Scaled.energy() / jetScaled.energy()));
    fts.fill("npfcand_erel", candP4Scaled.energy() / jetScaled.energy());
    fts.fill("npfcand_pt_log", std::log(candP4Scaled.pt()));

    fts.fill("npfcand_mask", 1);
    fts.fill("npfcand_pt_log_nopuppi", std::log(cand->pt() * k_scale));
    fts.fill("npfcand_e_log_nopuppi", std::log(std::sqrt(cand->massSqr() + std::pow(cand->p() * k_scale, 2))));
  }
}

void ParticleTransformerAK8TagInfoProducer::fillSVFeatures(DeepBoostedJetFeatures &fts, const reco::Jet &jet, const double &k_scale, const reco::Jet &jetScaled) {
  std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
  for (const auto &sv : *svs_) {
    if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
      jetSVs.push_back(&sv);
    }
  }
  // sort by dxy significance
  std::sort(jetSVs.begin(),
            jetSVs.end(),
            [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
              return sv_vertex_comparator(*sva, *svb, *pv_);
            });

  // reserve space
  for (const auto &name : sv_features_) {
    fts.reserve(name, jetSVs.size());
  }

  const float etasign = jet.eta() > 0 ? 1 : -1;

  for (const auto *sv : jetSVs) {

    /*
    // scaling SV p4 by iterating each daughter p4
    math::XYZTLorentzVector svP4;
    for (unsigned int i = 0; i < sv->numberOfDaughters(); i++) {
      const auto *dau = sv->daughter(i);
      svP4 += math::XYZTLorentzVector(dau->px() * k_scale, dau->py() * k_scale, dau->pz() * k_scale, std::sqrt(dau->massSqr() + std::pow(dau->p() * k_scale, 2)));
    }
    auto svScaled = new reco::VertexCompositePtrCandidate(*sv);
    svScaled->setP4(svP4);
    */

    // scaling SV p4 by maintaining its mass
    auto svScaled = new reco::VertexCompositePtrCandidate(*sv);
    svScaled->setP4(math::XYZTLorentzVector(sv->px() * k_scale, sv->py() * k_scale, sv->pz() * k_scale, std::sqrt(sv->massSqr() + std::pow(sv->p() * k_scale, 2))));

    // basic kinematics
    fts.fill("sv_mask", 1);

    fts.fill("sv_px", svScaled->px());
    fts.fill("sv_py", svScaled->py());
    fts.fill("sv_pz", svScaled->pz());
    fts.fill("sv_energy", svScaled->energy());

    fts.fill("sv_phirel", reco::deltaPhi(*svScaled, jetScaled));
    fts.fill("sv_etarel", etasign * (svScaled->eta() - jetScaled.eta()));
    fts.fill("sv_deltaR", reco::deltaR(*svScaled, jetScaled));
    fts.fill("sv_abseta", std::abs(svScaled->eta()));
    fts.fill("sv_mass", svScaled->mass());

    fts.fill("sv_ptrel_log", std::log(svScaled->pt() / jetScaled.pt()));
    fts.fill("sv_ptrel", svScaled->pt() / jetScaled.pt());
    fts.fill("sv_erel_log", std::log(svScaled->energy() / jetScaled.energy()));
    fts.fill("sv_erel", svScaled->energy() / jetScaled.energy());
    fts.fill("sv_pt_log", std::log(svScaled->pt()));
    fts.fill("sv_pt", svScaled->pt());

    // sv properties
    fts.fill("sv_ntracks", sv->numberOfDaughters());
    fts.fill("sv_normchi2", sv->vertexNormalizedChi2());

    const auto &dxy = vertexDxy(*sv, *pv_);
    fts.fill("sv_dxy", dxy.value());
    fts.fill("sv_dxysig", dxy.significance());

    const auto &d3d = vertexD3d(*sv, *pv_);
    fts.fill("sv_d3d", d3d.value());
    fts.fill("sv_d3dsig", d3d.significance());

    fts.fill("sv_costhetasvpv", (flip_ip_sign_ ? -1.f : 1.f) * vertexDdotP(*sv, *pv_));
  }
}

// define this as a plug-in
DEFINE_FWK_MODULE(ParticleTransformerAK8TagInfoProducer);
