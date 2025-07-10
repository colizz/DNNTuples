/*
 * PFCompleteFiller.cc
 *
 *  Created on: Sep 25, 2017
 *      Author: hqu
 */

#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DeepNTuples/Ntupler/interface/PFCompleteFiller.h"

namespace deepntuples {

void PFCompleteFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  transientTrackBuilderToken_ = cc.esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
  ltToken_ = cc.consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("losttracks"));
}

void PFCompleteFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
  iEvent.getByToken(ltToken_, LTs);
  builder_ = iSetup.getHandle(transientTrackBuilderToken_);
}

void PFCompleteFiller::book() {

  data.add<int>("n_cpfcands", 0);
  data.add<int>("n_npfcands", 0);
  data.add<int>("n_lts", 0);

  // ==================== charged pf candidates & lost tracks ====================
  
  // four-vector
  data.addMulti<float>("cpfcandlt_px");
  data.addMulti<float>("cpfcandlt_py");
  data.addMulti<float>("cpfcandlt_pz");
  data.addMulti<float>("cpfcandlt_energy");

  // no puppi scaled
  data.addMulti<float>("cpfcandlt_pt_nopuppi");
  data.addMulti<float>("cpfcandlt_pt_log_nopuppi");
  data.addMulti<float>("cpfcandlt_e_log_nopuppi");

  data.addMulti<float>("cpfcandlt_phirel");
  data.addMulti<float>("cpfcandlt_etarel");
  // data.addMulti<float>("cpfcandlt_deltaR");
  data.addMulti<float>("cpfcandlt_puppiw");
  data.addMulti<float>("cpfcandlt_abseta");

  data.addMulti<float>("cpfcandlt_drminsvin"); // restricted to within the jet cone

  // use uncorrected pT to order the two subjets
  data.addMulti<float>("cpfcandlt_dr_uncorrsj1");
  data.addMulti<float>("cpfcandlt_dr_uncorrsj2");

  data.addMulti<float>("cpfcandlt_charge");
  data.addMulti<float>("cpfcandlt_isMu");
  data.addMulti<float>("cpfcandlt_isEl");
  data.addMulti<float>("cpfcandlt_isChargedHad");
  data.addMulti<float>("cpfcandlt_isLostTrack");

  // for neutral
  data.addMulti<float>("cpfcandlt_hcalFrac");
  data.addMulti<float>("cpfcandlt_hcalFracCalib");

  // for charged
  data.addMulti<float>("cpfcandlt_VTX_ass");
  data.addMulti<float>("cpfcandlt_fromPV");
  data.addMulti<float>("cpfcandlt_lostInnerHits");
  data.addMulti<float>("cpfcandlt_trackHighPurity");

  // impact parameters
  data.addMulti<float>("cpfcandlt_dz");
  data.addMulti<float>("cpfcandlt_dzsig");
  data.addMulti<float>("cpfcandlt_dxy");
  data.addMulti<float>("cpfcandlt_dxysig");

  // track quality
  data.addMulti<float>("cpfcandlt_normchi2");
  data.addMulti<float>("cpfcandlt_quality");

  // track btag info
  // data.addMulti<float>("cpfcandlt_btagMomentum");
  // data.addMulti<float>("cpfcandlt_btagEta");
  data.addMulti<float>("cpfcandlt_btagEtaRel");
  data.addMulti<float>("cpfcandlt_btagPtRel");
  // data.addMulti<float>("cpfcandlt_btagPPar");
  // data.addMulti<float>("cpfcandlt_btagDeltaR");
  data.addMulti<float>("cpfcandlt_btagPtRatio");
  data.addMulti<float>("cpfcandlt_btagPParRatio");
  data.addMulti<float>("cpfcandlt_btagSip3dVal");
  data.addMulti<float>("cpfcandlt_btagSip3dSig");
  data.addMulti<float>("cpfcandlt_btagJetDistVal");
//  data.addMulti<float>("cpfcandlt_btagJetDistSig"); // always gives 0?

  // //Tracker per layer
  // //Pixel barrel 
  // data.addMulti<float>("cpfcandlt_nhitpixelBarrelLayer1");
  // data.addMulti<float>("cpfcandlt_nhitpixelBarrelLayer2");
  // data.addMulti<float>("cpfcandlt_nhitpixelBarrelLayer3");
  // data.addMulti<float>("cpfcandlt_nhitpixelBarrelLayer4");
  // //Pixel Endcap 
  // data.addMulti<float>("cpfcandlt_nhitpixelEndcapLayer1");
  // data.addMulti<float>("cpfcandlt_nhitpixelEndcapLayer2");
  // //Strip TIB
  // data.addMulti<float>("cpfcandlt_nhitstripTIBLayer1");
  // data.addMulti<float>("cpfcandlt_nhitstripTIBLayer2");
  // data.addMulti<float>("cpfcandlt_nhitstripTIBLayer3");
  // data.addMulti<float>("cpfcandlt_nhitstripTIBLayer4");
  // //Strip TID
  // data.addMulti<float>("cpfcandlt_nhitstripTIDLayer1");
  // data.addMulti<float>("cpfcandlt_nhitstripTIDLayer2");
  // data.addMulti<float>("cpfcandlt_nhitstripTIDLayer3");
  // //Strip TOB
  // data.addMulti<float>("cpfcandlt_nhitstripTOBLayer1");
  // data.addMulti<float>("cpfcandlt_nhitstripTOBLayer2");
  // data.addMulti<float>("cpfcandlt_nhitstripTOBLayer3");
  // data.addMulti<float>("cpfcandlt_nhitstripTOBLayer4");
  // data.addMulti<float>("cpfcandlt_nhitstripTOBLayer5");
  // data.addMulti<float>("cpfcandlt_nhitstripTOBLayer6");
  // //Strip TEC
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer1");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer2");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer3");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer4");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer5");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer6");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer7");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer8");
  // data.addMulti<float>("cpfcandlt_nhitstripTECLayer9");
  // //Tracker all layers together
  // //Valid hits
  // data.addMulti<float>("cpfcandlt_numberOfValidHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidTrackerHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidPixelHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidPixelBarrelHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidPixelEndcapHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidStripHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidStripTIBHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidStripTIDHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidStripTOBHits"); 
  // data.addMulti<float>("cpfcandlt_numberOfValidStripTECHits"); 
  //LayersWithMeasurement
  // data.addMulti<float>("cpfcandlt_trackerLayersWithMeasurementOld"); 
  // data.addMulti<float>("cpfcandlt_trackerLayersWithMeasurement"); 
  // data.addMulti<float>("cpfcandlt_pixelLayersWithMeasurementOld"); 
  // data.addMulti<float>("cpfcandlt_pixelLayersWithMeasurement"); 
  // data.addMulti<float>("cpfcandlt_stripLayersWithMeasurement"); 
  data.addMulti<float>("cpfcandlt_pixelBarrelLayersWithMeasurement"); 
  data.addMulti<float>("cpfcandlt_pixelEndcapLayersWithMeasurement"); 
  data.addMulti<float>("cpfcandlt_stripTIBLayersWithMeasurement"); 
  data.addMulti<float>("cpfcandlt_stripTIDLayersWithMeasurement"); 
  data.addMulti<float>("cpfcandlt_stripTOBLayersWithMeasurement"); 
  data.addMulti<float>("cpfcandlt_stripTECLayersWithMeasurement"); 
  // //Null
  // data.addMulti<float>("cpfcandlt_trackerLayersNull");
  // data.addMulti<float>("cpfcandlt_pixelLayersNull"); 
  // data.addMulti<float>("cpfcandlt_stripLayersNull"); 
  // data.addMulti<float>("cpfcandlt_pixelBarrelLayersNull");
  // data.addMulti<float>("cpfcandlt_pixelEndcapLayersNull");
  // data.addMulti<float>("cpfcandlt_stripTIBLayersNull");
  // data.addMulti<float>("cpfcandlt_stripTIDLayersNull");
  // data.addMulti<float>("cpfcandlt_stripTOBLayersNull");
  // data.addMulti<float>("cpfcandlt_stripTECLayersNull");
  
  // ==================== neutron pf candidates ====================
  data.addMulti<float>("npfcand_px");
  data.addMulti<float>("npfcand_py");
  data.addMulti<float>("npfcand_pz");
  data.addMulti<float>("npfcand_energy");

  // no puppi scaled
  data.addMulti<float>("npfcand_pt_nopuppi");
  data.addMulti<float>("npfcand_pt_log_nopuppi");
  data.addMulti<float>("npfcand_e_log_nopuppi");

  data.addMulti<float>("npfcand_phirel");
  data.addMulti<float>("npfcand_etarel");
  // data.addMulti<float>("npfcand_deltaR");
  data.addMulti<float>("npfcand_puppiw");
  data.addMulti<float>("npfcand_abseta");

  data.addMulti<float>("npfcand_drminsvin"); // restricted to within the jet cone

  // use uncorrected pT to order the two subjets
  data.addMulti<float>("npfcand_dr_uncorrsj1");
  data.addMulti<float>("npfcand_dr_uncorrsj2");

  data.addMulti<float>("npfcand_isGamma");
  data.addMulti<float>("npfcand_isNeutralHad");

  // for neutral
  data.addMulti<float>("npfcand_hcalFrac");
  data.addMulti<float>("npfcand_hcalFracCalib");

}

bool PFCompleteFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  std::vector<reco::CandidatePtr> cpfPtrs, npfPtrs;
  std::map<reco::CandidatePtr, bool> isLostTrackMap;
  const auto& pfCands = jet_helper.getJetConstituents();
  int n_cpfcands = 0, n_npfcands = 0, n_lts = 0;
  for (auto& cand : pfCands){
    if (cand->charge() != 0) {
      cpfPtrs.push_back(cand);
      isLostTrackMap[cand] = false;
      n_cpfcands++;
    }else {
      npfPtrs.push_back(cand);
      n_npfcands++;
    }
  }
  // lost tracks: fill to cpfcands
  for (size_t i = 0; i < LTs->size(); ++i) {
    auto cand = LTs->ptrAt(i);
    if (reco::deltaR(*cand, jet) < jetR_) {
      cpfPtrs.push_back(cand);
      isLostTrackMap[cand] = true;
      n_lts++;
    }
  }

  data.fill<int>("n_cpfcands", n_cpfcands);
  data.fill<int>("n_npfcands", n_npfcands);
  data.fill<int>("n_lts", n_lts);

  float etasign = jet.eta()>0 ? 1 : -1;

  // loop over all charged candidates and lost tracks
  for (const auto& cand : cpfPtrs){

    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));

    // basic kinematics, valid for both charged and neutral
    // four-vector
    data.fillMulti<float>("cpfcandlt_px", packed_cand->px());
    data.fillMulti<float>("cpfcandlt_py", packed_cand->py());
    data.fillMulti<float>("cpfcandlt_pz", packed_cand->pz());
    data.fillMulti<float>("cpfcandlt_energy", packed_cand->energy());

    // not puppi weighted
    data.fillMulti<float>("cpfcandlt_pt_nopuppi", packed_cand->pt());
    data.fillMulti<float>("cpfcandlt_pt_log_nopuppi", catchInfs(std::log(packed_cand->pt()), -99));
    data.fillMulti<float>("cpfcandlt_e_log_nopuppi", catchInfs(std::log(packed_cand->energy()), -99));

    data.fillMulti<float>("cpfcandlt_phirel", reco::deltaPhi(*packed_cand, jet));
    data.fillMulti<float>("cpfcandlt_etarel", etasign * (packed_cand->eta() - jet.eta()));
    // data.fillMulti<float>("cpfcandlt_deltaR", reco::deltaR(*packed_cand, jet));
    data.fillMulti<float>("cpfcandlt_abseta", std::abs(packed_cand->eta()));

    data.fillMulti<float>("cpfcandlt_puppiw", !isLostTrackMap[cand] ? jet_helper.getPuppiWeight(cand) : 0);

    double minDRin = 2.*jetR_;
    for (const auto &sv : *SVs){
      double dr = reco::deltaR(*packed_cand, sv);
      if (dr < minDRin && reco::deltaR(jet, sv) < jetR_) minDRin = dr;
    }
    data.fillMulti<float>("cpfcandlt_drminsvin", minDRin);

    // use uncorrected pT to order the two subjets
    {
      const auto& subjets = jet_helper.getUncorrSubJets();
      data.fillMulti<float>("cpfcandlt_dr_uncorrsj1", subjets.size()>0 ? reco::deltaR(*packed_cand, *subjets.at(0)) : -1);
      data.fillMulti<float>("cpfcandlt_dr_uncorrsj2", subjets.size()>1 ? reco::deltaR(*packed_cand, *subjets.at(1)) : -1);
    }

    data.fillMulti<float>("cpfcandlt_charge", packed_cand->charge());
    data.fillMulti<float>("cpfcandlt_isEl", std::abs(packed_cand->pdgId())==11);
    data.fillMulti<float>("cpfcandlt_isMu", std::abs(packed_cand->pdgId())==13);
    data.fillMulti<float>("cpfcandlt_isChargedHad", std::abs(packed_cand->pdgId())==211);
    data.fillMulti<float>("cpfcandlt_isLostTrack", isLostTrackMap[cand]);

    // for neutral
    float hcal_fraction = 0.;
    if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
      hcal_fraction = packed_cand->hcalFraction();
    } else if (packed_cand->isIsolatedChargedHadron()) {
      hcal_fraction = packed_cand->rawHcalFraction();
    }
    data.fillMulti<float>("cpfcandlt_hcalFrac", hcal_fraction);
    data.fillMulti<float>("cpfcandlt_hcalFracCalib", packed_cand->hcalFraction());

    // for charged
    data.fillMulti<float>("cpfcandlt_VTX_ass", packed_cand->pvAssociationQuality());
    data.fillMulti<float>("cpfcandlt_fromPV", packed_cand->fromPV());
    data.fillMulti<float>("cpfcandlt_lostInnerHits", packed_cand->lostInnerHits());
    data.fillMulti<float>("cpfcandlt_trackHighPurity", packed_cand->trackHighPurity());

    // impact parameters
    data.fillMulti<float>("cpfcandlt_dz", catchInfs(packed_cand->dz()));
    data.fillMulti<float>("cpfcandlt_dzsig", packed_cand->bestTrack() ? catchInfs(packed_cand->dz()/packed_cand->dzError()) : 0);
    data.fillMulti<float>("cpfcandlt_dxy", catchInfs(packed_cand->dxy()));
    data.fillMulti<float>("cpfcandlt_dxysig", packed_cand->bestTrack() ? catchInfs(packed_cand->dxy()/packed_cand->dxyError()) : 0);

    if (packed_cand->bestTrack()){
      const auto *trk = packed_cand->bestTrack();
      data.fillMulti<float>("cpfcandlt_normchi2", catchInfs(trk->normalizedChi2()));
      data.fillMulti<float>("cpfcandlt_quality", trk->qualityMask());
    }else{
      data.fillMulti<float>("cpfcandlt_normchi2", 999);
      data.fillMulti<float>("cpfcandlt_quality", 0);
    }

    // build track info map
    TrackInfoBuilder trkinfo;
    trkinfo.buildTrackInfo(builder_, *packed_cand, jet, vertices->at(0));

    // data.fillMulti<float>("cpfcandlt_btagMomentum", catchInfs(trkinfo.getTrackMomentum()));
    // data.fillMulti<float>("cpfcandlt_btagEta", catchInfs(trkinfo.getTrackEta()));
    data.fillMulti<float>("cpfcandlt_btagEtaRel", catchInfs(trkinfo.getTrackEtaRel()));
    data.fillMulti<float>("cpfcandlt_btagPtRel", catchInfs(trkinfo.getTrackPtRel()));
    // data.fillMulti<float>("cpfcandlt_btagPPar", catchInfs(trkinfo.getTrackPPar()));
    // data.fillMulti<float>("cpfcandlt_btagDeltaR", catchInfs(trkinfo.getTrackDeltaR()));
    data.fillMulti<float>("cpfcandlt_btagPtRatio", catchInfs(trkinfo.getTrackPtRatio()));
    data.fillMulti<float>("cpfcandlt_btagPParRatio", catchInfs(trkinfo.getTrackPParRatio()));
    data.fillMulti<float>("cpfcandlt_btagSip3dVal", catchInfs(trkinfo.getTrackSip3dVal()));
    data.fillMulti<float>("cpfcandlt_btagSip3dSig", catchInfs(trkinfo.getTrackSip3dSig()));
    data.fillMulti<float>("cpfcandlt_btagJetDistVal", catchInfs(trkinfo.getTrackJetDistVal()));
//    data.fillMulti<float>("cpfcandlt_btagJetDistSig", catchInfs(trkinfo.getTrackJetDistSig()));

    // pixel hit patterns
    const auto *trk = packed_cand->bestTrack();
    //get hit pattern information
    if(trk){
      const reco::HitPattern &p = trk->hitPattern();
      //Tracker per layer
      //Pixel barrel 
      int cpfcandlt_nhitpixelBarrelLayer1 = 0;
      int cpfcandlt_nhitpixelBarrelLayer2 = 0;
      int cpfcandlt_nhitpixelBarrelLayer3 = 0;
      int cpfcandlt_nhitpixelBarrelLayer4 = 0;
      //Pixel Endcap 
      int cpfcandlt_nhitpixelEndcapLayer1 = 0;
      int cpfcandlt_nhitpixelEndcapLayer2 = 0;
      //Strip TIB
      int cpfcandlt_nhitstripTIBLayer1 = 0;
      int cpfcandlt_nhitstripTIBLayer2 = 0;
      int cpfcandlt_nhitstripTIBLayer3 = 0;
      int cpfcandlt_nhitstripTIBLayer4 = 0;
      //Strip TID
      int cpfcandlt_nhitstripTIDLayer1 = 0;
      int cpfcandlt_nhitstripTIDLayer2 = 0;
      int cpfcandlt_nhitstripTIDLayer3 = 0;
      //Strip TOB
      int cpfcandlt_nhitstripTOBLayer1 = 0;
      int cpfcandlt_nhitstripTOBLayer2 = 0;
      int cpfcandlt_nhitstripTOBLayer3 = 0;
      int cpfcandlt_nhitstripTOBLayer4 = 0;
      int cpfcandlt_nhitstripTOBLayer5 = 0;
      int cpfcandlt_nhitstripTOBLayer6 = 0;
      //Strip TEC
      int cpfcandlt_nhitstripTECLayer1 = 0;
      int cpfcandlt_nhitstripTECLayer2 = 0;
      int cpfcandlt_nhitstripTECLayer3 = 0;
      int cpfcandlt_nhitstripTECLayer4 = 0;
      int cpfcandlt_nhitstripTECLayer5 = 0;
      int cpfcandlt_nhitstripTECLayer6 = 0;
      int cpfcandlt_nhitstripTECLayer7 = 0;
      int cpfcandlt_nhitstripTECLayer8 = 0;
      int cpfcandlt_nhitstripTECLayer9 = 0;
      // loop over the hits of the track.
      //const static unsigned short LayerOffset = 3;
      //const static unsigned short LayerMask = 0xF;
      for(int nh = 0; nh < p.numberOfAllHits(reco::HitPattern::TRACK_HITS); nh++){
        uint32_t hit = p.getHitPattern(reco::HitPattern::TRACK_HITS, nh);
        if(p.validHitFilter(hit)){// if the hit is valid
          //Pixel Barrel // it is in pixel barrel
          if(p.pixelBarrelHitFilter(hit)){
            //std::cout << "valid hit found in pixel Barrel layer " << p.getLayer(hit) << std::endl;
            //if(p.getLayer(hit)==1){
            // std::cout<< (hit >> LayerOffset) << " " << ((hit >> LayerOffset) & LayerMask) << std::endl;
            //}  
            if(p.getLayer(hit)==1) cpfcandlt_nhitpixelBarrelLayer1++;
            if(p.getLayer(hit)==2) cpfcandlt_nhitpixelBarrelLayer2++;
            if(p.getLayer(hit)==3) cpfcandlt_nhitpixelBarrelLayer3++;
            if(p.getLayer(hit)==4) cpfcandlt_nhitpixelBarrelLayer4++;
          } 
          //Pixel Endcap
          if(p.pixelEndcapHitFilter(hit)){
            //std::cout << "valid hit found in pixel Endcap layer " << p.getLayer(hit) << std::endl;
            if(p.getLayer(hit)==1) cpfcandlt_nhitpixelEndcapLayer1++;
            if(p.getLayer(hit)==2) cpfcandlt_nhitpixelEndcapLayer2++;
          } 
          //Strip TIB
          if(p.stripTIBHitFilter(hit)){
            //std::cout << "valid hit found in TIB layer " << p.getLayer(hit) << std::endl;
            if(p.getLayer(hit)==1) cpfcandlt_nhitstripTIBLayer1++;
            if(p.getLayer(hit)==2) cpfcandlt_nhitstripTIBLayer2++;
            if(p.getLayer(hit)==3) cpfcandlt_nhitstripTIBLayer3++;
            if(p.getLayer(hit)==4) cpfcandlt_nhitstripTIBLayer4++;
          } 
          //Strip TID
          if(p.stripTIDHitFilter(hit)){
            //std::cout << "valid hit found in TID layer " << p.getLayer(hit) << std::endl;
            if(p.getLayer(hit)==1) cpfcandlt_nhitstripTIDLayer1++;
            if(p.getLayer(hit)==2) cpfcandlt_nhitstripTIDLayer2++;
            if(p.getLayer(hit)==3) cpfcandlt_nhitstripTIDLayer3++;
          } 
          //Strip TOB
          if(p.stripTOBHitFilter(hit)){
            //std::cout << "valid hit found in TOB layer " << p.getLayer(hit) << std::endl;
            if(p.getLayer(hit)==1) cpfcandlt_nhitstripTOBLayer1++;
            if(p.getLayer(hit)==2) cpfcandlt_nhitstripTOBLayer2++;
            if(p.getLayer(hit)==3) cpfcandlt_nhitstripTOBLayer3++;
            if(p.getLayer(hit)==4) cpfcandlt_nhitstripTOBLayer4++;
            if(p.getLayer(hit)==5) cpfcandlt_nhitstripTOBLayer5++;
            if(p.getLayer(hit)==6) cpfcandlt_nhitstripTOBLayer6++;
          } 
          //Strip TEC
          if(p.stripTECHitFilter(hit)){
            //std::cout << "valid hit found in TEC layer " << p.getLayer(hit) << std::endl;
            if(p.getLayer(hit)==1) cpfcandlt_nhitstripTECLayer1++;
            if(p.getLayer(hit)==2) cpfcandlt_nhitstripTECLayer2++;
            if(p.getLayer(hit)==3) cpfcandlt_nhitstripTECLayer3++;
            if(p.getLayer(hit)==4) cpfcandlt_nhitstripTECLayer4++;
            if(p.getLayer(hit)==5) cpfcandlt_nhitstripTECLayer5++;
            if(p.getLayer(hit)==6) cpfcandlt_nhitstripTECLayer6++;
            if(p.getLayer(hit)==7) cpfcandlt_nhitstripTECLayer7++;
            if(p.getLayer(hit)==8) cpfcandlt_nhitstripTECLayer8++;
            if(p.getLayer(hit)==9) cpfcandlt_nhitstripTECLayer9++;
          } 
        }
        //// expert level: printout the hit in 11-bit binary format
        //std::cout << "hit in 11-bit binary format = ";
        //for (int j = 10; j >= 0; j--){
        //    int bit = (hit >> j) & 0x1;
        //    std::cout << bit;
        //}
        //std::cout << std::endl;
      }
      // //Pixel Barrel 
      // data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer1", (trk->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfs(cpfcandlt_nhitpixelBarrelLayer1, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer2", (trk->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfs(cpfcandlt_nhitpixelBarrelLayer2, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer3", (trk->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfs(cpfcandlt_nhitpixelBarrelLayer3, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer4", (trk->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfs(cpfcandlt_nhitpixelBarrelLayer4, 0) : 0);
      // //Pixel Endcap 
      // data.fillMulti<float>("cpfcandlt_nhitpixelEndcapLayer1", (trk->hitPattern().numberOfValidPixelEndcapHits()) ? catchInfs(cpfcandlt_nhitpixelEndcapLayer1, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitpixelEndcapLayer2", (trk->hitPattern().numberOfValidPixelEndcapHits()) ? catchInfs(cpfcandlt_nhitpixelEndcapLayer2, 0) : 0);
      // //Strip TIB
      // data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer1", (trk->hitPattern().numberOfValidStripTIBHits()) ? catchInfs(cpfcandlt_nhitstripTIBLayer1, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer2", (trk->hitPattern().numberOfValidStripTIBHits()) ? catchInfs(cpfcandlt_nhitstripTIBLayer2, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer3", (trk->hitPattern().numberOfValidStripTIBHits()) ? catchInfs(cpfcandlt_nhitstripTIBLayer3, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer4", (trk->hitPattern().numberOfValidStripTIBHits()) ? catchInfs(cpfcandlt_nhitstripTIBLayer4, 0) : 0);
      // //Strip TID
      // data.fillMulti<float>("cpfcandlt_nhitstripTIDLayer1", (trk->hitPattern().numberOfValidStripTIDHits()) ? catchInfs(cpfcandlt_nhitstripTIDLayer1, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTIDLayer2", (trk->hitPattern().numberOfValidStripTIDHits()) ? catchInfs(cpfcandlt_nhitstripTIDLayer2, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTIDLayer3", (trk->hitPattern().numberOfValidStripTIDHits()) ? catchInfs(cpfcandlt_nhitstripTIDLayer3, 0) : 0);
      // //Strip TOB
      // data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer1", (trk->hitPattern().numberOfValidStripTOBHits()) ? catchInfs(cpfcandlt_nhitstripTOBLayer1, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer2", (trk->hitPattern().numberOfValidStripTOBHits()) ? catchInfs(cpfcandlt_nhitstripTOBLayer2, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer3", (trk->hitPattern().numberOfValidStripTOBHits()) ? catchInfs(cpfcandlt_nhitstripTOBLayer3, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer4", (trk->hitPattern().numberOfValidStripTOBHits()) ? catchInfs(cpfcandlt_nhitstripTOBLayer4, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer5", (trk->hitPattern().numberOfValidStripTOBHits()) ? catchInfs(cpfcandlt_nhitstripTOBLayer5, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer6", (trk->hitPattern().numberOfValidStripTOBHits()) ? catchInfs(cpfcandlt_nhitstripTOBLayer6, 0) : 0);
      // //Strip TEC
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer1", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer1, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer2", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer2, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer3", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer3, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer4", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer4, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer5", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer5, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer6", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer6, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer7", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer7, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer8", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer8, 0) : 0);
      // data.fillMulti<float>("cpfcandlt_nhitstripTECLayer9", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(cpfcandlt_nhitstripTECLayer9, 0) : 0);
      // //Tracker all layers together   
      // //Valid hits
      // data.fillMulti<float>("cpfcandlt_numberOfValidHits", (trk->hitPattern().numberOfValidHits()) ? catchInfs(trk->hitPattern().numberOfValidHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidTrackerHits", (trk->hitPattern().numberOfValidTrackerHits()) ? catchInfs(trk->hitPattern().numberOfValidTrackerHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidPixelHits", (trk->hitPattern().numberOfValidPixelHits()) ? catchInfs(trk->hitPattern().numberOfValidPixelHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidPixelBarrelHits", (trk->hitPattern().numberOfValidPixelBarrelHits()) ? catchInfs(trk->hitPattern().numberOfValidPixelBarrelHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidPixelEndcapHits", (trk->hitPattern().numberOfValidPixelEndcapHits()) ? catchInfs(trk->hitPattern().numberOfValidPixelEndcapHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripHits", (trk->hitPattern().numberOfValidStripHits()) ? catchInfs(trk->hitPattern().numberOfValidStripHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTIBHits", (trk->hitPattern().numberOfValidStripTIBHits()) ? catchInfs(trk->hitPattern().numberOfValidStripTIBHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTIDHits", (trk->hitPattern().numberOfValidStripTIDHits()) ? catchInfs(trk->hitPattern().numberOfValidStripTIDHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTOBHits", (trk->hitPattern().numberOfValidStripTOBHits()) ? catchInfs(trk->hitPattern().numberOfValidStripTOBHits(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTECHits", (trk->hitPattern().numberOfValidStripTECHits()) ? catchInfs(trk->hitPattern().numberOfValidStripTECHits(), 0) : 0);
      //LayersWithMeasurement
      // data.fillMulti<float>("cpfcandlt_trackerLayersWithMeasurementOld", (trk->hitPattern().trackerLayersWithMeasurementOld()) ? catchInfs(trk->hitPattern().trackerLayersWithMeasurementOld(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_trackerLayersWithMeasurement", (trk->hitPattern().trackerLayersWithMeasurement()) ? catchInfs(trk->hitPattern().trackerLayersWithMeasurement(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_pixelLayersWithMeasurementOld", (trk->hitPattern().pixelLayersWithMeasurementOld()) ? catchInfs(trk->hitPattern().pixelLayersWithMeasurementOld(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_pixelLayersWithMeasurement", (trk->hitPattern().pixelLayersWithMeasurement()) ? catchInfs(trk->hitPattern().pixelLayersWithMeasurement(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_stripLayersWithMeasurement", (trk->hitPattern().stripLayersWithMeasurement()) ? catchInfs(trk->hitPattern().stripLayersWithMeasurement(), 0) : 0);
      data.fillMulti<float>("cpfcandlt_pixelBarrelLayersWithMeasurement", (trk->hitPattern().pixelBarrelLayersWithMeasurement()) ? catchInfs(trk->hitPattern().pixelBarrelLayersWithMeasurement(), 0) : 0);
      data.fillMulti<float>("cpfcandlt_pixelEndcapLayersWithMeasurement", (trk->hitPattern().pixelEndcapLayersWithMeasurement()) ? catchInfs(trk->hitPattern().pixelEndcapLayersWithMeasurement(), 0) : 0);
      data.fillMulti<float>("cpfcandlt_stripTIBLayersWithMeasurement", (trk->hitPattern().stripTIBLayersWithMeasurement()) ? catchInfs(trk->hitPattern().stripTIBLayersWithMeasurement(), 0) : 0);
      data.fillMulti<float>("cpfcandlt_stripTIDLayersWithMeasurement", (trk->hitPattern().stripTIDLayersWithMeasurement()) ? catchInfs(trk->hitPattern().stripTIDLayersWithMeasurement(), 0) : 0);
      data.fillMulti<float>("cpfcandlt_stripTOBLayersWithMeasurement", (trk->hitPattern().stripTOBLayersWithMeasurement()) ? catchInfs(trk->hitPattern().stripTOBLayersWithMeasurement(), 0) : 0);
      data.fillMulti<float>("cpfcandlt_stripTECLayersWithMeasurement", (trk->hitPattern().stripTECLayersWithMeasurement()) ? catchInfs(trk->hitPattern().stripTECLayersWithMeasurement(), 0) : 0);
      // //Null
      // data.fillMulti<float>("cpfcandlt_trackerLayersNull", (trk->hitPattern().trackerLayersNull()) ? catchInfs(trk->hitPattern().trackerLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_pixelLayersNull", (trk->hitPattern().pixelLayersNull()) ? catchInfs(trk->hitPattern().pixelLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_stripLayersNull", (trk->hitPattern().stripLayersNull()) ? catchInfs(trk->hitPattern().stripLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_pixelBarrelLayersNull", (trk->hitPattern().pixelBarrelLayersNull()) ? catchInfs(trk->hitPattern().pixelBarrelLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_pixelEndcapLayersNull", (trk->hitPattern().pixelEndcapLayersNull()) ? catchInfs(trk->hitPattern().pixelEndcapLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_stripTIBLayersNull", (trk->hitPattern().stripTIBLayersNull()) ? catchInfs(trk->hitPattern().stripTIBLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_stripTIDLayersNull", (trk->hitPattern().stripTIDLayersNull()) ? catchInfs(trk->hitPattern().stripTIDLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_stripTOBLayersNull", (trk->hitPattern().stripTOBLayersNull()) ? catchInfs(trk->hitPattern().stripTOBLayersNull(), 0) : 0);
      // data.fillMulti<float>("cpfcandlt_stripTECLayersNull", (trk->hitPattern().stripTECLayersNull()) ? catchInfs(trk->hitPattern().stripTECLayersNull(), 0) : 0);
    }else{
      //Tracker per layer
      //Pixel barrel 
      //data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer1", 0);
      //data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer2", 0);
      //data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer3", 0);
      //data.fillMulti<float>("cpfcandlt_nhitpixelBarrelLayer4", 0);
      //Pixel Endcap 
      //data.fillMulti<float>("cpfcandlt_nhitpixelEndcapLayer1", 0);
      //data.fillMulti<float>("cpfcandlt_nhitpixelEndcapLayer2", 0);
      //Strip TIB
      //data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer1", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer2", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer3", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTIBLayer4", 0);
      //Strip TID
      //data.fillMulti<float>("cpfcandlt_nhitstripTIDLayer1", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTIDLayer2", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTIDLayer3", 0);
      //Strip TOB
      //data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer1", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer2", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer3", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer4", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer5", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTOBLayer6", 0);
      //Strip TEC
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer1", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer2", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer3", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer4", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer5", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer6", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer7", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer8", 0);
      //data.fillMulti<float>("cpfcandlt_nhitstripTECLayer9", 0);
      //Tracker all layers together
      // //Valid hits
      // data.fillMulti<float>("cpfcandlt_numberOfValidHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidTrackerHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidPixelHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidPixelBarrelHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidPixelEndcapHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTIBHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTIDHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTOBHits", 0); 
      // data.fillMulti<float>("cpfcandlt_numberOfValidStripTECHits", 0); 
      //LayersWithMeasurement
      // data.fillMulti<float>("cpfcandlt_trackerLayersWithMeasurementOld", 0); 
      // data.fillMulti<float>("cpfcandlt_trackerLayersWithMeasurement", 0); 
      // data.fillMulti<float>("cpfcandlt_pixelLayersWithMeasurementOld", 0); 
      // data.fillMulti<float>("cpfcandlt_pixelLayersWithMeasurement", 0); 
      // data.fillMulti<float>("cpfcandlt_stripLayersWithMeasurement", 0); 
      data.fillMulti<float>("cpfcandlt_pixelBarrelLayersWithMeasurement", 0); 
      data.fillMulti<float>("cpfcandlt_pixelEndcapLayersWithMeasurement", 0); 
      data.fillMulti<float>("cpfcandlt_stripTIBLayersWithMeasurement", 0); 
      data.fillMulti<float>("cpfcandlt_stripTIDLayersWithMeasurement", 0); 
      data.fillMulti<float>("cpfcandlt_stripTOBLayersWithMeasurement", 0); 
      data.fillMulti<float>("cpfcandlt_stripTECLayersWithMeasurement", 0); 
      // //Null
      // data.fillMulti<float>("cpfcandlt_trackerLayersNull", 0);
      // data.fillMulti<float>("cpfcandlt_pixelLayersNull", 0); 
      // data.fillMulti<float>("cpfcandlt_stripLayersNull", 0); 
      // data.fillMulti<float>("cpfcandlt_pixelBarrelLayersNull", 0);
      // data.fillMulti<float>("cpfcandlt_pixelEndcapLayersNull", 0);
      // data.fillMulti<float>("cpfcandlt_stripTIBLayersNull", 0);
      // data.fillMulti<float>("cpfcandlt_stripTIDLayersNull", 0);
      // data.fillMulti<float>("cpfcandlt_stripTOBLayersNull", 0);
      // data.fillMulti<float>("cpfcandlt_stripTECLayersNull", 0);
    }
  }

  // neutral particles
  for (const auto& cand : npfPtrs){

    const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*cand));

    // basic kinematics, valid for neutral particles
    // four-vector
    data.fillMulti<float>("npfcand_px", packed_cand->px());
    data.fillMulti<float>("npfcand_py", packed_cand->py());
    data.fillMulti<float>("npfcand_pz", packed_cand->pz());
    data.fillMulti<float>("npfcand_energy", packed_cand->energy());

    // not puppi weighted
    data.fillMulti<float>("npfcand_pt_nopuppi", packed_cand->pt());
    data.fillMulti<float>("npfcand_pt_log_nopuppi", catchInfs(std::log(packed_cand->pt()), -99));
    data.fillMulti<float>("npfcand_e_log_nopuppi", catchInfs(std::log(packed_cand->energy()), -99));

    data.fillMulti<float>("npfcand_phirel", reco::deltaPhi(*packed_cand, jet));
    data.fillMulti<float>("npfcand_etarel", etasign * (packed_cand->eta() - jet.eta()));
    // data.fillMulti<float>("npfcand_deltaR", reco::deltaR(*packed_cand, jet));
    data.fillMulti<float>("npfcand_abseta", std::abs(packed_cand->eta()));

    data.fillMulti<float>("npfcand_puppiw", jet_helper.getPuppiWeight(cand));

    double minDRin = 2.*jetR_;
    for (const auto &sv : *SVs){
      double dr = reco::deltaR(*packed_cand, sv);
      if (dr < minDRin && reco::deltaR(jet, sv) < jetR_) minDRin = dr;
    }
    data.fillMulti<float>("npfcand_drminsvin", minDRin);

    // use uncorrected pT to order the two subjets
    {
      const auto& subjets = jet_helper.getUncorrSubJets();
      data.fillMulti<float>("npfcand_dr_uncorrsj1", subjets.size()>0 ? reco::deltaR(*packed_cand, *subjets.at(0)) : -1);
      data.fillMulti<float>("npfcand_dr_uncorrsj2", subjets.size()>1 ? reco::deltaR(*packed_cand, *subjets.at(1)) : -1);
    }

    data.fillMulti<float>("npfcand_isGamma", std::abs(packed_cand->pdgId())==22);
    data.fillMulti<float>("npfcand_isNeutralHad", std::abs(packed_cand->pdgId())==130);

    // for neutral
    float hcal_fraction = 0.;
    if (packed_cand->pdgId() == 1 || packed_cand->pdgId() == 130) {
      hcal_fraction = packed_cand->hcalFraction();
    } else if (packed_cand->isIsolatedChargedHadron()) {
      hcal_fraction = packed_cand->rawHcalFraction();
    }
    data.fillMulti<float>("npfcand_hcalFrac", hcal_fraction);
    data.fillMulti<float>("npfcand_hcalFracCalib", packed_cand->hcalFraction());

  }


  return true;
}

} /* namespace deepntuples */
