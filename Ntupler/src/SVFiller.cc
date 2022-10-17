/*
 * SVFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DeepNTuples/Ntupler/interface/SVFiller.h"

namespace deepntuples {

void SVFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void SVFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);
}

void SVFiller::book() {

  data.add<int>("n_sv", 0);
  data.add<float>("nsv", 0);

  // basic kinematics
  data.addMulti<float>("sv_ptrel");
  data.addMulti<float>("sv_erel");
  data.addMulti<float>("sv_phirel");
  data.addMulti<float>("sv_etarel");
  data.addMulti<float>("sv_deltaR");
  data.addMulti<float>("sv_pt");
  data.addMulti<float>("sv_abseta");
  data.addMulti<float>("sv_mass");

  data.addMulti<float>("sv_ptrel_log");
  data.addMulti<float>("sv_erel_log");
  data.addMulti<float>("sv_pt_log");
  data.addMulti<float>("sv_e_log");

  // sv properties
  data.addMulti<float>("sv_ntracks");
  data.addMulti<float>("sv_chi2");
  data.addMulti<float>("sv_ndf");
  data.addMulti<float>("sv_normchi2");
  data.addMulti<float>("sv_dxy");
  data.addMulti<float>("sv_dxyerr");
  data.addMulti<float>("sv_dxysig");
  data.addMulti<float>("sv_d3d");
  data.addMulti<float>("sv_d3derr");
  data.addMulti<float>("sv_d3dsig");
  data.addMulti<float>("sv_costhetasvpv");

  // four-vector
  data.addMulti<float>("sv_px");
  data.addMulti<float>("sv_py");
  data.addMulti<float>("sv_pz");
  data.addMulti<float>("sv_energy");

}

bool SVFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  std::vector<const reco::VertexCompositePtrCandidate*> jetSVs;
  for (const auto &sv : *SVs){
    if (reco::deltaR(sv, jet) < jetR_) {
      jetSVs.push_back(&sv);
    }
  }

  // sort by dxy significance
  const auto &pv = vertices->at(0);
  std::sort(jetSVs.begin(), jetSVs.end(), [&](const reco::VertexCompositePtrCandidate *sv1, const reco::VertexCompositePtrCandidate *sv2){
    return vertexDxy(*sv1, pv).significance() > vertexDxy(*sv2, pv).significance();
  });


  data.fill<int>("n_sv", jetSVs.size());
  data.fill<float>("nsv", jetSVs.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  for (const auto *sv : jetSVs){
    // basic kinematics
    data.fillMulti<float>("sv_ptrel", sv->pt() / jet.pt());
    data.fillMulti<float>("sv_erel", sv->energy() / jet.energy());
    data.fillMulti<float>("sv_phirel", reco::deltaPhi(*sv, jet));
    data.fillMulti<float>("sv_etarel", etasign * (sv->eta() - jet.eta()));
    data.fillMulti<float>("sv_deltaR", reco::deltaR(*sv, jet));
    data.fillMulti<float>("sv_pt", sv->pt());
    data.fillMulti<float>("sv_abseta", std::abs(sv->eta()));
    data.fillMulti<float>("sv_mass", sv->mass());

    data.fillMulti<float>("sv_ptrel_log", catchInfs(std::log(sv->pt()/jet.pt()), -99));
    data.fillMulti<float>("sv_erel_log", catchInfs(std::log(sv->energy()/jet.energy()), -99));
    data.fillMulti<float>("sv_pt_log", catchInfs(std::log(sv->pt()), -99));
    data.fillMulti<float>("sv_e_log", catchInfs(std::log(sv->energy()), -99));

    // sv properties
    data.fillMulti<float>("sv_ntracks", sv->numberOfDaughters());
    data.fillMulti<float>("sv_chi2", sv->vertexChi2());
    data.fillMulti<float>("sv_ndf", sv->vertexNdof());
    data.fillMulti<float>("sv_normchi2", catchInfs(sv->vertexNormalizedChi2()));

    const auto &dxy = vertexDxy(*sv, pv);
    data.fillMulti<float>("sv_dxy", dxy.value());
    data.fillMulti<float>("sv_dxyerr", dxy.error());
    data.fillMulti<float>("sv_dxysig", dxy.significance());

    const auto &d3d = vertexD3d(*sv, pv);
    data.fillMulti<float>("sv_d3d", d3d.value());
    data.fillMulti<float>("sv_d3derr", d3d.error());
    data.fillMulti<float>("sv_d3dsig", d3d.significance());
    data.fillMulti<float>("sv_costhetasvpv", vertexDdotP(*sv, pv));
    
    data.fillMulti<float>("sv_px", sv->px());
    data.fillMulti<float>("sv_py", sv->py());
    data.fillMulti<float>("sv_pz", sv->pz());
    data.fillMulti<float>("sv_energy", sv->energy());
  }

  return true;
}


Measurement1D SVFiller::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

Measurement1D SVFiller::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

float SVFiller::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}



} /* namespace deepntuples */
