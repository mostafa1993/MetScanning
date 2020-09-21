#include "BERegion.h"

#include <TVector2.h>


BERegion::BERegion(TTreeReader &reader, const std::string berCat)
  : Tree(reader),
    berCat_{berCat} {
  cfg_ = YAML::LoadFile("config/BER.yaml");
  minPfMet_ = Get<float>("min_pfmet");

  if (berCat == "BadChCand") {
    isBER = &BERegion::isBadChCandBER;
    minDPhi_ = Get<float>("min_dphi");
    pfCanMinPt_ = Get<float>("min_pfcandid_pt");
    pfCanPdgId_ = Get<int>("abs_pfcand_pdgid");
  }

  else if (berCat == "GlbSTightHalo") {
    isBER = &BERegion::isGlbSTHaloBER;
    maxPfMetPhi_ = Get<float>("max_abs_pfmet_phi");
    minPfMetPhi_ = Get<float>("min_abs_pfmet_phi");
    minJetPt_ = Get<float>("min_jet_pt");
    maxJetEta_ = Get<float>("max_abs_jet_eta");
    maxJetChef_ = Get<float>("max_jet_chef");
    jetPtGt25Num_ = Get<int>("num_jets_ptgt25");
    jetPtGt200Num_ = Get<int>("num_jets_ptgt200");
  }

  else if (berCat == "BadPFMuon") {
    isBER = &BERegion::isBadPFMuonBER;
    minDPhi_ = Get<float>("min_dphi");
    minMuPt_ = Get<float>("min_muon_pt");
    maxAddMuPt_ = Get<float>("max_additional_muon_pt");
    minDeltaR_ = Get<float>("min_delta_r");
    jetPtGt25Num_ = Get<int>("num_jets_ptgt25");
    jetPtGt200Num_ = Get<int>("num_jets_ptgt200");
    minMuNum_ = Get<int>("min_num_muon");
  }

  else if (berCat == "EcalBadCal_U") {
    isBER = &BERegion::isEcalBadCalBER;
    minDPhi_ = Get<float>("min_dphi");
    minJetPt_ = Get<float>("min_jet_pt");
    maxJetEta_ = Get<float>("max_abs_jet_eta");
    minJetEta_ = Get<float>("min_abs_jet_eta");
    minJetNeef_ = Get<float>("min_jet_neef");
    jetPtGt25Num_ = Get<int>("num_jets_ptgt25");
    jetPtGt200Num_ = Get<int>("num_jets_ptgt200");
  }

  else if (berCat == "HBHENoise") {
    isBER = &BERegion::isHBHENoiseBER;
    minDPhi_ = Get<float>("min_dphi");
    minJetPt_ = Get<float>("min_jet_pt");
    maxJetEta_ = Get<float>("max_abs_jet_eta");
    minJetNhef_ = Get<float>("min_jet_nhef");
    jetPtGt25Num_ = Get<int>("num_jets_ptgt25");
    jetPtGt200Num_ = Get<int>("num_jets_ptgt200");
  }

  else
    isBER = &BERegion::isOthersBER;

}


bool BERegion::operator()() {
  return (this->*isBER)();
}


bool BERegion::isBadChCandBER() {
  int i = 0;
  if (*pfmet > minPfMet_) {
    for (long unsigned int j = 0; j < (*pfCandPt).size(); j++)
      if ((*pfCandPt)[j] > pfCanMinPt_ and
          std::abs((*pfCandPdgId)[j]) == pfCanPdgId_ and
          std::fabs(
            TVector2::Phi_mpi_pi(*pfmetPhi-(*pfCandPhi)[j])) > minDPhi_)
        i++;
  }

  return (i >= 1);
}


bool BERegion::isGlbSTHaloBER() {
  if (*pfmet > minPfMet_ and
      (*jetPt).size() == 1 and
      (*numJetsPt25) == 1 and
      (std::fabs(TVector2::Phi_mpi_pi(*pfmetPhi)) < maxPfMetPhi_ or
       std::fabs(TVector2::Phi_mpi_pi(*pfmetPhi)) > minPfMetPhi_) and
      (*jetPt)[0] > minJetPt_ and
      std::fabs((*jetEta)[0]) < maxJetEta_ and
      (*jetCHEF)[0] < maxJetChef_)
    return true;

  return false;
}


bool BERegion::isBadPFMuonBER() {
  if (*pfmet > minPfMet_ and
      (*jetPt).size() == jetPtGt200Num_ and
      (*numJetsPt25) == jetPtGt25Num_ and
      (*muonPt).size() >= minMuNum_ and
      (*muonPt)[0] > minMuPt_ and
      std::fabs(
        TVector2::Phi_mpi_pi(*pfmetPhi - (*muonPhi)[0])) >= minDPhi_ and
      (((*muonPt).size() > minMuNum_) ? ((*muonPt)[1] < maxAddMuPt_) : true)
      ) {
    int k = 0;
    for (long unsigned int j = 0; j < (*jetPt).size(); j++) {
      // To be sure that the jets are not comming from muon jet cone
      // by requiring DR < 0.4
      float const deltaR = std::sqrt(
          std::pow((*muonEta)[0] - (*jetEta)[j], 2) +
          std::pow((*muonPhi)[0] - (*jetPhi)[j], 2));
      if (deltaR < minDeltaR_)
        continue;
      k++;
      break;
    }

    return (k == 0);
  }

  return false;
}


bool BERegion::isEcalBadCalBER() {
  if (*pfmet > minPfMet_ and
      (*jetPt).size() == jetPtGt200Num_ and
      (*numJetsPt25) == jetPtGt25Num_ and
      (*jetPt)[0] > minJetPt_ and
      std::fabs((*jetEta)[0]) > minJetEta_ and
      std::fabs((*jetEta)[0]) < maxJetEta_ and
      (*jetNEEF)[0] > minJetNeef_ and
      std::fabs(TVector2::Phi_mpi_pi(*pfmetPhi - (*jetPhi)[0])) > minDPhi_)
    return true;

  return false;
}


bool BERegion::isHBHENoiseBER() {
  if (*pfmet > minPfMet_ and
      (*jetPt).size() == jetPtGt200Num_ and
      (*numJetsPt25) == jetPtGt25Num_ and
      (*jetPt)[0] > minJetPt_ and
      std::fabs((*jetEta)[0]) < maxJetEta_ and
      (*jetNHEF)[0] > minJetNhef_ and
      std::fabs(TVector2::Phi_mpi_pi(*pfmetPhi - (*jetPhi)[0])) > minDPhi_)
    return true;

  return false;
}


bool BERegion::isOthersBER() {
  return (*pfmet > minPfMet_);
}


template <typename T>
T BERegion::Get(std::string const field) {
  if (cfg_[berCat_][field])
    return cfg_[berCat_][field].as<T>();
  else {
    std::ostringstream message;
    message << "BER configuration of " << berCat_
      << " doesn't contain [" << field << "] field.";
    throw std::runtime_error(message.str());
  }
}

