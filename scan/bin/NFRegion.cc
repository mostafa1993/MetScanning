#include "NFRegion.h"

#include <TLorentzVector.h>
#include <TVector2.h>


enum {pt, eta, phi, mass, ptr};

NFRegion::NFRegion(TTreeReader &reader, TString dataset, TString year,
    bool fillMore)
    : Tree(reader),
      fillMore_{fillMore} {

  if (dataset == "jetht") {
    isNFR = &NFRegion::isJetHtNFR;
    if (fillMore_) {
      ptDist_.emplace_back("jet0_pt", "jet0;pt;events", 1000, 0., 1000.);
      ptDist_.emplace_back("jet1_pt", "jet1;pt;events", 1000, 0., 1000.);
      etaDist_.emplace_back("jet0_eta", "jet0;eta;events", 100, -5, 5);
      etaDist_.emplace_back("jet1_eta", "jet1;eta;events", 100, -5, 5);
      phiDist_.emplace_back("jet0_phi", "jet0;phi;events", 100, -3.5, 3.5);
      phiDist_.emplace_back("jet1_phi", "jet1;phi;events", 100, -3.5, 3.5);
      ptRatioDist_ = TH1F("pt_ratio", "pt0/pt1;ratio;events", 100, 0.9, 1.3);
      massDist_ = TH1F("diJet_mass", "dijet mass;mass;events", 1000, 0, 5000);
    }
  }

  else if (dataset == "doublemuon") {
    isNFR = &NFRegion::isDMuonNFR;
    if (fillMore_) {
      ptDist_.emplace_back("muon0_pt", "muon0;pt;events", 1000, 0., 1000.);
      ptDist_.emplace_back("muon1_pt", "muon1;pt;events", 1000, 0., 1000.);
      etaDist_.emplace_back("muon0_eta", "muon0;eta;events", 100, -5, 5);
      etaDist_.emplace_back("muon1_eta", "muon1;eta;events", 100, -5, 5);
      phiDist_.emplace_back("muon0_phi", "muon0;phi;events", 90, -3.5, 3.5);
      phiDist_.emplace_back("muon1_phi", "muon1;phi;events", 90, -3.5, 3.5);
      massDist_ = TH1F("Z_mass", "", 1000, 0, 1000);
    }
  }

  else if (dataset == "egamma" or dataset == "singleelectron") {
    isNFR = &NFRegion::isEGammaNFR;
    if (fillMore_) {
      ptDist_.emplace_back("ele0_pt", "ele0;pt;events", 1000, 0., 1000.);
      ptDist_.emplace_back("ele1_pt", "ele1;pt;events", 1000, 0., 1000.);
      etaDist_.emplace_back("ele0_eta", "ele0;eta;events", 100, -5, 5);
      etaDist_.emplace_back("ele1_eta", "ele1;eta;events", 100, -5, 5);
      phiDist_.emplace_back("ele0_phi", "ele0;phi;events", 90, -3.5, 3.5);
      phiDist_.emplace_back("ele1_phi", "ele1;phi;events", 90, -3.5, 3.5);
      massDist_ = TH1F("Z_mass", "", 1000, 0, 1000);
    }
  }

  else if (dataset == "singlemuon") {
    isNFR = &NFRegion::isSMuonNFR;
    if (fillMore_) {
      ptDist_.emplace_back("muon0_pt", "muon0;pt;events", 1000, 0., 1000.);
      etaDist_.emplace_back("muon0_eta", "muon0;eta;events", 100, -5, 5);
      phiDist_.emplace_back("muon0_phi", "muon0;phi;events", 90, -3.5, 3.5);
      phiDist_.emplace_back("pfmetPhi", "PF MET;phi;events", 90, -3.5, 3.5);
      massDist_ = TH1F("mt", "mt of muon+pfmet;m_t;events", 1000, 0, 1000);
    }
  }

  else if (dataset == "mc") {
    isNFR = &NFRegion::isSimNFR;
  }


  if (year == "2018")
    passJetId = &NFRegion::passJetId_2018;
  else if (year == "2017")
    passJetId = &NFRegion::passJetId_2017;
  else if (year == "2016")
    passJetId = &NFRegion::passJetId_2016;
  
}


bool NFRegion::operator()() {
  favoredCandIdx_.clear();
  return (this->*isNFR)();
}


bool NFRegion::passJetId_2018(int i) {
  if (std::fabs((*jetEta)[i]) <= 2.6) {
    return ((*jetNHEF)[i] < 0.9 and (*jetNEEF)[i] < 0.9 and
        ((*jetChM)[i] + (*jetNM)[i]) > 1 and (*jetMuEF)[i] < 0.8 and
        (*jetCHEF)[i] > 0 and (*jetChM)[i] > 0 and (*jetCEEF)[i] < 0.8);
  }

  else if (std::fabs((*jetEta)[i]) <= 2.7) {
    return ((*jetNHEF)[i] < 0.9 and (*jetNEEF)[i] < 0.99 and
        (*jetMuEF)[i] < 0.8 and (*jetChM)[i] > 0 and (*jetCEEF)[i] < 0.8);
  }

  else if (std::fabs((*jetEta)[i]) <= 3.0) {
    return (*jetNEEF)[i] > 0.02 and (*jetNEEF)[i] < 0.99 and (*jetNM)[i] > 2;
  }

  else if (std::fabs((*jetEta)[i]) <= 5.0) {
    return (*jetNHEF)[i] > 0.2 and (*jetNEEF)[i] < 0.9 and (*jetNM)[i] > 10;
  }

  return false;
}


bool NFRegion::passJetId_2017(int i) {
  if (std::fabs((*jetEta)[i]) <= 2.7) {
    return ((*jetNHEF)[i] < 0.9 and (*jetNEEF)[i] < 0.9 and
        ((*jetChM)[i] + (*jetNM)[i]) > 1 and (*jetMuEF)[i] < 0.8 and
        // Applying additional constraints for abs(eta) <= 2.4
        ((std::fabs((*jetEta)[i]) <= 2.4) ?
        (*jetCHEF)[i] > 0 and (*jetChM)[i] > 0 and (*jetCEEF)[i] < 0.8 : true));
  }

  else if (std::fabs((*jetEta)[i]) <= 3.0) {
    return (*jetNEEF)[i] > 0.02 and (*jetNEEF)[i] < 0.99 and (*jetNM)[i] > 2;
  }

  else if (std::fabs((*jetEta)[i]) > 3.0) {
    return (*jetNHEF)[i] > 0.02 and (*jetNEEF)[i] < 0.9 and (*jetNM)[i] > 10;
  }

  return false;
}


bool NFRegion::passJetId_2016(int i) {
  if (std::fabs((*jetEta)[i]) <= 2.7) {
    return (*jetNHEF)[i] < 0.99 and
      (*jetNEEF)[i] < 0.99 and
      ((*jetChM)[i] + (*jetNM)[i]) > 1 and
      // Applyimg additional constraints for abs(eta) <= 2.4
      ((std::fabs((*jetEta)[i]) <= 2.4) ?
       ((*jetCHEF)[i] > 0 and (*jetChM)[i] > 0 and (*jetCEEF)[i] < 0.99)
       : true);
  }

  else if (std::fabs((*jetEta)[i]) <= 3.0) {
    return (*jetNEEF)[i] > 0.01 and
      (*jetNHEF)[i] < 0.98 and
      (*jetNM)[i] > 2;
  }

  else if (std::fabs((*jetEta)[i]) > 3.0) {
    return (*jetNEEF)[i] < 0.9 and (*jetNM)[i] > 10;
  }

  return false;
}


bool NFRegion::isJetHtNFR() {
  if (*pfmet > 100)
    return false;

  for (long unsigned int i = 0; i < (*jetPt).size(); i++) {
    if ((*jetPt)[i] < 200)
      continue;

    if ((this->*passJetId)(i))
      favoredCandIdx_.push_back(i);

    if (favoredCandIdx_.size() == 2)
      break;
  }

  if (favoredCandIdx_.size() >= 2) {
    short int const &O = favoredCandIdx_[0], &one = favoredCandIdx_[1];
    auto const deltaPhi = std::fabs(
        TVector2::Phi_mpi_pi((*jetPhi)[O] - (*jetPhi)[one]));
    if (deltaPhi >= 2.9) {
      auto const ptRatio = std::fabs((*jetPt)[O] / (*jetPt)[one]);
      if (ptRatio > 0.8 and ptRatio < 1.2) {
        if (not fillMore_)
          return true;

        TLorentzVector jet0P4(0, 0 ,0 ,0), jet1P4(0, 0, 0, 0);
        jet0P4.SetPtEtaPhiE(
            (*jetPt)[O], (*jetEta)[O], (*jetPhi)[O], (*jetEn)[O]);
        jet1P4.SetPtEtaPhiE(
            (*jetPt)[one], (*jetEta)[one], (*jetPhi)[one], (*jetEn)[one]);

        double const diJetMass = (jet0P4 + jet1P4).M();
        Fill(jetPt, pt);
        Fill(jetEta, eta);
        Fill(jetPhi, phi);
        Fill(diJetMass, mass);
        Fill(ptRatio, ptr);
        return true;
      }
    }
  }

  return false;
}


bool NFRegion::isDMuonNFR() {
  if (*pfmet > 50)
    return false;

  for (long unsigned int i = 0; i < (*isPFMuon).size(); i++)
    if ((*isPFMuon)[i] and
        (*isMediumMuon)[i] == 1 and
        (*muonIso)[i] < 0.2 and
        (*muonPt)[i] > 30)
      favoredCandIdx_.push_back(i);

  if (favoredCandIdx_.size() == 2) {
    short int const &O = favoredCandIdx_[0], &one = favoredCandIdx_[1];
    if ((*muonCharge)[O] * (*muonCharge)[one] != -1)
      return false;

    TLorentzVector mu0P4(0, 0 ,0 ,0), mu1P4(0, 0, 0, 0);

    mu0P4.SetPtEtaPhiE((*muonPt)[O], (*muonEta)[O],
        (*muonPhi)[O], (*muonEn)[O]);
    mu1P4.SetPtEtaPhiE((*muonPt)[one], (*muonEta)[one],
        (*muonPhi)[one], (*muonEn)[one]);

    double const invMass = (mu0P4 + mu1P4).M();
    if (invMass > 81 and invMass < 101) {
      if (not fillMore_)
        return true;

      Fill(muonPt, pt);
      Fill(muonEta, eta);
      Fill(muonPhi, phi);
      Fill(invMass, mass);
      return true;
    }
  }

  return false;
}


bool NFRegion::isEGammaNFR() {
  if (*pfmet > 50)
    return false;

  for (long unsigned int i = 0; i < (*elePt).size(); i++)
    if ((*elePt)[i] > 30.0 and (*isMediumEle)[i] == 1)
      favoredCandIdx_.push_back(i);

  if (favoredCandIdx_.size() == 2) {
    short int const &O = favoredCandIdx_[0], &one = favoredCandIdx_[1];
    if ((*eleCharge)[O] * (*eleCharge)[one] != -1)
      return false;

    TLorentzVector el0P4(0, 0 ,0 ,0), el1P4(0, 0, 0, 0);

    el0P4.SetPtEtaPhiE((*elePt)[O], (*eleEta)[O],
        (*elePhi)[O], (*eleEn)[O]);
    el1P4.SetPtEtaPhiE((*elePt)[one], (*eleEta)[one],
        (*elePhi)[one], (*eleEn)[one]);

    double const invMass = (el0P4 + el1P4).M();
    if (invMass > 81 and invMass < 101) {
      if (not fillMore_)
        return true;

      Fill(elePt, pt);
      Fill(eleEta, eta);
      Fill(elePhi, phi);
      Fill(invMass, mass);
      return true;
    }
  }

  return false;
}


bool NFRegion::isSMuonNFR() {
  if (*pfmet > 100.0)
    return false;

  TLorentzVector pfmetP4(0, 0 ,0 ,0);
  pfmetP4.SetPtEtaPhiE(*pfmet, 0, *pfmetPhi, *pfmet);
  for (long unsigned int i = 0; i < (*isPFMuon).size(); i++)
    if ((*isPFMuon)[i] and
        (*isMediumMuon)[i] == 1 and
        (*muonIso)[i] < 0.2 and
        (*muonPt)[i] > 30.0) {
      TLorentzVector muP4(0, 0 ,0 ,0);

      muP4.SetPtEtaPhiE(
          (*muonPt)[i], (*muonEta)[i], (*muonPhi)[i], (*muonEn)[i]);

      auto const mT = (muP4 + pfmetP4).Mt();
      if (mT <= 120.0) {
        favoredCandIdx_.push_back(i);
        if (not fillMore_)
          return true;

        Fill(muonPt, pt);
        Fill(muonEta, eta);
        Fill(muonPhi, phi);
        Fill(*pfmetPhi, phi);
        Fill(mT, mass);
        return true;
      }
    }

  return false;
}


bool NFRegion::isSimNFR() {
  if (*pfmet < 100)
    return true;

  return false;
}


void NFRegion::Fill(doubleVectorReader &var, int const &varType) {
  short int i = 0;
  switch (varType) {
    case pt:
      for (auto const &idx : favoredCandIdx_)
        ptDist_[i++].Fill((*var)[idx]);
      break;

    case eta:
      for (auto const &idx : favoredCandIdx_)
        etaDist_[i++].Fill((*var)[idx]);
      break;

    case phi:
      for (auto const &idx : favoredCandIdx_)
        phiDist_[i++].Fill((*var)[idx]);
      break;
  }
}


void NFRegion::Fill(double const &var, int const &varType) {
  switch (varType) {
    case mass:
      massDist_.Fill(var);
      break;

    case ptr:
      ptRatioDist_.Fill(var);
      break;
  }
}

