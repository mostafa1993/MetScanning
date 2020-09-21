#ifndef NFREGION_H
#define NFREGION_H

#include <TH1F.h>

#include "Tree.h"


class NFRegion final : private Tree {
 public:
   // Constructor
   NFRegion(TTreeReader &reader, TString dataset, TString year, bool fillMore);

   bool operator()();

 private:
   TString dataset_;

   std::vector<short int> favoredCandIdx_;

   // Histograms for cross-check plots.
   std::vector<TH1F> ptDist_, etaDist_, phiDist_;
   TH1F massDist_, ptRatioDist_;

   // Checks whether a given jet passes JetID for 2018 datasets
   bool passJetId_2018(int i);
   // Checks whether a given jet passes JetID for 2017 datasets
   bool passJetId_2017(int i);
   // Checks whether a given jet passes JetID for 2016 datasets
   bool passJetId_2016(int i);

   // A pointer to one of functions checking whether a given jet passes JetID
   bool (NFRegion::*passJetId)(int i);

   // Check whether the event is in Noise Free Region for JetHT dataset
   bool isJetHtNFR();

   // Check whether the event is in Noise Free Region for Double Muon dataset
   bool isDMuonNFR();

   // Check whether the event is in Noise Free Region for
   // EGamma or Single Electron dataset
   bool isEGammaNFR();

   // Check whether the event is in Noise Free Region for Single Muon dataset
   bool isSMuonNFR();

   // Check whether the event is in Noise Free Region for MC datasets
   bool isSimNFR();

   // A pointer to one of functions checking whether the event is in NFR or not
   bool (NFRegion::*isNFR)();

   void Fill(doubleVectorReader &var, int const &varType);

   void Fill(double const &var, int const &varType);

   // Indicates whether the cross-check plots should be filled or not.
   bool fillMore_;
};

#endif
