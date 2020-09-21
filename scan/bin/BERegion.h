#ifndef BEREGION_H
#define BEREGION_H

#include <yaml-cpp/yaml.h>

#include "Tree.h"


class BERegion final : private Tree {
 public:
   // Constructor
   BERegion(TTreeReader &reader, const std::string="Others");
   bool operator()();

 private:
   std::string berCat_;

   // BER configuration.
   YAML::Node cfg_;

   // All BER parameters
   float minPfMet_, minDPhi_, maxPfMetPhi_, minPfMetPhi_, minJetPt_, maxJetEta_,
         minJetEta_, maxJetChef_, minJetNeef_, minJetNhef_, minMuPt_,
         maxAddMuPt_, minDeltaR_, pfCanMinPt_;

   int pfCanPdgId_, jetPtGt25Num_;

   long unsigned int jetPtGt200Num_, minMuNum_; 

   // Check whether the event is in BER of BadChCand filters.
   bool isBadChCandBER();

   // Check whether the event is in BER of GlbSTHalo and GlbTHalo filters.
   bool isGlbSTHaloBER();

   // Check whether the event is in BER of BadPFMuon filters.
   bool isBadPFMuonBER();

   // Check whether the event is in BER of EcalBadCal filters.
   bool isEcalBadCalBER();

   // Check whether the event is in BER of HBHENoise and HBHEIsoNoise filters.
   bool isHBHENoiseBER();

   // Check whether the event is in BER of other filters.
   bool isOthersBER();

   // A pointer to one of functions checking whether the event is in BER or not.
   bool (BERegion::*isBER)();

   // Return the value of given field from configuration file.
   template <typename T>
   T Get(std::string const field);
};

#endif
