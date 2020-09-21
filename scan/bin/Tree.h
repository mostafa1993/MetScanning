#ifndef TREE_H
#define TREE_H

#include <TTreeReader.h>
#include <TTreeReaderValue.h>


typedef TTreeReaderValue<std::vector<double>> doubleVectorReader;
typedef TTreeReaderValue<std::vector<int>> intVectorReader;

class Tree {
 public:
   // Constructor
   Tree(TTreeReader &reader);

   TTreeReaderValue<ULong64_t> run, lumi, event;

   TTreeReaderValue<int> npvx, numJetsPt25;

   TTreeReaderValue<double> pfmet, pfmetPhi, puppimet;

   doubleVectorReader jetPt, jetEta, jetPhi, jetEn;

   doubleVectorReader jetCHEF, jetNEEF, jetNHEF, jetMuEF, jetCEEF;

   doubleVectorReader muonPt, muonEta, muonPhi, muonEn, muonIso;

   doubleVectorReader elePt, eleEta, elePhi, eleEn;

   doubleVectorReader pfCandPt, pfCandPhi;

   intVectorReader pfCandPdgId, jetChM, jetNM;

   intVectorReader isMediumMuon, muonCharge;

   intVectorReader isMediumEle,eleCharge;

   TTreeReaderValue<std::vector<bool>> isPFMuon;
};

#endif
