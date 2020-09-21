#ifndef METFILTER_H
#define METFILTER_H

#include <iostream>
#include <map>
#include <math.h>

#include <TEfficiency.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>


typedef TTreeReaderValue<std::vector<double>> doubleVectorReader;


class MetFilter {
 public:
   // COnstructor
   MetFilter(TString const &dataset, TTreeReader &reader, TString const &name,
       TString const &flag, bool const isRecmnded, int const &N_1Key,
       bool const isOldAnlz=false);
   // Destructor
   ~MetFilter() noexcept;

   // A method for filling filters efficiencies
   void Fill(double const &x, std::array<TString, 2> const &info,
       int const &passRecStatus);

   // A method for filling filters occupancy maps
   void Fill(doubleVectorReader &etaVec, doubleVectorReader &phiVec,
       TString const &region, int const &passRecStatus);

   // Setting directory of output root file
   void SetDirectory(TFile &outputRootFile);

   // Checkes whether the events id passed the filter or not.
   bool IsPassed();

   // Checks whether the filter is a recommended MET filters or not.
   bool IsRecommeded() const;

   // Retruns the filter name
   TString GetName() const;

   // Retruns the filter flag
   TString GetFlag() const;

   // A dictionary for different regions which the abbreviations are the keys.
   static const std::array<std::pair<TString, TString>, 3> regions;

   // Variables array
   static const std::array<TString, 4> vars;

   enum Filters{GoodVertices, GlbSTightHalo, HBHENoise, HBHEIsoNoise,
     EcalTP, EEBadSc, EcalBadCal_U, BadPFMuon, BadPFMuon_U, BadPFMuon_DzU,
     EcalDCBE, EcalDCBE_U, GlbTightHalo, BadChCand, BadChCand_U, Others};

 private:
   // Filter's name.
   TString const name_;

   // Filter's flag for accesing the filter's branch.
   TString const flag_;

   // Event's status w.r.t. the filter's respond.
   TTreeReaderValue<Bool_t> isPassed_;

   // Indicates whether the filters is a recommended MET filter or not.
   bool const isRecmnded_;

   // Indicates whether other filters are passed recommended filters.
   int passRecKey_;

   // Indicates whether the analysis is based on old strategy or not.
   // Old strategy means without applying NFR and BER constraints.
   bool const isOldAnlz_;

   std::vector<TEfficiency> efficiency_;

   std::vector<TH1F> numerator_, denominator_;

   std::vector<TH2F> failedOccupancyMap_, passedOccupancyMap_;

   std::map<TString, std::map<TString, unsigned short int>> effIndex_;

   std::map<TString, unsigned short int> occMapIndex_;
};

#endif
