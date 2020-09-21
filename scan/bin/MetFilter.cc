#include "MetFilter.h"


MetFilter::MetFilter(TString const &dataset, TTreeReader &reader,
    TString const &name, TString const &flag, bool const isRecmnded,
    int const &N_1Key, bool const isOldAnlz)
    : name_{name},
      flag_{flag},
      isPassed_{reader, flag},
      isRecmnded_{isRecmnded}, 
      passRecKey_{N_1Key}, 
      isOldAnlz_{isOldAnlz} {

  int min = 0, max = 5000, nbins = max - min;
  unsigned short int i = 0, j = 0;
  for(auto const &region : regions) {
    if ((isOldAnlz) ? region.first != "AR" : region.first == "AR")
      continue;
    TString tempTitle = name_ + ", " + dataset + ", " + region.second;
    TString occTitleF = "Jets Occupancy Map [pt > 200] for events failing ";
    TString occTitleP = "Jets Occupancy Map [pt > 200] for events passing ";
    TString fullNameF = "occMap_failed_" + name_ + "_" + region.first;
    TString fullNameP = "occMap_passed_" + name_ + "_" + region.first;
    TString titleF = occTitleF + tempTitle + ";#eta;#phi;Events";
    TString titleP = occTitleP + tempTitle + ";#eta;#phi;Events";
    failedOccupancyMap_.emplace_back(fullNameF, titleF,
        1000, -5, 5, 314, -3.1416, 3.1416);
    passedOccupancyMap_.emplace_back(fullNameP, titleP,
        1000, -5, 5, 314, -3.1416, 3.1416);
    occMapIndex_[region.first] = i++;

    for(auto const &var : vars) {
      TString tempName = name_+ "_" + region.first + "_" + var;
      TString fullName =  "eff_" + tempName;
      TString title = tempTitle + ";" + var +
        ((region.first != "BER") ? ";Efficiency" : ";NoiseRejectionFraction");

      auto fullNameDen =  "den_" + tempName;
      auto fullNameNum =  "num_" + tempName;
      auto numTitle = tempTitle + ", numerator" + ";" + var + ";Events";
      auto denTitle = tempTitle + ", denominator" + ";" + var + ";Events";

      efficiency_.emplace_back(fullName, title, nbins, min, max);
      denominator_.emplace_back(fullNameDen, denTitle, nbins, min, max);
      numerator_.emplace_back(fullNameNum, numTitle, nbins, min, max);
      effIndex_[region.first][var] = j++;
    }
  }
}


const std::array<std::pair<TString, TString>, 3> MetFilter::regions = {
  {{"NFR", "Noise Free Region"}, {"BER", "Background Enriched Region"},
    {"AR", "All Regions"}}};


const std::array<TString, 4> MetFilter::vars = {{"pfMet", "puppiMet",
  "leadingJetPt", "nPVx"}};


MetFilter::~MetFilter() noexcept {}


void MetFilter::Fill(double const &x, std::array<TString, 2> const &info,
    int const &passRecStatus) {
  if (passRecStatus == passRecKey_ or passRecStatus == 255 or
      passRecStatus == -1) {
    TString const &region = info[0];
    TString const &var = info[1];
    auto &i = effIndex_[region][var];
    efficiency_[i].Fill((region != "BER") ? *isPassed_ : not(*isPassed_), x);
    denominator_[i].Fill(x);
    if ((*isPassed_ and region != "BER") or
        (not(*isPassed_) and region == "BER"))
      numerator_[i].Fill(x);
  }
}


void MetFilter::Fill(
    doubleVectorReader &etaVec, doubleVectorReader &phiVec,
    TString const &region, int const &passRecStatus) {
  if (passRecStatus == passRecKey_ or passRecStatus == 255 or
      passRecStatus == -1) {
    auto &i = occMapIndex_[region];

    if (!(*isPassed_))
      for (long unsigned int j = 0; j < (*etaVec).size(); j++)
        failedOccupancyMap_[i].Fill((*etaVec)[j], (*phiVec)[j]);

    if (*isPassed_)
      for (long unsigned int j = 0; j < (*etaVec).size(); j++)
        passedOccupancyMap_[i].Fill((*etaVec)[j], (*phiVec)[j]);
  }
}


void MetFilter::SetDirectory(TFile &outputRootFile) {
  for (long unsigned int i = 0; i < efficiency_.size(); i++)
    efficiency_[i].SetDirectory(&outputRootFile);
}


bool MetFilter::IsPassed() {
  return *isPassed_;
}


bool MetFilter::IsRecommeded() const {
  return isRecmnded_;
}


TString MetFilter::GetName() const {
  return name_;
}


TString MetFilter::GetFlag() const {
  return flag_;
}

