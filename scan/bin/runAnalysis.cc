#include <bitset>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <forward_list>
#include <iostream>
#include <math.h>
#include <yaml-cpp/yaml.h>

#include <TChain.h>
#include <TVector2.h>

#include "BERegion.h"
#include "MetFilter.h"
#include "NFRegion.h"
#include "Tree.h"


namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  unsigned int step;
  int event_i = 0, maxEvents;
  std::string output, dataset, year, oldAnlz, verbosity;
  bool isVerbose, crossCheckPlots, passRecOption, isOldAnlz;
  std::vector<std::string> inputs;
  YAML::Node filtersConfig = YAML::LoadFile("config/met_filters.yaml");
  int const metFiltersNum = filtersConfig.size();
  // Stores which recommended MET filters are passed.
  int passRecStatus = 255;

  std::array<int, 16> statistics{};

  enum EffCat{nfr_pfmet, nfr_pupmet, nfr_lJpt, nfr_npvx,
    ber_pfmet, ber_pupmet, ber_lJpt, ber_npvx,
    ar_pfmet, ar_pupmet, ar_lJpt, ar_npvx};
  enum OccMapCat{nfr, ber, ar};

  std::map<unsigned short int, std::array<TString, 2>> effCat;
  std::map<unsigned short int, TString> occMap;

  // Managing program options.
  po::options_description optionsDescription{"Options"};
  optionsDescription.add_options()
    ("help,h", "Help screen")
    ("year", po::value<std::string>(&year)->required(), "Dataset year")
    ("output,o", po::value<std::string>(&output)->required(),
     "Output root file")
    ("dataset", po::value<std::string>(&dataset)->required(), "Dataset name")
    ("pass-rec", po::value<bool>(&passRecOption)->default_value(true), "Require"
     "events to pass the recommended MET filters except the one under study")
    ("max-events", po::value<int>(&maxEvents)->default_value(-1),
     "Maximum number of events to process.")
    ("step", po::value<unsigned int>(&step)->default_value(100000),
     "Progress steps (if --verbose[-v] is passed to the program)")
    ("cross-check", po::value<bool>(&crossCheckPlots)->default_value(false),
     "Plot cross-check histograms")
    ("old-analysis", "Enables the analysis without control regions Cuts.")
    ("verbose,v", "Show progress");

  po::options_description hiddenOptionsDescription;
  hiddenOptionsDescription.add_options()
    ("input-files", po::value<std::vector<std::string>>(), "");
  po::positional_options_description posOptionsDescription;
  posOptionsDescription.add("input-files", -1);

  po::options_description allOptionsDescription;
  allOptionsDescription.add(optionsDescription).add(hiddenOptionsDescription);

  po::variables_map options;
  po::store(
      po::command_line_parser(argc, argv).options(allOptionsDescription)
      .positional(posOptionsDescription).run(),
      options);

  if (options.count("help")) {
    std::cerr << "Usage:" << std::endl
      << " runAnalysis INPUT_FILE1 [INPUT_FILE2 [...]] [OPTIONS]" << std::endl;
    std::cerr << optionsDescription << std::endl;
    return EXIT_SUCCESS;
  }

  isOldAnlz = options.count("old-analysis");
  isVerbose = options.count("verbose");
  po::notify(options);

  // If no input file is provided throw an exception.
  if (options["input-files"].empty()) {
    std::ostringstream message;
    message << "No input file is provided. Please pass them as POSITIONAL "
      << "arguments.";
    throw std::runtime_error(message.str());
  }

  inputs = options["input-files"].as<std::vector<std::string>>();
  if( maxEvents == -1)
    maxEvents = INT_MAX;

  TString datasetYear = dataset + "_" + year;
  boost::algorithm::to_lower(dataset);
  unsigned short int I = 0, J = 0;
  for (auto const &region : MetFilter::regions) {
    occMap[I++] = region.first;
    for (auto const &var : MetFilter::vars)
      effCat[J++] = {{region.first, var}};
  }

  TFile outputRootFile(output.c_str(), "RECREATE");
  TChain chain("ntuplemakerminiaod/tree");
  for (auto const &file : inputs)
    chain.Add(file.c_str());

  if (isVerbose)
    std::cout << "Input files are added into the chain." << std::endl;

  TTreeReader reader(&chain);
  Tree tree(reader);
  NFRegion isNFR(reader, dataset, year, crossCheckPlots);
  BERegion isBadChCandBER(reader, "BadChCand"),
           isBadPFMuonBER(reader, "BadPFMuon"),
           isGlbSTHaloBER(reader, "GlbSTightHalo"),
           isEcalBadCalBER(reader, "EcalBadCal_U"),
           isHBHENoiseBER(reader, "HBHENoise"),
           isOthersBER(reader);

  std::map<unsigned int, BERegion*> isBER;
  isBER[MetFilter::GoodVertices] = &isOthersBER;
  isBER[MetFilter::GlbSTightHalo] = &isGlbSTHaloBER;
  isBER[MetFilter::HBHENoise] = &isHBHENoiseBER;
  isBER[MetFilter::HBHEIsoNoise] = &isHBHENoiseBER;
  isBER[MetFilter::EcalTP] = &isOthersBER;
  isBER[MetFilter::EEBadSc] = &isOthersBER;
  isBER[MetFilter::EcalBadCal_U] = &isEcalBadCalBER;
  isBER[MetFilter::BadPFMuon] = &isBadPFMuonBER;
  isBER[MetFilter::BadPFMuon_U] = &isBadPFMuonBER;
  isBER[MetFilter::BadPFMuon_DzU] = &isBadPFMuonBER;
  isBER[MetFilter::EcalDCBE] = &isOthersBER;
  isBER[MetFilter::EcalDCBE_U] = &isOthersBER;
  isBER[MetFilter::GlbTightHalo] = &isGlbSTHaloBER;
  isBER[MetFilter::BadChCand] = &isBadChCandBER;
  isBER[MetFilter::BadChCand_U] = &isBadChCandBER;

  // Uncomment below lines if you need to use run, lumi and event number.
  //auto run = tree.run;
  //auto lumi = tree.lumi;
  //auto event = tree.event;
  auto npvx = tree.npvx;
  auto pfmet = tree.pfmet;
  auto puppimet = tree.puppimet;
  auto jetPt = tree.jetPt;
  auto jetEta = tree.jetEta;
  auto jetPhi = tree.jetPhi;

  std::vector<MetFilter> metFilters;
  for (auto const &filter : filtersConfig) {
   metFilters.emplace_back(
       datasetYear, reader,
       filter["name"].as<std::string>(),
       filter["flag"].as<std::string>(),
       filter["isRecommended"].as<bool>(),
       filter["N-1key"].as<int>(),
       isOldAnlz);
  }

  for (auto &f : metFilters)
    f.SetDirectory(outputRootFile);

  // Main loop.
  while(reader.Next()) {
    event_i++;
    if (event_i > maxEvents)
      break;

    if (isVerbose and (event_i % step == 0))
      std::cout << event_i << " events are processed" << std::endl;

    // Applying N-1 performance analysis if the option is enabled.
    passRecStatus = -1;
    if (passRecOption) {
      passRecStatus = 255;
      for (int i = 0; i < metFiltersNum; i++)
        if (metFilters[i].IsRecommeded() and not metFilters[i].IsPassed())
          passRecStatus ^= 1 << i;
    }

    if (isOldAnlz) {
      for (int i = 0; i < metFiltersNum; i++) {
        metFilters[i].Fill(*pfmet, effCat[ar_pfmet], passRecStatus);
        metFilters[i].Fill(*puppimet, effCat[ar_pupmet], passRecStatus);
        metFilters[i].Fill(*npvx, effCat[ar_npvx], passRecStatus);
        if ((*jetPt).size() > 0)
          metFilters[i].Fill((*jetPt)[0], effCat[ar_lJpt], passRecStatus);
        metFilters[i].Fill(jetEta, jetPhi, occMap[ar], passRecStatus);
      }
      // As the old analysis option is enabled, don't go to the control regions.
      continue;
    }

    // Checking whether the event is in Noise Free Region w.r.t. the data-set.
    if (isNFR()) {
      statistics[0]++;
      for (int i = 0; i < metFiltersNum; i++) {
        metFilters[i].Fill(*pfmet, effCat[nfr_pfmet], passRecStatus);
        metFilters[i].Fill(*puppimet, effCat[nfr_pupmet], passRecStatus);
        metFilters[i].Fill(*npvx, effCat[nfr_npvx], passRecStatus);
        if ((*jetPt).size() > 0)
          metFilters[i].Fill((*jetPt)[0], effCat[nfr_lJpt], passRecStatus);
        metFilters[i].Fill(jetEta, jetPhi, occMap[nfr], passRecStatus);
      }
    }

    // Check whether the event is in BKG Enriched Region w.r.t. the filters.
      for (int i = 0; i < metFiltersNum; i++) {
        if (not (*isBER[i])())
          continue;
        statistics[i+1]++;
        metFilters[i].Fill(*pfmet, effCat[ber_pfmet], passRecStatus);
        metFilters[i].Fill(*puppimet, effCat[ber_pupmet], passRecStatus);
        metFilters[i].Fill(*npvx, effCat[ber_npvx], passRecStatus);
        if ((*jetPt).size() > 0)
          metFilters[i].Fill((*jetPt)[0], effCat[ber_lJpt], passRecStatus);
        metFilters[i].Fill(jetEta, jetPhi, occMap[ber], passRecStatus);
      }
  } // End of main loop.

  outputRootFile.cd();
  outputRootFile.Write();
  outputRootFile.Close();
  if (isVerbose)
    std::cout << event_i << " events are processed." << std::endl;

  if (isOldAnlz or not isVerbose)
    return 0;

  std::cout << std::endl << "------------------------------------" << std::endl;
  std::cout << "Statistics Summary:" << std::endl;
  std::string region = "";
  std::cout.fill(' ');
  for (long unsigned int i = 0; i < statistics.size(); i++) {
    if (i == 0)
      region = "NFR";
    else
      region = filtersConfig[i-1]["name"].as<std::string>() + " BER";
    std::cout << statistics[i] << "\t Events in " << region << std::endl;
  }

  return 0;
}

