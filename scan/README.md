# METScanning

This framework aims to centralize the MET scanners analysis in a more efficient framework.
For more details see the [Wiki page](https://gitlab.cern.ch/mmahdavi/metscanning/-/wikis/home).


# Mannual

1. Set the environment and bBuild the executable file if needed.
```
cmsenv
scram b -j $(nproc)
```
<br/>

2. Run the analysis
```
runAnalysis INPUT_FILE1 [INPUT_FILE2 [...]] --output output.root --dataset DATASET --year YYYY --pass-rec [ON/off] --cross-check [on/OFF]
```
***`--cross-check ON` option will allow you to have cross-check histograms in NFR***

- One can also run ```runAnalysis --help``` or ```runAnalysis -h``` for more detailed help message.
- One may need to modify the configuration files stored in [`config`](https://github.com/cms-met/MetScanning/tree/master/scan/config) directory based on her/his analysis, such as different years, scenarios and etc. 
- If one is interested in old-analysis strategy i.e. without BER and NFR control region, it is possible by passing the `--old-analysis` option to the program. Hence the obtained efficencies are with out any cut or ID applied.
<br/>

# Plotting Mannual
For running plot scripts, you need to logout and login again if you have set `CMSSW` environment i.e. `cmsenv`. Then you need to set the proper environment,
```
cd plot
. ./env
cd -
```
There are two different python scripts which provide `efficiencies` and `occupancy maps` plots such as [plot_efficencies](https://github.com/cms-met/MetScanning/tree/master/scan/plot/plot_efficiencies) and [plot_occupancy_maps](https://github.com/cms-met/MetScanning/tree/master/scan/plot/plot_occupancy_maps) respectively. They are mentioned below separately.

Before using the scripts, please make sure that the environment has been set up as mention above.

- **Occupancy maps:**
```
plot_occupancy_maps [list of output root-files separated by space] --dataset [list of datasets w.r.t provided output root-files separated by space] --year YYYY
```
<br/>

- **Efficiency plots:**

For full plots, i.e all filters plots individually and merged plots for each provided dataset, `--all` option is needed like below.
```
plot_efficencies [list of output root-files separated by space] --dataset [list of datasets w.r.t provided output root-files separated by space] --year 2018 --all
```

If one needs also to obtain the actual final merged plot, all datasets must be provided. In this case the `--all` option is optional, for example, below command will plot all kind of plots and the actual final merged plot as well.
```
plot_efficencies jetht_output.root singlemuon_output.root egamma_output.root --dataset jetht singlemuon egamma --year 2018 --all
```
<br/>

- **Plots Configurations:**

The `plot_efficencies` script uses `yaml` configuration files which contain `binning` settings according to regions (`BER` and `NFR`) variables and/or datasets,  for each filter separately.

The configuration files are provided in the [config](https://github.com/cms-met/MetScanning/tree/master/scan/config) directory. One may need to modify theses configuration files to obtian the proper plots.
<br/>
<br/>

- *One can also run* ```plot_efficencies [--help, -h]``` or ```plot_occupancy_maps [--help, -h]``` for more detailed help messages.
