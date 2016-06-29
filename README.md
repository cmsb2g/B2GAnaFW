# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a prroduction version of the B2G EDMNtuples. Updates since last production tag ``v7.6.x_v1.2``:
- Adding N-subjettiness for soft drop subjets (CHS and PUPPI jets) and subjet subjet indices. Needed for lepton subjet fraction.
- Bug fix in the jetToolBox

## Instructions

 * Make a new CMSSW area:
```
setenv SCRAM_ARCH slc6_amd64_gcc530 (or in bash: export SCRAM_ARCH=slc6_amd64_gcc530)
cmsrel CMSSW_8_0_10_patch2
cd CMSSW_8_0_10_patch2/src
cmsenv
```
 * Temporary checkouts:
```
```
 * Clone the github repository
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_8_0_X_V1
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_763
```
 * Compile
```
scram b -j 10
```

## Running

The python configuration file for cmsRun is B2GAnaFW/test/b2gedmntuples_cfg.py. It runs on the miniAOD data tier and produces an EDM-ntuple.

The configuration file contians a header explaining usage. Do
```
cd Analysis/B2GAnaFW/test
python b2gedmntuples_cfg.py 
```

for running instructions. 

## CRAB submission

The script test/submit_all.py can be used for mass submission of crab jobs. 

To run, prepare a text file CRAB/tosubmit.txt with dataset names of samples to submit.

Example usage: 

```
python submit_all.py -c b2gedmntuples_cfg.py -f <your_dataset_file> -s "T3_US_FNALLPC" -p "DataProcessing=MC_MiniAODv2_80X" -d B2GEDMNTuples_80x_V1p0 -o "/store/group/lpctlbsm/B2GAnaFW_80X_V1p0" -v RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0 -i 'Summer15_25nsV7*'
```

See all options with ```python submit_all.py --help```
