# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a prroduction version of the B2G EDMNtuples. Updates since last production tag ``v8.0.x_v1.2``:
- New JECs

## Instructions

 * Make a new CMSSW area:
```
setenv SCRAM_ARCH slc6_amd64_gcc530 (or in bash: export SCRAM_ARCH=slc6_amd64_gcc530)
cmsrel CMSSW_8_0_16
cd CMSSW_8_0_16/src
cmsenv
```
 * Temporary checkouts:
```
```
 * Clone the github repository
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b v8.0.x_v2.0
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
Prepare a text file CRAB/tosubmit.txt containing list of MiniAOD samples to process.

python submit_all.py -f CRAB/tosubmit.txt -i 'Spring16_25nsV6_MC.db' -s T2_CH_CERN -p "DataProcessing=MC_MiniAODv2_80X_reHLT" -o "/store/group/phys_b2g" -d B2GAnaFW -v v80x_v2p0

```

See all options with ```python submit_all.py --help```
