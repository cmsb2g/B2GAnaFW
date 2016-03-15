# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This version is used to produced B2G EDMNtuples for 76X samples, storing jet SV mass. Needs modified version of JetToolBox, see below:
N.B.: NOT stable version, testing phase.

## Instructions

 * Make a new CMSSW area:
```
cmsrel CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src
cmsenv
```
 * Clone the github repository
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_7_6_X_V1
git clone git@github.com:dmajumder/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_763
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

