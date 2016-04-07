# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a prroduction version of the B2G EDMNtuples. Updates since last production tag ``v7.6.x_v1.1'':
- No code change since tag ``v7.6.x_v1.1''
- Adding instructions for checking out extra CMSSW packages required to run BTagInfos

## Instructions

 * Make a new CMSSW area:
```
cmsrel CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src
cmsenv
```
 * Some temporary additional fixes for b-taging (only for CMSSW_7_6_3):
```
git cms-init
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw
git cms-merge-topic cms-btv-pog:fixTMVAEvaluatorMemoryProblem-from-CMSSW_7_6_3 
```
 * Clone the github repository
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b v7.6.x_v1.2
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

