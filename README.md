# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a prroduction version of the B2G EDMNtuples. Updates since last production tag ``v7.6.x_v1.2'':
- Adding N-subjettiness for soft drop subjets (CHS and PUPPI jets) and subjet subjet indices. Needed for lepton subjet fraction.
- Bug fix in the jetToolBox

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
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b v7.6.x_v2.0
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

python submit_all.py -f CRAB/tosubmit.txt -s T2_CH_CERN -i Fall15_25nsV2_MC.db -p "DataProcessing=MC25ns_MiniAOD_76X" -o "/store/group/phys_b2g/B2GAnaFW_76X_V2p0" -d  B2GEDMNTuples_76X_V2p0

See all options with 'python submit_all.py --help'
