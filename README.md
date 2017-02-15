# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a development branch of the B2G EDMNtuples to be used for 2016 re-reco (Runs2016B-G) and PromtReco (Run2016H) data, and the Summer16 MC

## Instructions

### Working release
 * Make a new CMSSW area:
```
setenv SCRAM_ARCH slc6_amd64_gcc530 ; ###chs/ tcsh 

export SCRAM_ARCH=slc6_amd64_gcc530 ; ### bash

cmsrel CMSSW_8_0_26_patch1

cd CMSSW_8_0_26_patch1/src

cmsenv
```
 * Mirror for github
```
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily  (this is not needed, it just boost your git clone)
git cms-init
```

### Temporary checkouts:
 * For MET significance in the data
```
git cms-merge-topic cms-met:METRecipe_8020
```
 * Apply latest Run II MET filters that are not in MINIAOD
```
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter
```
 * For running the new double b tagger training
```
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV4-WithWeightFiles-v1_from-CMSSW_8_0_21
```
 * Compile 
```
scram b -j 10
```

### Clone the github repositories
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b v8.0.x_v2.4
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X_V3
```
 * Compile (patience please!)
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

The script B2GAnaFW/test/submit_all.py can be used for mass submission of crab jobs. 

To run, prepare a text file CRAB/tosubmit.txt with dataset names of samples to submit.

 * Example usage: 

```
python submit_all.py -c b2gedmntuples_cfg.py -f CRAB/tosubmit.txt -s T2_CH_CERN -p "DataProcessing=MC_MiniAODv2_80X_Summer16" -o "/store/group/phys_b2g/" -v B2GAnaFW_80X_V2p4 -i 'JECs/*.db'
```

 * For data please note that there are different JECs for different run periods. Please refer to test/b2gedmntuples_cfg.py#L6-L13 for the full list of settings for the switch "DataProcessing". For instance, to run on Run2016G do
```
python submit_all.py -c b2gedmntuples_cfg.py -f CRAB/tosubmit.txt -s T2_CH_CERN -p "DataProcessing=Data_80X_Run2016G_23Sep2016" -o "/store/group/phys_b2g/" -v B2GAnaFW_80X_V2p4 -i 'JECs/*.db' -l "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt" usePrivateSQLite=True
```

Note that the ```-i``` option is not needed if the JECs are taken from the global tag, specified using the option "DataProcessing".

See all options with ```python submit_all.py --help```
