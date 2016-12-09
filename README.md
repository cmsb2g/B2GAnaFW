# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a development branch of the B2G EDMNtuples to be used for 2016 re-reco (Runs2016B-G) and PromtReco (Run2016H) data, and the Summer16 MC

## Instructions

 * Make a new CMSSW area:
```
setenv SCRAM_ARCH slc6_amd64_gcc530 
###OR in bash: export SCRAM_ARCH=slc6_amd64_gcc530)
cmsrel CMSSW_8_0_24_patch1
cd CMSSW_8_0_24_patch1/src
cmsenv
```

 * Mirror for github
```
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily  (this is not needed, it just boost your git clone)
git cms-init
```

 * Necessary for the running the new Electron IDs (cut-based and MVA)
```
git cms-merge-topic ikrav:egm_id_80X_v2

scram b -j 10

cd $CMSSW_BASE/external/slc6_amd64_gcc530

git clone git@github.com:ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data

cd data/RecoEgamma/ElectronIdentification/data

git checkout egm_id_80X_v1

cd $CMSSW_BASE/src
```

 * Apply latest Run II MET filters that are not in MINIAOD
```
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter
```

 * Temporary checkouts:
```
```

 * Clone the github repository
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_8_0_X_V2

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

Example usage: 

```
python submit_all.py -c b2gedmntuples_cfg.py -f CRAB/tosubmit.txt -s T2_CH_CERN -p "DataProcessing=MC_MiniAODv2_80X_reHLT" -o "/store/group/phys_b2g/" -v B2GAnaFW_80X_V2p1 -i 'JECs/*.db'
python submit_all.py -c b2gedmntuples_cfg.py -f CRAB/tosubmit.txt -s T2_CH_CERN -p "DataProcessing=Data_80X" -o "/store/group/phys_b2g/" -v B2GAnaFW_80X_V2p1 -i 'JECs/*.db' -l "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt"
```
Note that the ```-i``` option is not needed if the JECs are taken from the global tag, specified using the option "DataProcessing" (recommended).

See all options with ```python submit_all.py --help```
