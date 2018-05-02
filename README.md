# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a development branch of the B2G EDMNtuples to be used for 2017 data (Runs2017A-B) taken using `CMSSW_9_4_4`.

## Instructions

### Working release
 * Make a new CMSSW area:
```
setenv SCRAM_ARCH slc6_amd64_gcc630 ; ###chs/ tcsh 
export SCRAM_ARCH=slc6_amd64_gcc630 ; ### bash
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
```
 * Mirror for github
```
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily  (this is not needed, it just boost your git clone)
git cms-init
```

### Temporary checkouts:
```
```

### Clone the github repositories
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_9_4_X_V0
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_94X
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
python submit_all.py -c b2gedmntuples_cfg.py -f CRAB/tosubmit.txt -s T2_CH_CERN -p "DataProcessing=Data_92X_Run2017B" -o "/store/group/phys_b2g/" -v B2GAnaFW_92X 
```

Note that the ```-i``` option is not needed if the JECs are taken from the global tag, specified using the option "DataProcessing".

See all options with ```python submit_all.py --help```
