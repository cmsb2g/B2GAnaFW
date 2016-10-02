# B2GAnaFW

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

## Version

This is a development branch of the B2G EDMNtuples to be used for 2016 Data and Spring16 MC

## Instructions

 * Make a new CMSSW area:
```
setenv SCRAM_ARCH slc6_amd64_gcc530 (or in bash: export SCRAM_ARCH=slc6_amd64_gcc530)
cmsrel CMSSW_8_0_16
cd CMSSW_8_0_16/src
cmsenv

```
 * Mirror for github
```
setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily

git cms-init

```
 * Necessary for the VID tool and the EGamma Ids (and a temporary fix for loose cut based ID)
```
git cms-merge-topic ikrav:egm_id_80X_v1
sed -i "51s/0.0477/0.00477/" RecoEgamma/ElectronIdentification/python/Identification/cutBasedElectronID_Summer16_80X_V1_cff.py
```
 * Apply latest Run II MET filters that are not in MINIAOD
```
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
```
 * Temporary checkouts:
```
```
 * Clone the github repository
```
git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_8_0_X_V2
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
