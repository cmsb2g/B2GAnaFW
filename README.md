B2GAnaFW
========

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

Checkout Instructions
=====================

Make a new CMSSW area in 7_4_12_patch2

$ cmsrel CMSSW_7_4_12_patch2

$ cd CMSSW_7_4_12_patch2/src

$ cmsenv

Clone the github repository

This branch is used to run on MiniAODv2 reprocessed data: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Run2_Spring15MiniAODv2_re_miniAO

$ git clone https://github.com/cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_7_4_12

Compile

$ scram b -j 10

Running
=======

The python configuration file for cmsRun is B2GAnaFW/test/b2gedmntuples_cfg.py. It runs on the miniAOD data tier and produces an EDM-ntuple.

The configuration file contians a header explaining usage. Here is an example running on a Spring15 top sample:

$ cmsRun b2gedmntuples_cfg.py globalTag=MCRUN2_74_V9A isData=False LHE=1 DataProcessing=MC50ns
 

and here an example to run on data (PromptReco40ns):

$ cmsRun b2gedmntuples_cfg.py globalTag=74X_dataRun2_Prompt_v0 isData=Trued LHE=0 DataProcessing=PromptReco50ns
