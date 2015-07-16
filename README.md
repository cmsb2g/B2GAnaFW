B2GAnaFW
========

Analysis framework for Beyond Two Generations (B2G) Physics Analysis Group (PAG) of the Compact Muon Solenoid (CMS) Experiment

Checkout Instructions
=====================

Make a new CMSSW area in 7_4_7

$ cmsrel CMSSW_7_4_7_patch2

$ cd CMSSW_7_4_7_patch2/src

$ cmsenv

(check out any CMSSW packages needed)


Clone the github repository

$ git clone https://github.com/cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_7_4_X_V4

Compile

$ scram b -j 10

Running
=======

The python configuration file for cmsRun is B2GAnaFW/test/b2gedmntuples_cfg.py. It runs on the miniAOD data tier and produces an EDM-ntuple.

The configuration file contians a header explaining usage. Here is an example running on a phys14 top sample:

$ cmsRun b2gedmntuples_cfg.py LHE=0
