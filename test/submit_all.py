#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
import glob
from optparse import OptionParser

def make_list(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -p PARAMS -i INPUTS -d DIR -v VERSION -f FILE')

    parser = OptionParser(usage=usage)    
    parser.add_option("-c", "--config", dest="cfg", default="b2gedmntuples_cfg.py",
        help=("The crab script you want to submit "),
        metavar="CONFIG")
    parser.add_option("-i", "--inputFiles", 
        type='string',
        #action='callback',
        #callback=make_list,
        dest='inputFiles', 
        help=("Input files that need to be shipped with the job"),
        metavar="INPUTS")
    parser.add_option("-p", "--pyCfgParams",
        type='string',
        action='callback',
        callback=make_list,
        dest='pyCfgParams', 
        help=("input parameters for config file"),
        metavar="PARAMS")
    parser.add_option("-l", "--lumiMask", dest="lumiMask",
        default='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt',
        help=("The JSON file containing good lumi list"),
        metavar="LUMI")
    parser.add_option("-d", "--dir", dest="dir", default="B2GEDMNTuples",
        help=("The crab directory you want to use "),
        metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    parser.add_option("-v", "--version", dest="version", default="B2GAnaFW_80X_V1p2",
        help=("B2GAnaFW version"),
        metavar="VERSION")
    parser.add_option("-s", "--storageSite", dest="storageSite", #default="T3_US_FNALLPC",
        help=("B2GAnaFW version"),
        metavar="VERSION")
    parser.add_option("-o", "--outLFNDirBase", dest="outLFNDirBase", 
        help=("B2GAnaFW version"),
        metavar="LFN")
    (options, args) = parser.parse_args()


    if options.cfg == None or options.inputFiles == None or options.pyCfgParams == None or options.dir == None or options.datasets == None or options.version == None or options.storageSite == None:
        parser.error(usage)
    
    return options
    

def main():

    options = getOptions()

    from WMCore.Configuration import Configuration
    config = Configuration()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.section_("General")
    config.General.workArea = options.dir + '_' + options.version
    config.General.transferLogs = True

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.cfg
    inFiles = glob.glob( options.inputFiles )
    config.JobType.inputFiles = inFiles #options.inputFiles
    config.JobType.pyCfgParams = options.pyCfgParams
    
    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.splitting = ''
    config.Data.unitsPerJob = 1
    config.Data.ignoreLocality = False
    config.Data.publication = True    
    config.Data.publishDBS = 'phys03'
    config.Data.outputDatasetTag = options.version
    if options.outLFNDirBase and not options.outLFNDirBase.isspace(): 
      config.Data.outLFNDirBase = options.outLFNDirBase + options.dir + '_' + options.version
    
    config.section_("Site")
    config.Site.storageSite = options.storageSite

    print 'Using config ' + options.cfg
    print 'Writing to directory ' + options.dir
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute commend'
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        if (len(s)==0 or s[0][0]=='#'): continue
        s = ijob.rstrip()
        jobs.append( s )
        print '  --> added ' + s
        
    for ijob, job in enumerate(jobs) :

        ptbin = job.split('/')[1]
        cond = job.split('/')[2]
        datatier = job.split('/')[3]
        requestname = ptbin + '_' + cond
        if len(requestname) > 100: requestname = ''.join((requestname[:100-len(requestname)]).split('_')[:-1])
        config.General.requestName = requestname
        config.Data.inputDataset = job
        if datatier == 'MINIAODSIM': 
          config.Data.splitting = 'FileBased'
        elif datatier == 'MINIAOD': 
          config.Data.splitting = 'LumiBased'
          config.Data.lumiMask = options.lumiMask
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config
        try :
            from multiprocessing import Process
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
            #submit(config)
        except :
            print 'Not submitted.'



if __name__ == '__main__':
    main()            
