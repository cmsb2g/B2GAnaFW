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
    usage = ('usage: python submit_all.py -c CFG -p PYCFGPARAMS -o OUTPUTLFNDIRBASE -d DIR -v VERSION -f DATASETS -s STORAGESITE')

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
        #default='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt',
        default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt',
        help=("The JSON file containing good lumi list"),
        metavar="LUMI")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    parser.add_option("-v", "--version", dest="version", #default="B2GAnaFW_80X_V2p0",
        help=("B2GAnaFW version"),
        metavar="VERSION")
    parser.add_option("-s", "--storageSite", dest="storageSite", #default="T3_US_FNALLPC",
        help=("storage site"),
        metavar="VERSION")
    parser.add_option("-o", "--outLFNDirBase", dest="outLFNDirBase", 
        help=("EOS path for storage"),
        metavar="LFN")
    (options, args) = parser.parse_args()


    if options.cfg == None or options.pyCfgParams == None or options.datasets == None or options.version == None or options.storageSite == None or options.outLFNDirBase == None:
        parser.error(usage)
    
    return options
    

def main():

    options = getOptions()

    from CRABClient.UserUtilities import config
    config = config()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = options.version 
    config.General.transferLogs = False

    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.cfg
    if options.inputFiles != None:
      inFiles = glob.glob( options.inputFiles )
      config.JobType.inputFiles = inFiles #options.inputFiles
    config.JobType.pyCfgParams = options.pyCfgParams
    config.JobType.sendExternalFolder = True
    
    config.Data.inputDataset = None
    config.Data.splitting = ''
    config.Data.unitsPerJob = 1
    config.Data.ignoreLocality = False
    config.Data.publication = True    
    config.Data.publishDBS = 'phys03'
    config.Site.storageSite = options.storageSite

    print 'Using config ' + options.cfg
    print 'Writing to versionectory ' + options.version
    
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
	if options.outLFNDirBase and not options.outLFNDirBase.isspace(): 
          config.Data.outLFNDirBase = os.path.join(options.outLFNDirBase,options.version,ptbin,cond.split('-')[0])
        config.Data.outputDatasetTag = cond+'_'+options.version
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config
        try :
            from multiprocessing import Process
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
            submit(config)
        except :
            print 'Not submitted.'



if __name__ == '__main__':
    main()            
