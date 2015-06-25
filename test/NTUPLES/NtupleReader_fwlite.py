#! /usr/bin/env python

#CONFIGURATION

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--files', type='string', action='store',
                  dest='files',
                  help='Input files')

parser.add_option('--outname', type='string', action='store',
                  default='outplots.root',
                  dest='outname',
                  help='Name of output file')

parser.add_option('--verbose', action='store_true',
                  default=False,
                  dest='verbose',
                  help='Print debugging info')

parser.add_option('--maxevents', type='int', action='store',
                  default=-1,
                  dest='maxevents',
                  help='Number of events to run. -1 is all events')

parser.add_option('--maxjets', type='int', action='store',
                  default=999,
                  dest='maxjets',
                  help='Number of jets to plot. To plot all jets, set to a big number like 999')


parser.add_option('--bdisc', type='string', action='store',
                  default='combinedInclusiveSecondaryVertexV2BJetTags',
                  dest='bdisc',
                  help='Name of output file')


parser.add_option('--bDiscMin', type='float', action='store',
                  default=0.679,
                  dest='bDiscMin',
                  help='Minimum b discriminator')

parser.add_option('--minMuonPt', type='float', action='store',
                  default=30.,
                  dest='minMuonPt',
                  help='Minimum PT for muons')

parser.add_option('--maxMuonEta', type='float', action='store',
                  default=2.1,
                  dest='maxMuonEta',
                  help='Maximum muon pseudorapidity')

parser.add_option('--minElectronPt', type='float', action='store',
                  default=30.,
                  dest='minElectronPt',
                  help='Minimum PT for electrons')

parser.add_option('--maxElectronEta', type='float', action='store',
                  default=2.5,
                  dest='maxElectronEta',
                  help='Maximum electron pseudorapidity')


parser.add_option('--minAK4Pt', type='float', action='store',
                  default=30.,
                  dest='minAK4Pt',
                  help='Minimum PT for AK4 jets')

parser.add_option('--maxAK4Rapidity', type='float', action='store',
                  default=2.4,
                  dest='maxAK4Rapidity',
                  help='Maximum AK4 rapidity')

parser.add_option('--minAK8Pt', type='float', action='store',
                  default=400.,
                  dest='minAK8Pt',
                  help='Minimum PT for AK8 jets')

parser.add_option('--maxAK8Rapidity', type='float', action='store',
                  default=2.4,
                  dest='maxAK8Rapidity',
                  help='Maximum AK8 rapidity')


parser.add_option('--minMassCut', type='float', action='store',
                  default=50.,
                  dest='minMassCut',
                  help='Minimum Mass Pairing Cut for CMS Combined Tagger')

parser.add_option('--mAK8TrimmedCut', type='float', action='store',
                  default=100.,
                  dest='mAK8TrimmedCut',
                  help='Trimmed mass Cut for CMS Combined Tagger')

parser.add_option('--tau32Cut', type='float', action='store',
                  default=0.6,
                  dest='tau32Cut',
                  help='Tau3 / Tau2 n-subjettiness cut for CMS Combined Tagger')


(options, args) = parser.parse_args()
argv = []


#FWLITE STUFF

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
ROOT.gROOT.Macro("rootlogon.C")
from leptonic_nu_z_component import solve_nu_tmass, solve_nu
import copy


#muon labels
h_muPt = Handle("std::vector<float>")
l_muPt = ("muons" , "muPt")
h_muEta = Handle("std::vector<float>")
l_muEta = ("muons" , "muEta")
h_muPhi = Handle("std::vector<float>")
l_muPhi = ("muons" , "muPhi")
h_muTight = Handle("std::vector<float>")
l_muTight = ("muons" , "muIsTightMuon" )
h_muLoose = Handle("std::vector<float>")
l_muLoose = ("muons" , "muIsLooseMuon" )
h_muMass = Handle("std::vector<float>")
l_muMass = ("muons" , "muMass")
h_muDz = Handle("std::vector<float>")
l_muDz = ("muons", "muDz")

h_muKey = Handle("std::vector<std::vector<int> >")
l_muKey = ("muonKeys")

#electron label and handles
h_elPt = Handle("std::vector<float>")
l_elPt = ("electrons" , "elPt")
h_elEta = Handle("std::vector<float>")
l_elEta = ("electrons" , "elEta")
h_elPhi = Handle("std::vector<float>")
l_elPhi = ("electrons" , "elPhi")
h_elTight = Handle("std::vector<float>")
l_elTight = ("electrons" , "elisTight" )
h_elLoose = Handle("std::vector<float>")
l_elLoose = ("electrons" , "elisLoose" )
h_eldEtaIn = Handle("std::vector<float>")
l_eldEtaIn = ( "electrons" , "eldEtaIn" )
h_eldPhiIn = Handle("std::vector<float>")
l_eldPhiIn = ( "electrons" , "eldPhiIn" )
h_elHoE = Handle("std::vector<float>")
l_elHoE = ( "electrons" , "elHoE" )
h_elfull5x5siee = Handle("std::vector<float>")
l_elfull5x5siee = ( "electrons" , "elfull5x5siee")
h_elE = Handle("std::vector<float>")
l_elE = ( "electrons" , "elE" )
h_elD0 = Handle("std::vector<float>")
l_elD0 = ( "electrons" , "elD0" )
h_elDz = Handle("std::vector<float>")
l_elDz = ( "electrons" , "elDz")
h_elIso03 = Handle("std::vector<float>")
l_elIso03 = ( "electrons" , "elIso03" )
h_elisVeto = Handle("std::vector<float>")
l_elisVeto = ( "electrons" , "elisVeto" )
h_elhasMatchedConVeto = Handle("std::vector<float>")
l_elhasMatchedConVeto = ( "electrons" , "elhasMatchedConVeto" )
h_elooEmooP = Handle("std::vector<float>")
l_elooEmooP = ( "electrons" , "elooEmooP" )
h_elMass = Handle("std::vector<float>")
l_elMass = ( "electrons" , "elMass" )
h_elscEta = Handle("std::vector<float>")
l_elscEta = ( "electrons" , "elscEta" )

h_elKey = Handle("std::vector<std::vector<int> >")
l_elKey = ( "electronKeys" )

#AK4 Jet Label and Handles
h_jetsAK4Pt = Handle("std::vector<float>")
l_jetsAK4Pt = ("jetsAK4" , "jetAK4Pt") #
h_jetsAK4Eta = Handle("std::vector<float>")
l_jetsAK4Eta = ("jetsAK4" , "jetAK4Eta")
h_jetsAK4Phi = Handle("std::vector<float>")
l_jetsAK4Phi = ("jetsAK4" , "jetAK4Phi")
h_jetsAK4Mass = Handle("std::vector<float>")
l_jetsAK4Mass = ("jetsAK4" , "jetAK4Mass")
h_jetsAK4Energy = Handle("std::vector<float>")
l_jetsAK4Energy = ("jetsAK4" , "jetAK4E") #check! is this energy?
h_jetsAK4JEC = Handle("std::vector<float>")
l_jetsAK4JEC = ("jetsAK4" , "jetAK4jecFactor0") 
h_jetsAK4CSV = Handle("std::vector<float>")
l_jetsAK4CSV = ("jetsAK4" , "jetAK4CSV")
h_jetsAK4NumDaughters = Handle("std::vector<float>")
l_jetsAK4NumDaughters = ( "jetsAK4" , "jetAK4numberOfDaughters" )
h_jetsAK4Area = Handle("std::vector<float>")
l_jetsAK4Area = ( "jetsAK4" , "jetAK4jetArea" )

h_NPV = Handle("std::int")
l_NPV = ( "eventUserData" , "npv" )

h_jetsAK4Keys = Handle("std::vector<std::vector<int> >")
l_jetsAK4Keys = ( "jetKeysAK4" , "" )


h_jetsAK4nHadEnergy = Handle("std::vector<float>")
l_jetsAK4nHadEnergy = ("jetsAK4" , "jetAK4neutralHadronEnergy")
h_jetsAK4nEMEnergy = Handle("std::vector<float>")
l_jetsAK4nEMEnergy = ("jetsAK4" , "jetAK4neutralEmEnergy")
h_jetsAK4HFHadronEnergy = Handle("std::vector<float>")
l_jetsAK4HFHadronEnergy = ("jetsAK4" , "jetAK4HFHadronEnergy")
h_jetsAK4cHadEnergy = Handle("std::vector<float>")
l_jetsAK4cHadEnergy = ("jetsAK4" , "jetAK4chargedHadronEnergy")
h_jetsAK4cEMEnergy = Handle("std::vector<float>")
l_jetsAK4cEMEnergy = ("jetsAK4" , "jetAK4chargedEmEnergy")
h_jetsAK4numDaughters = Handle("std::vector<float>")
l_jetsAK4numDaughters = ("jetsAK4" , "jetAK4numberOfDaughters")
h_jetsAK4cMultip = Handle("std::vector<float>")
l_jetsAK4cMultip = ("jetsAK4" , "jetAK4chargedMultiplicity")
h_jetsAK4Y = Handle("std::vector<float>")
l_jetsAK4Y = ("jetsAK4" , "jetAK4Y")
h_jetsAK8Area = Handle("std::vector<float>")
l_jetsAK8Area = ( "jetsAK8" , "jetAK8jetArea" )


h_rho = Handle("double")
l_rho = ("fixedGridRhoFastjetAll", "")

#MET label and Handles
h_metPt = Handle("std::vector<float>")
l_metPt = ("met" , "metPt")
h_metPx = Handle("std::vector<float>")
l_metPx = ("met" , "metPx")
h_metPy = Handle("std::vector<float>")
l_metPy = ("met" , "metPy")
h_metPhi = Handle("std::vector<float>")
l_metPhi = ("met" , "metPhi")

#AK8 Jets label and Handles
h_jetsAK8Pt = Handle("std::vector<float>")
l_jetsAK8Pt = ("jetsAK8" , "jetAK8Pt") #
h_jetsAK8Eta = Handle("std::vector<float>")
l_jetsAK8Eta = ("jetsAK8" , "jetAK8Eta")
h_jetsAK8Phi = Handle("std::vector<float>")
l_jetsAK8Phi = ("jetsAK8" , "jetAK8Phi")
h_jetsAK8Mass = Handle("std::vector<float>")
l_jetsAK8Mass = ("jetsAK8" , "jetAK8Mass")
h_jetsAK8Energy = Handle("std::vector<float>")
l_jetsAK8Energy = ("jetsAK8" , "jetAK8E") #check! is this energy?
h_jetsAK8JEC = Handle("std::vector<float>")
l_jetsAK8JEC = ("jetsAK8" , "jetAK8jecFactor0")
h_jetsAK8Y = Handle("std::vector<float>")
l_jetsAK8Y = ("jetsAK8" , "jetAK8Y")
#h_jetsAK8CSV = Handle("std::vector<float>")
#l_jetsAK8CSV = ("jetsAK8" , "jetAK8CSV")


h_jetsAK8TrimMass = Handle("std::vector<float>")
l_jetsAK8TrimMass = ("jetsAK8", "jetAK8trimmedMass" )
h_jetsAK8PrunMass = Handle("std::vector<float>")
l_jetsAK8PrunMass = ("jetsAK8", "jetAK8prunedMass" )
h_jetsAK8FiltMass = Handle("std::vector<float>")
l_jetsAK8FiltMass = ("jetsAK8", "jetAK8filteredMass" )
h_jetsAK8Tau1 = Handle("std::vector<float>")
l_jetsAK8Tau1 = ("jetsAK8", "jetAK8tau1" )
h_jetsAK8Tau2 = Handle("std::vector<float>")
l_jetsAK8Tau2 = ("jetsAK8", "jetAK8tau2" )
h_jetsAK8Tau3 = Handle("std::vector<float>")
l_jetsAK8Tau3 = ("jetsAK8", "jetAK8tau3" )
h_jetsAK8nSubJets = Handle("std::vector<float>")
l_jetsAK8nSubJets = ("jetsAK8", "jetAK8nSubJets" )
h_jetsAK8minmass = Handle("std::vector<float>")
l_jetsAK8minmass = ("jetsAK8", "jetAK8minmass" )

h_jetsAK8VSubjetIndex0 = Handle("std::vector<float>")
l_jetsAK8VSubjetIndex0 = ("jetsAK8", "jetAK8vSubjetIndex0")
h_jetsAK8VSubjetIndex1 = Handle("std::vector<float>")
l_jetsAK8VSubjetIndex1 = ("jetsAK8", "jetAK8vSubjetIndex1")

h_jetsAK8TopSubjetIndex0 = Handle("std::vector<float>")
l_jetsAK8TopSubjetIndex0 = ("jetsAK8", "jetAK8topSubjetIndex0")
h_jetsAK8TopSubjetIndex1 = Handle("std::vector<float>")
l_jetsAK8TopSubjetIndex1 = ("jetsAK8", "jetAK8topSubjetIndex1")
h_jetsAK8TopSubjetIndex2 = Handle("std::vector<float>")
l_jetsAK8TopSubjetIndex2 = ("jetsAK8", "jetAK8topSubjetIndex2")
h_jetsAK8TopSubjetIndex3 = Handle("std::vector<float>")
l_jetsAK8TopSubjetIndex3 = ("jetsAK8", "jetAK8topSubjetIndex3")



h_subjetsAK8BDisc = Handle( "std::vector<float>")
l_subjetsAK8BDisc = ("subjetsCmsTopTag", "subjetCmsTopTagCSV")
h_subjetsAK8Pt = Handle( "std::vector<float>")
l_subjetsAK8Pt = ("subjetsCmsTopTag", "subjetCmsTopTagPt")
h_subjetsAK8Eta = Handle( "std::vector<float>")
l_subjetsAK8Eta = ("subjetsCmsTopTag", "subjetCmsTopTagEta")
h_subjetsAK8Phi = Handle( "std::vector<float>")
l_subjetsAK8Phi = ("subjetsCmsTopTag", "subjetCmsTopTagPhi")
h_subjetsAK8Mass = Handle( "std::vector<float>")
l_subjetsAK8Mass = ("subjetsCmsTopTag", "subjetCmsTopTagMass")



#HISTOGRAMS

f = ROOT.TFile(options.outname, "RECREATE")
f.cd()

h_mttbar = ROOT.TH1F("h_mttbar", ";m_{t#bar{t}} (GeV)", 200, 0, 6000)
h_mttbar_true = ROOT.TH1F("h_mttbar_true", "True m_{t#bar{t}};m_{t#bar{t}} (GeV)", 200, 0, 6000)

h_ptLep = ROOT.TH1F("h_ptLep", "Lepton p_{T};p_{T} (GeV)", 100, 0, 1000)
h_etaLep = ROOT.TH1F("h_etaLep", "Lepton #eta;p_{T} (GeV)#eta", 100, 0, ROOT.TMath.TwoPi() )
h_met = ROOT.TH1F("h_met", "Missing p_{T};p_{T} (GeV)", 100, 0, 1000)
h_ptRel = ROOT.TH1F("h_ptRel", "p_{T}^{REL};p_{T}^{REL} (GeV)", 100, 0, 100)
h_dRMin = ROOT.TH1F("h_dRMin", "#Delta R_{MIN};#Delta R_{MIN}", 100, 0, 5.0)
h_2DCut = ROOT.TH2F("h_2DCut", "2D Cut;#Delta R;p_{T}^{REL}", 20, 0, 5.0, 20, 0, 100 )

h_ptAK4 = ROOT.TH1F("h_ptAK4", "AK4 Jet p_{T};p_{T} (GeV)", 300, 0, 3000)
h_etaAK4 = ROOT.TH1F("h_etaAK4", "AK4 Jet #eta;#eta", 120, -6, 6)
h_yAK4 = ROOT.TH1F("h_yAK4", "AK4 Jet Rapidity;y", 120, -6, 6)
h_phiAK4 = ROOT.TH1F("h_phiAK4", "AK4 Jet #phi;#phi (radians)",100,-3.14, 3.14)#-ROOT.Math.Pi(),ROOT.Math.Pi())
h_mAK4 = ROOT.TH1F("h_mAK4", "AK4 Jet Mass;Mass (GeV)", 100, 0, 1000)
h_bdiscAK4 = ROOT.TH1F("h_bdiscAK4", "AK4 b discriminator;b discriminator", 100, 0, 1.0)

h_ptAK8 = ROOT.TH1F("h_ptAK8", "AK8 Jet p_{T};p_{T} (GeV)", 300, 0, 3000)
h_etaAK8 = ROOT.TH1F("h_etaAK8", "AK8 Jet #eta;#eta", 120, -6, 6)
h_yAK8 = ROOT.TH1F("h_yAK8", "AK8 Jet Rapidity;y", 120, -6, 6)
h_phiAK8 = ROOT.TH1F("h_phiAK8", "AK8 Jet #phi;#phi (radians)",100,-3.14, 3.14)#ROOT.Math.Pi(),ROOT.Math.Pi())
h_mAK8 = ROOT.TH1F("h_mAK8", "AK8 Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mprunedAK8 = ROOT.TH1F("h_mprunedAK8", "AK8 Pruned Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mfilteredAK8 = ROOT.TH1F("h_mfilteredAK8", "AK8 Filtered Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mtrimmedAK8 = ROOT.TH1F("h_mtrimmedAK8", "AK8 Trimmed Jet Mass;Mass (GeV)", 100, 0, 1000)
h_minmassAK8 = ROOT.TH1F("h_minmassAK8", "AK8 CMS Top Tagger Min Mass Paring;m_{min} (GeV)", 100, 0, 1000)
h_nsjAK8 = ROOT.TH1F("h_nsjAK8", "AK8 CMS Top Tagger N_{subjets};N_{subjets}", 5, 0, 5)
h_tau21AK8 = ROOT.TH1F("h_tau21AK8", "AK8 Jet #tau_{2} / #tau_{1};Mass#tau_{21}", 100, 0, 1.0)
h_tau32AK8 = ROOT.TH1F("h_tau32AK8", "AK8 Jet #tau_{3} / #tau_{2};Mass#tau_{32}", 100, 0, 1.0)



#JET CORRECTIONS

ROOT.gSystem.Load('libCondFormatsJetMETObjects')
#jecParStrAK4 = ROOT.std.string('JECs/PHYS14_25_V2_AK4PFchs.txt')
#jecUncAK4 = ROOT.JetCorrectionUncertainty( jecParStrAK4 )
#jecParStrAK8 = ROOT.std.string('JECs/PHYS14_25_V2_AK8PFchs.txt')
#jecUncAK8 = ROOT.JetCorrectionUncertainty( jecParStrAK8 )

print 'Getting L3 for AK4'
L3JetParAK4  = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L3Absolute_AK4PFchs.txt");
print 'Getting L2 for AK4'
L2JetParAK4  = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L2Relative_AK4PFchs.txt");
print 'Getting L1 for AK4'
L1JetParAK4  = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L1FastJet_AK4PFchs.txt");
# for data only :
#ResJetParAK4 = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L2L3Residual_AK4PFchs.txt");

print 'Getting L3 for AK8'
L3JetParAK8  = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L3Absolute_AK8PFchs.txt");
print 'Getting L2 for AK8'
L2JetParAK8  = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L2Relative_AK8PFchs.txt");
print 'Getting L1 for AK8'
L1JetParAK8  = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L1FastJet_AK8PFchs.txt");
# for data only :
#ResJetParAK8 = ROOT.JetCorrectorParameters("JECs/PHYS14_25_V2_L2L3Residual_AK8PFchs.txt"); 


#  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
vParJecAK4 = ROOT.vector('JetCorrectorParameters')()
vParJecAK4.push_back(L1JetParAK4)
vParJecAK4.push_back(L2JetParAK4)
vParJecAK4.push_back(L3JetParAK4)
# for data only :
#vParJecAK4.push_back(ResJetPar)

ak4JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK4)

vParJecAK8 = ROOT.vector('JetCorrectorParameters')()
vParJecAK8.push_back(L1JetParAK8)
vParJecAK8.push_back(L2JetParAK8)
vParJecAK8.push_back(L3JetParAK8)
# for data only :
#vParJecAK8.push_back(ResJetPar)

ak8JetCorrector = ROOT.FactorizedJetCorrector(vParJecAK8)


#EVENT LOOP

filelist = file( options.files )
filesraw = filelist.readlines()
files = []
nevents = 0
for ifile in filesraw :
    if len( ifile ) > 2 : 
        #s = 'root://cmsxrootd.fnal.gov/' + ifile.rstrip()
        s = ifile.rstrip()
        files.append( s )
        print 'Added ' + s


# loop over files
for ifile in files :
    print 'Processing file ' + ifile
    events = Events (ifile)
    if options.maxevents > 0 and nevents > options.maxevents :
        break

    # loop over events in this file
    i = 0
    for event in events:
        if options.maxevents > 0 and nevents > options.maxevents :
            break
        i += 1
        nevents += 1

        if nevents % 1000 == 0 : 
            print '    ---> Event ' + str(nevents)
        if options.verbose :
            print '==============================================='
            print '    ---> Event ' + str(nevents)


        #VERTEX SETS
        event.getByLabel( l_NPV, h_NPV )
        NPV = h_NPV.product()[0]
        if len(h_NPV.product()) == 0 :
            if options.verbose :
                print "Event has no good primary vertex."
            continue
            

        #RHO VALUE
        event.getByLabel( l_rho, h_rho )
        if len(h_rho.product()) == 0 :
            print "Event has no rho values."
            continue
        else:
            rho = h_rho.product()[0]
            if options.verbose :
                print 'rho = {0:6.2f}'.format( rho )

        #EVENT HANDLE FILLING

        event.getByLabel ( l_muPt, h_muPt )
        event.getByLabel ( l_muEta, h_muEta )
        event.getByLabel ( l_muPhi, h_muPhi )
        event.getByLabel ( l_muTight, h_muTight )
        event.getByLabel ( l_muLoose, h_muLoose )
        event.getByLabel ( l_muMass, h_muMass ) 
        event.getByLabel ( l_muDz, h_muDz )
        event.getByLabel ( l_muKey, h_muKey )

        #Muon Selection

        goodmuonPt = []
        goodmuonEta = []
        goodmuonPhi = []
        goodmuonMass = []
        goodmuonKeys = []

        #Use MuPt as iterater due to no definite value in ntuples
        if len(h_muPt.product()) > 0:
            muonPt = h_muPt.product()
            muonEta = h_muEta.product()
            muonPhi = h_muPhi.product()
            muonTight = h_muTight.product()
            muonLoose = h_muLoose.product()
            muonMass = h_muMass.product()
            muonDz = h_muDz.product()
            muKeys = h_muKey.product()
            for i in range(0,len(muonPt)):
                if muonPt[i] > options.minMuonPt and abs(muonEta[i]) < options.maxMuonEta and muonDz[i] < 5.0 and muonTight[i] :
                    goodmuonPt.append(muonPt[i])
                    goodmuonEta.append(muonEta[i])
                    goodmuonPhi.append(muonPhi[i])
                    goodmuonMass.append(muonMass[i])
                    goodmuonKeys.append(muKeys[i])
                    if options.verbose :
                        print "muon %2d: key %4d, pt %4.1f, eta %+5.3f phi %+5.3f dz(PV) %+5.3f, POG loose id %d, tight id %d." % ( i, muKeys[i][0], muonPt[i], muonEta[i],
                                                                                                                muonPhi[i], muonDz[i], muonLoose[i], muonTight[i])

        #Electron Selection
        event.getByLabel ( l_elPt, h_elPt )
        event.getByLabel ( l_elEta, h_elEta )
        event.getByLabel ( l_elPhi, h_elPhi )
        event.getByLabel ( l_elTight, h_elTight )
        event.getByLabel ( l_elLoose, h_elLoose )
        event.getByLabel ( l_eldEtaIn, h_eldEtaIn )
        event.getByLabel ( l_eldPhiIn, h_eldPhiIn )
        event.getByLabel ( l_elHoE, h_elHoE )
        event.getByLabel ( l_elfull5x5siee, h_elfull5x5siee )
        event.getByLabel ( l_elE, h_elE )
        event.getByLabel ( l_elD0, h_elD0)
        event.getByLabel ( l_elDz, h_elDz)
        event.getByLabel ( l_elIso03, h_elIso03)
        event.getByLabel ( l_elhasMatchedConVeto, h_elhasMatchedConVeto)
        event.getByLabel ( l_elooEmooP, h_elooEmooP)
        event.getByLabel ( l_elMass, h_elMass )
        # event.getByLabel ( l_isotropy, h_isotropy)
        event.getByLabel ( l_elscEta , h_elscEta )
        event.getByLabel ( l_elKey, h_elKey )
        
        goodelectronsPt = []
        goodelectronsEta = []
        goodelectronsPhi = []
        goodelectronsMass = []
        goodelectronKeys = []



        if len(h_elPt.product()) > 0:
            electronPt = h_elPt.product()
            electronEta = h_elEta.product()
            electronPhi = h_elPhi.product()
            electronTight = h_elTight.product()
            electronLoose = h_elLoose.product()
            electronecalEnergy = h_elE.product()
            electrondEtaIn = h_eldEtaIn.product()
            electrondPhiIn = h_eldPhiIn.product()
            electronHoE=h_elHoE.product()
            electronfullsiee=h_elfull5x5siee.product()
            electronooEmooP=h_elooEmooP.product()
            electronD0 = h_elD0.product()
            electronDz = h_elDz.product()
            electronabsiso = h_elIso03.product()
            electronMass = h_elMass.product()
            electronscEta = h_elscEta.product()
            elKeys = h_elKey.product()
            passConversionVeto = h_elhasMatchedConVeto.product()
            #for i in xrange( len(electronPt.size() ) ) :
            if len(electronPt) > 0 :
                for i in range(0,len(electronPt)):
                    iePt = electronPt[i]   ### << access like this
                    ieEta = electronEta[i]
                    iePhi = electronPhi[i]
                    ieEtaIn = electrondEtaIn[i]
                    iePhiIn = electronPhi[i]
                    ietight = electronTight[i]
                    ieloose = electronLoose[i]
                    ieEcal = electronecalEnergy[i]
                    ieooEmooP = electronooEmooP[i]
                    ieD0 = electronD0[i]
                    ieDz = electronDz[i]
                    ieMass = electronMass[i]
                    pfIso = electronabsiso[i]
                    ielscEta = electronscEta[i]
                    #electronpfIso = h_elIso03.product()
            
                    #for i in range(0,len(electronPt)):
                
                
                    if iePt < iePt and abs(ieEta) < options.maxElectronEta :
                        continue

                    ieHoE = electronHoE[i]
                
                    iefull = electronfullsiee[i]
                
                    #ieabsiso = electronabsiso[i]
                
                    ieabsiso = abs(pfIso)
                        
                    #absIso = electronabsiso[i]
                    #absIso = electronpfIso[i]
                    relIso = ieabsiso / iePt
                    iepass = passConversionVeto[i]
                                
                    goodElectron = electronTight[i] > 0
                    
                    if goodElectron == True :
                        
                        goodelectronsPt.append( iePt )
                        goodelectronsEta.append( ieEta )
                        goodelectronsPhi.append( iePhi )
                        goodelectronsMass.append( ieMass )
                        goodelectronKeys.append( elKeys[i] )
                        if options.verbose :
                            print "elec %2d: keys " % (i ),
                            for ikey in elKeys[i] :
                                print " %4d" % ikey ,
                            print ""
                            print "    ---> pt %4.1f, eta %+5.3f, phi %+5.3f" % ( electronPt[i], electronEta[i], electronPhi[i])
            
                
                    
                    
        if len(goodmuonPt) + len(goodelectronsPt) != 1 :
           continue
        elif len(goodmuonPt) > 0 :
            theLepton = ROOT.TLorentzVector()
            theLepton.SetPtEtaPhiM( goodmuonPt[0],
                                    goodmuonEta[0],
                                    goodmuonPhi[0],
                                    goodmuonMass[0] )
            theLeptonObjKeys = goodmuonKeys[0]
            
        else :
            theLepton = ROOT.TLorentzVector()
            theLepton.SetPtEtaPhiM( goodelectronsPt[0],
                                    goodelectronsEta[0],
                                    goodelectronsPhi[0],
                                    goodelectronsMass[0] )
            theLeptonObjKeys = goodelectronKeys[0]

        event.getByLabel ( l_jetsAK4Pt, h_jetsAK4Pt )
        event.getByLabel ( l_jetsAK4Eta, h_jetsAK4Eta )
        event.getByLabel ( l_jetsAK4Phi, h_jetsAK4Phi )
        event.getByLabel ( l_jetsAK4Mass, h_jetsAK4Mass )
        event.getByLabel ( l_jetsAK4Energy, h_jetsAK4Energy )
        event.getByLabel ( l_jetsAK4JEC, h_jetsAK4JEC )
        event.getByLabel ( l_jetsAK4CSV, h_jetsAK4CSV )
        event.getByLabel ( l_jetsAK4NumDaughters, h_jetsAK4NumDaughters)
        event.getByLabel ( l_jetsAK4Area, h_jetsAK4Area )

        event.getByLabel ( l_jetsAK4Keys, h_jetsAK4Keys )

        event.getByLabel ( l_jetsAK4nHadEnergy, h_jetsAK4nHadEnergy)
        event.getByLabel ( l_jetsAK4nEMEnergy, h_jetsAK4nEMEnergy )
        event.getByLabel ( l_jetsAK4cHadEnergy, h_jetsAK4cHadEnergy )
        event.getByLabel ( l_jetsAK4HFHadronEnergy, h_jetsAK4HFHadronEnergy )
        event.getByLabel ( l_jetsAK4cEMEnergy, h_jetsAK4cEMEnergy )
        event.getByLabel ( l_jetsAK4numDaughters, h_jetsAK4numDaughters )
        event.getByLabel ( l_jetsAK4cMultip, h_jetsAK4cMultip )
        event.getByLabel ( l_jetsAK4Y, h_jetsAK4Y )


        
        # These will hold all of the jets we need for the selection
        ak4JetsGood = []
        # For selecting leptons, look at 2-d cut of dRMin, ptRel of
        # lepton and nearest jet that has pt > 30 GeV
        dRMin = 9999.0
        inearestJet = -1    # Index of nearest jet
        nearestJetMass = None   # Nearest jet

        if len(h_jetsAK4Pt.product()) > 0:
            AK4Pt = h_jetsAK4Pt.product()
            AK4Eta = h_jetsAK4Eta.product()
            AK4Phi = h_jetsAK4Phi.product()
            AK4Mass = h_jetsAK4Mass.product()
            AK4Energy = h_jetsAK4Energy.product()
            AK4CSV = h_jetsAK4CSV.product()
            AK4NumDaughters = h_jetsAK4NumDaughters.product()
            AK4Area = h_jetsAK4Area.product()

            AK4Keys = h_jetsAK4Keys.product()

            if options.verbose :
                print '----------------'
                print 'N AK4 keys = ' + str( len(AK4Keys)) + ', N AK4Pt = ' + str(len(AK4Pt))

            AK4JEC = h_jetsAK4JEC.product()
            AK4nHadE = h_jetsAK4nHadEnergy.product()
            AK4nEME = h_jetsAK4nEMEnergy.product()
            AK4cHadE =  h_jetsAK4cHadEnergy.product()
            AK4HFHadE = h_jetsAK4HFHadronEnergy.product()
            AK4cEME =  h_jetsAK4cEMEnergy.product()
            AK4numDaughters = h_jetsAK4numDaughters.product()
            AK4cMultip =  h_jetsAK4cMultip.product()
            AK4Y =  h_jetsAK4Y.product()    
            

        for i in range(0,len(AK4Pt)):
            #get the jets transverse energy, Eta, Phi and Mass ( essentially mom-energy 4 vector)
            #v = TLorentzVector()
            jetP4Raw = ROOT.TLorentzVector()
            jetP4Raw.SetPtEtaPhiM( AK4Pt[i], AK4Eta[i], AK4Phi[i], AK4Mass[i])

            # Get correction applied to B2G ntuples
            AK4JECFromB2GAnaFW = AK4JEC[i]
            
            # Remove the old JEC's to get raw energy
            jetP4Raw *= AK4JEC[i] 
                        
            nhf = AK4nHadE[i] / jetP4Raw.E()
            nef = AK4nEME[i] / jetP4Raw.E()
            chf = AK4cHadE[i] / jetP4Raw.E()
            cef = AK4cEME[i] / jetP4Raw.E()
            nconstituents = AK4NumDaughters[i]
            nch = AK4cMultip[i] #jet.chargedMultiplicity()
            goodJet = \
              nhf < 0.99 and \
              nef < 0.99 and \
              chf > 0.00 and \
              cef < 0.99 and \
              nconstituents > 1 and \
              nch > 0

            if options.verbose :
                print '-----'
                print 'jet index = %2d, nhf = %6.2f, nef = %6.2f, chf = %6.2f, cef = %6.2f, nconstituents = %3d, nch = %3d' % (i, nhf, nef, chf, cef, nconstituents, nch)
                print '   keys : ', 
                for j in range( 0, len(AK4Keys[i]) ) :
                    print ' %d' % ( AK4Keys[i][j] ),
                print ''

            if not goodJet :
                if options.verbose : 
                    print '   bad jet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}'.format (
                        jetP4Raw.Perp(), jetP4Raw.Rapidity(), jetP4Raw.Phi(), jetP4Raw.M()
                        )
                continue
            if options.verbose :
                print '   raw jet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}'.format (
                    jetP4Raw.Perp(), jetP4Raw.Rapidity(), jetP4Raw.Phi(), jetP4Raw.M()
                    )

            cleaned = False
            if theLepton.DeltaR(jetP4Raw) < 0.4:
                # Check all daughters of jets close to the lepton
                pfcands = int(AK4NumDaughters[i])
                for j in range(0,pfcands) :                   
                    # If any of the jet daughters matches the good lepton, remove the lepton p4 from the jet p4
                    if AK4Keys[i][j] in theLeptonObjKeys : 
                        if options.verbose :
                            print '     -----> removing lepton, pt/eta/phi = {0:6.2f},{1:6.2f},{2:6.2f}'.format(
                                theLepton.Perp(), theLepton.Eta(), theLepton.Phi()
                                )
                        jetP4Raw -= theLepton
                        cleaned = True





            ak4JetCorrector.setJetEta( jetP4Raw.Eta() )
            ak4JetCorrector.setJetPt ( jetP4Raw.Perp() )
            ak4JetCorrector.setJetE  ( jetP4Raw.E() )
            ak4JetCorrector.setJetA  ( AK4Area[i] )
            ak4JetCorrector.setRho   ( rho )
            ak4JetCorrector.setNPV   ( NPV )
            newJEC = ak4JetCorrector.getCorrection()
            jetP4 = jetP4Raw*newJEC

            if jetP4.Perp() < options.minAK4Pt or abs(jetP4.Rapidity()) > options.maxAK4Rapidity :
                if options.verbose : 
                    print '   jet failed kinematic cuts'
                continue
            dR = jetP4.DeltaR(theLepton ) 
            ak4JetsGood.append(jetP4)
            if options.verbose :
                print '   corrjet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}, bdisc = {4:6.2f}'.format (
                    jetP4.Perp(), jetP4.Rapidity(), jetP4.Phi(), jetP4.M(), AK4CSV[i] )
            if dR < dRMin :
                inearestJet = i
                nearestJetP4 = jetP4
                dRMin = dR
                nearestJetbDiscrim = AK4CSV[i]

        if inearestJet < 0 :
            if options.verbose :
                print '   no nearest jet found, skipping'
            continue
        else :
            if options.verbose :
                print '>>>>>>>> nearest jet to lepton is ' + str( inearestJet )
                print '   corrjet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}, bdisc = {4:6.2f}'.format (
                    nearestJetP4.Perp(), nearestJetP4.Rapidity(), nearestJetP4.Phi(), nearestJetP4.M(), nearestJetbDiscrim )                 

        #Get MET HERE
        event.getByLabel ( l_metPt, h_metPt )
        event.getByLabel ( l_metPx, h_metPx )
        event.getByLabel ( l_metPy, h_metPy )
        event.getByLabel ( l_metPhi, h_metPhi )
        metPt = h_metPt.product()[0]
        metPx = h_metPx.product()[0]
        metPy = h_metPy.product()[0]
        metPhi = h_metPhi.product()[0]

        theLepJet = nearestJetP4
        theLepJetBDisc = nearestJetbDiscrim

        h_ptAK4.Fill( theLepJet.Perp() )
        h_etaAK4.Fill( theLepJet.Eta() )
        h_yAK4.Fill( theLepJet.Rapidity() )
        h_mAK4.Fill( theLepJet.M() )
        h_bdiscAK4.Fill( theLepJetBDisc )

        ptRel = theLepJet.Perp( theLepton.Vect() )
        h_ptLep.Fill(theLepton.Perp())
        h_etaLep.Fill(theLepton.Eta())
        h_met.Fill(metPt)
        h_ptRel.Fill( ptRel )
        h_dRMin.Fill( dRMin )
        h_2DCut.Fill( dRMin, ptRel )

        pass2D = ptRel > 20.0 or dRMin > 0.4
        if options.verbose :
            print '>>>>>>>>>>>>>>'
            print '2d cut : dRMin = {0:6.2f}, ptRel = {1:6.2f}'.format( dRMin, ptRel )
            print '>>>>>>>>>>>>>>'
        if pass2D == False :
            continue

        ############################################
        # Get the AK8 jet away from the lepton
        ############################################





        event.getByLabel ( l_jetsAK8Eta, h_jetsAK8Eta )
        event.getByLabel ( l_jetsAK8Pt, h_jetsAK8Pt )
        event.getByLabel ( l_jetsAK8Phi, h_jetsAK8Phi )
        event.getByLabel ( l_jetsAK8Mass, h_jetsAK8Mass )
        event.getByLabel ( l_jetsAK8Energy, h_jetsAK8Energy )
        event.getByLabel ( l_jetsAK8JEC, h_jetsAK8JEC )
        event.getByLabel ( l_jetsAK8Y, h_jetsAK8Y )
        event.getByLabel ( l_jetsAK8Area, h_jetsAK8Area )

        event.getByLabel ( l_jetsAK8TrimMass, h_jetsAK8TrimMass )
        event.getByLabel ( l_jetsAK8PrunMass, h_jetsAK8PrunMass )
        event.getByLabel ( l_jetsAK8FiltMass, h_jetsAK8FiltMass )
        event.getByLabel ( l_jetsAK8Tau1, h_jetsAK8Tau1 )
        event.getByLabel ( l_jetsAK8Tau2, h_jetsAK8Tau2 )
        event.getByLabel ( l_jetsAK8Tau3, h_jetsAK8Tau3 )
        event.getByLabel ( l_jetsAK8nSubJets, h_jetsAK8nSubJets )
        event.getByLabel ( l_jetsAK8minmass, h_jetsAK8minmass )

        event.getByLabel ( l_jetsAK8TopSubjetIndex0, h_jetsAK8TopSubjetIndex0 )
        event.getByLabel ( l_jetsAK8TopSubjetIndex1, h_jetsAK8TopSubjetIndex1 )
        event.getByLabel ( l_jetsAK8TopSubjetIndex2, h_jetsAK8TopSubjetIndex2 )
        event.getByLabel ( l_jetsAK8TopSubjetIndex3, h_jetsAK8TopSubjetIndex3 )

        event.getByLabel ( l_subjetsAK8BDisc, h_subjetsAK8BDisc)
        event.getByLabel ( l_subjetsAK8Pt, h_subjetsAK8Pt)
        event.getByLabel ( l_subjetsAK8Eta, h_subjetsAK8Eta)
        event.getByLabel ( l_subjetsAK8Phi, h_subjetsAK8Phi)
        
        ak8JetsGood = []
        ak8JetsGoodTrimMass = []
        ak8JetsGoodPrunMass = []
        ak8JetsGoodFiltMass = []
        ak8JetsGoodTau1 = []
        ak8JetsGoodTau2 = []
        ak8JetsGoodTau3 = []
        ak8JetsGoodNSubJets = []
        ak8JetsGoodMinMass = []
        ak8JetsGoodTopSubjetIndex0 = []
        ak8JetsGoodTopSubjetIndex1 = []
        ak8JetsGoodTopSubjetIndex2 = []
        ak8JetsGoodTopSubjetIndex3 = []


        if len( h_jetsAK8Pt.product()) > 0 : 
            AK8Pt = h_jetsAK8Pt.product()
            AK8Eta = h_jetsAK8Eta.product()
            AK8Phi = h_jetsAK8Phi.product()
            AK8Mass = h_jetsAK8Mass.product()
            AK8Energy = h_jetsAK8Energy.product()
            AK8Y = h_jetsAK8Y.product()

            AK8JEC = h_jetsAK8JEC.product()
            AK8Area = h_jetsAK8Area.product()


            AK8TrimmedM = h_jetsAK8TrimMass.product()
            AK8PrunedM = h_jetsAK8PrunMass.product()
            AK8FilteredM = h_jetsAK8FiltMass.product()
            AK8Tau1 = h_jetsAK8Tau1.product()
            AK8Tau2 = h_jetsAK8Tau2.product()
            AK8Tau3 = h_jetsAK8Tau3.product()
            AK8nSubJets = h_jetsAK8nSubJets.product()
            AK8minmass = h_jetsAK8minmass.product()
            AK8TopSubjetIndex0 = h_jetsAK8TopSubjetIndex0.product()
            AK8TopSubjetIndex1 = h_jetsAK8TopSubjetIndex1.product()
            AK8TopSubjetIndex2 = h_jetsAK8TopSubjetIndex2.product()
            AK8TopSubjetIndex3 = h_jetsAK8TopSubjetIndex3.product()


            
        ak8SubJetsBDisc = []
        ak8SubJetsPt = []
        ak8SubJetsEta = []
        ak8SubJetsPhi = []
        
        if len( h_subjetsAK8BDisc.product() ) > 0 : 
            AK8SubJetsBDisc = h_subjetsAK8BDisc.product()
            AK8SubJetsPt = h_subjetsAK8Pt.product()
            AK8SubJetsEta = h_subjetsAK8Eta.product()
            AK8SubJetsPhi = h_subjetsAK8Phi.product()
            
        
        for i in range(0,len(AK8Pt)):

            AK8JECFromB2GAnaFW = AK8JEC[i]   
            #AK8P4Raw = TLorentzVector()
            AK8P4Raw = ROOT.TLorentzVector()
            AK8P4Raw.SetPtEtaPhiM( AK8Pt[i] , AK8Eta[i], AK8Phi[i], AK8Mass[i])
            # Remove the old JEC's to get raw energy
            AK8P4Raw *= AK8JECFromB2GAnaFW 
         
            RawAK8Energy = AK8P4Raw.Energy()

            ak8JetCorrector.setJetEta( AK8P4Raw.Eta() )
            ak8JetCorrector.setJetPt ( AK8P4Raw.Perp() )
            ak8JetCorrector.setJetE  ( AK8P4Raw.E() )
            ak8JetCorrector.setJetA  ( AK8Area[i] )
            ak8JetCorrector.setRho   ( rho )
            ak8JetCorrector.setNPV   ( NPV )
            newJEC = ak8JetCorrector.getCorrection()
            AK8P4Corr = AK8P4Raw*newJEC

            
                        

            if AK8P4Raw.Perp() < options.minAK8Pt or abs(AK8P4Raw.Rapidity()) > options.maxAK8Rapidity :
                continue
            # Only keep AK8 jets "away" from the lepton, so we do not need
            # lepton-jet cleaning here. There's no double counting. 
            dR = jetP4.DeltaR(theLepton ) 
            if dR > ROOT.TMath.Pi()/2.0 :
                ak8JetsGood.append(AK8P4Corr)
                ak8JetsGoodTrimMass.append( AK8TrimmedM[i])
                ak8JetsGoodPrunMass.append( AK8PrunedM[i])
                ak8JetsGoodFiltMass.append( AK8FilteredM[i])
                ak8JetsGoodTau1.append( AK8Tau1[i])
                ak8JetsGoodTau2.append( AK8Tau2[i])
                ak8JetsGoodTau3.append( AK8Tau3[i])
                ak8JetsGoodNSubJets.append( AK8nSubJets[i])
                ak8JetsGoodMinMass.append( AK8minmass[i] )
                ak8JetsGoodTopSubjetIndex0.append( AK8TopSubjetIndex0[i] )
                ak8JetsGoodTopSubjetIndex1.append( AK8TopSubjetIndex1[i] )
                ak8JetsGoodTopSubjetIndex2.append( AK8TopSubjetIndex2[i] )
                ak8JetsGoodTopSubjetIndex3.append( AK8TopSubjetIndex3[i] )
                                                                                             


        #Tagging
        if len(ak8JetsGood) < 1 :
            if options.verbose :
                print 'Not enough AK8 jets, skipping'
            continue

            
        nttags = 0
        tJets = []



        for i in range(0,len(ak8JetsGood)): 
            if ak8JetsGood[i].Perp() < options.minAK8Pt :
                continue

            mAK8Pruned = AK8PrunedM[i] 
            mAK8Filtered = AK8FilteredM[i] 
            mAK8Trimmed = AK8TrimmedM[i]
            # Make sure there are top tags if we want to plot them
            minMass = AK8minmass[i]
            nsubjets = AK8nSubJets[i]
            tau1 = AK8Tau1[i]  
            tau2 = AK8Tau2[i] 
            tau3 = AK8Tau3[i]
            if tau1 > 0.0001 :
                tau21 = tau2 / tau1
                h_tau21AK8.Fill( tau21 )
            else :
                h_tau21AK8.Fill( -1.0 )
            if tau2 > 0.0001 :
                tau32 = tau3 / tau2
                h_tau32AK8.Fill( tau32 )
            else :
                h_tau32AK8.Fill( -1.0 )

            h_ptAK8.Fill( ak8JetsGood[i].Perp() )
            h_etaAK8.Fill( ak8JetsGood[i].Eta() )
            h_yAK8.Fill( ak8JetsGood[i].Rapidity() )
            h_mAK8.Fill( ak8JetsGood[i].M() )
            h_mprunedAK8.Fill( ak8JetsGoodPrunMass[i] )
            h_mfilteredAK8.Fill( ak8JetsGoodFiltMass[i] )
            h_mtrimmedAK8.Fill( ak8JetsGoodTrimMass[i] )
            h_minmassAK8.Fill( ak8JetsGoodMinMass[i] )
            h_nsjAK8.Fill( ak8JetsGoodNSubJets[i] )

            sbdisc0 = -1.
            sbdisc1 = -1.
            sbdisc2 = -1.
            sbdisc3 = -1.
            spt0 = -1.
            spt1 = -1.
            spt2 = -1.
            spt3 = -1.
            seta0 = -1.
            seta1 = -1.
            seta2 = -1.
            seta3 = -1.
            sphi0 = -1.
            sphi1 = -1.
            sphi2 = -1.
            sphi3 = -1.

            for check in [ak8JetsGoodTopSubjetIndex0[i],
                           ak8JetsGoodTopSubjetIndex1[i],
                           ak8JetsGoodTopSubjetIndex2[i],
                           ak8JetsGoodTopSubjetIndex3[i]] :
                if int(check) > len(AK8SubJetsPt) :
                    print '===================='
                    print ' Catastrophic failure. Index is out of range. Setup is completely wrong.'
                    print '===================='
                    exit(1)

            if ak8JetsGoodTopSubjetIndex0[i] > -1 :
                sbdisc0 = AK8SubJetsBDisc[ int(ak8JetsGoodTopSubjetIndex0[i]) ]
                spt0    = AK8SubJetsPt[ int(ak8JetsGoodTopSubjetIndex0[i]) ]
                seta0   = AK8SubJetsEta[ int(ak8JetsGoodTopSubjetIndex0[i]) ]
                sphi0   = AK8SubJetsPhi[ int(ak8JetsGoodTopSubjetIndex0[i]) ]

            if ak8JetsGoodTopSubjetIndex1[i] > -1 :
                sbdisc1 = AK8SubJetsBDisc[ int(ak8JetsGoodTopSubjetIndex1[i]) ]
                spt1    = AK8SubJetsPt[ int(ak8JetsGoodTopSubjetIndex1[i]) ]
                seta1   = AK8SubJetsEta[ int(ak8JetsGoodTopSubjetIndex1[i]) ]
                sphi1   = AK8SubJetsPhi[ int(ak8JetsGoodTopSubjetIndex1[i]) ]

            if ak8JetsGoodTopSubjetIndex2[i] > -1 :
                sbdisc2 = AK8SubJetsBDisc[ int(ak8JetsGoodTopSubjetIndex2[i]) ]
                spt2    = AK8SubJetsPt[ int(ak8JetsGoodTopSubjetIndex2[i]) ]
                seta2   = AK8SubJetsEta[ int(ak8JetsGoodTopSubjetIndex2[i]) ]
                sphi2   = AK8SubJetsPhi[ int(ak8JetsGoodTopSubjetIndex2[i]) ]

            if ak8JetsGoodTopSubjetIndex3[i] > -1 :
                sbdisc3 = AK8SubJetsBDisc[ int(ak8JetsGoodTopSubjetIndex3[i]) ]
                spt3    = AK8SubJetsPt[ int(ak8JetsGoodTopSubjetIndex3[i]) ]
                seta3   = AK8SubJetsEta[ int(ak8JetsGoodTopSubjetIndex3[i]) ]
                sphi3   = AK8SubJetsPhi[ int(ak8JetsGoodTopSubjetIndex3[i]) ]
                
                
            if options.verbose :

                print 'Top tag candidate : '
                print '   top jet pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, m = {3:6.2f}'.format (
                    ak8JetsGood[i].Perp(), ak8JetsGood[i].Rapidity(), ak8JetsGood[i].Phi(), ak8JetsGood[i].M() )                 

                
                print '   Subjet indices : %6d %6d %6d %6d' % (ak8JetsGoodTopSubjetIndex0[i],
                                                            ak8JetsGoodTopSubjetIndex1[i],
                                                            ak8JetsGoodTopSubjetIndex2[i],
                                                            ak8JetsGoodTopSubjetIndex3[i])

                for index in [ak8JetsGoodTopSubjetIndex0[i],
                              ak8JetsGoodTopSubjetIndex1[i],
                              ak8JetsGoodTopSubjetIndex2[i],
                              ak8JetsGoodTopSubjetIndex3[i] ] :
                    if index > len(AK8SubJetsBDisc) :
                        print "Problem! Subjet indices are wrong!"
                        exit(1)


                print '   -subjet0pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, B = {3:6.2f}'.format ( spt0,seta0,sphi0,sbdisc0 )
                print '   -subjet1pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, B = {3:6.2f}'.format ( spt1,seta1,sphi1,sbdisc1 )
                print '   -subjet2pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, B = {3:6.2f}'.format ( spt2,seta2,sphi2,sbdisc2 )
                print '   -subjet3pt = {0:6.2f}, y = {1:6.2f}, phi = {2:6.2f}, B = {3:6.2f}'.format ( spt3,seta3,sphi3,sbdisc3 )                 




                                
                print '  minMass = {0:6.2f}, trimmed mass = {1:6.2f}, tau32 = {2:6.2f}'.format(
                    minMass, mAK8Trimmed, tau32
                    ),
            if minMass > options.minMassCut and mAK8Trimmed > options.mAK8TrimmedCut and tau32 < options.tau32Cut :
                nttags += 1
                tJets.append( ak8JetsGood[i] )

                if options.verbose : 
                    print '  ------>Tagged jet!'
            else :
                if options.verbose : 
                    print ''

        #KINEMATICS
        if nttags == 0 :
            if options.verbose : 
                print 'No top tags'
        else :
            
            hadTopCandP4 = tJets[0]
            
            lepTopCandP4 = None
            # Check if the nearest jet to the lepton is b-tagged
            if theLepJetBDisc < options.bDiscMin :
                if options.verbose : 
                    print 'closest jet to lepton is not b-tagged'
            else  :

                if options.verbose :
                    print 'Event is fully tagged.'
                # Get the z-component of the lepton from the W mass constraint
                
                bJetCandP4 = theLepJet
                
                nuCandP4 = ROOT.TLorentzVector(metPx, metPy ,0.0, metPt)

                solution, nuz1, nuz2 = solve_nu( vlep=theLepton, vnu=nuCandP4 )
                # If there is at least one real solution, pick it up
                if solution :
                    if options.verbose : 
                        print '--- Have a solution --- '
                    nuCandP4.SetPz(nuz1)
                else :
                    if options.verbose : 
                        print '--- No solution for neutrino z ---'
                    nuCandP4.SetPz(nuz1.real)

                lepTopCandP4 = nuCandP4 + theLepton + bJetCandP4

                ttbarCand = hadTopCandP4 + lepTopCandP4
                h_mttbar.Fill( ttbarCand.M() )
        

#CLEANUP

f.cd()
f.Write()
f.Close()
