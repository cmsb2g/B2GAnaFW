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



parser.add_option('--minAK8Pt', type='float', action='store',
                  default=200.,
                  dest='minAK8Pt',
                  help='Minimum PT for AK8 jets')

parser.add_option('--maxAK8Rapidity', type='float', action='store',
                  default=2.4,
                  dest='maxAK8Rapidity',
                  help='Maximum AK8 rapidity')


(options, args) = parser.parse_args()
argv = []


#FWLITE STUFF

import ROOT
import sys
from DataFormats.FWLite import Events, Handle
ROOT.gROOT.Macro("rootlogon.C")
import copy
import random


h_NPV = Handle("std::int")
l_NPV = ( "eventUserData" , "npv" )

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
h_jetsAK8SoftDropMass = Handle("std::vector<float>")
l_jetsAK8SoftDropMass = ("jetsAK8", "jetAK8softDropMass" )
h_jetsAK8Tau1 = Handle("std::vector<float>")
l_jetsAK8Tau1 = ("jetsAK8", "jetAK8tau1" )
h_jetsAK8Tau2 = Handle("std::vector<float>")
l_jetsAK8Tau2 = ("jetsAK8", "jetAK8tau2" )
h_jetsAK8Tau3 = Handle("std::vector<float>")
l_jetsAK8Tau3 = ("jetsAK8", "jetAK8tau3" )


h_jetsAK8vSubjetIndex0 = Handle("std::vector<float>")
l_jetsAK8vSubjetIndex0 = ("jetsAK8", "jetAK8vSubjetIndex0")
h_jetsAK8vSubjetIndex1 = Handle("std::vector<float>")
l_jetsAK8vSubjetIndex1 = ("jetsAK8", "jetAK8vSubjetIndex1")

h_subjetsAK8Pt = Handle( "std::vector<float>")
l_subjetsAK8Pt = ("subjetsAK8", "subjetAK8Pt")
h_subjetsAK8Eta = Handle( "std::vector<float>")
l_subjetsAK8Eta = ("subjetsAK8", "subjetAK8Eta")
h_subjetsAK8Phi = Handle( "std::vector<float>")
l_subjetsAK8Phi = ("subjetsAK8", "subjetAK8Phi")
h_subjetsAK8Mass = Handle( "std::vector<float>")
l_subjetsAK8Mass = ("subjetsAK8", "subjetAK8Mass")


#HISTOGRAMS

f = ROOT.TFile(options.outname, "RECREATE")
f.cd()

h_mjj0 = ROOT.TH1F("h_mjj0", ";m_{jj} (GeV), Stage 0", 200, 0, 6000)
h_mjj1 = ROOT.TH1F("h_mjj1", ";m_{jj} (GeV), Stage 1", 200, 0, 6000)
h_mjj2 = ROOT.TH1F("h_mjj2", ";m_{jj} (GeV), Stage 2", 200, 0, 6000)
h_mjj3 = ROOT.TH1F("h_mjj3", ";m_{jj} (GeV), Stage 3", 200, 0, 6000)

h_rhoProbe0 = ROOT.TH2F("h_rhoProbe0", "AK8 SoftDrop Jet #rho, Stage 0; (m/p_{T}R)^{2}", 100, 0, 1.0, 25, 0, 500)
h_rhoProbe1 = ROOT.TH2F("h_rhoProbe1", "AK8 SoftDrop Jet #rho, Stage 1; (m/p_{T}R)^{2}", 100, 0, 1.0, 25, 0, 500)
h_rhoProbe2 = ROOT.TH2F("h_rhoProbe2", "AK8 SoftDrop Jet #rho, Stage 2; (m/p_{T}R)^{2}", 100, 0, 1.0, 25, 0, 500)
h_rhoProbe3 = ROOT.TH2F("h_rhoProbe3", "AK8 SoftDrop Jet #rho, Stage 3; (m/p_{T}R)^{2}", 100, 0, 1.0, 25, 0, 500)


h_ptAK8 = ROOT.TH1F("h_ptAK8", "AK8 Jet p_{T};p_{T} (GeV)", 300, 0, 3000)
h_etaAK8 = ROOT.TH1F("h_etaAK8", "AK8 Jet #eta;#eta", 120, -6, 6)
h_yAK8 = ROOT.TH1F("h_yAK8", "AK8 Jet Rapidity;y", 120, -6, 6)
h_phiAK8 = ROOT.TH1F("h_phiAK8", "AK8 Jet #phi;#phi (radians)",100,-3.14, 3.14)#ROOT.Math.Pi(),ROOT.Math.Pi())
h_mAK8 = ROOT.TH1F("h_mAK8", "AK8 Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mprunedAK8 = ROOT.TH1F("h_mprunedAK8", "AK8 Pruned Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mfilteredAK8 = ROOT.TH1F("h_mfilteredAK8", "AK8 Filtered Jet Mass;Mass (GeV)", 100, 0, 1000)
h_msoftdropAK8 = ROOT.TH1F("h_msoftdropAK8", "AK8 SoftDrop Jet Mass;Mass (GeV)", 100, 0, 1000)
h_mtrimmedAK8 = ROOT.TH1F("h_mtrimmedAK8", "AK8 Trimmed Jet Mass;Mass (GeV)", 100, 0, 1000)

h_rhosoftDropAK8 = ROOT.TH2F("h_rhosoftdropAK8", "AK8 SoftDrop Jet Rho; (m/p_{T}R)^{2}", 100, 0, 1.0, 25, 0, 500)


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
            


        ############################################
        # Get the AK8 jets
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
        event.getByLabel ( l_jetsAK8SoftDropMass, h_jetsAK8SoftDropMass )
        event.getByLabel ( l_jetsAK8Tau1, h_jetsAK8Tau1 )
        event.getByLabel ( l_jetsAK8Tau2, h_jetsAK8Tau2 )
        event.getByLabel ( l_jetsAK8Tau3, h_jetsAK8Tau3 )

        event.getByLabel ( l_jetsAK8vSubjetIndex0, h_jetsAK8vSubjetIndex0)
        event.getByLabel ( l_jetsAK8vSubjetIndex1, h_jetsAK8vSubjetIndex1)

        event.getByLabel ( l_subjetsAK8Pt, h_subjetsAK8Pt)
        event.getByLabel ( l_subjetsAK8Eta, h_subjetsAK8Eta)
        event.getByLabel ( l_subjetsAK8Phi, h_subjetsAK8Phi)
        event.getByLabel ( l_subjetsAK8Mass, h_subjetsAK8Mass)

        AK8Pt = h_jetsAK8Pt.product()
        AK8Eta = h_jetsAK8Eta.product()
        AK8Phi = h_jetsAK8Phi.product()
        AK8Mass = h_jetsAK8Mass.product()
        AK8Energy = h_jetsAK8Energy.product()
        AK8Y = h_jetsAK8Y.product()

        AK8JEC = h_jetsAK8JEC.product()
        AK8Area = h_jetsAK8Area.product()


        AK8vSubjetIndex0 = h_jetsAK8vSubjetIndex0.product()
        AK8vSubjetIndex1 = h_jetsAK8vSubjetIndex1.product()            
        AK8TrimmedM = h_jetsAK8TrimMass.product()
        AK8PrunedM = h_jetsAK8PrunMass.product()
        AK8FilteredM = h_jetsAK8FiltMass.product()
        AK8SoftDropM = h_jetsAK8SoftDropMass.product()
        AK8Tau1 = h_jetsAK8Tau1.product()
        AK8Tau2 = h_jetsAK8Tau2.product()
        AK8Tau3 = h_jetsAK8Tau3.product()



        AK8vSubJetsPt = h_subjetsAK8Pt.product()
        AK8vSubJetsEta = h_subjetsAK8Eta.product()
        AK8vSubJetsPhi = h_subjetsAK8Phi.product()
        AK8vSubJetsMass = h_subjetsAK8Mass.product()

        AK8Rho = []
        goodJetIndices = []


        if options.verbose : 
            print '------- Subjets'
            for isubjet in range(0,len(AK8vSubJetsPt)):
                print '%3d, pt=%6.2f eta=%6.2f phi=%6.2f m=%6.2f' % ( isubjet, AK8vSubJetsPt[isubjet], AK8vSubJetsEta[isubjet], AK8vSubJetsPhi[isubjet], AK8vSubJetsMass[isubjet])
            print '------- AK8 jets'


        # Make plots of all jets that pass the pt and eta cuts. 
        for i in range(0,len(AK8Pt)): 
            
            sp4_0 = None
            sp4_1 = None
            ival = int(AK8vSubjetIndex0[i])            
            if ival > -1 :
                spt0    = AK8vSubJetsPt[ival]
                seta0   = AK8vSubJetsEta[ival]
                sphi0   = AK8vSubJetsPhi[ival]
                sm0   = AK8vSubJetsMass[ival]
                sp4_0 = ROOT.TLorentzVector()
                sp4_0.SetPtEtaPhiM( spt0, seta0, sphi0, sm0 )
            ival = int(AK8vSubjetIndex1[i])
            if ival > -1 :
                spt1    = AK8vSubJetsPt[ival]
                seta1   = AK8vSubJetsEta[ival]
                sphi1   = AK8vSubJetsPhi[ival]
                sm1   = AK8vSubJetsMass[ival]
                sp4_1 = ROOT.TLorentzVector()
                sp4_1.SetPtEtaPhiM( spt1, seta1, sphi1, sm1 )


            if sp4_0 == None or sp4_1 == None :
                AK8Rho.append(0.0)
                continue 
            softdrop_p4 = sp4_0 + sp4_1
            jetR = 0.8
            jetrho = softdrop_p4.M() / (softdrop_p4.Perp() * jetR)
            jetrho *= jetrho


            AK8Rho.append( jetrho )
            
            if AK8Pt[i] < options.minAK8Pt or abs(AK8Eta[i]) > options.maxAK8Rapidity :
                continue


            goodJetIndices.append(i)
            for check in [AK8vSubjetIndex0[i],AK8vSubjetIndex1[i]] :
                if int(check) > len(AK8vSubJetsPt) :
                    print '===================='
                    print ' Catastrophic failure. Index is out of range. Setup is completely wrong.'
                    print '===================='
                    exit(1)



            if options.verbose : 
                print '%3d, pt=%6.2f eta=%6.2f phi=%6.2f m=%6.2f sdm=%6.2f, sj0=%3d sj1=%3d' % ( i, AK8Pt[i], AK8Eta[i], AK8Phi[i], AK8Mass[i], AK8SoftDropM[i], AK8vSubjetIndex0[i], AK8vSubjetIndex1[i] )
            mAK8Pruned = AK8PrunedM[i] 
            mAK8Filtered = AK8FilteredM[i] 
            mAK8Trimmed = AK8TrimmedM[i]
            mAK8SoftDrop = AK8SoftDropM[i]
            tau1 = AK8Tau1[i]  
            tau2 = AK8Tau2[i] 
            tau3 = AK8Tau3[i]



            
            h_ptAK8.Fill( AK8Pt[i] )
            h_etaAK8.Fill( AK8Eta[i] )
            h_yAK8.Fill( AK8Y[i] )
            h_mAK8.Fill( AK8Mass[i] )
            h_mprunedAK8.Fill( AK8PrunedM[i] )
            h_mfilteredAK8.Fill( AK8PrunedM[i] )
            h_mtrimmedAK8.Fill( AK8TrimmedM[i] )
            h_msoftdropAK8.Fill( AK8SoftDropM[i] )
            h_rhosoftDropAK8.Fill( jetrho, softdrop_p4.Perp() )
            if options.verbose :
                print '  << - good subjets, check m_softdrop = %6.2f' % (AK8SoftDropM[i])



        if len(goodJetIndices) < 2 :
            continue
        
        # Make plots of jets that pass the dijet selections
        jetP4_0 = ROOT.TLorentzVector()
        jetP4_0.SetPtEtaPhiM( AK8Pt[goodJetIndices[0]],
                              AK8Eta[goodJetIndices[0]],
                              AK8Phi[goodJetIndices[0]],
                              AK8Mass[goodJetIndices[0]] )

        jetP4_1 = ROOT.TLorentzVector()
        jetP4_1.SetPtEtaPhiM( AK8Pt[goodJetIndices[1]],
                              AK8Eta[goodJetIndices[1]],
                              AK8Phi[goodJetIndices[1]],
                              AK8Mass[goodJetIndices[1]] )

        tau1_0 = AK8Tau1[goodJetIndices[0]]
        tau2_0 = AK8Tau2[goodJetIndices[0]]
        tau3_0 = AK8Tau3[goodJetIndices[0]]

        tau1_1 = AK8Tau1[goodJetIndices[1]]
        tau2_1 = AK8Tau2[goodJetIndices[1]]
        tau3_1 = AK8Tau3[goodJetIndices[1]]

        rho_0 = AK8Rho[goodJetIndices[0]]
        rho_1 = AK8Rho[goodJetIndices[1]]

        tau21_0 = 0.0
        if tau1_0 > 0.01 :
            tau21_0 = tau2_0 / tau1_0
        tau21_1 = 0.0
        if tau1_1 > 0.01 :
            tau21_1 = tau2_1 / tau1_1        
            

        mjj = (jetP4_0 + jetP4_1).M()

        dijetsP4 = [ jetP4_0, jetP4_1 ]
        dijetsTau21 = [tau21_0, tau21_1 ]
        dijetsRho = [rho_0, rho_1]


        itag = random.randint(0,1)
        tagJetP4 = dijetsP4[itag]
        probeJetP4 = dijetsP4[1 - itag]
        tagJetTau21 = dijetsTau21[itag]
        probeJetTau21 = dijetsTau21[1 - itag]
        tagJetRho = dijetsRho[itag]
        probeJetRho = dijetsRho[1-itag]

        h_mjj0.Fill( mjj )
        h_rhoProbe0.Fill( probeJetRho, probeJetP4.Perp())

                
        if abs ( jetP4_0.Rapidity() - jetP4_1.Rapidity() ) > 1.0 :
            h_mjj1.Fill( mjj )
            h_rhoProbe1.Fill( probeJetRho, probeJetP4.Perp() )

            if tagJetP4.M() > 60. and tagJetTau21 < 0.6 :
                h_mjj2.Fill( mjj )
                h_rhoProbe2.Fill( probeJetRho, probeJetP4.Perp() )

                if probeJetP4.M() > 60. and probeJetTau21 < 0.6 :
                    h_mjj3.Fill( mjj )
                    h_rhoProbe3.Fill( probeJetRho, probeJetP4.Perp() )

                
        
        
#CLEANUP

f.cd()
f.Write()
f.Close()
