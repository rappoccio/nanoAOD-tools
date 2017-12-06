"""@TTbarResAnaHadronic Package to perform the data-driven mistag-rate-based ttbar hadronic analysis. 
"""

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSysColl import JetSysColl, JetSysObj

import random
import itertools
ROOT.gSystem.Load('libAnalysisPredictedDistribution')

"""@TTbarResAnaHadronic Package to perform the data-driven mistag-rate-based ttbar hadronic analysis. 

This module must be run twice: first to make the mistag rate in the "anti-tag and probe" selection,
and then applies that mistag rate to the single-tag selection. These are all done in bins of
b-tag categories (0, 1, >=2) and rapidity (|y| <= 1.0, |y| > 1.0).
The signal region is two top-tagged jets. 
The background estimate is the single-tag selection weighted by the mistag rate from the
"anti-tag and probe" region. 

The preselection is:
 - AK4-based HT > 1100 GeV (to be on the trigger plateau). 
 - >= 2 AK8 jets with pt > 400 and |eta| < 2.5, loose jet ID applied from matched AK4 jets

The 1-tag selection adds:
 - >=1 AK8 jet with top tagging applied to randomly-assigned tag jet. 

The anti-tag selection is disjoint from the 1-tag selection:
 - >=1 AK8 jet with top tagging VETO applied to randomly-assigned tag jet. 

The 2-tag selection is:
 - >=2 AK8 jets with top tagging applied to both leading jets. 

The ttbar candidate mass assumes the two leading top-tagged jets are the top quarks. 
"""
class TTbarResAnaHadronic(Module):
    def __init__(self, htCut=1100., minMSD=110., maxMSD=240., tau32Cut=0.6, ak8PtMin=400., bdisc=0.7, writePredDist=False ):
        """ Initialization for the module 
        """
        self.htCut = htCut
        self.minMSD = minMSD
        self.maxMSD = maxMSD
        self.tau32Cut = tau32Cut
        self.ak8PtMin = ak8PtMin
        self.bdisc = bdisc
        self.writePredDist = writePredDist
        self.writeHistFile = True

        self.systs = [
            'nom',
            'pu_up',  'pu_dn',
            'pdf_up', 'pdf_dn',
            'ps_up',  'ps_dn',
            'jec_up', 'jec_dn',
            'jer_up', 'jer_dn',
            'jms_up', 'jms_dn',
            'jmr_up', 'jmr_dn'
            ]
        self.btagcats = ["0b", "1b", "2b"]   # 0, 1, >=2 btags
        self.ycats = ['cen', 'fwd']          # Central and forward
        # Combine categories like "0bcen", "0bfwd", etc:
        self.anacats = [ b+y for b,y in itertools.product( self.btagcats, self.ycats) ]
        # Make string-based enumeration into aliases for faster processing speed
        self.systvals = []
        for isys,sys in enumerate(self.systs):
            setattr( self, sys, isys)
            self.systvals.append( getattr(self,sys) )
        self.anacatvals = []        
        for ianacat,anacat in enumerate(self.anacats):
            setattr( self, anacat, ianacat )
            self.anacatvals.append( getattr(self,anacat) )

    def formcat(self, nbtag, ycat ):
        """Form analysis category based on nbtag and ycat.

        """
        anacat = nbtag * len(self.ycats) + ycat
        return anacat
        
    def beginJob(self, histFile, histDirName):
        """Book control histograms and the predictions for the background.

        The background is data-driven and estimated by weighting the 1-tag region
        by the mistag rate to extrapolate to the 2-tag region. 
        """
        Module.beginJob(self, histFile, histDirName)
        self.addObjectList (self.systs, ROOT.TH1F('h_ak4ht',   'h_ak4ht',   25, 0, 2500) )
        self.addObjectList (self.systs, ROOT.TH1F('h_ak8pt',   'h_ak8pt',   25, 0, 2500) )
        self.addObjectList (self.systs, ROOT.TH1F('h_ak8msd',  'h_ak8msd',  25, 0, 500) )
        self.addObjectList (self.systs, ROOT.TH1F('h_ak8tau32','h_ak8tau32',25, 0, 1.0) )
        self.addObjectList (self.systs, ROOT.TH1F('h_ak8n3b1', 'h_ak8n3b1', 25, 0, 5.0) )
        self.addObjectList (self.systs, ROOT.TH1F('h_mttbar',  'h_mttbar',  25, 0, 5000) )

        
        if not self.writePredDist:
            self.predFile = ROOT.TFile( "ttbarreshad_predfile.root" )
            self.hpred = [ self.predFile.Get( "ttbarres/preddist" + str(ibtag) ) for ibtag in xrange(len(self.btagcatvals))]
            # PredictedDistribution needs to own this to ensure it doesn't go out of scope.
            for iana in xrange(len(self.anacatvals)) :
                ROOT.SetOwnership( self.hpred[iana], False )
                self.addObject( ROOT.PredictedDistribution(self.hpred[iana], "predJetP"+str(iana),        "Jet p_{T} (GeV), cat="+str(iana),   30, 0.0, 3000.) )
                self.addObject( ROOT.PredictedDistribution(self.hpred[iana], "predJetMTTBAR"+str(iana),   "M_{TTBAR} (GeV), cat="+str(iana),   50, 0.0, 5000.))
                self.addObject( ROOT.PredictedDistribution(self.hpred[iana], "predJetMTTBARMod"+str(iana),"M_{TTBAR} (GeV), cat="+str(iana),   50, 0.0, 5000.))
                self.addObject( ROOT.PredictedDistribution(self.hpred[iana], "predJetSDMass"+str(iana),    "Soft Drop Mass, cat="+str(iana),   50, 0.0, 250.))
                self.addObject( ROOT.PredictedDistribution(self.hpred[iana], "predJetMass"+str(iana),  "Ungroomed Jet Mass, cat="+str(iana),   50, 0.0, 250.))
                self.addObject( ROOT.PredictedDistribution(self.hpred[iana], "predJetMassMod"+str(iana), "Ungroomed Jet Mass, cat="+str(iana), 50, 0.0, 250.))
                self.addObject( ROOT.PredictedDistribution(self.hpred[iana], "predJetSDRho"+str(iana),        "Soft Drop Rho, cat="+str(iana), 50, 0.0, 1.))
            self.predJetP         = [ getattr( self, "predJetP"         + str (iana)) for iana in xrange(len(self.anacatvals))]
            self.predJetMTTBAR    = [ getattr( self, "predJetMTTBAR"    + str (iana)) for iana in xrange(len(self.anacatvals))]
            self.predJetMTTBARMod = [ getattr( self, "predJetMTTBARMod" + str (iana)) for iana in xrange(len(self.anacatvals))]
            self.predJetSDMass    = [ getattr( self, "predJetSDMass"    + str (iana)) for iana in xrange(len(self.anacatvals))]
            self.predJetMass      = [ getattr( self, "predJetMass"      + str (iana)) for iana in xrange(len(self.anacatvals))]
            self.predJetMassMod   = [ getattr( self, "predJetMassMod"   + str (iana)) for iana in xrange(len(self.anacatvals))]
            self.predJetSDRho     = [ getattr( self, "predJetSDRho"     + str (iana)) for iana in xrange(len(self.anacatvals))]

            # PredictedDistribution needs to own these also.
            for cat in self.anacatvals : 
                for hist in [
                        self.predJetP[cat], self.predJetMTTBAR[cat],
                        self.predJetMTTBARMod[cat],
                        self.predJetSDMass[cat], self.predJetMass[cat], self.predJetMassMod[cat],
                        self.predJetSDRho[cat] ] : 
                    ROOT.SetOwnership( hist, False )
        else:
            for iana in xrange(len(self.anacatvals)):
                self.addObject( ROOT.TH1D("preddist"+str(iana), "preddist"+str(iana), 25, 0, 2500) )
            self.preddist = [ getattr( self, "preddist"+str(iana)) for iana in xrange(len(self.anacatvals))]
            
    def endJob(self):
        """Calculate the correlated and uncorrelated errors.
        """
        if not self.writePredDist:
            for cat in self.anacatvals : 
                for hist in [
                        self.predJetP[cat], self.predJetMTTBAR[cat],
                        self.predJetMTTBARMod[cat],
                        self.predJetSDMass[cat], self.predJetMass[cat], self.predJetMassMod[cat],
                        self.predJetSDRho[cat] ] : 
                    hist.SetCalculatedErrors()
        Module.endJob(self)
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def passTopTag(self, jet):
        """Selects jets based on tau32 and soft drop mass. 
        """
        tau32 = jet.raw().tau3 / jet.raw().tau2 if jet.raw().tau2 > 0.0 else 0.
        passTau32 = tau32 < self.tau32Cut 
        passMSD = jet.msd() != None and self.minMSD < jet.msd() < self.maxMSD
        return passTau32 and passMSD


    
    def deriveWeights(self):
        """Derives all of the weights used that do not adjust the 4-vectors in the event.         
        """
        self.weightsDict = {}
        self.weightsDict[self.nom] = 1.0
        self.weightsDict[self.pu_up] = 1.0
        self.weightsDict[self.pu_dn] = 1.0
        self.weightsDict[self.pdf_up] = 1.0
        self.weightsDict[self.pdf_dn] = 1.0
        self.weightsDict[self.ps_up] = 1.0
        self.weightsDict[self.ps_dn] = 1.0
        return True
    
    def deriveJetSysts(self, jetSysCollAK4=None, jetSysCollAK8=None):
        """ Derive all of the jet-based kinematic systematic uncertainties. This includes JEC, JER, JMS, JMR.
        """
        if jetSysCollAK4 != None:
            for ijet, jet in enumerate(jetSysCollAK4.jets_raw()) :
                if ijet not in jetSysCollAK4[self.nom].keys() :
                    continue
                jetSysCollAK4[self.nom   ][ijet].p4().SetPtEtaPhiM( jet.pt_nom          , jet.eta, jet.phi, jet.mass_nom          )
                jetSysCollAK4[self.jer_up][ijet].p4().SetPtEtaPhiM( jet.pt_jerUp        , jet.eta, jet.phi, jet.mass_jerUp        )
                jetSysCollAK4[self.jer_dn][ijet].p4().SetPtEtaPhiM( jet.pt_jerDown      , jet.eta, jet.phi, jet.mass_jerDown      )
                jetSysCollAK4[self.jec_up][ijet].p4().SetPtEtaPhiM( jet.pt_jesTotalUp   , jet.eta, jet.phi, jet.mass_jesTotalUp   )
                jetSysCollAK4[self.jec_dn][ijet].p4().SetPtEtaPhiM( jet.pt_jesTotalDown , jet.eta, jet.phi, jet.mass_jesTotalDown )

        if jetSysCollAK8 != None:
            for ijet, jet in enumerate(jetSysCollAK8.jets_raw()) :
                if ijet not in jetSysCollAK8[self.nom].keys() :
                    continue
                jetSysCollAK8[self.nom   ][ijet].p4().SetPtEtaPhiM( jet.pt_nom          , jet.eta, jet.phi, jet.mass_nom          )
                jetSysCollAK8[self.jer_up][ijet].p4().SetPtEtaPhiM( jet.pt_jerUp        , jet.eta, jet.phi, jet.mass_jerUp        )
                jetSysCollAK8[self.jer_dn][ijet].p4().SetPtEtaPhiM( jet.pt_jerDown      , jet.eta, jet.phi, jet.mass_jerDown      )
                jetSysCollAK8[self.jec_up][ijet].p4().SetPtEtaPhiM( jet.pt_jesTotalUp   , jet.eta, jet.phi, jet.mass_jesTotalUp   )
                jetSysCollAK8[self.jec_dn][ijet].p4().SetPtEtaPhiM( jet.pt_jesTotalDown , jet.eta, jet.phi, jet.mass_jesTotalDown )
                jetSysCollAK8[self.jmr_up][ijet].p4().SetPtEtaPhiM( jet.pt_nom          , jet.eta, jet.phi, jet.mass_jmrUp        )
                jetSysCollAK8[self.jmr_dn][ijet].p4().SetPtEtaPhiM( jet.pt_nom          , jet.eta, jet.phi, jet.mass_jmrDown      )
                jetSysCollAK8[self.jms_up][ijet].p4().SetPtEtaPhiM( jet.pt_nom          , jet.eta, jet.phi, jet.mass_jmsUp        )
                jetSysCollAK8[self.jms_dn][ijet].p4().SetPtEtaPhiM( jet.pt_nom          , jet.eta, jet.phi, jet.mass_jmsDown      )
                
                jetSysCollAK8[self.nom   ][ijet].msd_ = jet.msoftdrop_nom
                jetSysCollAK8[self.jer_up][ijet].msd_ = jet.msoftdrop_jerUp                        
                jetSysCollAK8[self.jer_dn][ijet].msd_ = jet.msoftdrop_jerDown      
                jetSysCollAK8[self.jec_up][ijet].msd_ = jet.msoftdrop_jesTotalUp   
                jetSysCollAK8[self.jec_dn][ijet].msd_ = jet.msoftdrop_jesTotalDown 
                jetSysCollAK8[self.jmr_up][ijet].msd_ = jet.msoftdrop_jmrUp                        
                jetSysCollAK8[self.jmr_dn][ijet].msd_ = jet.msoftdrop_jmrDown      
                jetSysCollAK8[self.jms_up][ijet].msd_ = jet.msoftdrop_jmsUp        
                jetSysCollAK8[self.jms_dn][ijet].msd_ = jet.msoftdrop_jmsDown      
                
    def applyWeight(self, isys):
        """ Apply the weights for systematic "isys". If not found, use nominal. 
        """
        if isys not in self.weightsDict.keys():            
            weight = self.weightsDict[self.nom]
        else :
            weight = self.weightsDict[isys]
        return weight
        
    def analyze(self, event):
        """Perform either the anti-tag and probe (mistag estimate) or double tag (signal region) selection.
        """

        # Get the collections of AK4 and AK8 jets
        self.ak4JetsColl = Collection(event, "Jet")
        self.ak8JetsColl = Collection(event, "FatJet")

        # Select the jets that satisfy jet ID. 
        ak4Jets = [ x for x in self.ak4JetsColl if x.jetId > 0 ]
        ak8Jets = [ x for x in self.ak8JetsColl if x.jetId > 0 ]
        if len(ak8Jets) < 2 :
            return False
        
        # Make the systematic variations.
        jetSysCollAK4 = JetSysColl(self.ak4JetsColl, self.systvals, sel = lambda x : x.jetId > 0)
        jetSysCollAK8 = JetSysColl(self.ak8JetsColl, self.systvals, sel = lambda x : x.jetId > 0)
                
        # Derive the kinematic systematic effects. In this case,
        # jet-based systematic 4-vectors (AK4: JEC+JER, AK8:JEC+JER+JMS+JMR)
        self.deriveJetSysts(jetSysCollAK4, jetSysCollAK8)

        # Derive the weights to be used. 
        self.deriveWeights()
        
        # Loop over systematic uncertainties. These may change the kinematics,
        # or the weights. Both need to be adjusted.
        for isys,sys in enumerate(self.systs) :
            # Apply kinematic selection for this systematic. If no change, just use nominal. 
            ak4JetsSys = jetSysCollAK4[isys]
            ak8JetsSys = jetSysCollAK8[isys]
            
            # Adjust weight. If this systematic has no weight, just use nominal. 
            weight = self.applyWeight(isys)

            # Now get the AK4 and AK8 jets that pass the selection, and HT.
            # Don't copy the jet (expensive), copy the index (cheap)
            ak4JetsNdx = [i for i,x in ak4JetsSys.iteritems() if x.p4().Perp() > 20 and abs(x.p4().Eta())<2.5]
            ak8JetsNdx = [i for i,x in ak8JetsSys.iteritems() if x.p4().Perp() > self.ak8PtMin and abs(x.p4().Eta()) < 2.5]
            
            # Must have two AK8 jets that pass jet ID and kinematic cuts. 
            if len(ak8JetsNdx) < 2 :
                return False

            # Apply HT cut to ensure we are on the trigger plateau
            ht = sum( [ ak4JetsSys[j].p4().Perp() for j in ak4JetsNdx ] )
            self.h_ak4ht[isys].Fill( ht, weight )
            if ht < self.htCut :
                return False

            # Get a list of the jets that are top-tagged (ttag)
            isTagged = [ self.passTopTag(ak8JetsSys[x]) for x in ak8JetsNdx ]
            isTaggedDict = dict( zip( ak8JetsNdx ,isTagged) )

            

            # Make control plots
            for iak8Jet in ak8JetsNdx:
                jet = ak8JetsSys[iak8Jet]
                raw = jetSysCollAK8.jets_raw()[iak8Jet]
                self.h_ak8pt[isys].Fill( jet.p4().Perp(), weight )                
                self.h_ak8msd[isys].Fill( jet.msd(), weight )
                self.h_ak8tau32[isys].Fill( raw.tau3 / raw.tau2 if raw.tau2 > 0.0 else 0.0, weight )
                self.h_ak8n3b1[isys].Fill( raw.n3b1, weight )

            # Get a randomly assigned tag and probe jet from the leading two jets
            random.shuffle( ak8JetsNdx )
            iprobejet, itagjet = ak8JetsNdx[0:2]
            ttbarP4 =  ak8JetsSys[iprobejet].p4() + ak8JetsSys[itagjet].p4()

            # Find the analysis category: (0b,1b,2b) x (y<1,y>1)
            nbtag = min(2, sum( [ak8JetsSys[x].raw().maxCSVV2 > self.bdisc for x in [iprobejet,itagjet] ] ))
            yreg = 1 if abs( ak8JetsSys[iprobejet].p4().Rapidity() ) < 1.0 else 0
            anacat = self.formcat( nbtag, yreg )

            # Check if we have >=1 ttag
            if not self.passTopTag(ak8JetsSys[itagjet]) :
                # If we are in the signal selection, require at least 1 ttag.
                # Otherwise, we are vetoing the signal region
                # to derive the mistag weight (anti-tag and probe).
                if not self.writePredDist:
                    return False
                else:
                    # Here is the anti-ttag region selection to derive the mistag rate. 
                    if isys == self.nom:
                        self.preddist[anacat].Fill( ak8JetsSys[iprobejet].p4().P(), weight )

            # Here we have the actual signal region: 
            if not self.writePredDist:
                # Get the predicted background estimate
                if isys == self.nom : 
                    self.predJetP[anacat].Accumulate( ak8JetsSys[iprobejet].p4().P(), ak8JetsSys[iprobejet].p4().P(), isTaggedDict[iprobejet], weight )
                    self.predJetMTTBAR[anacat].Accumulate( ttbarP4.M(), ak8JetsSys[iprobejet].p4().P(), isTaggedDict[iprobejet], weight )
            # Now fill the double tagged histograms. 
            if isTaggedDict[iprobejet] : 
                self.h_mttbar[isys].Fill( ttbarP4.M(), weight )

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

ttbarreshad = lambda : TTbarResAnaHadronic() 
ttbarreshad_preddistwriter = lambda : TTbarResAnaHadronic(writePredDist=True)
