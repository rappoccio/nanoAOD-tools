import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection,Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import *

import random
import array

class ZPlusJetsXS(Module):
    def __init__(self ):
        self.writeHistFile = True
        self.verbose = False
    def beginJob(self, histFile, histDirName):
        Module.beginJob(self, histFile, histDirName)
        self.binsGen = array.array('d', [0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100]) 
        self.nGen = len(self.binsGen) - 1
        self.binsDet = array.array('d', [0, 0.5, 1, 3, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50, 75., 100.])
        self.nDet = len(self.binsDet) - 1
        self.nDetSD = 18
        self.nGenSD = 9

        self.minDPhiZJet = 1.57
        self.minZpt = 120.
        self.minJetPt = 220.
        
        self.addObject( ROOT.TH2D('h_response',     'h_response',   self.nDetSD, self.binsDet, self.nGenSD, self.binsGen) )
        self.addObject( ROOT.TH1D('h_reco',         'h_reco',       self.nDetSD, self.binsDet) )
        self.addObject( ROOT.TH1D('h_gen',          'h_gen',        self.nGenSD, self.binsGen) )
        self.addObject( ROOT.TH1D('h_fake',         'h_fake',       self.nDetSD, self.binsDet) )
        self.addObject( ROOT.TH1D('h_miss',         'h_miss',       self.nGenSD, self.binsGen) )
        
        self.addObject( ROOT.TH2D('h_response_u',   'h_response_u', self.nDet, self.binsDet, self.nGen, self.binsGen) )
        self.addObject( ROOT.TH1D('h_reco_u',       'h_reco_u',     self.nDet, self.binsDet) )
        self.addObject( ROOT.TH1D('h_gen_u',        'h_gen_u',      self.nGen, self.binsGen) )
        self.addObject( ROOT.TH1D('h_fake_u',       'h_fake_u',     self.nDet, self.binsDet) )
        self.addObject( ROOT.TH1D('h_miss_u',       'h_miss_u',     self.nGen, self.binsGen) )
        self.addObject( ROOT.TH1D('h_zpt',          'h_zpt',        100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_zmass',        'h_zmass',      100, 50, 150 ) )
        self.addObject( ROOT.TH1D('h_genjetpt',     'h_genjetpt',   100, 0, 500 ) )
        self.addObject( ROOT.TH1D('h_recojetpt',    'h_recojetpt',  100, 0, 500 ) )

        self.addObject( ROOT.TH1D('h_drGenReco',    'h_drGenReco',   40, 0, 0.8) )
        self.addObject( ROOT.TH1D('h_drGenGroomed', 'h_drGenGroomed',40, 0, 0.8) )
                            
    def endJob(self):
        Module.endJob(self)
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    def getSubjets(self, p4, subjets, dRmax=0.8):
        ret = []
        for subjet in subjets :
            if p4.DeltaR(subjet.p4()) < dRmax and len(ret) < 2 :
                ret.append(subjet.p4())
        return ret

    def printP4( self, c ):
        if hasattr( c, "p4"):
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.p4().Perp(), c.p4().Eta(), c.p4().Phi(), c.p4().M() )
        else :
            s = ' %6.2f %5.2f %5.2f %6.2f ' % ( c.Perp(), c.Eta(), c.Phi(), c.M() )
        return s
    def printCollection(self,coll):
        for ic,c in enumerate(coll):
            s = self.printP4( c )
            print ' %3d : %s' % ( ic, s )
            
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        weight = 1.0

        
        isMC = event.run == 1
        if self.verbose:
            print '------------------------ ', event.event

        if isMC:
            ###### Get gen Z candidate #######
            genleptons = Collection(event, "GenDressedLepton")

            if len(genleptons) < 2 :
                return False
            if abs(genleptons[0].pdgId) != 13 :
                return False
            if self.verbose :
                print '----'
                print 'Gen leptons:'
                self.printCollection( genleptons )
            Zboson = genleptons[0].p4() + genleptons[1].p4()
            if Zboson.Perp() < self.minZpt * 0.9 :
                return False
            if self.verbose:
                print '-----'
                print 'Gen Z:'
                print self.printP4( Zboson )

            ###### Get list of gen jets #######
            # List of gen jets:
            allgenjets = list(Collection(event, "GenJetAK8"))
            if self.verbose:
                print '-----'
                print 'all genjets:'
                self.printCollection( allgenjets )
            genjets = [ x for x in allgenjets if x.p4().Perp() > self.minJetPt * 0.8 and x.p4().DeltaPhi( Zboson ) > self.minDPhiZJet ]
            # List of gen subjets (no direct link from Genjet):
            gensubjets = list(Collection(event, "SubGenJetAK8"))
            # Dictionary to hold ungroomed-->groomed for gen
            genjetsGroomed = {}
            # Get the groomed gen jets
            for igen,gen in enumerate(genjets):
                gensubjetsMatched = self.getSubjets( p4=gen.p4(),subjets=gensubjets, dRmax=0.8)
                for isub,sub in enumerate(gensubjetsMatched) : 
                    self.h_drGenGroomed.Fill( gen.p4().DeltaR( sub ) )
                genjetsGroomed[gen] = sum( gensubjetsMatched, ROOT.TLorentzVector() ) if len(gensubjetsMatched) > 0 else None
                
            if self.verbose:
                print '----'
                print 'opposite-Z genjets:'
                for genjet in genjets:
                    sdmassgen = genjetsGroomed[genjet].M() if genjet in genjetsGroomed else -1.0
                    print '         : %s %6.2f' % ( self.printP4(genjet), sdmassgen )            
            

            
        ###### Get reco Z candidate #######
        # List of reco muons
        allmuons = Collection(event, "Muon")
        # Select reco muons:
        muons = [ x for x in allmuons if x.tightId ]
        if len(muons) < 2 :
            return False
        Zcand = muons[0].p4() + muons[1].p4()
        if Zcand.Perp() < self.minZpt or Zcand.M() < 50. or Zcand.M() > 150. :
            return False
        self.h_zpt.Fill( Zcand.Perp() )
        self.h_zmass.Fill( Zcand.M() )
        if self.verbose:
            print '-----'
            print ' recoZ:', self.printP4( Zcand )
        
        ###### Get list of reco jets #######
        # List of reco jets:
        allrecojets = list(Collection(event, "FatJet"))
        if self.verbose:
            print '----'
            print 'all recojets:'
            self.printCollection( allrecojets )
        recojets = [ x for x in allrecojets if x.p4().Perp() > self.minJetPt and x.p4().DeltaPhi( Zcand ) > self.minDPhiZJet ]
        if isMC == False:
            genjets = [None] * len(recojets)
        # List of reco subjets:
        recosubjets = list(Collection(event,"SubJet"))
        # Dictionary to hold reco--> gen matching
        recoToGen = matchObjectCollection( recojets, genjets, dRmax=0.05 )
        # Dictionary to hold ungroomed-->groomed for reco
        recojetsGroomed = {}        
        # Get the groomed reco jets
        for ireco,reco in enumerate(recojets):
            if reco.subJetIdx1 >= 0 and reco.subJetIdx2 >= 0 :
                recojetsGroomed[reco] = recosubjets[reco.subJetIdx1].p4() + recosubjets[reco.subJetIdx2].p4()
            elif reco.subJetIdx1 >= 0 :
                recojetsGroomed[reco] = recosubjets[reco.subJetIdx1].p4()
            else :
                recojetsGroomed[reco] = None

        if self.verbose:
            print '----'
            print 'opposite-Z recojets:'
            for recojet in recojets:
                sdmassreco = recojetsGroomed[recojet].M() if recojet in recojetsGroomed and recojetsGroomed[recojet] != None else -1.0
                print '         : %s %6.2f' % ( self.printP4( recojet),  sdmassreco )            

                
        # Loop over the reco,gen pairs.
        # Check if there are reco and gen SD jets
        # If both reco+gen: "fill"
        # If only reco: "fake"
        # (See below for "misses")
        for reco,gen in recoToGen.iteritems():
            recoSD = recojetsGroomed[reco]
            if reco == None :
                continue
            # Always fill the ungroomed det
            self.h_reco_u.Fill( reco.p4().M() )
            if recoSD != None :
                # Fill the groomed det if available
                self.h_reco.Fill( recoSD.M() )

            # Now check ungroomed gen
            genSDVal = None
            if gen != None:

                # Ungroomed gen OK, fill ungroomed response and truth
                self.h_response_u.Fill( reco.p4().M(), gen.p4().M() )
                self.h_gen_u.Fill( gen.p4().M() )
                self.h_genjetpt.Fill( gen.p4().Perp() )
                self.h_recojetpt.Fill( reco.p4().Perp() )
                self.h_drGenReco.Fill( reco.p4().DeltaR(gen.p4()) )

                genSD = genjetsGroomed[gen]
                if recoSD != None and genSD != None:
                    # Groomed gen OK, fill groomed response and truth
                    if self.verbose : 
                        print ' reco: %s %8.4f, gen : %s %8.4f ' % (
                            self.printP4(reco), recoSD.M(), 
                            self.printP4(gen), genSD.M()
                            )
                    genSDVal = genSD.M()
                    self.h_response.Fill( recoSD.M(), genSD.M() )
                    self.h_gen.Fill( genSD.M() )
            else :
                # No ungroomed jet at all. Ungroomed fake.
                self.h_fake_u.Fill( reco.p4().M() )
                # Here we have a groomed det, but no groomed gen. Groomed fake. 
                if genSDVal == None and recoSD != None :
                    self.h_fake.Fill( recoSD.M() )
        # Now loop over gen jets. If not in reco-->gen list,
        # then we have a "miss"
        for igen,gen in enumerate(genjets):
            if gen != None and gen not in recoToGen.values() :
                #print 'Missed!'
                #print 'Zboson: ', self.printP4( Zboson )
                #print 'Zcand : ', self.printP4( Zcand )
                #print 'GenJet: ', self.printP4( gen )
                #print 'Reco Jets:',
                #self.printCollection( recojets )
                #print ''
                #for ireco,igen in recoToGen :
                #    if ireco != None : 
                #        print self.printP4( ireco )
                #    if igen != None:
                #        print self.printP4( igen )

                # Ungroomed miss: 
                self.h_response_u.Fill( -1.0, gen.p4().M() )
                self.h_gen_u.Fill(gen.p4().M())
                self.h_miss_u.Fill(gen.p4().M())
                
                genSD = genjetsGroomed[gen]
                # Groomed miss: check if there is a groomed gen.
                # If there isn't, it gets skipped. 
                if genSD == None :
                    continue
                self.h_response.Fill( -1.0, genSD.M() )
                self.h_gen.Fill( genSD.M() )
                self.h_miss.Fill( genSD.M() )

        return True
# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

zplusjetsxs = lambda : ZPlusJetsXS() 
