#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.analysis.smp.xs.ZPlusJetsXS import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties import *
#from PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *

if False :
    with open('zjets_files.txt') as f:
        files = f.readlines()
    files = [x.strip() for x in files] 


files = [ 'DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root']
    
import random
random.seed(12345)
p1=PostProcessor(".",files,'nMuon >= 2 && nFatJet + nGenJetAK8 >= 1',"keep_and_drop.txt",[ZPlusJetsXS()],provenance=False, noOut=True,
                     histFileName='zplusjetsxs_hists.root', histDirName='zjets', postfix='zjets')
p1.run()
