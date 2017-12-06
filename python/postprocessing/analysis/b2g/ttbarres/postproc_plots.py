#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from PhysicsTools.NanoAODTools.postprocessing.analysis.b2g.ttbarres.TTbarResAnaHadronic import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetSmearer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.ak8JetIDProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.ak8SubjetVariableProducer import *
from PhysicsTools.NanoAODTools.postprocessing.examples.puWeightProducer import *


files=[
    "test94X_NANO_14_Skim.root",
    ]

import random
random.seed(12345)

p1=PostProcessor(".",files,'','',[ttbarreshad_preddistwriter()],provenance=False, noOut=True, histFileName='ttbarreshad_predfile.root', histDirName='ttbarres', postfix='predwrite')




#p0=PostProcessor(".",files,'FatJet_pt > 400.',"keep_and_drop.txt",[],outputbranchsel="output_keep_and_drop.txt",histFileName='ttbarreshad_predfile.root', histDirName='ttbarres', postfix='predwrite')

#p0=PostProcessor(".",files,'FatJet_pt > 400.',"keep_and_drop.txt",[jetmetUncertaintiesAK4Puppi(), jetmetUncertaintiesAK8Puppi(), ak8JetID(),ak8SubjetVariables(),ttbarreshad_preddistwriter()],outputbranchsel="output_keep_and_drop.txt", histFileName='ttbarreshad_predfile.root', histDirName='ttbarres', postfix='predwrite')


#p2=PostProcessor(".",['test94X_NANO_addPU.root'],'',"keep_and_drop.txt",[TTbarResAnaHadronic()],provenance=False, noOut=True,histFileName='hists.root', histDirName='ttbarres', postfix='predread')

p1.run()
#p2.run()
