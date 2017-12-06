import ROOT

f = ROOT.TFile("zplusjetsxs_hists.root")

isGroomed = False
if not isGroomed :
    postfix = '_u'
else :
    postfix = ''

h_response = f.Get('zjets/h_response' + postfix)
h_reco = f.Get('zjets/h_reco' + postfix)
h_gen = f.Get('zjets/h_gen' + postfix)
h_fake = f.Get('zjets/h_fake' + postfix)
h_miss = f.Get('zjets/h_miss' + postfix)

h_reco.Add( h_fake, -1.0 )

tunfolder = ROOT.TUnfoldDensity(h_response,ROOT.TUnfold.kHistMapOutputVert,ROOT.TUnfold.kRegModeNone, ROOT.TUnfold.kEConstraintNone, ROOT.TUnfoldDensity.kDensityModeBinWidth) 
tunfolder.SetInput( h_reco )

#scanResult = ROOT.TSpline3()
#tunfolder.ScanTau(1000,0.00001,0.00005,scanResult,ROOT.TUnfoldDensity.kEScanTauRhoAvg)


tunfolder.DoUnfold(0.0)
h_unfolded = tunfolder.GetOutput("unfolded")


for ihist in [ h_gen, h_unfolded ] :
    for x in xrange(ihist.GetNbinsX()+1):
        ihist.SetBinContent(x, ihist.GetBinContent(x)/ihist.GetXaxis().GetBinWidth(x))


h_reco.SetLineColor(ROOT.kBlack)
h_gen.SetLineColor(ROOT.kRed)
h_unfolded.SetMarkerStyle(20)
h_unfolded.SetMarkerColor(ROOT.kRed)

c1 = ROOT.TCanvas("c1", "c1")
hs = ROOT.THStack("hs", "hs")
hs.Add( h_gen, "hist" )
hs.Add( h_reco, "hist" )
hs.Add( h_unfolded, "e" )
hs.Draw("nostack")

c2 = ROOT.TCanvas("c2", "c2")
h_response.Draw("colz")
