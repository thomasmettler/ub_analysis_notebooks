
import global

data = ROOT.TChain()
ext = ROOT.TChain()
dirt = ROOT.TChain()
overlay = ROOT.TChain()
scale = {data:1.0,ext:1.0,overlay:1.0,dirt:1.0}
overlay_signals = {}
c1 = ROOT.TCanvas("c1","c1",1600,1200)
outputdir_png = ''
outputdir_root = ''

