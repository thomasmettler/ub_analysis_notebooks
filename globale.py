
import ROOT

data = ROOT.TChain()
ext = ROOT.TChain()
dirt = ROOT.TChain()
overlay = ROOT.TChain()
detsys = ROOT.TChain()
data_out = ROOT.TChain()
ext_out = ROOT.TChain()
dirt_out = ROOT.TChain()
overlay_out = ROOT.TChain()
detsys_out = ROOT.TChain()
scale = {data:1.0,ext:1.0,overlay:1.0,dirt:1.0,detsys:0}
scale_out = {data_out:1.0,ext_out:1.0,overlay_out:1.0,dirt_out:1.0,detsys_out:1.0}
tot_num_fidVol = 0.0
overlay_signals = {}
sample = {}
sample_out = {}
name = {}
name_out = {}
c1 = ROOT.TCanvas("c1","c1",1600,1200)
outputdir_png = ''
outputdir_root = ''
outputdir_pdf = ''

