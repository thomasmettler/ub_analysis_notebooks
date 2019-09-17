import uproot
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
import globale

# active volume : 
# lower = [-1.55, -115.53, 0.1]
# upper = [254.8, 117.47, 1036.9]

def printGlobal():
    print globale.data
    print globale.ext
    print globale.dirt
    print globale.overlay
    print globale.scale
    print globale.tot_num_fidVol
    print globale.overlay_signals
    print globale.sample
    print globale.outputdir_png
    print globale.outputdir_root
    return


def loadGlobal(data_in,ext_in,dirt_in,overlay_in,scale_in,tot_num_fidVol_in,overlay_signals_in,sample_in,name_in,outputdir_png_in,outputdir_root_in):
    globale.data = data_in
    globale.ext = ext_in
    globale.dirt = dirt_in
    globale.overlay = overlay_in
    globale.scale = scale_in
    globale.tot_num_fidVol = tot_num_fidVol_in
    globale.overlay_signals = overlay_signals_in
    globale.sample = sample_in
    globale.name = name_in
    globale.outputdir_png = outputdir_png_in
    globale.outputdir_root = outputdir_root_in
    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    return

def prepareOutput(outputdir):
    outputdir_png = outputdir+'/png/'
    outputdir_root = outputdir+'/root/'
    try:
        os.stat(outputdir)
    except:
        os.mkdir(outputdir)
    try:
        os.stat(outputdir_png)
    except:
        os.mkdir(outputdir_png)
    try:
        os.stat(outputdir_root)
    except:
        os.mkdir(outputdir_root)
    return outputdir_png, outputdir_root

def openTrees(input_dir, file_data, file_ext, file_dirt, file_overlay, tree_name):
    data = ROOT.TChain(tree_name + "/event","event")
    data.Add( input_dir + file_data)
    ext = ROOT.TChain(tree_name + "/event","event")
    ext.Add( input_dir + file_ext)
    dirt = ROOT.TChain(tree_name + "/event","event")
    dirt.Add( input_dir + file_dirt)
    overlay = ROOT.TChain(tree_name + "/event","event")
    overlay.Add( input_dir + file_overlay)
    return data, ext, dirt, overlay

def printNumberOfEntries(data,ext,dirt,overlay):
    print "Overlay: Number of Entries:\t", overlay.GetEntries()
    print "Data: Number of Entries:\t", data.GetEntries()
    print "Ext: Number of Entries:\t\t", ext.GetEntries()
    print "Dirt: Number of Entries:\t", dirt.GetEntries()
    print ''
    return

def calculateScale(data_trigger, ext_trigger, pot_data, pot_dirt, pot_overlay):
    scale_data = 1.0
    scale_ext = data_trigger/ext_trigger
    scale_overlay = pot_data/pot_overlay
    scale_dirt = pot_data/pot_dirt
    return scale_data, scale_ext, scale_dirt, scale_overlay

def getPOT(input_dir, file_name, tree_name):
    tree_pot = ROOT.TChain(tree_name + "/pottree","pottree")
    tree_pot.Add( input_dir + file_name)
    h_x = ROOT.TH1F("h_x","Pot",1000,0,1e22)
    tree_pot.Draw('pot>>h_x','1')
    mean_pot = h_x.GetMean()
    entries_pot = h_x.GetEntries()
    pot = mean_pot*entries_pot
    #print file_name,' POT: ', pot #, entries_pot, mean_pot
    return pot

def printEff(cut):
    
    print 'Rejection on each sample:'
    for x in globale.sample:
        print 'Keep of %7s'%globale.name[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon')) + '\t reject:\t{0:0.2f}%'.format(100.0 - x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon'))
    print ''    
    print 'Efficiency:\t{0:0.2f}%'.format( globale.overlay.GetEntries(cut+' && numu_signal')*100.0/globale.overlay.GetEntries('numu_signal'))
    print 'Purity:\t\t{0:0.2f}%'.format(  globale.overlay.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/getTotNum(cut))
    print ''
    for x in globale.sample:
        print 'Amount of %7s'%globale.name[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*globale.scale[x]*100.0/getTotNum(cut))
    print ''
    for x in globale.overlay_signals:
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(globale.overlay.GetEntries(cut+' && '+x)*globale.scale[globale.overlay]*100.0/getTotNum(cut))+'\t({0:0.2f}%)'.format(globale.overlay.GetEntries(cut + ' && ' + x)*100.0/(globale.overlay.GetEntries(cut+' && fidVol')))
    return

def getTotNum(cut):
    num_fidVol = 0.0
    for x in globale.sample:
        if x != globale.data:
            num_fidVol = num_fidVol+x.GetEntries(cut + ' && fidVol')*globale.scale[x]
    return num_fidVol


def make_stacked_histo_plot(cut,variable,title,xstart,xend,xbins,file_name):
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {}
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    variable = variable
    globale.data.Draw(variable+'>>h_data_func',cut,'')
    globale.ext.Draw(variable+'>>h_ext_func',cut,'')
    globale.dirt.Draw(variable+'>>h_dirt_func',cut,'')
    #overlay.Draw('TimFla>>h_overlay',cut,'')
    #h_data_func.SetXTitle(title)
    h_data_func.SetYTitle("Entries per bin")
    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay.Draw(variable+'>>'+histo,cut+x,'')
    #legend = ROOT.TLegend(0.15,0.65,0.3,0.9) # LEGEND LEFT
    legend = ROOT.TLegend(0.75,0.65,0.9,0.9); #LEGEND RIGHT
    legend.AddEntry(h_data_func,"data","lep");
    legend.AddEntry(h_ext_func,"ext","f");
    legend.AddEntry(h_dirt_func,"dirt","f");
    for x in globale.overlay_signals:
        legend.AddEntry(h_overlay_func[x],x,"f");

    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(3)
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetFillColor((2*i+11))
        h_overlay_func[x].SetLineColor((2*i+11))
        h_overlay_func[x].SetLineColor(1)
        hs.Add(h_overlay_func[x])
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    h_data_func.SetMinimum(0)
    h_data_func.Draw('E')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E same')

    c1.cd()
    pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.33)
    pad2.SetGrid(1)
    pad2.Draw()
    pad2.cd()
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func)
    h_div_func.SetMinimum(0.5)
    h_div_func.SetMaximum(1.5)
    h_div_func.SetTitle('')
    h_div_func.SetYTitle('Data/(Ext+MC)')
    h_div_func.GetYaxis().SetTitleSize(0.08)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.06)
    h_div_func.SetXTitle(title)
    h_div_func.GetXaxis().SetTitleSize(0.1)
    h_div_func.GetXaxis().SetTitleOffset(0.4)
    h_div_func.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    #img = mpimg.imread(globale.outputdir_png+ file_name + ".png")
    #plt.rcParams['figure.figsize'] = [16, 12]
    #plt.imshow(img)
    #plt.show()
    return

def make_stacked_histo_plot_crtcorr(cut,variable,title,xstart,xend,xbins,file_name):
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {}
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    variable = variable
    globale.data.Draw(variable+'>>h_data_func',cut,'')
    globale.ext.Draw(variable+'-3.57+3.195>>h_ext_func',cut,'')
    globale.dirt.Draw(variable+'>>h_dirt_func',cut,'')
    #overlay.Draw('TimFla>>h_overlay',cut,'')
    #h_data_func.SetXTitle(title)
    h_data_func.SetYTitle("Entries per bin")
    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay.Draw(variable+'>>'+histo,cut+x,'')
    #legend = ROOT.TLegend(0.15,0.65,0.3,0.9) # LEGEND LEFT
    legend = ROOT.TLegend(0.75,0.65,0.9,0.9); #LEGEND RIGHT
    legend.AddEntry(h_data_func,"data","lep");
    legend.AddEntry(h_ext_func,"ext","f");
    legend.AddEntry(h_dirt_func,"dirt","f");
    for x in globale.overlay_signals:
        legend.AddEntry(h_overlay_func[x],x,"f");

    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(3)
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetFillColor((2*i+11))
        h_overlay_func[x].SetLineColor((2*i+11))
        h_overlay_func[x].SetLineColor(1)
        hs.Add(h_overlay_func[x])
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    h_data_func.Draw('E')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E same')

    c1.cd()
    pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.33)
    pad2.SetGrid(1)
    pad2.Draw()
    pad2.cd()
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func)
    h_div_func.SetMinimum(0.5)
    h_div_func.SetMaximum(1.5)
    h_div_func.SetTitle('')
    h_div_func.SetYTitle('Data/(Ext+MC)')
    h_div_func.GetYaxis().SetTitleSize(0.08)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.06)
    h_div_func.SetXTitle(title)
    h_div_func.GetXaxis().SetTitleSize(0.1)
    h_div_func.GetXaxis().SetTitleOffset(0.4)
    h_div_func.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    #img = mpimg.imread(globale.outputdir_png+ file_name + ".png")
    #plt.rcParams['figure.figsize'] = [16, 12]
    #plt.imshow(img)
    #plt.show()
    return

def make_stacked_histo_plot_flashcorr(cut,variable,title,xstart,xend,xbins,file_name):
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {}
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    variable = variable
    globale.data.Draw(variable+'>>h_data_func',cut,'')
    globale.ext.Draw(variable+'-3.57+3.195>>h_ext_func',cut,'')
    globale.dirt.Draw(variable+'-3.57+3.195>>h_dirt_func',cut,'')
    #overlay.Draw('TimFla>>h_overlay',cut,'')
    #h_data_func.SetXTitle(title)
    h_data_func.SetYTitle("Entries per bin")
    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay.Draw(variable+'-3.57+3.195>>'+histo,cut+x,'')
    #legend = ROOT.TLegend(0.15,0.65,0.3,0.9) # LEGEND LEFT
    legend = ROOT.TLegend(0.75,0.65,0.9,0.9); #LEGEND RIGHT
    legend.AddEntry(h_data_func,"data","lep");
    legend.AddEntry(h_ext_func,"ext","f");
    legend.AddEntry(h_dirt_func,"dirt","f");
    for x in globale.overlay_signals:
        legend.AddEntry(h_overlay_func[x],x,"f");

    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(3)
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetFillColor((2*i+11))
        h_overlay_func[x].SetLineColor((2*i+11))
        h_overlay_func[x].SetLineColor(1)
        hs.Add(h_overlay_func[x])
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    h_data_func.Draw('E')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E same')

    c1.cd()
    pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.33)
    pad2.SetGrid(1)
    pad2.Draw()
    pad2.cd()
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func)
    h_div_func.SetMinimum(0.5)
    h_div_func.SetMaximum(1.5)
    h_div_func.SetTitle('')
    h_div_func.SetYTitle('Data/(Ext+MC)')
    h_div_func.GetYaxis().SetTitleSize(0.08)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.06)
    h_div_func.SetXTitle(title)
    h_div_func.GetXaxis().SetTitleSize(0.1)
    h_div_func.GetXaxis().SetTitleOffset(0.4)
    h_div_func.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    #img = mpimg.imread(globale.outputdir_png+ file_name + ".png")
    #plt.rcParams['figure.figsize'] = [16, 12]
    #plt.imshow(img)
    #plt.show()
    return
def make_stacked_histo_plot_crtt0corr(cut,variable,title,xstart,xend,xbins,file_name):
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {}
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    variable = variable
    globale.data.Draw(variable+'-(69000-crt_trig_corr_med)/1000>>h_data_func',cut,'')
    globale.ext.Draw(variable+'-(69000-crt_trig_corr_med)/1000-3.57+3.195>>h_ext_func',cut,'')
    globale.dirt.Draw(variable+'>>h_dirt_func',cut,'')
    #overlay.Draw('TimFla>>h_overlay',cut,'')
    #h_data_func.SetXTitle(title)
    h_data_func.SetYTitle("Entries per bin")
    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay.Draw(variable+'>>'+histo,cut+x,'')
    #legend = ROOT.TLegend(0.15,0.65,0.3,0.9) # LEGEND LEFT
    legend = ROOT.TLegend(0.75,0.65,0.9,0.9); #LEGEND RIGHT
    legend.AddEntry(h_data_func,"data","lep");
    legend.AddEntry(h_ext_func,"ext","f");
    legend.AddEntry(h_dirt_func,"dirt","f");
    for x in globale.overlay_signals:
        legend.AddEntry(h_overlay_func[x],x,"f");

    hs = ROOT.THStack("hs","");
    h_ext_func.SetFillColor(2)
    h_ext_func.SetLineColor(1)
    h_dirt_func.SetFillColor(42)
    h_dirt_func.SetLineColor(1)
    h_data_func.SetLineWidth(3)
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetFillColor((2*i+11))
        h_overlay_func[x].SetLineColor((2*i+11))
        h_overlay_func[x].SetLineColor(1)
        hs.Add(h_overlay_func[x])
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    c1.cd()
    pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
    pad1.SetGrid(1)
    pad1.Draw()
    pad1.cd()
    h_data_func.Draw('E')
    legend.Draw();
    hs.Draw('same hist')
    h_data_func.Draw('E same')

    c1.cd()
    pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.33)
    pad2.SetGrid(1)
    pad2.Draw()
    pad2.cd()
    h_tot_func = h_ext_func.Clone()
    h_div_func = h_data_func.Clone()
    h_tot_func.Add(h_dirt_func)
    for i,x in enumerate(globale.overlay_signals):
        h_tot_func.Add(h_overlay_func[x])
    h_div_func.Divide(h_tot_func)
    h_div_func.SetMinimum(0.5)
    h_div_func.SetMaximum(1.5)
    h_div_func.SetTitle('')
    h_div_func.SetYTitle('Data/(Ext+MC)')
    h_div_func.GetYaxis().SetTitleSize(0.08)
    h_div_func.GetYaxis().SetTitleOffset(0.3)
    h_div_func.GetYaxis().SetLabelSize(0.06)
    h_div_func.SetXTitle(title)
    h_div_func.GetXaxis().SetTitleSize(0.1)
    h_div_func.GetXaxis().SetTitleOffset(0.4)
    h_div_func.Draw()
    c1.Draw()
    c1.SaveAs(globale.outputdir_png+ file_name + ".png")
    c1.SaveAs(globale.outputdir_root+ file_name + ".root")
    #img = mpimg.imread(globale.outputdir_png+ file_name + ".png")
    #plt.rcParams['figure.figsize'] = [16, 12]
    #plt.imshow(img)
    #plt.show()
    return


'''

xstart = 2
xend = 6
xbins = 50
title = 'Flash time [us]'

h_data = ROOT.TH1F("h_data",title,xbins,xstart,xend)
h_ext = ROOT.TH1F("h_ext",title,xbins,xstart,xend)
h_dirt = ROOT.TH1F("h_dirt",title,xbins,xstart,xend)
h_overlay = {}
for x in overlay_signals:
    h_overlay[x] = ROOT.TH1F(x,title,xbins,xstart,xend)

cut = 'fidVol'
data.Draw('TimFla>>h_data',cut)
ext.Draw('TimFla-3.57+3.195>>h_ext',cut,'')
dirt.Draw('TimFla-3.57+3.195>>h_dirt',cut,'')
#overlay.Draw('TimFla-3.57+3.195>>h_overlay',cut,'')
#h_data.SetXTitle("Flash time [us]")
h_data.SetYTitle("Entries per bin")
cut = cut +' && '
for x in overlay_signals:
    histo = x
    overlay.Draw('TimFla-3.57+3.195>>'+histo,cut+x,'')
    
legend = ROOT.TLegend(0.15,0.7,0.48,0.9);
legend.AddEntry(h_data,"data","lep");
legend.AddEntry(h_ext,"ext","f");
legend.AddEntry(h_dirt,"dirt","f");
for x in overlay_signals:
    legend.AddEntry(h_overlay[x],x,"f");

hs = ROOT.THStack("hs","");
h_ext.SetFillColor(2)
h_ext.SetLineColor(1)
h_dirt.SetFillColor(42)
h_dirt.SetLineColor(1)
h_data.SetLineWidth(3)
h_ext.Scale(scale[ext])
h_dirt.Scale(scale[dirt])
hs.Add(h_ext)
hs.Add(h_dirt)
for i,x in enumerate(overlay_signals):
    h_overlay[x].Scale(scale[overlay])
    h_overlay[x].SetFillColor((2*i+11))
    h_overlay[x].SetLineColor((2*i+11))
    h_overlay[x].SetLineColor(1)
    hs.Add(h_overlay[x])
c1.cd()
pad1 = ROOT.TPad('pad1','pad1',0,0.3,1,1)
pad1.SetGrid(1)
pad1.Draw()
pad1.cd()
#c_all.cd(1)
h_data.Draw('E')
legend.Draw();
hs.Draw('same hist')
h_data.Draw('E same')
c1.cd()
pad2 = ROOT.TPad('pad2','pad2',0,0.05,1,0.33)
pad2.SetGrid(1)
pad2.Draw()
pad2.cd()
h_tot = h_ext.Clone()
h_div = h_data.Clone()
h_tot.Add(h_dirt)
for i,x in enumerate(overlay_signals):
    h_tot.Add(h_overlay[x])
h_div.Divide(h_tot)
h_div.SetMinimum(0.5)
h_div.SetMaximum(1.5)
h_div.SetTitle('')
h_div.SetYTitle('Data/(Ext+MC)')
h_div.GetYaxis().SetTitleSize(0.08)
h_div.GetYaxis().SetTitleOffset(0.3)
h_div.GetYaxis().SetLabelSize(0.06)
h_div.SetXTitle("Flash time [us]")
h_div.GetXaxis().SetTitleSize(0.1)
h_div.GetXaxis().SetTitleOffset(0.45)
h_div.Draw()
c1.Draw()
#c_all.Draw()
c1.SaveAs(outputdir_png+ "Flash.png")
c1.SaveAs(outputdir_root+ "Flash.root")

'''