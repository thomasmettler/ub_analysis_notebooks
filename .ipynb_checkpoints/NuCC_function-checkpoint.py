import uproot
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
import os
import ROOT
import time
import math
import imp
globale = imp.load_source('globale','/home/tmettler/Desktop/uBoone/do_plots/globale.py')
#import globale

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


def loadGlobal(data_in,ext_in,dirt_in,overlay_in,data_in_out,ext_in_out,dirt_in_out,overlay_in_out,scale_in,scale_in_out,tot_num_fidVol_in,overlay_signals_in,sample_in,sample_in_out,name_in, name_in_out,outputdir_png_in,outputdir_root_in):
    globale.data = data_in
    globale.ext = ext_in
    globale.dirt = dirt_in
    globale.overlay = overlay_in
    globale.data_out = data_in_out
    globale.ext_out = ext_in_out
    globale.dirt_out = dirt_in_out
    globale.overlay_out = overlay_in_out
    globale.scale = scale_in
    globale.scale_out = scale_in_out
    globale.tot_num_fidVol = tot_num_fidVol_in
    globale.overlay_signals = overlay_signals_in
    globale.sample = sample_in
    globale.sample_out = sample_in_out
    globale.name = name_in
    globale.name_out = name_in_out
    globale.outputdir_png = outputdir_png_in
    globale.outputdir_root = outputdir_root_in
    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)
    return

def loadGlobal_detsys(data_in, ext_in,dirt_in, overlay_in, detsys_in, data_in_out, ext_in_out, dirt_in_out, overlay_in_out, detsys_in_out, scale_in,scale_in_out, tot_num_fidVol_in, overlay_signals_in, sample_in, sample_in_out,name_in, name_in_out, outputdir_png_in, outputdir_root_in):
    globale.data = data_in 
    globale.ext = ext_in
    globale.dirt = dirt_in
    globale.overlay = overlay_in
    globale.detsys = detsys_in
    globale.data_out = data_in_out
    globale.ext_out = ext_in_out
    globale.dirt_out = dirt_in_out
    globale.overlay_out = overlay_in_out
    globale.detsys_out = detsys_in_out
    globale.scale = scale_in
    globale.scale_out = scale_in_out
    globale.tot_num_fidVol = tot_num_fidVol_in
    globale.overlay_signals = overlay_signals_in
    globale.sample = sample_in
    globale.sample_out = sample_in_out
    globale.name = name_in
    globale.name_out = name_in_out
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

def openTreesOut(input_dir, file_data, file_ext, file_dirt, file_overlay, tree_name):
    data = ROOT.TChain(tree_name,"t_out")
    data.Add( input_dir + file_data)
    ext = ROOT.TChain(tree_name,"t_out")
    ext.Add( input_dir + file_ext)
    dirt = ROOT.TChain(tree_name,"t_out")
    dirt.Add( input_dir + file_dirt)
    overlay = ROOT.TChain(tree_name,"t_out")
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
    h_x = ROOT.TH1F("h_x","Pot",1000,0,1e17)
    tree_pot.Draw('pot>>h_x','1')
    mean_pot = h_x.GetMean()
    entries_pot = h_x.GetEntries()
    pot = mean_pot*entries_pot
    #pot = h_x.Integral('width')
    #print file_name,' POT: ', pot #, entries_pot, mean_pot
    return pot

def printEff(cut):
    print "Efficiency for cut: ", cut
    print 'Rejection on each sample:'
    for x in globale.sample:
        if x != globale.data:
            print 'Keep of %7s'%globale.name[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon')) + '\t reject:\t{0:0.2f}%'.format(100.0 - x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon'))
    print ''    
    efficiency_func = globale.overlay.GetEntries(cut+' && numu_signal')*100.0/globale.overlay.GetEntries('numu_true')
    print 'Efficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    purity_func = globale.overlay.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/getTotNum(cut)
    print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    print 'Purity*Efficiency\t{0:0.2f}'.format(purity_func*efficiency_func/100)
    print ''
    for x in globale.sample:
        if x != globale.data:
            print 'Amount of %7s'%globale.name[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*globale.scale[x]*100.0/getTotNum(cut))
    print ''
    for x in globale.overlay_signals:
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(globale.overlay.GetEntries(cut+' && '+x)*globale.scale[globale.overlay]*100.0/getTotNum(cut))+'\t({0:0.2f}%)'.format(globale.overlay.GetEntries(cut + ' && ' + x)*100.0/(globale.overlay.GetEntries(cut+' && fidVol')))
    return

def printEff_out(cut):
    print "Efficiency for cut: ", cut
    print 'Rejection on each sample:'
    for x in globale.sample_out:
        if x != globale.data_out:
            print 'Keep of %7s'%globale.name_out[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon')) + '\t reject:\t{0:0.2f}%'.format(100.0 - x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon'))
    print ''    
    efficiency_func = globale.overlay_out.GetEntries(cut+' && numu_signal')*100.0/globale.overlay_out.GetEntries('numu_true')
    print 'Efficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    purity_func = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/getTotNum_out(cut)
    print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    print 'Purity*Efficiency\t{0:0.2f}'.format(purity_func*efficiency_func/100)
    print ''
    for x in globale.sample_out:
        if x != globale.data_out:
            print 'Amount of %7s'%globale.name_out[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*globale.scale_out[x]*100.0/getTotNum_out(cut))
    print ''
    for x in globale.overlay_signals:
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(globale.overlay_out.GetEntries(cut+' && '+x)*globale.scale[globale.overlay]*100.0/getTotNum_out(cut))+'\t({0:0.2f}%)'.format(globale.overlay_out.GetEntries(cut + ' && ' + x)*100.0/(globale.overlay_out.GetEntries(cut+' && fidVol')))
    return

def printEffonly_out(cut, data_sample):
    if data_sample == "overlay_out":
        this_sample =  globale.overlay_out
    if data_sample == 'ext_out':
        this_sample =  globale.ext_out
    if data_sample == 'dirt_out':
        this_sample =  globale.dirt_out
    if data_sample == 'data_out':
        this_sample =  globale.data_out
    print 'Rejection on each sample:'
    #for x in globale.sample_out:
    #    if x != globale.data_out:
    #        print 'Keep of %7s'%globale.name_out[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon')) + '\t reject:\t{0:0.2f}%'.format(100.0 - x.GetEntries(cut+' && fidVol')*100.0/x.GetEntries('fidVol && muon'))
    print 'Keep of %7s'%data_sample+':\t{0:0.2f}%'.format( this_sample.GetEntries(cut+' && fidVol')*100.0/this_sample.GetEntries('fidVol && muon')) + '\t reject:\t{0:0.2f}%'.format(100.0 - this_sample.GetEntries(cut+' && fidVol')*100.0/this_sample.GetEntries('fidVol && muon'))
    #print ''    
    #efficiency_func = this_sample.GetEntries(cut+' && numu_signal')*100.0/this_sample.GetEntries('numu_true')
    #print 'Efficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    #purity_func = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/getTotNum_out(cut)
    #print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    #print 'Purity*Efficiency\t{0:0.2f}'.format(purity_func*efficiency_func/100)
    #print ''
    #for x in globale.sample_out:
    #    if x != globale.data_out:
    #        print 'Amount of %7s'%globale.name_out[x]+':\t{0:0.2f}%'.format( x.GetEntries(cut+' && fidVol')*globale.scale_out[x]*100.0/getTotNum_out(cut))
    #print ''
    if data_sample == 'overlay_out':
        efficiency_func = this_sample.GetEntries(cut+' && numu_signal')*100.0/this_sample.GetEntries('numu_true')
        print 'Efficiency:\t\t{0:0.2f}%'.format(efficiency_func)
        for x in globale.overlay_signals:
            print 'Signal definition= %12s'%x+'\t({0:0.2f}%)'.format(globale.overlay_out.GetEntries(cut + ' && ' + x)*100.0/(globale.overlay_out.GetEntries(cut+' && fidVol')))
    return


def printEff_sample(cut, sample_cut):
    print "Efficiency for cut: ", cut
    print 'Rejection on each sample:'
    for x in globale.sample:
        if x != globale.data:
            print 'Keep of %7s'%globale.name[x]+':\t{0:0.2f}%'.format( x.GetEntries(sample_cut+cut+' && fidVol')*100.0/x.GetEntries(sample_cut+'fidVol && muon')) + '\t reject:\t{0:0.2f}%'.format(100.0 - x.GetEntries(sample_cut+cut+' && fidVol')*100.0/x.GetEntries(sample_cut+'fidVol && muon'))
    print ''    
    efficiency_func = globale.overlay.GetEntries(sample_cut+cut+' && numu_signal')*100.0/globale.overlay.GetEntries(sample_cut+'numu_true')
    print 'Efficiency:\t\t{0:0.2f}%'.format(efficiency_func)
    purity_func = globale.overlay.GetEntries(sample_cut+cut+' && numu_signal')*globale.scale[globale.overlay]*100/getTotNum(sample_cut+cut)
    print 'Purity:\t\t\t{0:0.2f}%'.format(purity_func)
    print 'Purity*Efficiency\t{0:0.2f}'.format(purity_func*efficiency_func/100)
    print ''
    for x in globale.sample:
        if x != globale.data:
            print 'Amount of %7s'%globale.name[x]+':\t{0:0.2f}%'.format( x.GetEntries(sample_cut+cut+' && fidVol')*globale.scale[x]*100.0/getTotNum(sample_cut+cut))
    print ''
    for x in globale.overlay_signals:
        print 'Signal definition= %12s'%x+': {0:0.2f}%'.format(globale.overlay.GetEntries(sample_cut+cut+' && '+x)*globale.scale[globale.overlay]*100.0/getTotNum(sample_cut+cut))+'\t({0:0.2f}%)'.format(globale.overlay.GetEntries(sample_cut+cut + ' && ' + x)*100.0/(globale.overlay.GetEntries(sample_cut+cut+' && fidVol')))
    return

def printNumber(cut):
    print "Efficiency for cut: ", cut
    print 'Rejection on each sample:'  
    n_signal = globale.overlay.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]
    n_signal_tot =getTotNum(cut)
    n_signal_data = globale.data.GetEntries(cut)*globale.scale[globale.data]
    n_background = globale.overlay.GetEntries(cut+' && !numu_signal')*globale.scale[globale.overlay]+\
    globale.ext.GetEntries(cut)*globale.scale[globale.ext]+\
    globale.dirt.GetEntries(cut)*globale.scale[globale.dirt]
    
    print 'Number of true signal:\t\t{0:0.2f}'.format(n_signal)
    print 'Number of total events:\t\t{0:0.2f}'.format(n_signal_tot)
    print 'Number of data events:\t\t{0:0.2f}'.format(n_signal_data)
    print 'Number of bk events:\t\t{0:0.2f}'.format(n_background)
    
    return

def getTotNum(cut):
    num_fidVol = 0.0
    for x in globale.sample:
        if x != globale.data:
            num_fidVol = num_fidVol+x.GetEntries(cut + ' && fidVol')*globale.scale[x]
    return num_fidVol

def getTotNum_out(cut):
    num_fidVol = 0.0
    for x in globale.sample_out:
        if x != globale.data_out:
            num_fidVol = num_fidVol+x.GetEntries(cut + ' && fidVol')*globale.scale_out[x]
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
    h_data_func.SetLineWidth(1)
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
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

def make_stacked_histo_plot_out(cut,variable,title,xstart,xend,xbins,file_name):
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {}
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    variable = variable
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'>>h_ext_func',cut,'')
    globale.dirt_out.Draw(variable+'>>h_dirt_func',cut,'')
    #overlay.Draw('TimFla>>h_overlay',cut,'')
    #h_data_func.SetXTitle(title)
    h_data_func.SetYTitle("Entries per bin")
    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,cut+x,'')
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
    h_data_func.SetLineWidth(1)
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
    h_ext_func.Sumw2()	
    h_ext_func.Scale(globale.scale[globale.ext])
    h_dirt_func.Sumw2()	
    h_dirt_func.Scale(globale.scale[globale.dirt])
    hs.Add(h_ext_func)
    hs.Add(h_dirt_func)
    mc_events = h_dirt_func.GetEntries()*globale.scale[globale.dirt]
    for i,x in enumerate(globale.overlay_signals):
        h_overlay_func[x].Sumw2()
        h_overlay_func[x].Scale(globale.scale[globale.overlay])
        h_overlay_func[x].SetFillColor((2*i+11))
        h_overlay_func[x].SetLineColor((2*i+11))
        h_overlay_func[x].SetLineColor(1)
        print 'MC entries, scaled: ', h_overlay_func[x].GetEntries()*globale.scale[globale.overlay]
        mc_events = mc_events+h_overlay_func[x].GetEntries()*globale.scale[globale.overlay]
        hs.Add(h_overlay_func[x])
    data_events = h_data_func.GetEntries()*globale.scale[globale.data]
    ext_events = h_ext_func.GetEntries()*globale.scale[globale.ext]
    normalization = (data_events- ext_events)/(mc_events+0.0000001)
    print 'Normalization (data-ext)/mc = ', normalization
    print 'Normalization (data)/(mc +ext) = ', (data_events)/(mc_events+ext_events+0.0000001)
    
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
    h_data_func.SetLineWidth(1)
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
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

def make_stacked_histo_plot_crtcorr_out(cut,variable,title,xstart,xend,xbins,file_name):
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {}
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    variable = variable
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'-3.57+3.195>>h_ext_func',cut,'')
    globale.dirt_out.Draw(variable+'>>h_dirt_func',cut,'')
    #overlay.Draw('TimFla>>h_overlay',cut,'')
    #h_data_func.SetXTitle(title)
    h_data_func.SetYTitle("Entries per bin")
    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay_out.Draw(variable+'>>'+histo,cut+x,'')
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
    h_data_func.SetLineWidth(1)
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
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
    h_data_func.SetLineWidth(1)
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
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

def make_stacked_histo_plot_flashcorr_out(cut,variable,title,xstart,xend,xbins,file_name):
    h_data_func = ROOT.TH1F("h_data_func",title,xbins,xstart,xend)
    h_ext_func = ROOT.TH1F("h_ext_func",title,xbins,xstart,xend)
    h_dirt_func = ROOT.TH1F("h_dirt_func",title,xbins,xstart,xend)
    h_overlay_func = {}
    for x in globale.overlay_signals:
        h_overlay_func[x] = ROOT.TH1F(x,title,xbins,xstart,xend)
    variable = variable
    globale.data_out.Draw(variable+'>>h_data_func',cut,'')
    globale.ext_out.Draw(variable+'-3.57+3.195>>h_ext_func',cut,'')
    globale.dirt_out.Draw(variable+'-3.57+3.195>>h_dirt_func',cut,'')
    #overlay.Draw('TimFla>>h_overlay',cut,'')
    #h_data_func.SetXTitle(title)
    h_data_func.SetYTitle("Entries per bin")
    cut = cut +' && '
    for x in globale.overlay_signals:
        histo = x
        globale.overlay_out.Draw(variable+'-3.57+3.195>>'+histo,cut+x,'')
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
    h_data_func.SetLineWidth(1)
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
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
    h_data_func.SetLineWidth(1)
    h_data_func.Sumw2()	
    h_data_func.Scale(globale.scale[globale.data])
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

def plot_eff_out(nenner_cut, zahler_cut, cut, name, title):
    efficiency = globale.overlay_out.GetEntries(zahler_cut)*100.0/(globale.overlay_out.GetEntries(nenner_cut)+0.001)
    purity = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/(getTotNum_out(cut)+0.001)

    print 'Efficiency:\t{0:0.2f}%'.format(efficiency)
    print 'Purity:\t\t{0:0.2f}%'.format(  purity)

    title = title+' Eff/Pur={0:0.2f}%'.format( efficiency)+'/{0:0.2f}%'.format(purity)

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)

    xstart = 0
    xend = 2
    xbin = 100
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Energy>>h_init_eff_energy',zahler_cut)
    globale.overlay_out.Draw('MCNu_Energy>>h_init_eff_energy_1',nenner_cut)
    h_init_eff_energy_1.Sumw2()
    h_init_eff_energy.Divide(h_init_eff_energy_1)
    h_init_eff_energy.Scale(100)
    h_init_eff_energy.SetMaximum(100)
    h_init_eff_energy.SetMinimum(0)
    h_init_eff_energy.SetXTitle("Truth neutrino energy [GeV]")
    h_init_eff_energy.SetYTitle("Signal efficiency [%]")
    h_init_eff_energy.Draw("E")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_energy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_energy"+name+".root")


    xstart = -1
    xend = 1
    xbin = 100
    h_init_eff_theta = ROOT.TH1F("h_init_eff_theta",title,xbin,xstart,xend)
    h_init_eff_theta_1 = ROOT.TH1F("h_init_eff_theta_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta',zahler_cut)
    globale.overlay_out.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta_1',nenner_cut)
    h_init_eff_theta_1.Sumw2()
    h_init_eff_theta.Divide(h_init_eff_theta_1)
    h_init_eff_theta.Scale(100)
    h_init_eff_theta.SetMaximum(100)
    h_init_eff_theta.SetMinimum(0)
    h_init_eff_theta.SetXTitle("Truth lepton cos(theta)")
    h_init_eff_theta.SetYTitle("Signal efficiency [%]")
    h_init_eff_theta.Draw("E")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_theta"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_theta"+name+".root")

    xstart = -3.14159
    xend = 3.14159
    xbin = 100
    h_init_eff_phi = ROOT.TH1F("h_init_eff_phi",title,xbin,xstart,xend)
    h_init_eff_phi_1 = ROOT.TH1F("h_init_eff_phi_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('TrackPhi>>h_init_eff_phi',zahler_cut+' && muon')
    globale.overlay_out.Draw('TrackPhi>>h_init_eff_phi_1',nenner_cut+'&& muon')
    h_init_eff_phi_1.Sumw2()
    h_init_eff_phi.Divide(h_init_eff_phi_1)
    h_init_eff_phi.Scale(100)
    h_init_eff_phi.SetMaximum(100)
    h_init_eff_phi.SetMinimum(0)
    h_init_eff_phi.SetXTitle("Reco lepton phi [pi]")
    h_init_eff_phi.SetYTitle("Signal efficiency [%]")
    h_init_eff_phi.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_phi"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_phi"+name+".root")

    xstart = -50
    xend = 300
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Vx>>h_init_eff_V',zahler_cut)
    globale.overlay_out.Draw('MCNu_Vx>>h_init_eff_V_1',nenner_cut)
    h_init_eff_V_1.Sumw2()
    h_init_eff_V.Divide(h_init_eff_V_1)
    h_init_eff_V.Scale(100)
    h_init_eff_V.SetMaximum(100)
    h_init_eff_V.SetMinimum(0)
    h_init_eff_V.SetXTitle("Truth nu vertex X [cm]")
    h_init_eff_V.SetYTitle("Signal efficiency [%]")
    h_init_eff_V.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vx"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vx"+name+".root")

    xstart = -150
    xend = 150
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Vy>>h_init_eff_V',zahler_cut)
    globale.overlay_out.Draw('MCNu_Vy>>h_init_eff_V_1',nenner_cut)
    h_init_eff_V_1.Sumw2()
    h_init_eff_V.Divide(h_init_eff_V_1)
    h_init_eff_V.Scale(100)
    h_init_eff_V.SetMaximum(100)
    h_init_eff_V.SetMinimum(0)
    h_init_eff_V.SetXTitle("Truth nu vertex Y [cm]")
    h_init_eff_V.SetYTitle("Signal efficiency [%]")
    h_init_eff_V.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vy"+name+".root")

    xstart = -50
    xend = 1050
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Vz>>h_init_eff_V',zahler_cut)
    globale.overlay_out.Draw('MCNu_Vz>>h_init_eff_V_1',nenner_cut)
    h_init_eff_V_1.Sumw2()
    h_init_eff_V.Divide(h_init_eff_V_1)
    h_init_eff_V.Scale(100)
    h_init_eff_V.SetMaximum(100)
    h_init_eff_V.SetMinimum(0)
    h_init_eff_V.SetXTitle("Truth nu vertex Z [cm]")
    h_init_eff_V.SetYTitle("Signal efficiency [%]")
    h_init_eff_V.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vz"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vz"+name+".root")
    
def plot_eff_outV2(nenner_cut, zahler_cut, cut, name, title):
    #efficiency = globale.overlay_out.GetEntries(zahler_cut)*100.0/globale.overlay_out.GetEntries(nenner_cut)
    efficiency = globale.overlay_out.GetEntries(cut+' && numu_signal')*100.0/(globale.overlay_out.GetEntries('numu_true')+0.001)
    #purity = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/getTotNum_out(cut)
    purity = globale.overlay_out.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/(getTotNum_out(cut)+0.001)

    print 'Efficiency:\t{0:0.2f}%'.format(efficiency)
    print 'Purity:\t\t{0:0.2f}%'.format(  purity)

    title = title+' Eff/Pur={0:0.2f}%'.format( efficiency)+'/{0:0.2f}%'.format(purity)

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)

    xstart = 0
    xend = 2
    xbin = 100
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Energy>>h_init_eff_energy',zahler_cut)
    globale.overlay_out.Draw('MCNu_Energy>>h_init_eff_energy_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_energy,h_init_eff_energy_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")
    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth neutrino energy [GeV]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    #graph1 = graph.Clone()
    #graph.SetFillColor(2);
    #graph.SetFillStyle(1001);
    #graph.Draw("a4");
    #graph1.Draw("same")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_energy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_energy"+name+".root")


    xstart = -1
    xend = 1
    xbin = 100
    h_init_eff_theta = ROOT.TH1F("h_init_eff_theta",title,xbin,xstart,xend)
    h_init_eff_theta_1 = ROOT.TH1F("h_init_eff_theta_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta',zahler_cut)
    globale.overlay_out.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_theta,h_init_eff_theta_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth lepton cos(theta)")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_theta"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_theta"+name+".root")

    xstart = -3.14159
    xend = 3.14159
    xbin = 100
    h_init_eff_phi = ROOT.TH1F("h_init_eff_phi",title,xbin,xstart,xend)
    h_init_eff_phi_1 = ROOT.TH1F("h_init_eff_phi_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('TrackPhi>>h_init_eff_phi',zahler_cut+' && muon')
    globale.overlay_out.Draw('TrackPhi>>h_init_eff_phi_1',nenner_cut+' && muon')
    eff =  ROOT.TEfficiency(h_init_eff_phi,h_init_eff_phi_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Reco lepton phi [pi]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_phi"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_phi"+name+".root")

    xstart = -50
    xend = 300
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Vx>>h_init_eff_V',zahler_cut)
    globale.overlay_out.Draw('MCNu_Vx>>h_init_eff_V_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth neutrino vertex X [cm]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vx"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vx"+name+".root")

    xstart = -150
    xend = 150
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Vy>>h_init_eff_V',zahler_cut)
    globale.overlay_out.Draw('MCNu_Vy>>h_init_eff_V_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth neutrino vertex Y [cm]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vy"+name+".root")

    xstart = -50
    xend = 1050
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay_out.Draw('MCNu_Vz>>h_init_eff_V',zahler_cut)
    globale.overlay_out.Draw('MCNu_Vz>>h_init_eff_V_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle('Truth neutrino vertex Z [cm]')
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vz"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vz"+name+".root")

def plot_eff(nenner_cut, zahler_cut, cut, name, title):
    #efficiency = globale.overlay.GetEntries(zahler_cut)*100.0/globale.overlay.GetEntries(nenner_cut)
    efficiency = globale.overlay.GetEntries(cut+' && numu_signal')*100.0/(globale.overlay.GetEntries('numu_true')+0.001)
    purity = globale.overlay.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/(getTotNum(cut)+0.001)

    print 'Efficiency:\t{0:0.2f}%'.format(efficiency)
    print 'Purity:\t\t{0:0.2f}%'.format(  purity)

    title = title+' Eff/Pur={0:0.2f}%'.format( efficiency)+'/{0:0.2f}%'.format(purity)

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)

    xstart = 0
    xend = 2
    xbin = 100
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Energy>>h_init_eff_energy',zahler_cut)
    globale.overlay.Draw('MCNu_Energy>>h_init_eff_energy_1',nenner_cut)
    h_init_eff_energy_1.Sumw2()
    h_init_eff_energy.Divide(h_init_eff_energy_1)
    h_init_eff_energy.Scale(100)
    h_init_eff_energy.SetMaximum(100)
    h_init_eff_energy.SetMinimum(0)
    h_init_eff_energy.SetXTitle("Truth neutrino energy [GeV]")
    h_init_eff_energy.SetYTitle("Signal efficiency [%]")
    h_init_eff_energy.Draw("E")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_energy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_energy"+name+".root")


    xstart = -1
    xend = 1
    xbin = 100
    h_init_eff_theta = ROOT.TH1F("h_init_eff_theta",title,xbin,xstart,xend)
    h_init_eff_theta_1 = ROOT.TH1F("h_init_eff_theta_1",title,xbin,xstart,xend)
    globale.overlay.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta',zahler_cut)
    globale.overlay.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta_1',nenner_cut)
    h_init_eff_theta_1.Sumw2()
    h_init_eff_theta.Divide(h_init_eff_theta_1)
    h_init_eff_theta.Scale(100)
    h_init_eff_theta.SetMaximum(100)
    h_init_eff_theta.SetMinimum(0)
    h_init_eff_theta.SetXTitle("Truth lepton cos(theta)")
    h_init_eff_theta.SetYTitle("Signal efficiency [%]")
    h_init_eff_theta.Draw("E")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_theta"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_theta"+name+".root")

    xstart = -3.14159
    xend = 3.14159
    xbin = 100
    h_init_eff_phi = ROOT.TH1F("h_init_eff_phi",title,xbin,xstart,xend)
    h_init_eff_phi_1 = ROOT.TH1F("h_init_eff_phi_1",title,xbin,xstart,xend)
    globale.overlay.Draw('TrackPhi>>h_init_eff_phi',zahler_cut)
    globale.overlay.Draw('TrackPhi>>h_init_eff_phi_1',nenner_cut)
    h_init_eff_phi_1.Sumw2()
    h_init_eff_phi.Divide(h_init_eff_phi_1)
    h_init_eff_phi.Scale(100)
    h_init_eff_phi.SetMaximum(100)
    h_init_eff_phi.SetMinimum(0)
    h_init_eff_phi.SetXTitle("Reco lepton phi [pi]")
    h_init_eff_phi.SetYTitle("Signal efficiency [%]")
    h_init_eff_phi.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_phi"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_phi"+name+".root")

    xstart = -50
    xend = 300
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Vx>>h_init_eff_V',zahler_cut)
    globale.overlay.Draw('MCNu_Vx>>h_init_eff_V_1',nenner_cut)
    h_init_eff_V_1.Sumw2()
    h_init_eff_V.Divide(h_init_eff_V_1)
    h_init_eff_V.Scale(100)
    h_init_eff_V.SetMaximum(100)
    h_init_eff_V.SetMinimum(0)
    h_init_eff_V.SetXTitle("Truth nu vertex X [cm]")
    h_init_eff_V.SetYTitle("Signal efficiency [%]")
    h_init_eff_V.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vx"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vx"+name+".root")

    xstart = -150
    xend = 150
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Vy>>h_init_eff_V',zahler_cut)
    globale.overlay.Draw('MCNu_Vy>>h_init_eff_V_1',nenner_cut)
    h_init_eff_V_1.Sumw2()
    h_init_eff_V.Divide(h_init_eff_V_1)
    h_init_eff_V.Scale(100)
    h_init_eff_V.SetMaximum(100)
    h_init_eff_V.SetMinimum(0)
    h_init_eff_V.SetXTitle("Truth nu vertex Y [cm]")
    h_init_eff_V.SetYTitle("Signal efficiency [%]")
    h_init_eff_V.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vy"+name+".root")

    xstart = -50
    xend = 1050
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Vz>>h_init_eff_V',zahler_cut)
    globale.overlay.Draw('MCNu_Vz>>h_init_eff_V_1',nenner_cut)
    h_init_eff_V_1.Sumw2()
    h_init_eff_V.Divide(h_init_eff_V_1)
    h_init_eff_V.Scale(100)
    h_init_eff_V.SetMaximum(100)
    h_init_eff_V.SetMinimum(0)
    h_init_eff_V.SetXTitle("Truth nu vertex Z [cm]")
    h_init_eff_V.SetYTitle("Signal efficiency [%]")
    h_init_eff_V.Draw("e")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vz"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vz"+name+".root")
    
def plot_effV2(nenner_cut, zahler_cut, cut, name, title):
    #efficiency = globale.overlay.GetEntries(zahler_cut)*100.0/globale.overlay.GetEntries(nenner_cut)
    efficiency = globale.overlay.GetEntries(cut+' && numu_signal')*100.0/(globale.overlay.GetEntries('numu_true')+0.001)
    purity = globale.overlay.GetEntries(cut+' && numu_signal')*globale.scale[globale.overlay]*100/(getTotNum(cut)+0.001)

    print 'Efficiency:\t{0:0.2f}%'.format(efficiency)
    print 'Purity:\t\t{0:0.2f}%'.format(  purity)

    title = title+' Eff/Pur={0:0.2f}%'.format( efficiency)+'/{0:0.2f}%'.format(purity)

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","c1",1600,1200)
    c1.SetGrid(1)
    c1.SetLeftMargin(0.14)
    c1.SetRightMargin(0.18)
    c1.SetBottomMargin(0.14)

    xstart = 0
    xend = 2
    xbin = 100
    h_init_eff_energy = ROOT.TH1F("h_init_eff_energy",title,xbin,xstart,xend)
    h_init_eff_energy_1 = ROOT.TH1F("h_init_eff_energy_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Energy>>h_init_eff_energy',zahler_cut)
    globale.overlay.Draw('MCNu_Energy>>h_init_eff_energy_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_energy,h_init_eff_energy_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")
    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth neutrino energy [GeV]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    #graph1 = graph.Clone()
    #graph.SetFillColor(2);
    #graph.SetFillStyle(1001);
    #graph.Draw("a4");
    #graph1.Draw("same")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_energy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_energy"+name+".root")


    xstart = -1
    xend = 1
    xbin = 100
    h_init_eff_theta = ROOT.TH1F("h_init_eff_theta",title,xbin,xstart,xend)
    h_init_eff_theta_1 = ROOT.TH1F("h_init_eff_theta_1",title,xbin,xstart,xend)
    globale.overlay.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta',zahler_cut)
    globale.overlay.Draw('cos(MCNu_leptonTheta)>>h_init_eff_theta_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_theta,h_init_eff_theta_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth lepton cos(theta)")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_theta"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_theta"+name+".root")

    xstart = -3.14159
    xend = 3.14159
    xbin = 100
    h_init_eff_phi = ROOT.TH1F("h_init_eff_phi",title,xbin,xstart,xend)
    h_init_eff_phi_1 = ROOT.TH1F("h_init_eff_phi_1",title,xbin,xstart,xend)
    globale.overlay.Draw('TrackPhi>>h_init_eff_phi',zahler_cut)
    globale.overlay.Draw('TrackPhi>>h_init_eff_phi_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_phi,h_init_eff_phi_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Reco lepton phi")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_phi"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_phi"+name+".root")

    xstart = -50
    xend = 300
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Vx>>h_init_eff_V',zahler_cut)
    globale.overlay.Draw('MCNu_Vx>>h_init_eff_V_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth nu vertex X [cm]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vx"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vx"+name+".root")

    xstart = -150
    xend = 150
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Vy>>h_init_eff_V',zahler_cut)
    globale.overlay.Draw('MCNu_Vy>>h_init_eff_V_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth nu vertex Y [cm]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vy"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vy"+name+".root")

    xstart = -50
    xend = 1050
    xbin = 100
    h_init_eff_V = ROOT.TH1F("h_init_eff_V",title,xbin,xstart,xend)
    h_init_eff_V_1 = ROOT.TH1F("h_init_eff_V_1",title,xbin,xstart,xend)
    globale.overlay.Draw('MCNu_Vz>>h_init_eff_V',zahler_cut)
    globale.overlay.Draw('MCNu_Vz>>h_init_eff_V_1',nenner_cut)
    eff =  ROOT.TEfficiency(h_init_eff_V,h_init_eff_V_1)
    eff.SetStatisticOption(ROOT.TEfficiency.kFCP)#;  // to set option for errors (see ref doc)
    eff.SetConfidenceLevel(0.68)
    eff.SetTitle(title)
    eff.Draw("AP")

    ROOT.gPad.Update()
    graph = eff.GetPaintedGraph()
    graph.SetMinimum(0)
    graph.SetMaximum(1)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("Truth nu vertex Z [cm]")
    graph.GetYaxis().SetTitle("Signal efficiency")
    graph.Draw("AP")
    c1.Draw()
    c1.SaveAs(globale.outputdir_png + "h_eff_Vz"+name+".png")
    c1.SaveAs(globale.outputdir_root + "h_eff_Vz"+name+".root")
    
def pdg_content_out(cut):
    xstart = -3000
    xend = 3000
    xbin = 6000
    h_1d = ROOT.TH1F("h_1d","Track PDG code",xbin, xstart, xend)
    globale.overlay_out.Draw('MCle_PDG>>h_1d',cut,'')
    h_1d.SetXTitle("Track PDG code")
    h_1d.SetYTitle("number of entries")
    tot = h_1d.GetEntries()
    particle = [(0,0)]
    #num_particle = []
    print 'Total number of entries: ',h_1d.GetEntries()
    for i in range(xbin):
        num = h_1d.GetBinContent(i)
        if num!=0 :
            particle.append((num,i+xstart-1))
            #num_particle.append(num)
            #print 'PDG: ',i+xstart-1,': ',num
    #h_1d.Draw()
    #c1.SetLogy(0)
    #c1.Draw()
    particle.sort() 
    for x in reversed(range(len(particle)-1)): 
        print 'PDG: ', particle[x+1][1],'\t=\t{0:0.1f}%'.format(particle[x+1][0]*100/tot), ',\terror: {0:0.1f}%'.format(math.sqrt(particle[x+1][0])*100/tot),',\tnumber: ',particle[x+1][0]
        
def pdg_content(cut):
    xstart = -3000
    xend = 3000
    xbin = 6000
    h_1d = ROOT.TH1F("h_1d","Track PDG code",xbin, xstart, xend)
    globale.overlay.Draw('MCle_PDG>>h_1d',cut,'')
    h_1d.SetXTitle("Track PDG code")
    h_1d.SetYTitle("number of entries")
    tot = h_1d.GetEntries()
    particle = [(0,0)]
    #num_particle = []
    print 'Total number of entries: ',h_1d.GetEntries()
    for i in range(xbin):
        num = h_1d.GetBinContent(i)
        if num!=0 :
            particle.append((num,i+xstart-1))
            #num_particle.append(num)
            #print 'PDG: ',i+xstart-1,': ',num
    #h_1d.Draw()
    #c1.SetLogy(0)
    #c1.Draw()
    particle.sort() 
    for x in reversed(range(len(particle)-1)): 
        print 'PDG: ', particle[x+1][1],'\t=\t{0:0.1f}%'.format(particle[x+1][0]*100/tot), ',\terror: {0:0.1f}%'.format(math.sqrt(particle[x+1][0])*100/tot),',\tnumber: ',particle[x+1][0]
        
def pdg_content33(cut):
    xstart = -3000
    xend = 3000
    xbin = 6000
    h_1d = ROOT.TH1F("h_1d","Track PDG code",xbin, xstart, xend)
    globale.overlay.Draw('MCTrackPDG>>h_1d',cut,'')
    h_1d.SetXTitle("Track PDG code")
    h_1d.SetYTitle("number of entries")
    tot = h_1d.GetEntries()
    particle = [(0,0)]
    #num_particle = []
    print 'Total number of entries: ',h_1d.GetEntries()
    for i in range(xbin):
        num = h_1d.GetBinContent(i)
        if num!=0 :
            particle.append((num,i+xstart-1))
            #num_particle.append(num)
            #print 'PDG: ',i+xstart-1,': ',num
    #h_1d.Draw()
    #c1.SetLogy(0)
    #c1.Draw()
    particle.sort() 
    for x in reversed(range(len(particle)-1)): 
        print 'PDG: ', particle[x+1][1],'\t=\t{0:0.1f}%'.format(particle[x+1][0]*100/tot), ',\terror: {0:0.1f}%'.format(math.sqrt(particle[x+1][0])*100/tot),',\tnumber: ',particle[x+1][0]

def pdg_content33_out(cut):
    xstart = -3000
    xend = 3000
    xbin = 6000
    h_1d = ROOT.TH1F("h_1d","Track PDG code",xbin, xstart, xend)
    globale.overlay_out.Draw('MCTrackPDG>>h_1d',cut,'')
    h_1d.SetXTitle("Track PDG code")
    h_1d.SetYTitle("number of entries")
    tot = h_1d.GetEntries()
    particle = [(0,0)]
    #num_particle = []
    print 'Total number of entries: ',h_1d.GetEntries()
    for i in range(xbin):
        num = h_1d.GetBinContent(i)
        if num!=0 :
            particle.append((num,i+xstart-1))
            #num_particle.append(num)
            #print 'PDG: ',i+xstart-1,': ',num
    #h_1d.Draw()
    #c1.SetLogy(0)
    #c1.Draw()
    particle.sort() 
    for x in reversed(range(len(particle)-1)): 
        print 'PDG: ', particle[x+1][1],'\t=\t{0:0.1f}%'.format(particle[x+1][0]*100/tot), ',\terror: {0:0.1f}%'.format(math.sqrt(particle[x+1][0])*100/tot),',\tnumber: ',particle[x+1][0]
        
def entry33_out(cut,max_entry):
    xstart = 0
    xend = max_entry
    xbin = max_entry
    h_1d = ROOT.TH1F("h_1d","Track PDG code",xbin, xstart, xend)
    globale.overlay_out.Draw('event_counter>>h_1d',cut,'')
    #h_1d.SetXTitle("Track PDG code")
    #h_1d.SetYTitle("number of entries")
    tot = h_1d.GetEntries()
    particle = [(0,0)]
    #num_particle = []
    print 'Total number of entries: ',h_1d.GetEntries()
    for i in range(xbin):
        num = h_1d.GetBinContent(i)
        if num!=0 :
            particle.append((num,i+xstart-1))
            #num_particle.append(num)
            #print 'PDG: ',i+xstart-1,': ',num
    #h_1d.Draw()
    #c1.SetLogy(0)
    #c1.Draw()
    particle.sort() 
    for x in reversed(range(len(particle)-1)): 
        print 'Event counter: ', particle[x+1][1]-1 #,'\t=\t{0:0.1f}%'.format(particle[x+1][0]*100/tot), ',\terror: {0:0.1f}%'.format(math.sqrt(particle[x+1][0])*100/tot),',\tnumber: ',particle[x+1][0]