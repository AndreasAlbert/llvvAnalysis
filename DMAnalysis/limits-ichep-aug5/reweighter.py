#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import glob
import json
import os
import re
import subprocess
import sys

# Won't work for TTreeFormula :<
class ReweightFcn:
    def __init__(self, numerator, denominator, name):
        self.function = ROOT.TF1("reweight"+name, self, 0, 1.e4, 1)
        self.function.SetParName(0, 'genMet')
        self.hist = numerator.Clone("reweight_hist")
        self.hist.Divide(denominator)

    def __call__(self, x, par):
        genmet = x[0]
        res = self.hist.GetBinContent(self.hist.FindBin(genmet))
        if res == 0.:
            res = 1.
        return res
    
    def get_TF1(self):
        return self.function


dmfiles = glob.glob("../run_fullana_20160726_1056/MC13TeV_DM_*.root")
fn2point = re.compile(r".*/MC13TeV_DM_(.)_Mx\-(\d+)_Mv\-(\d+).root")
dmpoints = dict(map(lambda fn: (fn2point.match(fn).groups(), fn), dmfiles))

genwfile = ROOT.TFile("../data/weights/genmet_acc15_simplifiedModels.root")
weightRE = re.compile(r"monoz_genmet_acc15_cV(.)_cA(.)_gDM1_gQ(.*)_Mx(\d*)_Mmed(\d*)")
def weight2point(name):
    m = weightRE.match(name)
    if not m:
        return
    (isV, isA, gQ, mx, mmed) = m.groups()
    if gQ != '1':
        return
    tup1 = 'V'
    if isA == '1':
        tup1 = 'A'
    return (tup1, mx, mmed)
allpoints = map(weight2point, map(lambda k: k.GetName(), genwfile.GetListOfKeys()))
allpoints = filter(lambda p: p is not None, allpoints)
def point2genWeight(point):
    isV = point[0] == 'V'
    isA = point[0] == 'A'
    return "monoz_genmet_acc15_cV%d_cA%d_gDM1_gQ1_Mx%s_Mmed%s" % (isV, isA, point[1], point[2])

point2title = lambda p: "DM(%s)M%s(%s)gQ(1.00)" % (p[1], p[0], p[2])
point2process = lambda p: "dm%sm%s%s" % (p[1], p[0].lower(), p[2])

def chooseReweight(point, n=1):
    distance = lambda p: (float(p[1])-float(point[1]))**2+(float(p[2])-float(point[2]))**2
    rightProcess = (p for p in dmpoints if p[0]==point[0])
    if n>1:
        is5000 = (p for p in rightProcess if p[2]=='5000')
        options = sorted(is5000, key=distance)
        if point == options[0]:
            options.pop(0)
        return options[:min(n,len(options))]
    lowerMmed = (p for p in rightProcess if int(p[2])>int(point[2]))
    if int(point[2]) < 3*int(point[1]):
        lowerMchi = (p for p in lowerMmed if int(p[1])>int(point[1]))
        options = sorted(lowerMchi, key=distance)
    else:
        options = sorted(lowerMmed, key=distance)
    if len(options) == 0:
        print "Point has no reweight options with mass conditions, trying closest"
        rightProcess = (p for p in dmpoints if p[0]==point[0])
        options = sorted(rightProcess, key=distance)
    return options[min(0,len(options)-1)]

def reweightString(variable, point, reweightPoint):
    numerator = genwfile.Get(point2genWeight(point))
    if not numerator:
        print "Could not find weights for point", point2title(point)
    denominator = genwfile.Get(point2genWeight(reweightPoint))
    if not denominator:
        print "Could not find weights for point", point2title(reweightPoint)
    # New plan, let acceptance change!
    # numerator.SetNormFactor(1.)
    # denominator.SetNormFactor(1.)
    ratioHist = numerator.Clone("reweight_hist")
    ratioHist.Divide(denominator)
    ternaries = []
    for i in range(1, ratioHist.GetNbinsX()+1):
        binEdge = ratioHist.GetXaxis().GetBinUpEdge(i)
        value = ratioHist.GetBinContent(i)
        if value == 0.:
            value = 1.
        ternaries.append("(%s<%s)?%f" % (variable, binEdge, value))
    ternaries.append("1.")
    return "(" + ":".join(ternaries) + ")"

with open("sample_xs.json") as fin:
    sample_xs = json.load(fin)
fiducial_xs = {}
with open("xs_fiducial.txt") as fin:
    for line in fin:
        print line.strip().split(' ')
        (point, xs) = line.strip().split(' ')
        fiducial_xs[point] = float(xs)
luminosity = 12900.

systematics = {
    '{process}_CMS_res_jUp' : '_jerup',
    '{process}_CMS_res_jDown' : '_jerdown',
    '{process}_CMS_scale_jUp' : '_jesup',
    '{process}_CMS_scale_jDown' : '_jesdown',
    '{process}_CMS_scale_{lep}Up' : '_lesup',
    '{process}_CMS_scale_{lep}Down' : '_lesdown',
    '{process}_CMS_zllwimps_puUp' : '_puup',
    '{process}_CMS_zllwimps_puDown' : '_pudown',
    '{process}_CMS_eff_bUp' : '_btagup',
    '{process}_CMS_eff_bDown' : '_btagdown',
    '{process}_pdf_qqbarUp' : '_pdfup',
    '{process}_pdf_qqbarDown' : '_pdfdown',
    '{process}_QCDscale_VDMUp' : '_qcdscaleup',
    '{process}_QCDscale_VDMDown' : '_qcdscaledown'
}
variations = ['', '_jerup', '_jerdown', '_jesup', '_jesdown', '_umetup', '_umetdown', '_lesup', '_lesdown', '_puup', '_pudown', '_btagup', '_btagdown', '_pdfup', '_pdfdown', '_qcdscaleup', '_qcdscaledown']
channelNames = ['', 'mumulesq1jets', 'eelesq1jets']
def drawShapes(point, templateHist, outputDir, reweightPoint=None):
    if type(reweightPoint) is list:
        dmfiles = []
        dmtrees = []
        reweights = []
        for rwPoint in reweightPoint:
            dmfile = ROOT.TFile(dmpoints[rwPoint])
            if not dmfile:
                continue
            dmtree = dmfile.Get("dmReweight")
            if not dmtree:
                continue
            dmfiles.append(dmfile)
            dmtrees.append(dmtree)
            reweights.append("*"+reweightString("genmet", point, rwPoint))
    elif reweightPoint:
        dmfiles = [ROOT.TFile(dmpoints[reweightPoint])]
        dmtrees = map(lambda dmfile: dmfile.Get("dmReweight"), dmfiles)
        reweights = ["*"+reweightString("genmet", point, reweightPoint)]
        reweightPoint = [reweightPoint]
    else:
        dmfiles = [ROOT.TFile(dmpoints[point])]
        dmtrees = map(lambda dmfile: dmfile.Get("dmReweight"), dmfiles)
        reweights = [""]
    normHists = map(lambda dmfile: dmfile.Get("all_cutflow"), dmfiles)
    eventNorms = []
    if reweightPoint is not None:
        for normHist, rwPoint in zip(normHists, reweightPoint):
            weight = luminosity/(normHist.GetBinContent(1)*normHist.GetBinContent(3)/normHist.GetBinContent(2))
            weight *= fiducial_xs[point2genWeight(point).replace('monoz_genmet_acc15_','')] * sample_xs[point2title(rwPoint)] / fiducial_xs[point2genWeight(rwPoint).replace('monoz_genmet_acc15_','')]
            eventNorms.append(weight)
    else:
        eventNorms = map(lambda normHist: sample_xs[point2title(point)]*luminosity/(normHist.GetBinContent(1)*normHist.GetBinContent(3)/normHist.GetBinContent(2)), normHists)
    puUpNorms   = map(lambda normHist: normHist.GetBinContent(3)/normHist.GetBinContent(5), normHists)
    puDownNorms = map(lambda normHist: normHist.GetBinContent(3)/normHist.GetBinContent(4), normHists)
    tmpHist = templateHist.Clone("temp")
    tmpHist.SetTitle(point2title(point))
    def writeHist(drawCmd, cutCmd, outputHistName, pu=0):
        for i, dmfile in enumerate(dmfiles):
            dmfile.cd()
            tmpHist.SetDirectory(dmfile)
            plus = '+'
            if i == 0:
                plus = ''
            norm = eventNorms[i]
            if pu==1:
                norm = norm * puUpNorms[i]
            elif pu==-1:
                norm = norm * puDownNorms[i]
            dmtrees[i].Draw(drawCmd+">>%stemp" % plus, cutCmd+reweights[i]+"*%f" % norm, "goff")
        tmpHist.Scale(1./len(dmfiles))
        outputHist = tmpHist.Clone(outputHistName)
        outputDir.cd()
        outputHist.Write()
        return outputHist

    channel = outputDir.GetName()
    evcat = channelNames.index(channel)
    processName = point2process(point)

    nominal = writeHist("pfmet[0]", "(evcat==%d && nJets<2)*weight[0]" % (evcat), processName)
    outputDir.cd()

    nonzeroBins = filter(lambda b: nominal.GetBinContent(b)>0, range(1,nominal.GetNbinsX()+1))
    zeroBins = filter(lambda b: nominal.GetBinContent(b)==0, range(1,nominal.GetNbinsX()+1))
    for i in nonzeroBins:
        nameStatUp = '{process}_CMS_zllwimps_stat_{channel}_{process}_13TeV_Bin{bin}Up'.format(process=processName, channel=channelNames[evcat], bin=i)
        statUp = nominal.Clone(nameStatUp)
        statUp.SetBinContent(i, statUp.GetBinContent(i)+statUp.GetBinError(i))
        statUp.Write()
        nameStatDown = '{process}_CMS_zllwimps_stat_{channel}_{process}_13TeV_Bin{bin}Down'.format(process=processName, channel=channelNames[evcat], bin=i)
        statDown = nominal.Clone(nameStatDown)
        statDown.SetBinContent(i, statDown.GetBinContent(i)-statDown.GetBinError(i))
        statDown.Write()

    shapes_to_remove = []
    for i in zeroBins:
        shapes_to_remove.append('CMS_zllwimps_stat_{channel}_{process}_13TeV_Bin{bin}'.format(process=processName, channel=channelNames[evcat], bin=i))

    for systName, variation in systematics.iteritems():
        name = systName.format(process=processName, lep=channelNames[evcat][0])
        ivar = variations.index(variation)
        pu = 0
        if variation == '_puup':
            pu = 1
        elif variation == '_pudown':
            pu = -1
        writeHist("pfmet[%d]" % ivar, "(evcat==%d && nJets<2)*weight[%d]" % (evcat, ivar), name, pu)

    return (nominal.Integral(), shapes_to_remove)


template_shapes = ROOT.TFile("template/zllwimps_50_13TeV.root")
with open("template/card_combined.dat") as fin:
    template_card = fin.read()

def getLimit(point, nReweight=1, noCombine=False):
    output_dir = point2process(point)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_shapes = ROOT.TFile(output_dir+"/shapes.root", "recreate")
    output_card = open(output_dir+"/card.dat", "w")
    if nReweight == 0:
        reweightPoint = None
    else:
        reweightPoint = chooseReweight(point, nReweight)
        print "Reweighting ", repr(point)
        print "Using points: ", repr(reweightPoint)
    channels = ['eelesq1jets', 'mumulesq1jets']
    cardInfo = {'signalName': point2process(point), 'filename':"shapes.root"}
    shapes_to_remove = []
    for channel in channels:
        inHists = map(lambda k: (k.GetName(), k.ReadObj()), template_shapes.Get(channel).GetListOfKeys())
        copyable = filter(lambda (n, h): 'dm1ma500' not in n, inHists)
        # To match key name to hist name since runLimits writes hists with different key names :<
        map(lambda (n,h): h.SetName(n), copyable)
        outDir = output_shapes.mkdir(channel)
        outDir.cd()
        map(lambda (n,h): h.Write(), copyable)
        (rate, shapes) = drawShapes(point, copyable[0][1], outDir, reweightPoint)
        shapes_to_remove.extend(shapes)
        if channel == 'eelesq1jets':
            cardInfo['eeRate'] = rate
        elif channel == 'mumulesq1jets':
            cardInfo['mmRate'] = rate

    cardInfo['kmax'] = 81 - len(shapes_to_remove)
    for line in template_card.format(**cardInfo).split("\n"):
        shape = line.split(' ')[0]
        if shape in shapes_to_remove:
            continue
        output_card.write(line+"\n")
    output_card.close()
    output_shapes.Close()

    if noCombine:
        return

    combine = subprocess.check_output(["combine", "-M", "Asymptotic", "--cl", "0.95", output_dir+"/card.dat"], stderr=subprocess.STDOUT)
    with open(output_dir+"/combine.log", "w") as fout:
        fout.write(combine)

    limits = {}
    for line in combine.split("\n"):
        if 'Observed' in line:
            limits['obs'] = float(line.split('<')[-1])
        elif 'Expected' in line:
            (percent, value) = line.replace('Expected','').split('%: r <')
            limits[percent.strip()] = float(value)
    return limits

getLimit(('V', '1', '20'), 1)
infos = []
fullsimResults = {('A', '10', '10'): {'obs': 0.2716, '16.0': 0.1638, '84.0': 0.3348, '50.0': 0.2334, '2.5': 0.1199, '97.5': 0.4622}, ('V', '10', '20'): {'obs': 0.1192, '16.0': 0.0665, '84.0': 0.1389, '50.0': 0.0952, '2.5': 0.0487, '97.5': 0.1936}, ('A', '50', '50'): {'obs': 0.7987, '16.0': 0.4985, '84.0': 1.0397, '50.0': 0.7168, '2.5': 0.3654, '97.5': 1.4538}, ('V', '10', '10'): {'obs': 0.183, '16.0': 0.1035, '84.0': 0.2134, '50.0': 0.1479, '2.5': 0.0754, '97.5': 0.2956}, ('V', '1', '300'): {'obs': 0.374, '16.0': 0.2357, '84.0': 0.4915, '50.0': 0.3389, '2.5': 0.1727, '97.5': 0.6915}, ('A', '500', '2000'): {'obs': 87.7863, '16.0': 44.9093, '84.0': 106.0646, '50.0': 68.25, '2.5': 31.5923, '97.5': 159.1797}, ('V', '50', '200'): {'obs': 0.2439, '16.0': 0.1487, '84.0': 0.3085, '50.0': 0.2139, '2.5': 0.109, '97.5': 0.43}, ('V', '1000', '1995'): {'obs': 205.0778, '16.0': 103.4501, '84.0': 247.7663, '50.0': 157.8125, '2.5': 72.4335, '97.5': 376.5195}, ('A', '150', '200'): {'obs': 3.5761, '16.0': 2.1346, '84.0': 4.5728, '50.0': 3.1016, '2.5': 1.5447, '97.5': 6.5276}, ('V', '10', '100'): {'obs': 0.139, '16.0': 0.0826, '84.0': 0.1704, '50.0': 0.1182, '2.5': 0.0605, '97.5': 0.2361}, ('A', '1', '300'): {'obs': 0.3701, '16.0': 0.22, '84.0': 0.4629, '50.0': 0.3174, '2.5': 0.1606, '97.5': 0.6532}, ('A', '1', '100'): {'obs': 0.4597, '16.0': 0.2874, '84.0': 0.5888, '50.0': 0.4082, '2.5': 0.2113, '97.5': 0.8156}, ('A', '1', '2000'): {'obs': 42.0493, '16.0': 22.0442, '84.0': 51.6007, '50.0': 33.375, '2.5': 15.5793, '97.5': 76.7395}, ('V', '50', '50'): {'obs': 0.4729, '16.0': 0.2999, '84.0': 0.6198, '50.0': 0.4297, '2.5': 0.2207, '97.5': 0.8639}, ('A', '1', '500'): {'obs': 0.8664, '16.0': 0.4987, '84.0': 1.0683, '50.0': 0.7246, '2.5': 0.3609, '97.5': 1.525}, ('V', '500', '2000'): {'obs': 51.8878, '16.0': 26.7318, '84.0': 63.4576, '50.0': 40.625, '2.5': 18.8049, '97.5': 94.9042}, ('V', '500', '10'): {'obs': 45.4192, '16.0': 24.4453, '84.0': 55.659, '50.0': 36.375, '2.5': 17.406, '97.5': 82.4911}, ('V', '150', '10'): {'obs': 2.5403, '16.0': 1.4771, '84.0': 3.1504, '50.0': 2.1484, '2.5': 1.0742, '97.5': 4.485}, ('V', '1', '100'): {'obs': 0.1418, '16.0': 0.0822, '84.0': 0.1688, '50.0': 0.1177, '2.5': 0.0602, '97.5': 0.2345}, ('V', '1', '10'): {'obs': 0.0659, '16.0': 0.0379, '84.0': 0.0835, '50.0': 0.0562, '2.5': 0.0276, '97.5': 0.1193}, ('V', '1', '500'): {'obs': 0.9576, '16.0': 0.5531, '84.0': 1.187, '50.0': 0.8008, '2.5': 0.402, '97.5': 1.6892}, ('V', '150', '500'): {'obs': 1.0555, '16.0': 0.621, '84.0': 1.3232, '50.0': 0.9023, '2.5': 0.4494, '97.5': 1.8837}, ('A', '1', '10'): {'obs': 0.0657, '16.0': 0.0407, '84.0': 0.0854, '50.0': 0.0581, '2.5': 0.0295, '97.5': 0.12}, ('V', '1', '20'): {'obs': 0.0605, '16.0': 0.0441, '84.0': 0.0919, '50.0': 0.063, '2.5': 0.032, '97.5': 0.1289}, ('V', '500', '995'): {'obs': 12.2644, '16.0': 6.5292, '84.0': 14.8887, '50.0': 9.7812, '2.5': 4.6423, '97.5': 22.0259}, ('A', '10', '100'): {'obs': 0.1572, '16.0': 0.087, '84.0': 0.1796, '50.0': 0.1245, '2.5': 0.0637, '97.5': 0.2488}, ('V', '1', '2000'): {'obs': 41.8005, '16.0': 21.7377, '84.0': 50.181, '50.0': 32.625, '2.5': 15.4204, '97.5': 74.5029}, ('A', '50', '200'): {'obs': 0.2992, '16.0': 0.1981, '84.0': 0.4093, '50.0': 0.2822, '2.5': 0.145, '97.5': 0.5689}, ('V', '150', '200'): {'obs': 1.5954, '16.0': 0.9124, '84.0': 1.9303, '50.0': 1.3164, '2.5': 0.6659, '97.5': 2.7319}, ('A', '50', '95'): {'obs': 0.5169, '16.0': 0.328, '84.0': 0.679, '50.0': 0.4707, '2.5': 0.2409, '97.5': 0.9464}, ('V', '1', '200'): {'obs': 0.2526, '16.0': 0.149, '84.0': 0.3085, '50.0': 0.2139, '2.5': 0.1094, '97.5': 0.4273}, ('V', '150', '295'): {'obs': 0.8029, '16.0': 0.4616, '84.0': 0.9766, '50.0': 0.666, '2.5': 0.3369, '97.5': 1.3821}, ('A', '10', '20'): {'obs': 0.1614, '16.0': 0.1071, '84.0': 0.2217, '50.0': 0.1528, '2.5': 0.0782, '97.5': 0.3081}, ('V', '50', '95'): {'obs': 0.2576, '16.0': 0.1554, '84.0': 0.3212, '50.0': 0.2227, '2.5': 0.1144, '97.5': 0.4449}, ('V', '1', '50'): {'obs': 0.1012, '16.0': 0.058, '84.0': 0.119, '50.0': 0.0825, '2.5': 0.0429, '97.5': 0.1649}, ('A', '1', '200'): {'obs': 0.2478, '16.0': 0.151, '84.0': 0.313, '50.0': 0.2158, '2.5': 0.1109, '97.5': 0.4377}, ('A', '10', '50'): {'obs': 0.0944, '16.0': 0.0658, '84.0': 0.1359, '50.0': 0.0942, '2.5': 0.0482, '97.5': 0.1883}, ('A', '150', '10'): {'obs': 4.8185, '16.0': 2.7807, '84.0': 6.0451, '50.0': 4.0781, '2.5': 2.0152, '97.5': 8.652}, ('V', '50', '300'): {'obs': 0.3616, '16.0': 0.22, '84.0': 0.4629, '50.0': 0.3174, '2.5': 0.1606, '97.5': 0.6532}, ('A', '50', '10'): {'obs': 0.8955, '16.0': 0.5515, '84.0': 1.1438, '50.0': 0.793, '2.5': 0.4042, '97.5': 1.5944}, ('A', '1', '20'): {'obs': 0.085, '16.0': 0.0535, '84.0': 0.1119, '50.0': 0.0771, '2.5': 0.0389, '97.5': 0.1555}, ('A', '150', '295'): {'obs': 2.1144, '16.0': 1.2356, '84.0': 2.6519, '50.0': 1.7891, '2.5': 0.898, '97.5': 3.752}, ('A', '500', '10'): {'obs': 136.1339, '16.0': 71.4775, '84.0': 163.3428, '50.0': 106.75, '2.5': 50.6646, '97.5': 242.0871}, ('A', '500', '500'): {'obs': 104.8858, '16.0': 54.882, '84.0': 127.0868, '50.0': 82.625, '2.5': 39.0532, '97.5': 188.6837}, ('V', '50', '10'): {'obs': 0.5951, '16.0': 0.3552, '84.0': 0.7394, '50.0': 0.5098, '2.5': 0.2609, '97.5': 1.0276}, ('A', '50', '300'): {'obs': 0.446, '16.0': 0.2565, '84.0': 0.5398, '50.0': 0.3701, '2.5': 0.1872, '97.5': 0.7571}, ('V', '1', '1000'): {'obs': 4.4249, '16.0': 2.4029, '84.0': 5.3943, '50.0': 3.5625, '2.5': 1.7186, '97.5': 7.88}, ('A', '150', '500'): {'obs': 1.513, '16.0': 0.9087, '84.0': 1.9466, '50.0': 1.3203, '2.5': 0.6576, '97.5': 2.7788}, ('A', '500', '995'): {'obs': 44.5251, '16.0': 23.5189, '84.0': 54.0263, '50.0': 35.125, '2.5': 16.6707, '97.5': 79.7996}}
for point, fullsim in fullsimResults.iteritems():
    reweight = getLimit(point, 1)
    if '50.0' not in reweight:
        continue
    diff = (reweight['50.0']-fullsim['50.0'])*2/(fullsim['84.0']-fullsim['16.0'])
    print diff
    infos.append((point, fullsim, reweight))

print repr(infos)
h = ROOT.TH1F("sigmas2", ";#Delta 95% CL Limit [#sigma];Points", 15, -3, 3)
sigma = lambda p: (p[2]['50.0']-p[1]['50.0'])*2/(p[1]['84.0']-p[1]['16.0'])
map(lambda x: h.Fill(sigma(x)), infos)
h.Draw()
ROOT.gPad.Print("~/reweight.root")
