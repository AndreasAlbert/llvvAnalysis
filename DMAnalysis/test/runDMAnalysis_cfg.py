import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os
import json

options = VarParsing('analysis')
options.parseArguments()
# Currently analyzer only supports one root file at a time :<
inputFile = options.inputFiles[0]

with open(os.path.expandvars("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/sample_13TeV_25ns.json")) as fin:
    procList = json.load(fin)['proc']
thisDataset = {}
for proc in procList :
    for data in proc['data']:
        if data['dtag'] in inputFile:
            thisDataset = data
            break

isMC = False
if 'MC13TeV' in inputFile:
    isMC = True

doWIMPreweighting = False
genwimpweights = []
if 'doWIMPreweighting' in thisDataset.keys():
    doWIMPreweighting = True
    genwimpweights.append(thisDataset['genwimpweights'])

process = cms.Process('magic')

process.config = cms.PSet(
    input = cms.string(inputFile),
    outdir = cms.string("."),
    doWIMPreweighting = cms.bool(doWIMPreweighting),
    MCweights = cms.vstring([]),
    genwimpweights = cms.vstring(genwimpweights),
    WIMPreweighting_DM_V_Mx = cms.string("root://eoscms//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_ReReco_17Feb2016/MC13TeV_DM_V_Mx1Mv200.root"),
    WIMPreweighting_DM_A_Mx = cms.string("root://eoscms//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_ReReco_17Feb2016/MC13TeV_DM_A_Mx1Mv200.root"),
    usemetNoHF = cms.bool(False),
    pdfInput = cms.string("root://eoscms//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_25ns_PDFs/"),
    PU_Central = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/weight_PU_Central_0_50.root"),
    PU_Up      = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/weight_PU_Up_0_50.root"),
    PU_Down    = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/weight_PU_Down_0_50.root"),
    isMC = cms.bool(isMC),
    mctruthmode=cms.int32(0),
    runSystematics = cms.bool(True),
    runOptimization = cms.bool(False),
    BtagEffFiles = cms.vstring("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/BtagEff_13TeV_10March2016.root"),
    BtagSF = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/CSVv2_76x.csv"),
    MuonTrigEffSF = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/muon_trigeff_sf_eta_eta.root"),
    ElectronTrigEffSF = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/electron_trigeff_sf_abseta_abseta.root"),
    ElectronMediumWPSF = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/CutBasedID_MediumWP_76X_18Feb.txt_SF2D.root"),
    ElectronVetoWPSF   = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/CutBasedID_VetoWP_76X_18Feb.txt_SF2D.root"),
    ElectronRECOSF     = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/eleRECO.txt.egamma_SF2D.root"),
    UnpartQCDscaleWeights = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/UnpartWeights_QCDscale.root"),
    ADDQCDscaleWeights = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/ADDWeights_QCDscale.root"),
    evStart = cms.int32(0),
    evEnd = cms.int32(-1),
    dirName = cms.string("mainAnalyzer/data"),
    jesUncFileName = cms.string("${CMSSW_BASE}/src/llvvAnalysis/DMAnalysis/data/weights/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt"),
)