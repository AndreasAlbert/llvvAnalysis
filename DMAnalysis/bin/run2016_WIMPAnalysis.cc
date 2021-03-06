#include <iostream>

//#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "llvvAnalysis/DMAnalysis/interface/MacroUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/DataEvtSummaryHandler.h"
#include "llvvAnalysis/DMAnalysis/interface/DMPhysicsEvent.h"
#include "llvvAnalysis/DMAnalysis/interface/SmartSelectionMonitor.h"
#include "llvvAnalysis/DMAnalysis/interface/METUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/BTagUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/WIMPReweighting.h"
#include "llvvAnalysis/DMAnalysis/interface/EventCategory.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"


#include "llvvAnalysis/DMAnalysis/interface/PDFInfo.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
#include "llvvAnalysis/DMAnalysis/interface/BTagCalibrationStandalone.h"

#include "llvvAnalysis/DMAnalysis/interface/LeptonEfficiencySF.h"

using namespace std;


namespace LHAPDF {
void initPDFSet(int nset, const std::string& filename, int member=0);
int numberPDF(int nset);
void usePDFMember(int nset, int member);
double xfx(int nset, double x, double Q, int fl);
double getXmin(int nset, int member);
double getXmax(int nset, int member);
double getQ2min(int nset, int member);
double getQ2max(int nset, int member);
void extrapolate(bool extrapolate=true);
}

struct stPDFval {
    stPDFval() {}
    stPDFval(const stPDFval& arg) :
        qscale(arg.qscale),
        x1(arg.x1),
        x2(arg.x2),
        id1(arg.id1),
        id2(arg.id2) {
    }

    double qscale;
    double x1;
    double x2;
    int id1;
    int id2;
};

struct DMTree {
  int evcat;
  int nJets;
  float genmet;
  float weight[20];
  float pfmet[20];
};

//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation76X
const float CSVLooseWP = 0.460;
const float CSVMediumWP = 0.800;
const float CSVTightWP = 0.935;

void dumpList( std::vector<LorentzVector> particles, std::string tag ) {
    std::cout << std::endl  << "--- DUMP " << tag << std::endl;
    for (auto const & part : particles ) {
        std::cout << part.pt() << " " << part.phi() << " " << part.eta() << std::endl;
    }
}
void dumpList( PhysicsObjectLeptonCollection particles, std::string tag ) {
    std::cout << std::endl  << "--- DUMP " << tag << std::endl;
    for (auto const & part : particles ) {
        std::cout << part.pt() << " " << part.phi() << " " << part.eta() << std::endl;
    }
}
void dumpList( PhysicsObjectJetCollection particles, std::string tag ) {
    std::cout << std::endl  << "--- DUMP " << tag << std::endl;
    for (auto const & part : particles ) {
        std::cout << part.pt() << " " << part.phi() << " " << part.eta() << " " << part.isPFLoose << std::endl;
    }
}
int main(int argc, char* argv[])
{
    //##################################################################################
    //##########################    GLOBAL INITIALIZATION     ##########################
    //##################################################################################

    // check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }

    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    //AutoLibraryLoader::enable();
    FWLiteEnabler::enable();

    // configure the process
    std::shared_ptr<edm::ParameterSet> config = edm::readConfig(argv[1], argc, argv);
    const edm::ParameterSet &runProcess = config->getParameter<edm::ParameterSet>("config");

    bool isMC       = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");
    bool doWIMPreweighting = runProcess.getParameter<bool>("doWIMPreweighting");
    bool usemetNoHF = runProcess.getParameter<bool>("usemetNoHF");

    TString url=runProcess.getParameter<std::string>("input");
    TString outFileUrl(gSystem->BaseName(url));
    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir=runProcess.getParameter<std::string>("outdir");
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);

    TString outTxtUrl_final= outUrl + "/" + outFileUrl + "_FinalList.txt";
    FILE* outTxtFile_final = NULL;
    outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl_final.Data());
    fprintf(outTxtFile_final,"run lumi event\n");

    //save control plots to file
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results saved in %s\n", outUrl.Data());
    //save all to the file
    TFile *ofile=TFile::Open(outUrl, "recreate");

    bool isSingleMuPD(!isMC && url.Contains("SingleMuon"));
    bool isDoubleMuPD(!isMC && url.Contains("DoubleMuon"));
    bool isSingleElePD(!isMC && url.Contains("SingleElectron"));
    bool isDoubleElePD(!isMC && url.Contains("DoubleEG"));
    bool isMuEGPD(!isMC && url.Contains("MuonEG"));
    if (!isMC && !(isSingleMuPD||isSingleElePD||isDoubleMuPD||isDoubleElePD||isMuEGPD) )
      cout << "WARNING: Data sample comes from an unrecognized primary dataset.  Please check filenames!" << endl;

    bool isMC_ZZ2L2Nu  = isMC && ( string(url.Data()).find("TeV_ZZTo2L2Nu")  != string::npos);
    bool isMC_ZZTo4L   = isMC && ( string(url.Data()).find("TeV_ZZTo4L")  != string::npos);
    bool isMC_ZZTo2L2Q = isMC && ( string(url.Data()).find("TeV_ZZTo2L2Q")  != string::npos);

    bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZamcatnloFXFX")  != string::npos
                              || string(url.Data()).find("MC13TeV_WZpowheg")  != string::npos );

    bool isMC_VVV = isMC && ( string(url.Data()).find("MC13TeV_WZZ")  != string::npos
                              || string(url.Data()).find("MC13TeV_WWZ")  != string::npos
                              || string(url.Data()).find("MC13TeV_ZZZ")  != string::npos );

    bool isMCBkg_runPDFQCDscale = (isMC_ZZ2L2Nu || isMC_ZZTo4L || isMC_ZZTo2L2Q ||
                                   isMC_WZ || isMC_VVV);

    bool isMC_ttbar = isMC && (string(url.Data()).find("TeV_TT")  != string::npos);
    bool isMC_stop  = isMC && (string(url.Data()).find("TeV_SingleT")  != string::npos);
    bool isMC_WIMP  = isMC && (string(url.Data()).find("TeV_DM_V_Mx") != string::npos
                               || string(url.Data()).find("TeV_DM_A_Mx") != string::npos
                               || string(url.Data()).find("TeV_EWKDM_S_Mx") != string::npos
                               || string(url.Data()).find("TeV_EWKDM_P_Mx") != string::npos);
    bool isMC_ADD  = isMC && (string(url.Data()).find("TeV_ADD_D") != string::npos);
    bool isMC_Unpart = isMC && (string(url.Data()).find("TeV_Unpart") != string::npos);


    bool isSignal = (isMC_WIMP || isMC_ADD || isMC_Unpart);

    // Save DM reweight tree for faster interpolation
    TTree * dmReweightTree{nullptr};
    DMTree dmtree;
    if (isSignal && !doWIMPreweighting) {
      dmReweightTree = new TTree("dmReweight", "DM model reweighting tree");
      dmReweightTree->Branch("evcat", &dmtree.evcat);
      dmReweightTree->Branch("nJets", &dmtree.nJets);
      dmReweightTree->Branch("genmet", &dmtree.genmet);
      dmReweightTree->Branch("pfmet", &dmtree.pfmet, "pfmet[20]/F");
      dmReweightTree->Branch("weight", &dmtree.weight, "weight[20]/F");
    }


    BTagUtils myBtagUtils(runProcess);


    // https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
    // setup calibration readers
    TString BtagSF = runProcess.getParameter<std::string>("BtagSF");
    gSystem->ExpandPathName(BtagSF);
    cout << "Loading btag sacle factor: " << BtagSF << endl;
    BTagCalibration btagcalib("csvv2", BtagSF.Data());
    BTagCalibrationReader btag_reader(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "central");
    BTagCalibrationReader btag_reader_up(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "up");  // sys up
    BTagCalibrationReader btag_reader_down(&btagcalib, BTagEntry::OP_MEDIUM, "mujets", "down");  // sys down



    WIMPReweighting myWIMPweights;

    //systematics
    bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
    std::vector<TString> varNames(1,"");
    if(runSystematics) {
        cout << "Systematics will be computed for this analysis: " << endl;
        varNames.push_back("_jerup"); 	//1
        varNames.push_back("_jerdown"); //2
        varNames.push_back("_jesup"); 	//3
        varNames.push_back("_jesdown"); //4
        varNames.push_back("_umetup"); 	//5
        varNames.push_back("_umetdown");//6
        varNames.push_back("_lesup"); 	//7
        varNames.push_back("_lesdown"); //8
        varNames.push_back("_puup"); 	//9
        varNames.push_back("_pudown"); 	//10
        varNames.push_back("_btagup"); 	//11
        varNames.push_back("_btagdown");//12
        if(isMCBkg_runPDFQCDscale) {
            varNames.push_back("_pdfup");
            varNames.push_back("_pdfdown");
            varNames.push_back("_qcdscaleup");
            varNames.push_back("_qcdscaledown");

        }
        if(isSignal) {
            varNames.push_back("_pdfup");
            varNames.push_back("_pdfdown");
            varNames.push_back("_qcdscaleup");
            varNames.push_back("_qcdscaledown");
//            varNames.push_back("_pdfacceptup");
//            varNames.push_back("_pdfacceptdown");
//            varNames.push_back("_qcdscaleacceptup");
//            varNames.push_back("_qcdscaleacceptdown");
        }
        if(isMC_ZZ2L2Nu) {
            varNames.push_back("_qqZZewkup");
            varNames.push_back("_qqZZewkdown");
        }

        for(size_t sys=1; sys<varNames.size(); sys++) {
            cout << varNames[sys] << endl;
        }
    }
    size_t nvarsToInclude=varNames.size();


    //tree info
    int evStart     = runProcess.getParameter<int>("evStart");
    int evEnd       = runProcess.getParameter<int>("evEnd");
    TString dirname = runProcess.getParameter<std::string>("dirName");

    //jet energy scale uncertainties
    TString uncFile = runProcess.getParameter<std::string>("jesUncFileName");
    gSystem->ExpandPathName(uncFile);
    cout << "Loading jet energy scale uncertainties: " << uncFile << endl;
    JetCorrectionUncertainty jecUnc(uncFile.Data());

    //pdf info
    PDFInfo *mPDFInfo=0;
    if(isMC_Unpart || isMC_ADD) {
        TString pdfUrl = runProcess.getParameter<std::string>("pdfInput");
        std::string Url = runProcess.getParameter<std::string>("input");
        std::size_t found = Url.find_last_of("/\\");
        pdfUrl += '/';
        pdfUrl += Url.substr(found+1);
        pdfUrl.ReplaceAll(".root","_pdf.root");

        mPDFInfo=new PDFInfo(pdfUrl,"NNPDF30_lo_as_0130.LHgrid");//cteq66.LHgrid");
        cout << "Readout " << mPDFInfo->numberPDFs() << " pdf variations: " << pdfUrl << endl;
    }


    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;


    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("scalefactors", ";;Sum Weights", 10,0,10) );
    h->GetXaxis()->SetBinLabel(1,"Lep1 ID");
    h->GetXaxis()->SetBinLabel(2,"Lep2 ID");
    h->GetXaxis()->SetBinLabel(3,"Trigger");
    h->GetXaxis()->SetBinLabel(4,"Product");
    h=(TH1F*) mon.addHistogram( new TH1F ("scalefactors_final", ";;Sum Weights", 10,0,10) );
    h->GetXaxis()->SetBinLabel(1,"Lep1 ID");
    h->GetXaxis()->SetBinLabel(2,"Lep2 ID");
    h->GetXaxis()->SetBinLabel(3,"Trigger");
    h->GetXaxis()->SetBinLabel(4,"Product");

    int nbin = 1;
    h=(TH1F*) mon.addHistogram( new TH1F ("sync_cutflow", ";;Events", 12,0,12) );
    h->GetXaxis()->SetBinLabel(nbin++,"Trigger && 2 leptons");
    h->GetXaxis()->SetBinLabel(nbin++,"No third lepton");
    h->GetXaxis()->SetBinLabel(nbin++,"b veto");
    h->GetXaxis()->SetBinLabel(nbin++,"#tau veto");
    h->GetXaxis()->SetBinLabel(nbin++,"#lesq 1 jets"); 
    h->GetXaxis()->SetBinLabel(nbin++,"76.1876<#it{m}_{ll}<101.1876");
    h->GetXaxis()->SetBinLabel(nbin++,"#it{p}_{T}^{ll}>60");
    h->GetXaxis()->SetBinLabel(nbin++,"MET>100");
    h->GetXaxis()->SetBinLabel(nbin++,"#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss})>2.8");
    h->GetXaxis()->SetBinLabel(nbin++,"|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T}<0.4");
    h->GetXaxis()->SetBinLabel(nbin++,"#Delta#it{#phi}(jet,E_{T}^{miss})>0.5");
    h->GetXaxis()->SetBinLabel(nbin++,"|u_{ll}/p_{T}^{ll}|<1");


    h=(TH1F*) mon.addHistogram( new TH1F ("DMAcceptance", ";;Events", 3,0,3) );
    h->GetXaxis()->SetBinLabel(1,"Trigger && 2 Tight Leptons");
    h->GetXaxis()->SetBinLabel(2,"Trigger && 2 Tight Leptons && #geq 1 Leptons in MEX/1");
    h->GetXaxis()->SetBinLabel(3,"Trigger && 2 Tight Leptons && =2 Leptons in MEX/1");

    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;
    //for WIMP reweighting norms
    TH1F* sumWIMPweights  = (TH1F*) mon.addHistogram( new TH1F ("sumWIMPweights"    , "sumWIMPweights"    ,2,0,2) ) ;

    mon.addHistogram( new TH1F( "nvtx_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "zpt_raw",      ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "mt_raw",       ";#it{m}_{T} [GeV];Events", 100,0,2000) );
    double MTBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000,2000};
    const int nBinsMT = sizeof(MTBins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "mt2_raw",       ";#it{m}_{T} [GeV];Events", nBinsMT,MTBins) );
    mon.addHistogram( new TH1F( "zmass_raw",    ";#it{m}_{ll} [GeV];Events", 100,40,250) );

    mon.addHistogram( new TH2F( "ptlep1vs2_raw",";#it{p}_{T}^{l1} [GeV];#it{p}_{T}^{l2} [GeV];Events",250,0,500, 250,0,500) );

    mon.addHistogram( new TH1F( "leadlep_pt_raw", ";Leading lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "leadlep_eta_raw",";Leading lepton #eta^{l};Events", 52,-2.6,2.6) );
    mon.addHistogram( new TH1F( "trailep_pt_raw", ";Trailing lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "trailep_eta_raw",";Trailing lepton #eta^{l};Events", 52,-2.6,2.6) );


    mon.addHistogram( new TH1F( "jet_pt_raw", ";all jet #it{p}_{T}^{j};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "jet_eta_raw",";all jet #eta^{j};Events", 52,-2.6,2.6) );


    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nleptons_raw", ";Lepton multiplicity;Events", 3,2,5) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="=";
        label += (ibin+1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h2 = (TH1F *)mon.addHistogram( new TH1F("njets_raw",  ";Jet multiplicity (#it{p}_{T}>30 GeV);Events",5,0,5) );
    for(int ibin=1; ibin<=h2->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h2->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h2->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h3 = (TH1F *)mon.addHistogram( new TH1F("nbjets_raw",    ";b-tag Jet multiplicity;Events",5,0,5) );
    for(int ibin=1; ibin<=h3->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h3->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h3->GetXaxis()->SetBinLabel(ibin,label);
    }

   //// SYNCHRONIZATION
    mon.addHistogram( new TH1F( "sync_nlep_minus",  ";Veto lepton multiplicity (#it{p}_{T}>30 GeV);Events",5,0,5) );
    mon.addHistogram( new TH1F( "sync_njet_minus",  ";Jet multiplicity (#it{p}_{T}>30 GeV);Events",5,0,5) );
    mon.addHistogram( new TH1F( "sync_zmass_minus",    ";#it{m}_{ll} [GeV];Events", 100,40,250) );
    mon.addHistogram( new TH1F( "sync_ptz_minus",      ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );

    mon.addHistogram( new TH1F( "leadlep_pt_met_minus", ";Leading lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "leadlep_eta_met_minus",";Leading lepton #eta^{l};Events", 52,-2.6,2.6) );
    mon.addHistogram( new TH1F( "trailep_pt_met_minus", ";Leading lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "trailep_eta_met_minus",";Leading lepton #eta^{l};Events", 52,-2.6,2.6) );

    mon.addHistogram( new TH1F( "sync_met_minus",        ";E_{T}^{miss} [GeV];Events", 50,0,500 ) );
    mon.addHistogram( new TH1F( "sync_dPhiZMET_minus",   ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "sync_balance_minus",    ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events",10,0,2.0) );
    mon.addHistogram( new TH1F( "sync_dPhiJetMET_minus", ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "sync_response_minus",   ";Response(#it{l^{+}l^{-}},E_{T}^{miss});Events", 20,-1,1 ) );

    double MET2Bins[]= {0,100,180,260,340,420,500,580,800,1200};
    const int nBinsMET2 = sizeof(MET2Bins)/sizeof(double) - 1;

    double MT2Bins[]= {0,100,200,300,400,500,600,700,800,1000,1200};
    const int nBinsMT2 = sizeof(MT2Bins)/sizeof(double) - 1;

    //// Track MET
    mon.addHistogram( new TH1F( "dPhiTrackMET",";#Delta#it{#phi}(E_{T}^{miss}(trk),E_{T}^{miss});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH2F( "dPhiTrackMET_v_pfmet",";E_{T}^{miss} [GeV];#Delta#it{#phi}(E_{T}^{miss}(trk),E_{T}^{miss});Events",50,0,500, 10,0,TMath::Pi()) );

 /*
    // preselection plots
    double METBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};
    const int nBinsMET = sizeof(METBins)/sizeof(double) - 1;
    double METBins2[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000};
    const int nBinsMET2 = sizeof(METBins2)/sizeof(double) - 1;
    //double METBins3[]= {0,40,50,60,120,260,1000};
//    double METBins3[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000};
//    const int nBinsMET3 = sizeof(METBins3)/sizeof(double) - 1;
    double MET2Bins[]= {0,80,160,240,320,400,480,560,640,800,1200};
    const int xnBinsMET2 = sizeof(MET2Bins)/sizeof(double) - 1;
    double MT2Bins[]= {0,100,200,300,400,500,600,700,800,1000,1200};
    const int xnBinsMT2 = sizeof(MT2Bins)/sizeof(double) - 1;
*/
    mon.addHistogram( new TH1F( "pfmet2_presel",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, MET2Bins));
    mon.addHistogram( new TH1F( "dphiZMET_presel",   ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "dphiJetMET_presel",   ";#Delta#it{#phi}(#it{jet},E_{T}^{miss});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "balancedif_presel", ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events", 5,0,1.0) );
/*
    //adding N-1 plots
    mon.addHistogram( new TH1F( "pfmet_nm1",       ";E_{T}^{miss} [GeV];Events / 80 GeV", 15,0,1200));
    mon.addHistogram( new TH1F( "pfmet2_nm1",      ";E_{T}^{miss} [GeV];Events / 1 GeV", xnBinsMET2,MET2Bins));
    mon.addHistogram( new TH1F( "mt_nm1",          ";#it{m}_{T} [GeV];Events / 100 GeV", 12,0,1200));
    mon.addHistogram( new TH1F( "mt2_nm1",         ";#it{m}_{T} [GeV];Events / 1 GeV", xnBinsMT2,MT2Bins));
    mon.addHistogram( new TH1F( "pfmet_met80",       ";E_{T}^{miss} [GeV];Events / 80 GeV", 15,0,1200));
    mon.addHistogram( new TH1F( "pfmet2_met80",      ";E_{T}^{miss} [GeV];Events / 1 GeV", xnBinsMET2,MET2Bins));
    mon.addHistogram( new TH1F( "mt_met80",          ";#it{m}_{T} [GeV];Events / 100 GeV", 12,0,1200));
    mon.addHistogram( new TH1F( "mt2_met80",         ";#it{m}_{T} [GeV];Events / 1 GeV", xnBinsMT2,MT2Bins));
    mon.addHistogram( new TH1F( "pfmet_met140",       ";E_{T}^{miss} [GeV];Events / 80 GeV", 15,0,1200));
    mon.addHistogram( new TH1F( "pfmet2_met140",      ";E_{T}^{miss} [GeV];Events / 1 GeV", xnBinsMET2,MET2Bins));
    mon.addHistogram( new TH1F( "mt_met140",          ";#it{m}_{T} [GeV];Events / 100 GeV", 12,0,1200));
    mon.addHistogram( new TH1F( "mt2_met140",         ";#it{m}_{T} [GeV];Events / 1 GeV", xnBinsMT2,MT2Bins));
    //MET X-Y shift correction
    mon.addHistogram( new TH2F( "pfmetx_vs_nvtx",";Vertices;E_{X}^{miss} [GeV];Events",50,0,50, 200,-75,75) );
    mon.addHistogram( new TH2F( "pfmety_vs_nvtx",";Vertices;E_{Y}^{miss} [GeV];Events",50,0,50, 200,-75,75) );
    mon.addHistogram( new TH1F( "pfmetphi_wocorr",";#it{#phi}(E_{T}^{miss});Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "pfmetphi_wicorr",";#it{#phi}(E_{T}^{miss});Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "pfmet_wicorr",      ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET, METBins));
    mon.addHistogram( new TH1F( "pfmet2_wicorr",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, METBins2));
    // generator level plots
    mon.addHistogram( new TH1F( "pileup", ";pileup;Events", 50,0,50) );
    mon.addHistogram( new TH1F( "met_Gen", ";#it{p}_{T}(#bar{#chi}#chi) [GeV];Events", nBinsMET, METBins) );
    mon.addHistogram( new TH1F( "met2_Gen", ";#it{p}_{T}(#bar{#chi}#chi) [GeV];Events", 500, 0,1000) );
    mon.addHistogram( new TH1F( "zpt_Gen", ";#it{p}_{T}(Z) [GeV];Events", 800,0,800) );
    mon.addHistogram( new TH1F( "dphi_Gen", ";#Delta#phi(Z,#bar{#chi}#chi) [rad];Events", 100,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "zmass_Gen", ";#it{m}_{ll} [GeV] [GeV];Events", 250,0,250) );
    mon.addHistogram( new TH2F( "ptlep1vs2_Gen",";#it{p}_{T}^{l1} [GeV];#it{p}_{T}^{l2} [GeV];Events",250,0,500, 250,0,500) );
    h=(TH1F *)mon.addHistogram( new TH1F ("acceptance", ";;Events", 2,0,2) );
    h->GetXaxis()->SetBinLabel(1,"Gen");
    h->GetXaxis()->SetBinLabel(2,"Gen Acc");
*/
    // btaging efficiency
    std::vector<TString> CSVkey;
    CSVkey.push_back("CSVL");
    CSVkey.push_back("CSVM");
    CSVkey.push_back("CSVT");
/*
    double JetPTBins[]= {20,40,60,80,100,120,150,300,500};
    const int nJetPTBins = sizeof(JetPTBins)/sizeof(double) - 1;
    for(size_t csvtag=0; csvtag<CSVkey.size(); csvtag++) {
        mon.addHistogram( new TH1F( TString("beff_Denom_")+CSVkey[csvtag],      "; Jet #it{p}_{T} [GeV];Events", nJetPTBins,JetPTBins) );
        mon.addHistogram( new TH1F( TString("ceff_Denom_")+CSVkey[csvtag],      "; Jet #it{p}_{T} [GeV];Events", nJetPTBins,JetPTBins) );
        mon.addHistogram( new TH1F( TString("udsgeff_Denom_")+CSVkey[csvtag],   "; Jet #it{p}_{T} [GeV];Events", nJetPTBins,JetPTBins) );
        mon.addHistogram( new TH1F( TString("beff_Num_")+CSVkey[csvtag],        "; Jet #it{p}_{T} [GeV];Events", nJetPTBins,JetPTBins) );
        mon.addHistogram( new TH1F( TString("ceff_Num_")+CSVkey[csvtag],        "; Jet #it{p}_{T} [GeV];Events", nJetPTBins,JetPTBins) );
        mon.addHistogram( new TH1F( TString("udsgeff_Num_")+CSVkey[csvtag],     "; Jet #it{p}_{T} [GeV];Events", nJetPTBins,JetPTBins) );
    }
*/

    mon.addHistogram( new TH1F( "pfmet2_final",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, MET2Bins));

    //#################################################
    //############# CONTROL PLOTS #####################
    //#################################################
    // WW control plots, for k-method (for emu channel)
    mon.addHistogram( new TH1F( "zpt_wwctrl_raw",   ";#it{p}_{T}^{ll} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "zmass_wwctrl_raw", ";#it{m}_{ll} [GeV];Events", 100,20,300) );
    mon.addHistogram( new TH1F( "pfmet_wwctrl_raw", ";E_{T}^{miss} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "mt_wwctrl_raw",";#it{m}_{T}(#it{ll}, E_{T}^{miss}) [GeV];Events", 50,0,300) );


    //WZ control
    mon.addHistogram( new TH1F( "pfmet_WZctrl",         ";E_{T}^{miss} [GeV];Events", 20,0,500));
    mon.addHistogram( new TH1F( "balancedif_WZctrl",    ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events", 20,0,1.0) );
    mon.addHistogram( new TH1F( "DphiZMET_WZctrl",      ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 20,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "zpt_WZctrl",           ";#it{p}_{T}^{ll} [GeV];Events", 20,0,300) );

    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET40",    ";Fake E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET40",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET50",    ";Fake E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET50",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET60",    ";Fake E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET60",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET70",    ";Fake E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET70",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET80",    ";Fake E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET80",       ";#it{m}_{T} [GeV];Events", 20,0,800) );



    //##################################################################################
    //########################## STUFF FOR CUTS OPTIMIZATION  ##########################
    //##################################################################################
    //optimization
    std::vector<double> optim_Cuts1_MET;
    std::vector<double> optim_Cuts1_Balance;
    std::vector<double> optim_Cuts1_DphiZMET;

    bool runOptimization = runProcess.getParameter<bool>("runOptimization");
    if(runOptimization) {
        // for optimization
        cout << "Optimization will be performed for this analysis" << endl;
        for(double met=100; met<=150; met+=10) {
            for(double balance=0.4; balance<=0.4; balance+=0.05) {
                for(double dphi=2.8; dphi<=2.8; dphi+=0.1) {
                    optim_Cuts1_MET     .push_back(met);
                    optim_Cuts1_Balance .push_back(balance);
                    optim_Cuts1_DphiZMET.push_back(dphi);
                }
            }
        }
    }

    size_t nOptims = optim_Cuts1_MET.size();


    //make it as a TProfile so hadd does not change the value
    TProfile* Hoptim_cuts1_MET      = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_MET",";cut index;met",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_Balance  = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_Balance",";cut index;balance",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_DphiZMET = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_DphiZMET",";cut index;dphi",nOptims,0,nOptims) );

    for(size_t index=0; index<nOptims; index++) {
        Hoptim_cuts1_MET        ->Fill(index, optim_Cuts1_MET[index]);
        Hoptim_cuts1_Balance    ->Fill(index, optim_Cuts1_Balance[index]);
        Hoptim_cuts1_DphiZMET   ->Fill(index, optim_Cuts1_DphiZMET[index]);
    }

    TH1F* Hoptim_systs  =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;


    //for extrapolation of DY process based on MET
    mon.addHistogram( new TH2F ("pfmet_minus_shapes",";cut index; E_{T}^{miss} [GeV];#Events ",nOptims,0,nOptims, 100,0,500) );



    for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
        Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);

        //1D shapes for limit setting
        mon.addHistogram( new TH2F (TString("mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",nOptims,0,nOptims,  12,0,1200) );
        mon.addHistogram( new TH2F (TString("mt2_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events",nOptims,0,nOptims, nBinsMT2,MT2Bins) );

        mon.addHistogram( new TH2F (TString("met_shapes")+varNames[ivar],";cut index; E_{T}^{miss} [GeV];#Events (/80GeV)",nOptims,0,nOptims, 15,0,1200) );
        mon.addHistogram( new TH2F (TString("met2_shapes")+varNames[ivar],";cut index; E_{T}^{miss} [GeV];#Events",nOptims,0,nOptims, nBinsMET2,MET2Bins) );

        //2D shapes for limit setting
        //
        //
    }

    //##################################################################################
    //#############         GET READY FOR THE EVENT LOOP           #####################
    //##################################################################################

    //open the file and get events tree
    DataEvtSummaryHandler summaryHandler_;
    if(doWIMPreweighting) {
        // Only for simplified models
        if(url.Contains("TeV_DM_")) {
          bool isreweighted = myWIMPweights.Init(runProcess, url);
          if(!isreweighted) {
            cerr << " *** WARNING: WIMP re-weighting initialization failed! ***" << endl;
          }
        }

        if(url.Contains("K1_0.1_K2_1")) url.ReplaceAll("K1_0.1_K2_1","K1_1_K2_1");
        if(url.Contains("K1_0.2_K2_1")) url.ReplaceAll("K1_0.2_K2_1","K1_1_K2_1");
        if(url.Contains("K1_0.3_K2_1")) url.ReplaceAll("K1_0.3_K2_1","K1_1_K2_1");
        if(url.Contains("K1_0.5_K2_1")) url.ReplaceAll("K1_0.5_K2_1","K1_1_K2_1");
        if(url.Contains("K1_2_K2_1"))   url.ReplaceAll("K1_2_K2_1","K1_1_K2_1");
        if(url.Contains("K1_3_K2_1"))   url.ReplaceAll("K1_3_K2_1","K1_1_K2_1");
        if(url.Contains("K1_5_K2_1"))   url.ReplaceAll("K1_5_K2_1","K1_1_K2_1");
        if(url.Contains("K1_10_K2_1"))  url.ReplaceAll("K1_10_K2_1","K1_1_K2_1");
    }
    TFile *file = TFile::Open(url);
    printf("Looping on %s\n",url.Data());
    if(file==0) return -1;
    if(file->IsZombie()) return -1;
    if( !summaryHandler_.attachToTree( (TTree *)file->Get(dirname) ) ) {
        file->Close();
        return -1;
    }


    //check run range to compute scale factor (if not all entries are used)
    const Int_t totalEntries= summaryHandler_.getEntries();
    float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
    if(evEnd<0 || evEnd>summaryHandler_.getEntries() ) evEnd=totalEntries;
    if(evStart > evEnd ) {
        file->Close();
        return -1;
    }

    //MC normalization (to 1/pb)
    float cnorm=1.0;
    if(isMC) {
        //TH1F* cutflowH = (TH1F *) file->Get("mainAnalyzer/llvv/nevents");
        //if(cutflowH) cnorm=cutflowH->GetBinContent(1);
        TH1F* posH = (TH1F *) file->Get("mainAnalyzer/llvv/n_posevents");
        TH1F* negH = (TH1F *) file->Get("mainAnalyzer/llvv/n_negevents");
        if(posH && negH) cnorm = posH->GetBinContent(1) - negH->GetBinContent(1);
        if(rescaleFactor>0) cnorm /= rescaleFactor;
        printf("cnorm = %f\n",cnorm);
    }
    Hcutflow->SetBinContent(1,cnorm);



    //pileup weighting
    TString PU_Central = runProcess.getParameter<std::string>("PU_Central");
    gSystem->ExpandPathName(PU_Central);
    cout << "Loading PU weights Central: " << PU_Central << endl;
    TFile *PU_Central_File = TFile::Open(PU_Central);
    TH1F* weight_pileup_Central = (TH1F *) PU_Central_File->Get("pileup");

    TString PU_Up = runProcess.getParameter<std::string>("PU_Up");
    gSystem->ExpandPathName(PU_Up);
    cout << "Loading PU weights Up: " << PU_Up << endl;
    TFile *PU_Up_File = TFile::Open(PU_Up);
    TH1F* weight_pileup_Up = (TH1F *) PU_Up_File->Get("pileup");

    TString PU_Down = runProcess.getParameter<std::string>("PU_Down");
    gSystem->ExpandPathName(PU_Down);
    cout << "Loading PU weights Down: " << PU_Down << endl;
    TFile *PU_Down_File = TFile::Open(PU_Down);
    TH1F* weight_pileup_Down = (TH1F *) PU_Down_File->Get("pileup");



    // muon trigger efficiency SF
    TString MuonTrigEffSF_ = runProcess.getParameter<std::string>("MuonTrigEffSF");
    gSystem->ExpandPathName(MuonTrigEffSF_);
    cout << "Loading Muon Trigger Eff SF: " << MuonTrigEffSF_ << endl;
    TFile *MuonTrigEffSF_File = TFile::Open(MuonTrigEffSF_);
    TH2F* h_MuonTrigEffSF_0 = (TH2F *) MuonTrigEffSF_File->Get("h_dimuon_eff_0");
    TH2F* h_MuonTrigEffSF_1 = (TH2F *) MuonTrigEffSF_File->Get("h_dimuon_eff_1");
    TH2F* h_MuonTrigEffSF_2 = (TH2F *) MuonTrigEffSF_File->Get("h_dimuon_eff_2");
    TH2F* h_MuonTrigEffSF_3 = (TH2F *) MuonTrigEffSF_File->Get("h_dimuon_eff_3");
    // electron trigger efficiency SF
    TString ElectronTrigEffSF_ = runProcess.getParameter<std::string>("ElectronTrigEffSF");
    gSystem->ExpandPathName(ElectronTrigEffSF_);
    cout << "Loading Electron Trigger Eff SF: " << ElectronTrigEffSF_ << endl;
    TFile *ElectronTrigEffSF_File = TFile::Open(ElectronTrigEffSF_);
    TH2F* h_ElectronTrigEffSF_0 = (TH2F *) ElectronTrigEffSF_File->Get("h_dielectron_eff_0");
    TH2F* h_ElectronTrigEffSF_1 = (TH2F *) ElectronTrigEffSF_File->Get("h_dielectron_eff_1");
    TH2F* h_ElectronTrigEffSF_2 = (TH2F *) ElectronTrigEffSF_File->Get("h_dielectron_eff_2");
    TH2F* h_ElectronTrigEffSF_3 = (TH2F *) ElectronTrigEffSF_File->Get("h_dielectron_eff_3");

    //Electron ID RECO SF
    TString ElectronMediumWPSF_ = runProcess.getParameter<std::string>("ElectronMediumWPSF");
    gSystem->ExpandPathName(ElectronMediumWPSF_);
    cout << "Loading Electron MediumWP SF: " << ElectronMediumWPSF_ << endl;
    TFile *ElectronMediumWPSF_File = TFile::Open(ElectronMediumWPSF_);
    TH2F* h_ElectronMediumWPSF = (TH2F *) ElectronMediumWPSF_File->Get("scalefactors_Medium_Electron"); // POG name: EGamma_SF2D

    TString ElectronRECOSF_ = runProcess.getParameter<std::string>("ElectronRECOSF");
    gSystem->ExpandPathName(ElectronRECOSF_);
    cout << "Loading Electron RECO SF: " << ElectronRECOSF_ << endl;
    TFile *ElectronRECOSF_File = TFile::Open(ElectronRECOSF_);
    TH2F* h_ElectronRECOSF = (TH2F *) ElectronRECOSF_File->Get("EGamma_SF2D");

    // Muon ID SF
    TString MuonTightWPSF_ = runProcess.getParameter<std::string>("MuonTightWPSF");
    gSystem->ExpandPathName(MuonTightWPSF_);
    cout << "Loading Muon TightWP SF: " << MuonTightWPSF_ << endl;
    TFile *MuonTightWPSF_File = TFile::Open(MuonTightWPSF_);
    TH2F* h_MuonTightWPSF = (TH2F *) MuonTightWPSF_File->Get("scalefactors_Tight_Muon");


    //
    TString UnpartQCDscaleWeights_ = runProcess.getParameter<std::string>("UnpartQCDscaleWeights");
    gSystem->ExpandPathName(UnpartQCDscaleWeights_);
    cout << "Loading UnpartQCDscaleWeights: " << UnpartQCDscaleWeights_ << endl;
    TFile *UnpartQCDscaleWeights_File = TFile::Open(UnpartQCDscaleWeights_);


    TString ADDQCDscaleWeights_ = runProcess.getParameter<std::string>("ADDQCDscaleWeights");
    gSystem->ExpandPathName(ADDQCDscaleWeights_);
    cout << "Loading ADDQCDscaleWeights: " << ADDQCDscaleWeights_ << endl;
    TFile *ADDQCDscaleWeights_File = TFile::Open(ADDQCDscaleWeights_);


    // event categorizer
    EventCategory eventCategoryInst(1);   //jet(0,1,>=2) binning

    //####################################################################################################################
    //###########################################           EVENT LOOP         ###########################################
    //####################################################################################################################


    // loop on all the events
    int treeStep = (evEnd-evStart)/50;
    if(treeStep==0)treeStep=1;
    printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    printf("Scanning the ntuple :");

    for( int iev=evStart; iev<evEnd; iev++) {
        if((iev-evStart)%treeStep==0) {
            printf("#");
            fflush(stdout);
        }

        //##############################################   EVENT LOOP STARTS   ##############################################
        //load the event content from tree
        summaryHandler_.getEntry(iev);
        DataEvtSummary_t &ev=summaryHandler_.getEvent();

        //prepare the tag's vectors for histo filling
        std::vector<TString> tags(1,"all");

        //genWeight
        float genWeight = 1.0;
        if(isMC && ev.genWeight<0) genWeight = -1.0;

        //systematical weight
        float weight = 1.0;
        if(isMC) weight *= genWeight;
        //only take up and down from pileup effect
        double TotalWeight_plus = 1.0;
        double TotalWeight_minus = 1.0;

        if(isMC) mon.fillHisto("pileup", "all", ev.ngenTruepu, 1.0);

        if(isMC) {
            weight 		*= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Central);
            TotalWeight_plus 	*= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Up);
            TotalWeight_minus 	*= getSFfrom1DHist(ev.ngenTruepu, weight_pileup_Down);
        }


        Hcutflow->Fill(1,genWeight);
        Hcutflow->Fill(2,weight);
        Hcutflow->Fill(3,weight*TotalWeight_minus);
        Hcutflow->Fill(4,weight*TotalWeight_plus);


        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        // FIXME need to have a function: loop all leptons, find a Z candidate,
        // can have input, ev.mn, ev.en
        // assign ee,mm,emu channel
        // check if channel name is consistent with trigger
        // store dilepton candidate in lep1 lep2, and other leptons in 3rdleps


        bool hasMMtrigger = ev.triggerType & 0x1;
        bool hasMtrigger  = (ev.triggerType >> 1 ) & 0x1;
        bool hasEEtrigger = (ev.triggerType >> 2 ) & 0x1;
        bool hasEtrigger  = (ev.triggerType >> 3 ) & 0x1;
        bool hasEMtrigger = (ev.triggerType >> 4 ) & 0x1;



        //#########################################################################
        //####################  Generator Level Reweighting  ######################
        //#########################################################################


        //for Wimps
        double weight_QCDscale_muR1_muF2(0.);
        double weight_QCDscale_muR1_muF0p5(0.);
        double weight_QCDscale_muR2_muF1(0.);
        double weight_QCDscale_muR2_muF2(0.);
        double weight_QCDscale_muR0p5_muF1(0.);
        double weight_QCDscale_muR0p5_muF0p5(0.);

        if(isMC_WIMP || isMC_ADD || isMC_Unpart || isMC_ZZ2L2Nu) {
            if(phys.genleptons.size()!=2) continue;
            if(phys.genGravitons.size()!=1 && phys.genWIMPs.size()!=2 && phys.genneutrinos.size()!=2) continue;

            LorentzVector genmet(0,0,0,0);
            if(phys.genWIMPs.size()==2) 	  genmet = phys.genWIMPs[0]+phys.genWIMPs[1];
            else if(phys.genGravitons.size()==1)  genmet = phys.genGravitons[0];
            else if(phys.genneutrinos.size()==2)  genmet = phys.genneutrinos[0]+phys.genneutrinos[1];

            if ( dmReweightTree != nullptr ) dmtree.genmet = genmet.pt();

            LorentzVector dilep = phys.genleptons[0]+phys.genleptons[1];
            double dphizmet = fabs(deltaPhi(dilep.phi(),genmet.phi()));

            // WIMP interpolation reweighting
            sumWIMPweights->Fill(0., 1.);
            if(doWIMPreweighting) {
                double wimpWeight = 1.;
                if(url.Contains("TeV_DM_")) wimpWeight *= myWIMPweights.get1DWeights(genmet.pt(),"genmet_acc_simplmod");
                if(url.Contains("TeV_EWKDM_S_Mx")) wimpWeight *= myWIMPweights.get1DWeights(genmet.pt(),"pt_chichi");
                sumWIMPweights->Fill(1., wimpWeight);
                weight *= wimpWeight;
            } else {
                sumWIMPweights->Fill(1., 1.);
            }

            mon.fillHisto("met_Gen", tags, genmet.pt(), weight, true);
            mon.fillHisto("met2_Gen", tags, genmet.pt(), weight);
            mon.fillHisto("zpt_Gen", tags, dilep.pt(), weight);
            mon.fillHisto("dphi_Gen", tags, dphizmet, weight);
            mon.fillHisto("zmass_Gen", tags, dilep.mass(), weight);
            if(phys.genleptons[0].pt() > phys.genleptons[1].pt()) mon.fillHisto("ptlep1vs2_Gen", tags, phys.genleptons[0].pt(), phys.genleptons[1].pt(), weight);
            else mon.fillHisto("ptlep1vs2_Gen", tags, phys.genleptons[1].pt(), phys.genleptons[0].pt(), weight);


            bool isInGenAcceptance = (phys.genleptons[0].pt()>20 && phys.genleptons[1].pt()>20);
            isInGenAcceptance &= ((fabs(dilep.mass()-91)<15) && (genmet.pt()>80) && (dilep.pt()>50) && (dphizmet>2.7));
            isInGenAcceptance &= (fabs(genmet.pt()-dilep.pt())/dilep.pt() < 0.2);

            mon.fillHisto("acceptance",tags,0,1.0);
            if(isInGenAcceptance) mon.fillHisto("acceptance",tags,1,1.0);

            if(isMC_Unpart || isMC_ADD) {

                TString tag_unpartadd="";
                if     (string(url.Data()).find("TeV_Unpart_1p01") != string::npos) tag_unpartadd = "Unpart_dU_1.01";
                else if(string(url.Data()).find("TeV_Unpart_1p02") != string::npos) tag_unpartadd = "Unpart_dU_1.02";
                else if(string(url.Data()).find("TeV_Unpart_1p04") != string::npos) tag_unpartadd = "Unpart_dU_1.04";
                else if(string(url.Data()).find("TeV_Unpart_1p06") != string::npos) tag_unpartadd = "Unpart_dU_1.06";
                else if(string(url.Data()).find("TeV_Unpart_1p09") != string::npos) tag_unpartadd = "Unpart_dU_1.09";
                else if(string(url.Data()).find("TeV_Unpart_1p10") != string::npos) tag_unpartadd = "Unpart_dU_1.1";
                else if(string(url.Data()).find("TeV_Unpart_1p20") != string::npos) tag_unpartadd = "Unpart_dU_1.2";
                else if(string(url.Data()).find("TeV_Unpart_1p30") != string::npos) tag_unpartadd = "Unpart_dU_1.3";
                else if(string(url.Data()).find("TeV_Unpart_1p40") != string::npos) tag_unpartadd = "Unpart_dU_1.4";
                else if(string(url.Data()).find("TeV_Unpart_1p50") != string::npos) tag_unpartadd = "Unpart_dU_1.5";
                else if(string(url.Data()).find("TeV_Unpart_1p60") != string::npos) tag_unpartadd = "Unpart_dU_1.6";
                else if(string(url.Data()).find("TeV_Unpart_1p70") != string::npos) tag_unpartadd = "Unpart_dU_1.7";
                else if(string(url.Data()).find("TeV_Unpart_1p80") != string::npos) tag_unpartadd = "Unpart_dU_1.8";
                else if(string(url.Data()).find("TeV_Unpart_1p90") != string::npos) tag_unpartadd = "Unpart_dU_1.9";
                else if(string(url.Data()).find("TeV_Unpart_2p00") != string::npos) tag_unpartadd = "Unpart_dU_2.0";
                else if(string(url.Data()).find("TeV_Unpart_2p20") != string::npos) tag_unpartadd = "Unpart_dU_2.2";
                else if(string(url.Data()).find("TeV_ADD_D2") != string::npos)  tag_unpartadd = "ADD_d_2";
                else if(string(url.Data()).find("TeV_ADD_D3") != string::npos)  tag_unpartadd = "ADD_d_3";
                else if(string(url.Data()).find("TeV_ADD_D4") != string::npos)  tag_unpartadd = "ADD_d_4";
                else if(string(url.Data()).find("TeV_ADD_D5") != string::npos)  tag_unpartadd = "ADD_d_5";
                else if(string(url.Data()).find("TeV_ADD_D6") != string::npos)  tag_unpartadd = "ADD_d_6";
                else if(string(url.Data()).find("TeV_ADD_D7") != string::npos)  tag_unpartadd = "ADD_d_7";
                else if(string(url.Data()).find("TeV_ADD_D8") != string::npos)  tag_unpartadd = "ADD_d_8";

                TH1F* h_muR1_muF2 = NULL;
                TH1F* h_muR1_muF0p5 = NULL;
                TH1F* h_muR2_muF1 = NULL;
                TH1F* h_muR2_muF2 = NULL;
                TH1F* h_muR0p5_muF1 = NULL;
                TH1F* h_muR0p5_muF0p5 = NULL;

                if(isMC_Unpart) {
                    h_muR1_muF2             = (TH1F *) UnpartQCDscaleWeights_File->Get(tag_unpartadd+"_muR_1.0_muF_2.0");
                    h_muR1_muF0p5           = (TH1F *) UnpartQCDscaleWeights_File->Get(tag_unpartadd+"_muR_1.0_muF_0.5");
                    h_muR2_muF1             = (TH1F *) UnpartQCDscaleWeights_File->Get(tag_unpartadd+"_muR_2.0_muF_1.0");
                    h_muR2_muF2             = (TH1F *) UnpartQCDscaleWeights_File->Get(tag_unpartadd+"_muR_2.0_muF_2.0");
                    h_muR0p5_muF1           = (TH1F *) UnpartQCDscaleWeights_File->Get(tag_unpartadd+"_muR_0.5_muF_1.0");
                    h_muR0p5_muF0p5         = (TH1F *) UnpartQCDscaleWeights_File->Get(tag_unpartadd+"_muR_0.5_muF_0.5");
                } else if(isMC_ADD) {
                    h_muR1_muF2             = (TH1F *) ADDQCDscaleWeights_File->Get(tag_unpartadd+"_muR_1.0_muF_2.0");
                    h_muR1_muF0p5           = (TH1F *) ADDQCDscaleWeights_File->Get(tag_unpartadd+"_muR_1.0_muF_0.5");
                    h_muR2_muF1             = (TH1F *) ADDQCDscaleWeights_File->Get(tag_unpartadd+"_muR_2.0_muF_1.0");
                    h_muR2_muF2             = (TH1F *) ADDQCDscaleWeights_File->Get(tag_unpartadd+"_muR_2.0_muF_2.0");
                    h_muR0p5_muF1           = (TH1F *) ADDQCDscaleWeights_File->Get(tag_unpartadd+"_muR_0.5_muF_1.0");
                    h_muR0p5_muF0p5         = (TH1F *) ADDQCDscaleWeights_File->Get(tag_unpartadd+"_muR_0.5_muF_0.5");
                }

                weight_QCDscale_muR1_muF2    = getSFfrom1DHist(genmet.pt(),h_muR1_muF2);
                weight_QCDscale_muR1_muF0p5  = getSFfrom1DHist(genmet.pt(),h_muR1_muF0p5);
                weight_QCDscale_muR2_muF1    = getSFfrom1DHist(genmet.pt(),h_muR2_muF1);
                weight_QCDscale_muR2_muF2    = getSFfrom1DHist(genmet.pt(),h_muR2_muF2);
                weight_QCDscale_muR0p5_muF1  = getSFfrom1DHist(genmet.pt(),h_muR0p5_muF1);
                weight_QCDscale_muR0p5_muF0p5 = getSFfrom1DHist(genmet.pt(),h_muR0p5_muF0p5);

                delete h_muR1_muF2;
                delete h_muR1_muF0p5;
                delete h_muR2_muF1;
                delete h_muR2_muF2;
                delete h_muR0p5_muF1;
                delete h_muR0p5_muF0p5;


            }

        }


        double weight_ewkup(1.0);
        double weight_ewkdown(1.0);

        if(isMC_ZZ2L2Nu) {
            if(phys.genleptons.size()!=2) continue;
            if(phys.genneutrinos.size()!=2) continue;
            double pt_dilep = (phys.genleptons[0]+phys.genleptons[1]).pt();
            double pt_zvv = (phys.genneutrinos[0]+phys.genneutrinos[1]).pt();
            double trailing_pt = pt_dilep < pt_zvv ? pt_dilep : pt_zvv;
            // Variable to assess the hadronic activity in the event
            float rhoZZ = (phys.genleptons[0]+phys.genleptons[1]+phys.genneutrinos[0]+phys.genneutrinos[1]).pt();
            rhoZZ /= (phys.genleptons[0].pt() + phys.genleptons[1].pt() + phys.genneutrinos[0].pt() + phys.genneutrinos[1].pt());

            // ewk correction ZZ
            double qqZZ_EWKNLO = getNLOEWKZZWeight(trailing_pt);

            // NNLO/NLO k-factor
            double dPhiZZ = deltaPhi((phys.genleptons[0]+phys.genleptons[1]).Phi(), (phys.genneutrinos[0]+phys.genneutrinos[1]).Phi());
            double qqZZ_NNLO = kfactor_qqZZ_qcd_dPhi(dPhiZZ);

            //double qqZZwgt = (weight * qqZZ_EWKNLO * qqZZ_NNLO); //qqZZ
            weight *= qqZZ_EWKNLO * qqZZ_NNLO;

            weight_ewkup   = 1 /  (rhoZZ<0.3 ? (1. - (1.6 - 1.)*(1. - qqZZ_EWKNLO)) : qqZZ_EWKNLO);
            weight_ewkdown = 1 *  (rhoZZ<0.3 ? (1. - (1.6 - 1.)*(1. - qqZZ_EWKNLO)) : qqZZ_EWKNLO);
        }

        // WZ NNLO QCD k-factor
        if (isMC_WZ) { weight *= 1.1 ; }

        //#########################################################################
        //#####################      Objects Selection       ######################
        //#########################################################################

        //
        // MET ANALYSIS
        //
        //apply Jet Energy Resolution corrections to jets (and compute associated variations on the MET variable)
        std::vector<PhysicsObjectJetCollection> variedJets;
        LorentzVectorCollection variedMET;
        METUtils::computeVariation(phys.jets, phys.leptons, (usemetNoHF ? phys.metNoHF : phys.met), variedJets, variedMET, &jecUnc);

        LorentzVector metP4=variedMET[0];
        LorentzVector trkMETP4 = phys.trkMET;
        double dPhiTrkMET = fabs(deltaPhi(trkMETP4.phi(), metP4.phi()));

        PhysicsObjectJetCollection &corrJets = variedJets[0];

        //
        // LEPTON ANALYSIS
        //

        // looping leptons (electrons + muons)
        int nGoodLeptons(0);
        std::vector<std::pair<int,LorentzVector> > goodLeptons;
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            int lepid = phys.leptons[ilep].id;
            if(lep.pt()<20) continue;
            if(abs(lepid)==13 && fabs(lep.eta())> 2.4) continue;
            if(abs(lepid)==11 && fabs(lep.eta())> 2.5) continue;

            bool hasTightIdandIso(true);
            if(abs(lepid)==13) { //muon
                hasTightIdandIso &= phys.leptons[ilep].isTightMu;
                hasTightIdandIso &= ( fabs(phys.leptons[ilep].mn_dZ) < 0.1 );
                hasTightIdandIso &= ( fabs(phys.leptons[ilep].mn_d0) < 0.02 );

                //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2?sortcol=1;table=7;up=0#Muon_Isolation
                hasTightIdandIso &= ( phys.leptons[ilep].m_pfRelIsoDbeta() < 0.15 );
            } else if(abs(lepid)==11) { //electron
                hasTightIdandIso &= phys.leptons[ilep].isElpassMedium;
            } else continue;


            if(!hasTightIdandIso) continue;
            nGoodLeptons++;
            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(lepid,lep);
            goodLeptons.push_back(goodlep);

        }

        if(nGoodLeptons<2) continue; // 2 tight leptons

        float _MASSDIF_(999.);
        int id1(0),id2(0);
        LorentzVector lep1(0,0,0,0),lep2(0,0,0,0);
        for(size_t ilep=0; ilep<goodLeptons.size(); ilep++) {
            int id1_ = goodLeptons[ilep].first;
            LorentzVector lep1_ = goodLeptons[ilep].second;

            for(size_t jlep=ilep+1; jlep<goodLeptons.size(); jlep++) {
                int id2_ = goodLeptons[jlep].first;
                LorentzVector lep2_ = goodLeptons[jlep].second;
                if(id1_*id2_>0) continue; // opposite charge

                LorentzVector dilepton=lep1_+lep2_;
                double massdif = fabs(dilepton.mass()-91.);
                if(massdif < _MASSDIF_) {
                    _MASSDIF_ = massdif;
                    lep1.SetPxPyPzE(lep1_.px(),lep1_.py(),lep1_.pz(),lep1_.energy());
                    lep2.SetPxPyPzE(lep2_.px(),lep2_.py(),lep2_.pz(),lep2_.energy());
                    id1 = id1_;
                    id2 = id2_;
                }
            }
        }






        if(id1*id2==0) continue;
        LorentzVector zll(lep1+lep2);
        bool passZmass( ( 76.1876 < zll.mass() ) && ( zll.mass() < 101.1876 ) );
        bool passZpt(zll.pt()>60);


        TString tag_cat;
        int evcat = getDileptonId(abs(id1),abs(id2));
        switch(evcat) {
        case MUMU :
            tag_cat = "mumu";
            break;
        case EE   :
            tag_cat = "ee";
            break;
        case EMU  :
            tag_cat = "emu";
            break;
        default   :
            continue;
        }
        // ID + ISO scale factors
        // Need to implement variations for errors (unused for now)
        if(isMC) {
            // Muon ID SF
            if(abs(id1)==13) weight *= getSFfrom2DHist( fabs(lep1.eta()), lep1.pt(), h_MuonTightWPSF );
            if(abs(id2)==13) weight *= getSFfrom2DHist( fabs(lep2.eta()), lep2.pt(), h_MuonTightWPSF );

            //electron ID SF
            if(abs(id1)==11) weight *= getSFfrom2DHist( fabs(lep1.eta()), lep1.pt(), h_ElectronMediumWPSF );
            if(abs(id2)==11) weight *= getSFfrom2DHist( fabs(lep2.eta()), lep2.pt(), h_ElectronMediumWPSF );
/*
            //electron RECO SF
            if(abs(id1)==11) weight *= getSFfrom2DHist( lep1.eta(), lep1.pt(), h_ElectronRECOSF );
            if(abs(id2)==11) weight *= getSFfrom2DHist( lep2.eta(), lep2.pt(), h_ElectronRECOSF );
*/
        }

        // Trigger Efficiency Scale Factors
        // Pre-ICHEP 2016, the efficiency in MC is 1
        // Thus, the weight is simply the efficiency in data
        // The efficiency is binned in |eta| of tag and probe,
        // for now we just take the leading lepton as tag
        double triggerSF = weight;
        if(isMC) {
            double eta_tag   = lep1.pt() > lep2.pt() ? lep1.eta() : lep2.eta();
            double eta_probe = lep1.pt() > lep2.pt() ? lep2.eta() : lep1.eta();
            if(evcat==MUMU) {
                if (lep1.pt() < 40 && lep2.pt() < 40) 
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_MuonTrigEffSF_0 );
                else if (lep1.pt() < 40 && lep2.pt() >= 40) 
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_MuonTrigEffSF_1 );
                else if (lep1.pt() >= 40 && lep2.pt() < 40) 
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_MuonTrigEffSF_2 );
                else if (lep1.pt() >= 40 && lep2.pt() >= 40) 
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_MuonTrigEffSF_3 );
            } else if(evcat==EE) {
                if (lep1.pt() < 40 && lep2.pt() < 40)
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_ElectronTrigEffSF_0 );
                else if (lep1.pt() < 40 && lep2.pt() >= 40) 
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_ElectronTrigEffSF_1 );
                else if (lep1.pt() >= 40 && lep2.pt() < 40) 
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_ElectronTrigEffSF_2 );
                else if (lep1.pt() >= 40 && lep2.pt() >= 40)
                    weight *= getSFfrom2DHist( fabs(eta_probe), fabs(eta_tag), h_ElectronTrigEffSF_3 );
            }
        }
        triggerSF = weight/triggerSF;

        //split inclusive DY sample into DYToLL and DYToTauTau
        if(isMC && mctruthmode==1) {
            //if(phys.genleptons.size()!=2) continue;
            if(phys.genleptons.size()==2 && isDYToTauTau(phys.genleptons[0].id, phys.genleptons[1].id) ) continue;
        }

        if(isMC && mctruthmode==2) {
            if(phys.genleptons.size()!=2) continue;
            if(!isDYToTauTau(phys.genleptons[0].id, phys.genleptons[1].id) ) continue;
        }


        // Minimum unprescaled double-electron trigger will likely be Ele23_Ele12
        // So offline leading electron >25 to stay out of the strong inefficiency region
        double minLeadingElectronPt = 25;
        if ( evcat==EE && (
              (lep1.pt() > lep2.pt() && lep1.pt() < minLeadingElectronPt)
              || (lep1.pt() < lep2.pt() && lep2.pt() < minLeadingElectronPt)
            ) )
        {
          continue;
        }


        bool hasTrigger(false);

        if(!isMC) {
            // Trigger requirements and duplicate removal (same event could come from different datasets)
            if(evcat==EE) {
              // Seed triggers and their datasets
              if( hasEtrigger && isSingleElePD ) hasTrigger = true;
              if( hasEEtrigger && isDoubleElePD ) hasTrigger = true;
              // Deduplicate: Prefer DoubleEG for double triggers
              if( isSingleElePD && hasEEtrigger ) hasTrigger = false;
            }
            else if(evcat==MUMU) {
              // Seed triggers and their datasets
              if( hasMtrigger && isSingleMuPD ) hasTrigger = true;
              if( hasMMtrigger && isDoubleMuPD ) hasTrigger = true;
              // Deduplicate: Prefer DoubleMu for double triggers
              if( isSingleMuPD && hasMMtrigger ) hasTrigger = false;
            }
            else if(evcat==EMU) {
              // Seed triggers and their datasets
              // We allow emu to be seeded by single lepton triggers as well
              if( hasEtrigger && isSingleElePD ) hasTrigger = true;
              if( hasMtrigger && isSingleMuPD ) hasTrigger = true;
              if( hasEMtrigger && isMuEGPD ) hasTrigger = true;
              // Deduplicate: Prefer MuEG for cross-triggers
              if( (isSingleElePD||isSingleMuPD) && hasEMtrigger ) hasTrigger = false;
              // Deduplicate: Prefer SingleMu for emu events without emu trigger
              if(isSingleElePD && hasEtrigger && hasMtrigger) hasTrigger = false;
            }
        } else {
            // No trigger requirements for 80X ICHEP MC
            // When trigger is simulated, copy logic block above
            // but without the seed dataset or preference conditions
            hasTrigger=true;
        }
        if ( !hasTrigger ) continue;

        tags.push_back(tag_cat); //add ee, mumu, emu category

        // pielup reweightiing
        mon.fillHisto("nvtx_raw",   tags, phys.nvtx,      1.0);
        mon.fillHisto("nvtxwgt_raw",tags, phys.nvtx,      weight);




/*
        if(isMC_WIMP) {
            mon.fillHisto("DMAcceptance",tags,0,1);
            if( (fabs(lep1.eta())>1.6 && fabs(lep1.eta())<2.4) || (fabs(lep2.eta())>1.6 && fabs(lep2.eta())<2.4) ) mon.fillHisto("DMAcceptance",tags,1,1);
            if( (fabs(lep1.eta())>1.6 && fabs(lep1.eta())<2.4) && (fabs(lep2.eta())>1.6 && fabs(lep2.eta())<2.4) ) mon.fillHisto("DMAcceptance",tags,2,1);
        }
*/
        mon.fillHisto("nleptons_raw",tags, nGoodLeptons, weight);


        if(lep1.pt()>lep2.pt()) {
            mon.fillHisto("leadlep_pt_raw",   tags, lep1.pt(), weight);
            mon.fillHisto("leadlep_eta_raw",  tags, lep1.eta(), weight);
            mon.fillHisto("trailep_pt_raw",   tags, lep2.pt(), weight);
            mon.fillHisto("trailep_eta_raw",  tags, lep2.eta(), weight);
        } else {
            mon.fillHisto("leadlep_pt_raw",   tags, lep2.pt(), weight);
            mon.fillHisto("leadlep_eta_raw",  tags, lep2.eta(), weight);
            mon.fillHisto("trailep_pt_raw",   tags, lep1.pt(), weight);
            mon.fillHisto("trailep_eta_raw",  tags, lep1.eta(), weight);
        }


        //
        // 3rd LEPTON ANALYSIS
        //

        //loop over all lepton again, check deltaR with dilepton,
        bool pass3dLeptonVeto(true);
        bool passTauVeto(true);
        bool hasTight3dLepton(false);
        int n3rdLeptons(0), nTight3rdLeptons(0), nTaus(0);
        vector<LorentzVector> allLeptons;
        allLeptons.push_back(lep1);
        allLeptons.push_back(lep2);
        std::vector<LorentzVector> extraTight10Leptons;
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            int lepid = phys.leptons[ilep].id;
            //muon veto
            if(abs(lepid)==13 && fabs(lep.eta()) > 2.4) continue;
            //electront veto
            if(abs(lepid)==11 && fabs(lep.eta()) > 2.5) continue;
            //tau veto
            if(abs(lepid)==15 && fabs(lep.eta())> 2.3) continue;

            bool isMatched(false);
            if( abs(lepid)==15 ) {
                isMatched |= (deltaR(lep1,lep) < 0.4);
                isMatched |= (deltaR(lep2,lep) < 0.4);
            } else {
                isMatched |= (deltaR(lep1,lep) < 0.01);
                isMatched |= (deltaR(lep2,lep) < 0.01);
            }
            if(isMatched) continue;

            bool hasLooseIdandIso(true);
            bool hasTightIdandIso(true);
            if(abs(lepid)==13) { //muon
                hasLooseIdandIso &= ( phys.leptons[ilep].isLooseMu && phys.leptons[ilep].m_pfRelIsoDbeta()<0.25 && phys.leptons[ilep].pt()>10 );
                //
                hasTightIdandIso &= phys.leptons[ilep].pt()>10;
                hasTightIdandIso &= phys.leptons[ilep].isTightMu;
                hasTightIdandIso &= fabs(phys.leptons[ilep].mn_dZ) < 0.1;
                hasTightIdandIso &= fabs(phys.leptons[ilep].mn_d0) < 0.02;
                hasTightIdandIso &= phys.leptons[ilep].m_pfRelIsoDbeta() < 0.15;
            } else if(abs(lepid)==11) { //electron
                hasLooseIdandIso &= ( phys.leptons[ilep].isElpassVeto && phys.leptons[ilep].pt()>10 );
                hasLooseIdandIso |= ( phys.leptons[ilep].isElpassMedium && phys.leptons[ilep].pt()>10 );
                //
                hasTightIdandIso &= ( phys.leptons[ilep].isElpassMedium && phys.leptons[ilep].pt()>10 );
            } else if(abs(lepid)==15) { //tau
                hasLooseIdandIso &= ( phys.leptons[ilep].isTauDM && ( phys.leptons[ilep].ta_combIsoDBeta3Hits < 5 ) && phys.leptons[ilep].pt()>18 );
                //
                hasTightIdandIso &= ( phys.leptons[ilep].isTauDM && phys.leptons[ilep].ta_IsTightIso && phys.leptons[ilep].pt()>20 );
            } else continue;



            if(!hasLooseIdandIso) continue;
            if( abs(lepid)==15 ) nTaus++;
            else n3rdLeptons++;

            allLeptons.push_back(lep);

            if(hasTightIdandIso) {
                extraTight10Leptons.push_back(lep);
                nTight3rdLeptons++;
            }
        }

        pass3dLeptonVeto=(n3rdLeptons==0);
        passTauVeto=(nTaus==0);
        hasTight3dLepton=(nTight3rdLeptons==1);



        //
        //JET AND BTAGING ANALYSIS
        //
        PhysicsObjectJetCollection GoodIdJets;
        bool passBveto(true);
        int nJetsGood30(0);
        int nCSVLtags(0),nCSVMtags(0),nCSVTtags(0);
        double BTagWeights(1.0);
        double dphiJetMET(999);
        for(size_t ijet=0; ijet<corrJets.size(); ijet++) {

            if(corrJets[ijet].pt()<20) continue;
            if(fabs(corrJets[ijet].eta())>5.0) continue;

            //jet ID
            if(!corrJets[ijet].isPFLoose) continue;
            //if(corrJets[ijet].pumva<0.5) continue;

            //check overlaps with selected leptons
            double minDR(999.);
            for(vector<LorentzVector>::iterator lIt = allLeptons.begin(); lIt != allLeptons.end(); lIt++) {
                double dR = deltaR( corrJets[ijet], *lIt );
                if(dR > minDR) continue;
                minDR = dR;
            }
            if(minDR < 0.4) continue;

            GoodIdJets.push_back(corrJets[ijet]);
            if(corrJets[ijet].pt()>30) nJetsGood30++;


            if(corrJets[ijet].pt()>30) {
               dphiJetMET = std::min(fabs(deltaPhi(corrJets[ijet].phi(),metP4.phi())), dphiJetMET );
            }
            //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
            if(corrJets[ijet].pt()>20 && fabs(corrJets[ijet].eta())<2.4)  {

                nCSVLtags += (corrJets[ijet].btag0>CSVLooseWP);
                nCSVMtags += (corrJets[ijet].btag0>CSVMediumWP);
                nCSVTtags += (corrJets[ijet].btag0>CSVTightWP);

                if(!isMC) continue;
                bool isCSVtagged(corrJets[ijet].btag0>CSVMediumWP);


                if(abs(corrJets[ijet].flavid)==5) {
                    double btag_sf = btag_reader.eval( BTagEntry::FLAV_B, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
                    BTagWeights *= myBtagUtils.getNewBTagWeight(isCSVtagged, corrJets[ijet].pt(), btag_sf, "CSVM","CSVM/b_eff");
                } else if(abs(corrJets[ijet].flavid)==4) {
                    double btag_sf = btag_reader.eval( BTagEntry::FLAV_C, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
                    BTagWeights *= myBtagUtils.getNewBTagWeight(isCSVtagged, corrJets[ijet].pt(), btag_sf, "CSVM","CSVM/c_eff");
                } else {
                }


                // begin btagging efficiency
                if(isMC_ttbar||isMC_stop) {
                    int flavid = abs(corrJets[ijet].flavid);
                    for(size_t csvtag=0; csvtag<CSVkey.size(); csvtag++) {
                        bool isBTag(false);
                        if     (CSVkey[csvtag]=="CSVL" && (corrJets[ijet].btag0>CSVLooseWP)) isBTag = true;
                        else if(CSVkey[csvtag]=="CSVM" && (corrJets[ijet].btag0>CSVMediumWP)) isBTag = true;
                        else if(CSVkey[csvtag]=="CSVT" && (corrJets[ijet].btag0>CSVTightWP)) isBTag = true;

                        if(flavid==5) {
                            mon.fillHisto(TString("beff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                            if(isBTag) mon.fillHisto(TString("beff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        } else if(flavid==4) {
                            mon.fillHisto(TString("ceff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                            if(isBTag) mon.fillHisto(TString("ceff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        } else {
                            mon.fillHisto(TString("udsgeff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                            if(isBTag) mon.fillHisto(TString("udsgeff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        }

                    }
                } //end btagging efficiency
            }


        }

        //using CSV Medium WP
        passBveto=(nCSVMtags==0);

        for(size_t ij=0; ij<GoodIdJets.size(); ij++) {
            mon.fillHisto("jet_pt_raw",   tags, GoodIdJets[ij].pt(),weight);
            mon.fillHisto("jet_eta_raw",  tags, GoodIdJets[ij].eta(),weight);
        }

        double dphiZMET=fabs(deltaPhi(zll.phi(),metP4.phi()));
        bool passDphiZMETcut(dphiZMET>2.8);

        TVector2 dil2(zll.px(),zll.py());
        TVector2 met2(metP4.px(),metP4.py());
        double axialmet=dil2*met2;
        axialmet /= -zll.pt();

        //missing ET
        bool passMETcut=(metP4.pt()>100);
        bool passMETcut120=(metP4.pt()>120);


        bool passDphiJetMETCut = dphiJetMET > 0.5;

        //missing ET balance
        double balanceDif = fabs(1-metP4.pt()/zll.pt());
        bool passBalanceCut= balanceDif < 0.4;

        //transverse mass
        double MT_massless = METUtils::transverseMass(zll,metP4,false);

        double response = METUtils::response(zll,metP4);
        bool passResponseCut = (response>-1 && response<1);

        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################


        //event category
        int eventSubCat  = eventCategoryInst.Get(phys,&GoodIdJets);
        TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);

        tags.push_back(tag_cat+tag_subcat); // add jet binning category
        if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back(tag_cat+"lesq1jets");
        if(tag_cat=="mumu" || tag_cat=="ee") {
            tags.push_back("ll"+tag_subcat);
            if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back("lllesq1jets");
        }

        double lep1SF{1.}, lep2SF{1.};
        if(abs(id1)==13) lep1SF = getSFfrom2DHist( fabs(lep1.eta()), lep1.pt(), h_MuonTightWPSF );
        if(abs(id1)==11) lep1SF = getSFfrom2DHist( fabs(lep1.eta()), lep1.pt(), h_ElectronMediumWPSF );
        if(abs(id2)==13) lep2SF = getSFfrom2DHist( fabs(lep2.eta()), lep2.pt(), h_MuonTightWPSF );
        if(abs(id2)==11) lep2SF = getSFfrom2DHist( fabs(lep2.eta()), lep2.pt(), h_ElectronMediumWPSF );
        mon.fillHisto("scalefactors", tags, 0, lep1SF);
        mon.fillHisto("scalefactors", tags, 1, lep2SF);
        mon.fillHisto("scalefactors", tags, 2, triggerSF);
        mon.fillHisto("scalefactors", tags, 3, lep1SF*lep2SF*triggerSF);

        //apply weights
        if(isMC) weight *= BTagWeights;

    // Blinding
//    if( isMC or ( metP4.pt() < 100 ) or ( tag_cat == "emu" )  ) {
        //// Cut Flow synchronization
        int ncut=0;
        mon.fillHisto( "sync_nlep_minus", tags, n3rdLeptons, weight );
        mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
        if( pass3dLeptonVeto ) {
            mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
            if ( passBveto ) {
                mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                if( passTauVeto ) {
                    mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                    mon.fillHisto( "sync_njet_minus", tags, nJetsGood30, weight );
                    if( nJetsGood30 < 2 ) {
                        mon.fillHisto( "sync_zmass_minus", tags, zll.mass(), weight );
                        mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                        if( passZmass ) {
                            mon.fillHisto( "sync_ptz_minus", tags, zll.pt(), weight );
                            mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                            mon.fillHisto( "dPhiTrackMET", tags,  dPhiTrkMET, weight );
                            mon.fillHisto( "dPhiTrackMET_v_pfmet", tags, metP4.pt(),  dPhiTrkMET, weight );
                            if( passZpt ) {
                                mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                                mon.fillHisto( "sync_met_minus",        tags, metP4.pt(), weight );

                                mon.fillHisto( "pfmet2_presel",     tags, metP4.pt(),weight);
                                mon.fillHisto( "dphiZMET_presel",   tags,dphiZMET,weight);
                                mon.fillHisto( "dphiJetMET_presel", tags,dphiJetMET,weight);
                                mon.fillHisto( "balancedif_presel", tags,balanceDif,weight);
                                if(lep1.pt()>lep2.pt()) {
                                    mon.fillHisto("leadlep_pt_met_minus",   tags, lep1.pt(), weight);
                                    mon.fillHisto("leadlep_eta_met_minus",  tags, lep1.eta(), weight);
                                    mon.fillHisto("trailep_pt_met_minus",   tags, lep2.pt(), weight);
                                    mon.fillHisto("trailep_eta_met_minus",  tags, lep2.eta(), weight);
                                } else {
                                    mon.fillHisto("leadlep_pt_met_minus",   tags, lep2.pt(), weight);
                                    mon.fillHisto("leadlep_eta_met_minus",  tags, lep2.eta(), weight);
                                    mon.fillHisto("trailep_pt_met_minus",   tags, lep1.pt(), weight);
                                    mon.fillHisto("trailep_eta_met_minus",  tags, lep1.eta(), weight);
                                }
                                if( passMETcut ) {
                                    mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                                    mon.fillHisto( "sync_dPhiZMET_minus",   tags, dphiZMET, weight );
                                    if( passDphiZMETcut ) {
                                        mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                                        mon.fillHisto( "sync_balance_minus",    tags, balanceDif, weight );
                                        if(passBalanceCut) {
                                            mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                                            mon.fillHisto( "sync_dPhiJetMET_minus", tags, dphiJetMET, weight );
                                            if(passDphiJetMETCut) {
                                                mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                                                mon.fillHisto( "sync_response_minus",   tags, response, weight );
                                                if(passResponseCut) {
                                                    mon.fillHisto( "sync_cutflow",  tags, ncut++, weight);
                                                    mon.fillHisto( "pfmet2_final",tags, metP4.pt(), weight, true);
                                                    if(!isMC) fprintf(outTxtFile_final,"%d | %d | %d | pfmet: %f | mt: %f | mass: %f |jpt: %f | Cat: %s \n",ev.run,ev.lumi
,ev.event,metP4.pt(), MT_massless,zll.mass(),corrJets[0].pt(),(const char*) tag_cat);
        mon.fillHisto("scalefactors_final", tags, 0, lep1SF);
        mon.fillHisto("scalefactors_final", tags, 1, lep2SF);
        mon.fillHisto("scalefactors_final", tags, 2, triggerSF);
        mon.fillHisto("scalefactors_final", tags, 3, lep1SF*lep2SF*triggerSF);

                                                } // passResponseCut
                                            } // passDphiJetMET
                                        } // passBalanceCut
                                    } // passDphiZMET
                                } // passMET
                            } //passZpt
                        } // passZmass
                    } // nJetsGood30
                }// passTauVeto
            }// passBveto
        } // pass3dLeptonVeto


        mon.fillHisto("zpt_raw"                         ,tags, zll.pt(),   weight);
        mon.fillHisto("pfmet_raw"                       ,tags, metP4.pt(), weight);

        if(!isMC) {
            if(MT_massless<300) {
                mon.fillHisto("mt_raw"                          ,tags, MT_massless, weight);
                mon.fillHisto("mt2_raw"                          ,tags, MT_massless, weight,true);
            }
        } else {
            mon.fillHisto("mt_raw"                          ,tags, MT_massless, weight);
            mon.fillHisto("mt2_raw"                          ,tags, MT_massless, weight,true);
        }
        mon.fillHisto("zmass_raw"                       ,tags, zll.mass(), weight);
        mon.fillHisto("njets_raw"                       ,tags, nJetsGood30, weight);
        mon.fillHisto("nbjets_raw"                      ,tags, nCSVLtags, weight);
        if(lep1.pt()>lep2.pt()) mon.fillHisto("ptlep1vs2_raw"                   ,tags, lep1.pt(), lep2.pt(), weight);
        else 			mon.fillHisto("ptlep1vs2_raw"                   ,tags, lep2.pt(), lep1.pt(), weight);


        // WW/ttbar/Wt/tautau control
        mon.fillHisto("zpt_wwctrl_raw"                      ,tags, zll.pt(),   weight);
        mon.fillHisto("zmass_wwctrl_raw"                    ,tags, zll.mass(), weight);
        mon.fillHisto("pfmet_wwctrl_raw"                    ,tags, metP4.pt(), weight);
        mon.fillHisto("mt_wwctrl_raw"              	    ,tags, MT_massless, weight);

        //
        // WZ control region
        //
        if(hasTight3dLepton && passZmass && passBveto) {
            mon.fillHisto("pfmet_WZctrl"                ,tags, metP4.pt(), weight);
            mon.fillHisto("balancedif_WZctrl"           ,tags, balanceDif, weight);
            mon.fillHisto("DphiZMET_WZctrl"             ,tags, dphiZMET, weight);
            mon.fillHisto("zpt_WZctrl"                  ,tags, zll.pt(),   weight);

            if(metP4.pt()>40) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET40"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET40"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>50) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET50"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET50"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>60) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET60"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET60"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>70) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET70"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET70"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>80) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET80"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET80"     ,tags, fakeMT, weight);
            }


        }
    

        //##############################################################################
        //### HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //##############################################################################

        // Clear dmtree before variations loop
        // Since met can shift or be excluded in a variation, make default -1 to make cut clear
        if ( dmReweightTree != nullptr ) {
          dmtree.evcat = evcat;
          dmtree.nJets = nJetsGood30;
          for(size_t i=0; i<20; ++i) {
            dmtree.pfmet[i] = -1.;
            dmtree.weight[i] = 0.;
          }
        }

        //Fill histogram for posterior optimization, or for control regions
        for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
            float iweight = weight;                                               //nominal

            //pileup
            if(varNames[ivar]=="_puup")        iweight *=TotalWeight_plus;        //pu up
            if(varNames[ivar]=="_pudown")      iweight *=TotalWeight_minus;       //pu down

            //PDF variation: Acceptance+Normalization
            if( (isSignal || isMCBkg_runPDFQCDscale) && (varNames[ivar]=="_pdfup" || varNames[ivar]=="_pdfdown") ) {

                if(isMC_Unpart || isMC_ADD) {
                    //for Pythia8 based signal samples, PDF weights do not exist in miniAOD
                    if(mPDFInfo) {
                        float PDFWeight_plus(1.0), PDFWeight_down(1.0);
                        std::vector<float> wgts=mPDFInfo->getWeights(iev);
                        for(size_t ipw=0; ipw<wgts.size(); ipw++) {
                            PDFWeight_plus = TMath::Max(PDFWeight_plus,wgts[ipw]);
                            PDFWeight_down = TMath::Min(PDFWeight_down,wgts[ipw]);
                        }
                        if(varNames[ivar]=="_pdfup")    iweight *= PDFWeight_plus;
                        else if(varNames[ivar]=="_pdfdown")  iweight *= PDFWeight_down;
                    }
                } else {
                    // for POWHEG, MADGRAPH based samples
                    // retrive PDFweights from ntuple
                    TH1F *pdf_h = new TH1F();
                    for(int npdf=0; npdf<ev.npdfs; npdf++) {
                        pdf_h->Fill(ev.pdfWeights[npdf]);
                    }
                    double pdfError = pdf_h->GetRMS();
                    delete pdf_h;
                    //cout << "pdfError: " << pdfError << endl;

                    // retrive alphaSweights from ntuple
                    TH1F *alphaS_h = new TH1F();
                    for(int nalphaS=0; nalphaS<ev.nalphaS; nalphaS++) {
                        alphaS_h->Fill(ev.alphaSWeights[nalphaS]);
                    }
                    double alphaSError = alphaS_h->GetRMS();
                    delete alphaS_h;
                    //cout << "alphaSError: " << alphaSError << endl;
                    double PDFalphaSWeight = sqrt(pdfError*pdfError + alphaSError*alphaSError);

                    if(varNames[ivar]=="_pdfup")    iweight *= (1.+PDFalphaSWeight);
                    else if(varNames[ivar]=="_pdfdown")  iweight *= (1.-PDFalphaSWeight);
                }
            }

            //QCDscale: Acceptance+Normalization
            if( (isMCBkg_runPDFQCDscale || isSignal) && (varNames[ivar]=="_qcdscaleup" || varNames[ivar]=="_qcdscaledown") ) {
                float QCDscaleWeight_plus(1.0), QCDscaleWeight_down(1.0);

                std::vector<float> QCDscaleWgts;
                QCDscaleWgts.clear();
                if(isMC_Unpart || isMC_ADD) {
                    QCDscaleWgts.push_back(weight_QCDscale_muR1_muF2*genWeight);
                    QCDscaleWgts.push_back(weight_QCDscale_muR1_muF0p5*genWeight);
                    QCDscaleWgts.push_back(weight_QCDscale_muR2_muF1*genWeight);
                    QCDscaleWgts.push_back(weight_QCDscale_muR2_muF2*genWeight);
                    QCDscaleWgts.push_back(weight_QCDscale_muR0p5_muF1*genWeight);
                    QCDscaleWgts.push_back(weight_QCDscale_muR0p5_muF0p5*genWeight);
                } else {
                    QCDscaleWgts.push_back(ev.weight_QCDscale_muR1_muF2);
                    QCDscaleWgts.push_back(ev.weight_QCDscale_muR1_muF0p5);
                    QCDscaleWgts.push_back(ev.weight_QCDscale_muR2_muF1);
                    QCDscaleWgts.push_back(ev.weight_QCDscale_muR2_muF2);
                    QCDscaleWgts.push_back(ev.weight_QCDscale_muR0p5_muF1);
                    QCDscaleWgts.push_back(ev.weight_QCDscale_muR0p5_muF0p5);
                }
                for(size_t ipw=0; ipw<QCDscaleWgts.size(); ipw++) {
                    if(ipw==0) {
                        // need to add genWeight to remove the negative sign from negative event
                        // the negative sign is already assigned when using iweight = weight
                        QCDscaleWeight_plus = QCDscaleWgts[ipw]/genWeight;
                        QCDscaleWeight_down = QCDscaleWgts[ipw]/genWeight;
                    } else {
                        QCDscaleWeight_plus = TMath::Max(QCDscaleWeight_plus,QCDscaleWgts[ipw]/genWeight);
                        QCDscaleWeight_down = TMath::Min(QCDscaleWeight_down,QCDscaleWgts[ipw]/genWeight);
                    }
                    //cout << "QCDscaleWgts[" << ipw << "]: " << QCDscaleWgts[ipw] << endl;
                }
                if(varNames[ivar]=="_qcdscaleup")    	   iweight *= QCDscaleWeight_plus;
                else if(varNames[ivar]=="_qcdscaledown")   iweight *= QCDscaleWeight_down;
            }

            if( isMC_ZZ2L2Nu && (varNames[ivar]=="_qqZZewkup" || varNames[ivar]=="_qqZZewkdown") ) {
                if(varNames[ivar]=="_qqZZewkup")        iweight *= weight_ewkup;
                if(varNames[ivar]=="_qqZZewkdown")      iweight *= weight_ewkdown;
            }




            //##############################################
            // recompute MET/MT if JES/JER was varied
            //##############################################
            //LorentzVector vMET = variedMET[ivar>8 ? 0 : ivar];
            //PhysicsObjectJetCollection &vJets = ( ivar<=4 ? variedJets[ivar] : variedJets[0] );


            LorentzVector vMET = variedMET[0];
            if(varNames[ivar]=="_jerup" || varNames[ivar]=="_jerdown" || varNames[ivar]=="_jesup" || varNames[ivar]=="_jesdown" ||
                    varNames[ivar]=="_umetup" || varNames[ivar]=="_umetdown" || varNames[ivar]=="_lesup" || varNames[ivar]=="_lesdown") {
                vMET = variedMET[ivar];
            }

            PhysicsObjectJetCollection &vJets = variedJets[0];
            if(varNames[ivar]=="_jerup" || varNames[ivar]=="_jerdown" || varNames[ivar]=="_jesup" || varNames[ivar]=="_jesdown") {
                vJets = variedJets[ivar];
            }

            double LocalDphiJetMET(999.);
            bool passLocalBveto(true);
            for(size_t ijet=0; ijet<vJets.size(); ijet++) {

                if(vJets[ijet].pt()<20) continue;
                if(fabs(vJets[ijet].eta())>5.0) continue;

                //jet ID
                if(!vJets[ijet].isPFLoose) continue;
                //if(vJets[ijet].pumva<0.5) continue;

                //check overlaps with selected leptons
                double minDR(999.);
                for(vector<LorentzVector>::iterator lIt = allLeptons.begin(); lIt != allLeptons.end(); lIt++) {
                    double dR = deltaR( vJets[ijet], *lIt );
                    if(dR > minDR) continue;
                    minDR = dR;
                }
                if(minDR < 0.4) continue;

                // Separation between jet and MET
                if(vJets[ijet].pt() > 30 ) {
                    LocalDphiJetMET = std::min( LocalDphiJetMET, fabs(deltaPhi(vJets[ijet].phi(),vMET.phi())));
                }

                if(vJets[ijet].pt()>20 && fabs(vJets[ijet].eta())<2.4) {

                    passLocalBveto &= (vJets[ijet].btag0<CSVMediumWP);
                    bool isLocalCSVtagged(vJets[ijet].btag0>CSVMediumWP);

                    double BTagWeights_Up=1., BTagWeights_Down=1.;
                    if(abs(vJets[ijet].flavid)==5) {
                        double BTagSF_Up   = btag_reader_up.eval( BTagEntry::FLAV_B, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                        double BTagSF_Down = btag_reader_down.eval( BTagEntry::FLAV_B, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                        BTagWeights_Up   *= myBtagUtils.getNewBTagWeight(isLocalCSVtagged,vJets[ijet].pt(),BTagSF_Up, "CSVM","CSVM/b_eff");
                        BTagWeights_Down *= myBtagUtils.getNewBTagWeight(isLocalCSVtagged,vJets[ijet].pt(),BTagSF_Down,"CSVM","CSVM/b_eff");

                    } else if(abs(vJets[ijet].flavid)==4) {
                        double BTagSF_Up   = btag_reader_up.eval( BTagEntry::FLAV_C, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                        double BTagSF_Down = btag_reader_down.eval( BTagEntry::FLAV_C, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                        BTagWeights_Up 	   *= myBtagUtils.getNewBTagWeight( isLocalCSVtagged,vJets[ijet].pt(),BTagSF_Up,"CSVM","CSVM/c_eff" );
                        BTagWeights_Down   *= myBtagUtils.getNewBTagWeight( isLocalCSVtagged,vJets[ijet].pt(),BTagSF_Down,"CSVM","CSVM/c_eff" );

                    } else {
                    }
                    if(varNames[ivar]=="_btagup") 	iweight *= BTagWeights_Up;
                    if(varNames[ivar]=="_btagdown")	iweight *= BTagWeights_Down;
                }

            }


            double responseVal = METUtils::response(zll,vMET);
            bool passResponseValCut = (fabs(responseVal)<1);

            bool passLocalDphiJetMETCut = LocalDphiJetMET > 0.5;
            bool passBaseSelection( passZmass && passZpt && pass3dLeptonVeto && passTauVeto && passLocalBveto && passResponseValCut && passLocalDphiJetMETCut );

            double mt_massless = METUtils::transverseMass(zll,vMET,false); //massless mt
            double LocalDphiZMET=fabs(deltaPhi(zll.phi(),vMET.phi()));


            //############
            //optimization
            //############
            for(unsigned int index=0; index<nOptims; index++) {
                double minMET = optim_Cuts1_MET[index];
                double minBalance = optim_Cuts1_Balance[index];
                double minDphi = optim_Cuts1_DphiZMET[index];

                bool passLocalMETcut(vMET.pt()>minMET);
                bool passLocalBalanceCut=(vMET.pt()/zll.pt()>(1.-minBalance) && vMET.pt()/zll.pt()<(1.+minBalance));
                bool passLocalDphiZMETcut(LocalDphiZMET>minDphi);

                bool passOptimSelection(passBaseSelection && passLocalMETcut && passLocalBalanceCut && passLocalDphiZMETcut);


                //for extrapolation of DY process based on MET
                if(ivar==0 && passBaseSelection && passLocalBalanceCut && passLocalDphiZMETcut ) {
                    mon.fillHisto("pfmet_minus_shapes",tags,index, vMET.pt(), iweight);
                }

                // fill shapes for limit setting
                if( passOptimSelection ) {
                    mon.fillHisto(TString("mt_shapes")+varNames[ivar],tags,index, mt_massless, iweight);
                    mon.fillHisto(TString("mt2_shapes")+varNames[ivar],tags,index, mt_massless, iweight);

                    mon.fillHisto(TString("met_shapes")+varNames[ivar],tags,index, vMET.pt(), iweight);
                    mon.fillHisto(TString("met2_shapes")+varNames[ivar],tags,index, vMET.pt(), iweight);
                    bool ismmlesq1 = false;
                    for(auto& tag : tags) if (tag=="mumulesq1jets") ismmlesq1 = true;
                    if ( ivar == 0 && index == 0 && ismmlesq1 && dmtree.nJets >1 ) {
                      std::cout << "Filling met2shapes " << vMET.pt() << " wgt " << iweight << std::endl;
                      std::cout << "nJets = " << dmtree.nJets <<std::endl;
                    }

                    if (index==0 && dmReweightTree != nullptr ) {
                        dmtree.pfmet[ivar] = vMET.pt();
                        dmtree.weight[ivar] = iweight;
                    }
                }


            }//all optimization END
        }//Systematic variation END
        
        if ( dmReweightTree != nullptr ) dmReweightTree->Fill();
    } // loop on all events END




    printf("\n");
    file->Close();

    //##############################################
    //########     SAVING HISTO TO FILE     ########
    //##############################################

    ofile->cd();
    mon.Write();
    if ( dmReweightTree != nullptr ) dmReweightTree->Write();
    ofile->Close();

    PU_Central_File->Close();
    PU_Up_File->Close();
    PU_Down_File->Close();

    MuonTrigEffSF_File->Close();
    //ElectronTrigEffSF_File->Close();
    ElectronMediumWPSF_File->Close();
    ElectronRECOSF_File->Close();

    UnpartQCDscaleWeights_File->Close();
    ADDQCDscaleWeights_File->Close();

    if(outTxtFile_final)fclose(outTxtFile_final);
}
