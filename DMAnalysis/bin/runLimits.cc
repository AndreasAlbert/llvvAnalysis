#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include <sstream>
#include "llvvAnalysis/DMAnalysis/src/tdrstyle.C"
#include "llvvAnalysis/DMAnalysis/src/JSONWrapper.cc"
#include "llvvAnalysis/DMAnalysis/interface/setStyle.h"
#include "llvvAnalysis/DMAnalysis/interface/MacroUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/plotter.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>

using namespace std;

double NonResonnantSyst = 0.25;

double WWtopSyst_ee0jet = 0.;
double WWtopSyst_ee1jet = 0.;
double WWtopSyst_eelesq1jet = 0.;
double WWtopSyst_mm0jet = 0.;
double WWtopSyst_mm1jet = 0.;
double WWtopSyst_mmlesq1jet = 0.;
double WWtopSyst_ll0jet = 0.;
double WWtopSyst_ll1jet = 0.;
double WWtopSyst_lllesq1jet = 0.;

double WjetsSyst_ee = 0.;
double WjetsSyst_mm = 0.;


struct YIELDS_T {
    double WZ;
    double ZZ;
    double Zjets;
    double WWtop;
    double Wjets;
    double Data;
    double Sig;
    double totBkg;

    double WZ_StatErr;
    double ZZ_StatErr;
    double Zjets_StatErr;
    double WWtop_StatErr;
    double Wjets_StatErr;
    double Data_StatErr;
    double Sig_StatErr;
    double totBkg_StatErr;

    YIELDS_T():WZ(0.),ZZ(0.),Zjets(0.),WWtop(0.),Wjets(0.),Data(0.),Sig(0.),totBkg(0.),
        WZ_StatErr(0.),ZZ_StatErr(0.),Zjets_StatErr(0.),WWtop_StatErr(0.),
        Wjets_StatErr(0.),Data_StatErr(0.),Sig_StatErr(0.),totBkg_StatErr(0.) { }
};





// Map for non-shape systematics (lnN)
std::map<TString, double> normSysts;
void initNormalizationSysts();

//wrapper for a projected shape for a given set of cuts
class Shape_t {
public:
    TH1* data, *totalBckg;
    std::vector<TH1 *> bckg, signal, bckgInSignal;
    //the key corresponds to the proc name
    //the key is the name of the variation: e.g. jesup, jesdown, etc.
    std::map<TString,std::vector<std::pair<TString, TH1*> > > bckgVars, signalVars, bckgInSignalVars;

    std::map<TString, double> xsections;
    std::map<TString, double> BRs;

    Shape_t() {}
    ~Shape_t() {}


    void clear() {
        std::cout<<"shape is destructed...";
        if(data)delete data;
        if(totalBckg)totalBckg->Delete();
        for(unsigned int i=0; i<bckg.  size(); i++) {
            delete bckg  [i];
        }
        bckg  .clear();
        for(unsigned int i=0; i<signal.size(); i++) {
            delete signal[i];
        }
        signal.clear();
        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it=bckgVars  .begin(); it!=bckgVars  .end(); it++) {
            for(unsigned int i=0; i<(*it).second.size(); i++) {
                delete (*it).second[i].second;
            }
        }
        bckgVars  .clear();
        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it=signalVars.begin(); it!=signalVars.end(); it++) {
            for(unsigned int i=0; i<(*it).second.size(); i++) {
                delete (*it).second[i].second;
            }
        }
        signalVars.clear();
        std::cout<<"done\n";

    }


};

typedef std::pair<TString,TString> RateKey_t;
struct DataCardInputs {
    TString shapesFile;
    std::vector<TString> ch;
    std::vector<TString> procs;
    std::map<RateKey_t, Double_t> obs;
    std::map<RateKey_t, Double_t> rates;
    std::map<TString, std::map<RateKey_t,Double_t> > systs;
    int nsignalproc;
};


void printHelp();
Shape_t getShapeFromFile(TFile* inF, TString ch, TString shapeName, int cutBin,JSONWrapper::Object &Root,double minCut=0, double maxCut=9999, bool onlyData=false);

void getYieldsFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName, bool isdataBlinded);




void convertHistosForLimits_core(DataCardInputs& dci, TString& proc, TString& bin, TString& ch, std::vector<TString>& systs, std::vector<TH1*>& hshapes);
DataCardInputs convertHistosForLimits(Int_t mass,TString histo="finalmt",TString url="plotter.root",TString Json="");
std::vector<TString> buildDataCard(Int_t mass, TString histo="finalmt", TString url="plotter.root",TString Json="");
void doBackgroundSubtraction(std::vector<TString>& signal_channels,TString control_channel,map<TString, Shape_t> &allShapes, TString mainHisto, TString sideBandHisto, TString url, JSONWrapper::Object &Root, bool isMCclosureTest);
void dodataDrivenWWtW(std::vector<TString>& signal_channels,TString control_channel,map<TString, Shape_t> &allShapes, TString mainHisto,bool isMCclosureTest);
void doWjetsBackground(std::vector<TString>& signal_channels,map<TString, Shape_t> &allShapes, TString mainHisto);
void doDYextrapolation(std::vector<TString>& signal_channels,map<TString, Shape_t> &allShapes,TString mainHisto,TString DY_EXTRAPOL_SHAPES,float DY_EXTRAPOL,bool isleftCR);
void doQCDBackground(std::vector<TString>& signal_channels,map<TString, Shape_t> &allShapes, TString mainHisto);
void doWZSubtraction(std::vector<TString>& signal_channels,TString control_channel,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto);
void BlindData(std::vector<TString>& signal_channels, map<TString, Shape_t>& allShapes, TString mainHisto, bool addSignal);



bool subNRB2011 = false;
bool subNRB2012 = false;
bool MCclosureTest = false;

bool mergeWWandZZ = false;
bool skipWW = true;
std::vector<TString> Channels;
std::vector<TString> AnalysisBins;
bool fast = false;
bool subDY = false;
bool subWZ = false;
double DDRescale = 1.0;
double MCRescale = 1.0;
bool blindData = false;
bool blindWithSignal = false;
TString DYFile ="";
TString inFileUrl(""),jsonFile(""), histo("");
TString postfix="";
TString systpostfix="";
double shapeMin = 0;
double shapeMax = 9999;
double shapeMinVBF = 0;
double shapeMaxVBF = 9999;

float DYMET_EXTRAPOL = 9999;
float DYRESP_EXTRAPOL = 9999;
float DYDPHI_EXTRAPOL = 9999;

int indexvbf = -1;
int indexcut   = -1, indexcutL=-1, indexcutR=-1;
int mass=-1, massL=-1, massR=-1;
int MV=-1, MA=-1;
float K1=-1, K2=-1;
bool runSystematics = false;
bool shape = false;
float sysSherpa=1.;

TString signal="";
void initNormalizationSysts()
{
    normSysts["lumi_13TeV"] = 0.062;
    normSysts["CMS_eff_e"] = 0.03;
    normSysts["CMS_eff_m"] = 0.04;

    //unparticle
    normSysts["QCDscale_UnPart1p01"]=1.027561608;
    normSysts["QCDscale_UnPart1p02"]=1.027560594;
    normSysts["QCDscale_UnPart1p04"]=1.027107579;
    normSysts["QCDscale_UnPart1p06"]=1.026726974;
    normSysts["QCDscale_UnPart1p09"]=1.026031746;
    normSysts["QCDscale_UnPart1p10"]=1.025737817;
    normSysts["QCDscale_UnPart1p20"]=1.023378111;
    normSysts["QCDscale_UnPart1p30"]=1.019920319;
    normSysts["QCDscale_UnPart1p40"]=1.016645015;
    normSysts["QCDscale_UnPart1p50"]=1.012751678;
    normSysts["QCDscale_UnPart1p60"]=1.009233738;
    normSysts["QCDscale_UnPart1p70"]=1.004749134;
    normSysts["QCDscale_UnPart1p80"]=1.001779359;
    normSysts["QCDscale_UnPart1p90"]=1.006003132;
    normSysts["QCDscale_UnPart2p00"]=1.008754209;
    normSysts["QCDscale_UnPart2p20"]=1.017029328;

}

void printHelp()
{
    printf("Options\n");
    printf("--in        --> input file with from plotter\n");
    printf("--json      --> json file with the sample descriptor\n");
    printf("--histo     --> name of histogram to be used\n");
    printf("--shapeMin  --> left cut to apply on the shape histogram\n");
    printf("--shapeMax  --> right cut to apply on the shape histogram\n");
    printf("--shapeMinVBF  --> left cut to apply on the shape histogram for Vbf bin\n");
    printf("--shapeMaxVBF  --> right cut to apply on the shape histogram for Vbf bin\n");
    printf("--indexvbf  --> index of selection to be used for the vbf bin (if unspecified same as --index)\n");
    printf("--index     --> index of selection to be used (Xbin in histogram to be used)\n");
    printf("--indexL    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
    printf("--indexR    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
    printf("--syst      --> use this flag if you want to run systematics, default is no systematics\n");
    printf("--shape     --> use this flag if you want to run shapeBased analysis, default is cut&count\n");
    printf("--subNRB    --> use this flag if you want to subtract non-resonant-backgounds similarly to what was done in 2011 (will also remove H->WW)\n");
    printf("--subNRB12  --> use this flag if you want to subtract non-resonant-backgounds using a new technique that keep H->WW\n");
    printf("--subDY     --> histogram that contains the Z+Jets background estimated from Gamma+Jets)\n");
    printf("--subWZ     --> use this flag if you want to subtract WZ background by the 3rd lepton SB)\n");
    printf("--DDRescale --> factor to be used in order to multiply/rescale datadriven estimations\n");
    printf("--closure   --> use this flag if you want to perform a MC closure test (use only MC simulation)\n");
    printf("--bins      --> list of bins to be used (they must be comma separated without space)\n");
    printf("--HWW       --> use this flag to consider HWW signal)\n");
    printf("--skipGGH   --> use this flag to skip GGH signal)\n");
    printf("--skipQQH   --> use this flag to skip GGH signal)\n");
    printf("--blind     --> use this flag to replace observed data by total predicted background)\n");
    printf("--blindWithSignal --> use this flag to replace observed data by total predicted background+signal)\n");
    printf("--fast      --> use this flag to only do assymptotic prediction (very fast but inaccurate))\n");
    printf("--postfix    --> use this to specify a postfix that will be added to the process names)\n");
    printf("--systpostfix    --> use this to specify a syst postfix that will be added to the process names)\n");
    printf("--MCRescale    --> use this to specify a syst postfix that will be added to the process names)\n");
    printf("--interf     --> use this to rescale xsection according to WW interferences)\n");
}

//
int main(int argc, char* argv[])
{
    setTDRStyle();
    gStyle->SetPadTopMargin   (0.06);
    gStyle->SetPadBottomMargin(0.12);
    //gStyle->SetPadRightMargin (0.16);
    gStyle->SetPadRightMargin (0.06);
    gStyle->SetPadLeftMargin  (0.14);
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.45);
    gStyle->SetPalette(1);
    gStyle->SetNdivisions(505);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);


    //get input arguments
    for(int i=1; i<argc; i++) {
        string arg(argv[i]);
        if(arg.find("--help")          !=string::npos) {
            printHelp();
            return -1;
        } else if(arg.find("--subNRB12") !=string::npos) {
            subNRB2012=true;
            skipWW=false;
            printf("subNRB2012 = True\n");
        } else if(arg.find("--subNRB")   !=string::npos) {
            subNRB2011=true;
            skipWW=true;
            printf("subNRB2011 = True\n");
        } else if(arg.find("--subDY")    !=string::npos) {
            subDY=true;
            DYFile=argv[i+1];
            i++;
            printf("Z+Jets will be replaced by %s\n",DYFile.Data());
        } else if(arg.find("--subWZ")    !=string::npos) {
            subWZ=true;
            printf("WZ will be estimated from 3rd lepton SB\n");
        } else if(arg.find("--DDRescale")!=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&DDRescale);
            i++;
        } else if(arg.find("--MCRescale")!=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&MCRescale);
            i++;
        } else if(arg.find("--metEXTRAPOL") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%f",&DYMET_EXTRAPOL);
            i++;
            printf("DYMET_EXTRAPOL = %f\n", DYMET_EXTRAPOL);
        } else if(arg.find("--respEXTRAPOL") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%f",&DYRESP_EXTRAPOL);
            i++;
            printf("DYRESP_EXTRAPOL = %f\n", DYRESP_EXTRAPOL);
        } else if(arg.find("--dphiEXTRAPOL") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%f",&DYDPHI_EXTRAPOL);
            i++;
            printf("DYDPHI_EXTRAPOL = %f\n", DYDPHI_EXTRAPOL);
        } else if(arg.find("--HWW")      !=string::npos) {
            skipWW=false;
            printf("HWW = True\n");
        } else if(arg.find("--blindWithSignal")  !=string::npos) {
            blindData=true;
            blindWithSignal=true;
            printf("blindData = True; blindWithSignal = True\n");
        } else if(arg.find("--blind")    !=string::npos) {
            blindData=true;
            printf("blindData = True\n");
        } else if(arg.find("--closure")  !=string::npos) {
            MCclosureTest=true;
            printf("MCclosureTest = True\n");
        } else if(arg.find("--shapeMinVBF") !=string::npos && i+1<argc) {
            sscanf(argv[i+1],"%lf",&shapeMinVBF);
            i++;
            printf("Min cut on shape for VBF = %f\n", shapeMinVBF);
        } else if(arg.find("--shapeMaxVBF") !=string::npos && i+1<argc) {
            sscanf(argv[i+1],"%lf",&shapeMaxVBF);
            i++;
            printf("Max cut on shape for VBF = %f\n", shapeMaxVBF);
        } else if(arg.find("--shapeMin") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMin);
            i++;
            printf("Min cut on shape = %f\n", shapeMin);
        } else if(arg.find("--shapeMax") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMax);
            i++;
            printf("Max cut on shape = %f\n", shapeMax);
        } else if(arg.find("--indexvbf") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexvbf);
            i++;
            printf("indexVBF = %i\n", indexvbf);
        } else if(arg.find("--indexL")   !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcutL);
            i++;
            printf("indexL = %i\n", indexcutL);
        } else if(arg.find("--indexR")   !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcutR);
            i++;
            printf("indexR = %i\n", indexcutR);
        } else if(arg.find("--index")    !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcut);
            i++;
            printf("index = %i\n", indexcut);
        } else if(arg.find("--in")       !=string::npos && i+1<argc)  {
            inFileUrl = argv[i+1];
            i++;
            printf("in = %s\n", inFileUrl.Data());
        } else if(arg.find("--json")     !=string::npos && i+1<argc)  {
            jsonFile  = argv[i+1];
            i++;
            printf("json = %s\n", jsonFile.Data());
        } else if(arg.find("--histo")    !=string::npos && i+1<argc)  {
            histo     = argv[i+1];
            i++;
            printf("histo = %s\n", histo.Data());
        } else if(arg.find("--bins")     !=string::npos && i+1<argc)  {
            char* pch = strtok(argv[i+1],",");
            printf("bins are : ");
            while (pch!=NULL) {
                printf(" %s ",pch);
                AnalysisBins.push_back(pch);
                pch = strtok(NULL,",");
            }
            printf("\n");
            i++;
        } else if(arg.find("--channels") !=string::npos && i+1<argc)  {
            char* pch = strtok(argv[i+1],",");
            printf("channels are : ");
            while (pch!=NULL) {
                printf(" %s ",pch);
                Channels.push_back(pch);
                pch = strtok(NULL,",");
            }
            printf("\n");
            i++;
        } else if(arg.find("--fast")     !=string::npos) {
            fast=true;
            printf("fast = True\n");
        } else if(arg.find("--postfix")  !=string::npos && i+1<argc)  {
            postfix = argv[i+1];
            systpostfix = argv[i+1];
            i++;
            printf("postfix '%s' will be used\n", postfix.Data());
        } else if(arg.find("--systpostfix") !=string::npos && i+1<argc)  {
            systpostfix = argv[i+1];
            i++;
            printf("systpostfix '%s' will be used\n", systpostfix.Data());
        } else if(arg.find("--syst")     !=string::npos) {
            runSystematics=true;
            printf("syst = True\n");
        } else if(arg.find("--shape")    !=string::npos) {
            shape=true;
            printf("shapeBased = True\n");
        } else if(arg.find("--aTGC_Syst")!=string::npos) {
            sscanf(argv[i+1],"%f",&sysSherpa);
            i++;
            printf("Additional systematic on the shape setted: %f\n", sysSherpa);
        } else if(arg.find("--signal")!=string::npos) {
            signal = argv[i+1];
            i++;
            printf("Chosen signal: %s\n", signal.Data());
        }
    }

    if(jsonFile.IsNull() || inFileUrl.IsNull() || histo.IsNull() || indexcut == -1 || signal.IsNull()) {
        printHelp();
        return -1;
    }
    if(AnalysisBins.size()==0)AnalysisBins.push_back("");
    if(Channels.size()==0) {
        Channels.push_back("ee");
        Channels.push_back("mumu");
        Channels.push_back("ll"); //RENJIE, add all
    }

    initNormalizationSysts();

    //build the datacard for this mass point
    std::vector<TString> dcUrls = buildDataCard(mass,histo,inFileUrl, jsonFile);

}




//
Shape_t getShapeFromFile(TFile* inF, TString ch, TString shapeName, int cutBin, JSONWrapper::Object &Root, double minCut, double maxCut, bool onlyData)
{
    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

    Shape_t shape;
    shape.totalBckg=NULL;
    shape.data=NULL;

    std::vector<TString> BackgroundsInSignal;

    //iterate over the processes required
    std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
    for(unsigned int i=0; i<Process.size(); i++) {
        TString procCtr("");
        procCtr+=i;
        TString proc=(Process[i])["tag"].toString();

        bool isData(Process[i]["isdata"].toBool());
        if(onlyData && !isData)continue; //just here to speedup the NRB prediction

        bool isSignal(Process[i].isTag("issignal") && Process[i]["issignal"].toBool());

        // Only keep the one signal explicitly specified on the command line
        if( isSignal and not proc.Contains(signal) ) continue;
        if(Process[i]["spimpose"].toBool() && (proc.Contains("ggH") || proc.Contains("qqH")))isSignal=true;
        int color(1);
        if(Process[i].isTag("color" ) ) color  = (int)Process[i]["color" ].toInt();
        int lcolor(color);
        if(Process[i].isTag("lcolor") ) lcolor = (int)Process[i]["lcolor"].toInt();
        int mcolor(color);
        if(Process[i].isTag("mcolor") ) mcolor = (int)Process[i]["mcolor"].toInt();
        int lwidth(1);
        if(Process[i].isTag("lwidth") ) lwidth = (int)Process[i]["lwidth"].toInt();
        int lstyle(1);
        if(Process[i].isTag("lstyle") ) lstyle = (int)Process[i]["lstyle"].toInt();
        int fill(1001);
        if(Process[i].isTag("fill"  ) ) fill   = (int)Process[i]["fill"  ].toInt();
        int marker(20);
        if(Process[i].isTag("marker") ) marker = (int)Process[i]["marker"].toInt();

        TDirectory *pdir = (TDirectory *)inF->Get(proc);
        if(pdir==0) {
            printf("Skip Proc=%s because its directory is missing in root file\n", proc.Data());
            continue;
        }
        TH1* syst = (TH1*)pdir->Get("optim_systs");
        if(syst==NULL) {
            cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> NULL" << endl;
            syst =new TH1F("optim_systs","optim_systs",1,0,1);
            syst->GetXaxis()->SetBinLabel(1,"");
        }
        for(int ivar = 1; ivar<=syst->GetNbinsX(); ivar++) {
            TH1D* hshape   = NULL;

            TString varName = syst->GetXaxis()->GetBinLabel(ivar);
            TString histoName = ch+"_"+shapeName+varName ;
            TH2* hshape2D = (TH2*)pdir->Get(histoName );
            if(!hshape2D) {
                hshape2D = (TH2*)pdir->Get(shapeName+varName);
                if(hshape2D)hshape2D->Reset();
            }

            if(hshape2D) {
                histoName.ReplaceAll(ch,ch+"_proj"+procCtr);
                hshape   = hshape2D->ProjectionY(histoName,cutBin,cutBin);
                if(hshape->Integral()<=0 && varName=="" && !isData) {
                    hshape->Reset();
                    hshape->SetBinContent(1, 1E-10);
                }

                if(isnan((float)hshape->Integral())) {
                    hshape->Reset();
                }
                hshape->SetDirectory(0);
                hshape->SetTitle(proc);
                fixExtremities(hshape,true,true);
                hshape->SetFillColor(color);
                hshape->SetLineColor(lcolor);
                hshape->SetMarkerColor(mcolor);
                hshape->SetFillStyle(fill);
                hshape->SetLineWidth(lwidth);
                hshape->SetMarkerStyle(marker);
                hshape->SetLineStyle(lstyle);
            } else {
                if(!histoName.Contains("FR_WjetCtrl") && !histoName.Contains("FR_QCDCtrl")
                        && !histoName.Contains("minus_shapes") && !histoName.Contains("_NRBctrl") && !histoName.Contains("_NRBsyst"))
                    printf("Histo %s does not exist for syst:%s\n", histoName.Data(), varName.Data());
                continue;
            }


            //if current shape is the one to cut on, then apply the cuts
            if(shapeName == histo) {
                for(int x=0; x<=hshape->GetXaxis()->GetNbins()+1; x++) {
                    if(hshape->GetXaxis()->GetBinCenter(x)<=minCut || hshape->GetXaxis()->GetBinCenter(x)>=maxCut) {
                        hshape->SetBinContent(x,0);
                        hshape->SetBinError(x,0);
                    }
                }
                hshape->GetYaxis()->SetTitle("Entries (/25GeV)");
            }

            hshape->Scale(MCRescale);

            //save in structure
            if(isData) {
                if(varName=="")  shape.data=hshape;
                else continue;
            } else if(isSignal) {

                if(varName=="") {
                    if(Process[i]["data"].daughters()[0].isTag("xsec"))shape.xsections[proc] = Process[i]["data"].daughters()[0]["xsec"].toDouble();
                    if(Process[i]["data"].daughters()[0].isTag("br")) {
                        std::vector<JSONWrapper::Object> BRs = Process[i]["data"].daughters()[0]["br"].daughters();
                        double totalBR=1.0;
                        for(size_t ipbr=0; ipbr<BRs.size(); ipbr++) {
                            totalBR*=BRs[ipbr].toDouble();
                        }
                        shape.BRs[proc] = totalBR;
                    }

                    int procIndex = -1;
                    for(unsigned int i=0; i<shape.signal.size(); i++) {
                        if(string(proc.Data())==shape.signal[i]->GetTitle() ) {
                            procIndex=i;
                            break;
                        }
                    }
                    if(procIndex>=0) shape.signal[procIndex]->Add(hshape);
                    else             {
                        hshape->SetTitle(proc);
                        shape.signal.push_back(hshape);
                    }

                    //printf("Adding signal %s\n",proc.Data());
                } else {
                    std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it = shape.signalVars.find(proc);

                    bool newVar = true;
                    if(it!=shape.signalVars.end()) {
                        for(unsigned int i=0; i<it->second.size(); i++) {
                            if( string(it->second[i].first.Data()) == varName ) {
                                it->second[i].second->Add(hshape);
                                newVar=false;
                                break;
                            }
                        }
                    }

                    if(newVar) {
                        shape.signalVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
                    }
                }
            } else {
                if(varName=="")  shape.bckg.push_back(hshape);
                else             shape.bckgVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
            }
        }
        //delete syst;
    }

    //compute the total
    for(size_t i=0; i<shape.bckg.size(); i++) {
        if(i==0) {
            shape.totalBckg = (TH1 *)shape.bckg[i]->Clone(ch+"_"+shapeName+"_total");
            shape.totalBckg->SetDirectory(0);
        } else     {
            shape.totalBckg->Add(shape.bckg[i]);
        }
    }

    if(MCclosureTest) {
        if(shape.totalBckg) {
            if(!shape.data) {
                shape.data=(TH1F*)shape.totalBckg->Clone("data");
                shape.data->SetDirectory(0);
                shape.data->SetTitle("data");
            } else {
                shape.data->Reset();
                shape.data->Add(shape.totalBckg, 1);
            }
        }
    }

    //all done
    return shape;
}


//
void getYieldsFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName, bool isdataBlinded)
{
    cout << endl;
    cout << "########################## getYieldsFromShape ##########################" << endl;

    TString massStr("");
    if(mass>0)massStr += mass;

    TString MVStr("");
    if(MV>0)MVStr += MV;

    TString MAStr("");
    if(MA>0)MAStr += MA;

    TString K1Str("");
    std::ostringstream K1out;
    K1out << std::setprecision(2) << K1;
    std::string K1str = K1out.str();
    if(K1>0) K1Str = K1str;

    TString K2Str("");
    std::ostringstream K2out;
    K2out << std::setprecision(2) << K2;
    std::string K2str = K2out.str();
    if(K2>0) K2Str = K2str;

    YIELDS_T ee0jet_Yields;
    YIELDS_T ee1jet_Yields;
    YIELDS_T eelesq1jet_Yields;
    YIELDS_T mm0jet_Yields;
    YIELDS_T mm1jet_Yields;
    YIELDS_T mmlesq1jet_Yields;

    TH1* h;
    Double_t valerr, val;// syst;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t ich=0; ich<ch.size(); ich++) {
            YIELDS_T fortableYields;
            TString YieldsLabel = ch[ich]+AnalysisBins[b];
            cout << "Channels: " << YieldsLabel << endl;

            //nbckgs
            size_t nbckg=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.bckg.size();
            double sum_allbkgs(0.);
            double err_allbkgs(0.);

            for(size_t ibckg=0; ibckg<nbckg; ibckg++) {
                TH1* h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.bckg[ibckg];
                TString procTitle(h->GetTitle());
                cout << "backgrounds:  " << procTitle << endl;

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);

                if(procTitle.Contains("ZZ#rightarrow 2l2#nu")) {
                    fortableYields.ZZ = val;
                    fortableYields.ZZ_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("WZ#rightarrow 3l#nu")) {
                    fortableYields.WZ = val;
                    fortableYields.WZ_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("Z+jets (data)")) {
                    fortableYields.Zjets = val;
                    fortableYields.Zjets_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("Top/WW/Ztautau (data)")) {
                    fortableYields.WWtop = val;
                    fortableYields.WWtop_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("Wjets (data)")) {
                    fortableYields.Wjets = val;
                    fortableYields.Wjets_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                }
            } //nbckgs END!

            fortableYields.totBkg = sum_allbkgs;
            fortableYields.totBkg_StatErr = sqrt(err_allbkgs);

            //signal
            size_t nsig=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal.size();
            for(size_t isig=0; isig<nsig; isig++) {
                h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal[isig];
                TString procTitle(h->GetTitle());
                procTitle.ReplaceAll("#","\\");

                cout << "Signals >>>>>>> " << procTitle << endl;

                if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("ZZ"))procTitle = "ggH("+massStr+")";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("ZZ"))procTitle = "qqH("+massStr+")";
                else if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("WW"))procTitle = "ggH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("WW"))procTitle = "qqH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("ZH")                             )procTitle = "ZH("+massStr+")2lMET";

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);

                fortableYields.Sig = val;
                fortableYields.Sig_StatErr = valerr;

            }// signal END!


            //data
            h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.data;
            val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);

            fortableYields.Data = val;
            fortableYields.Data_StatErr = valerr;

            if(YieldsLabel.Contains("eeeq0jets")) 	ee0jet_Yields = fortableYields;
            if(YieldsLabel.Contains("eeeq1jets")) 	ee1jet_Yields = fortableYields;
            if(YieldsLabel.Contains("eelesq1jets")) 	eelesq1jet_Yields = fortableYields;
            if(YieldsLabel.Contains("mumueq0jets")) 	mm0jet_Yields = fortableYields;
            if(YieldsLabel.Contains("mumueq1jets")) 	mm1jet_Yields = fortableYields;
            if(YieldsLabel.Contains("mumulesq1jets")) 	mmlesq1jet_Yields = fortableYields;

        }//end channel
    }// end analysis bin

    //print out yields table

    FILE* pFile = fopen("Yields.tex","w");
    fprintf(pFile,"\n");
    fprintf(pFile,"\\begin{table}[h!]\n");
    fprintf(pFile,"\\begin{center}\n");
    fprintf(pFile,"\\begin{tabular}{c|cc|cc}\n");
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"Process & $ee$ & $\\mu\\mu$ & $ee$ & $\\mu\\mu$ \\\\\n");
    fprintf(pFile,"        & \\multicolumn{2}{c|}{0 jet selection} & \\multicolumn{2}{c}{1 jet selection} \\\\\n");
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"%40s & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f \\\\\n"
            ,"Signal"
            ,ee0jet_Yields.Sig, ee0jet_Yields.Sig_StatErr
            ,mm0jet_Yields.Sig, mm0jet_Yields.Sig_StatErr
            ,ee1jet_Yields.Sig, ee1jet_Yields.Sig_StatErr
            ,mm1jet_Yields.Sig, mm1jet_Yields.Sig_StatErr
           );

    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$Z/\\gamma^*\\rightarrow\\ell^+\\ell^-$"
            ,ee0jet_Yields.Zjets, ee0jet_Yields.Zjets_StatErr
            ,mm0jet_Yields.Zjets, mm0jet_Yields.Zjets_StatErr
            ,ee1jet_Yields.Zjets, ee1jet_Yields.Zjets_StatErr
            ,mm1jet_Yields.Zjets, mm1jet_Yields.Zjets_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$WZ\\rightarrow 3\\ell\\nu$"
            ,ee0jet_Yields.WZ, ee0jet_Yields.WZ_StatErr
            ,mm0jet_Yields.WZ, mm0jet_Yields.WZ_StatErr
            ,ee1jet_Yields.WZ, ee1jet_Yields.WZ_StatErr
            ,mm1jet_Yields.WZ, mm1jet_Yields.WZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$ZZ\\rightarrow 2\\ell2\\nu$"
            ,ee0jet_Yields.ZZ, ee0jet_Yields.ZZ_StatErr
            ,mm0jet_Yields.ZZ, mm0jet_Yields.ZZ_StatErr
            ,ee1jet_Yields.ZZ, ee1jet_Yields.ZZ_StatErr
            ,mm1jet_Yields.ZZ, mm1jet_Yields.ZZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"Top/WW/$Z\\to\\tau^+\\tau^-$"
            ,ee0jet_Yields.WWtop, ee0jet_Yields.WWtop_StatErr
            ,mm0jet_Yields.WWtop, mm0jet_Yields.WWtop_StatErr
            ,ee1jet_Yields.WWtop, ee1jet_Yields.WWtop_StatErr
            ,mm1jet_Yields.WWtop, mm1jet_Yields.WWtop_StatErr
           );


    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"W+jets"
            ,ee0jet_Yields.Wjets, ee0jet_Yields.Wjets_StatErr
            ,mm0jet_Yields.Wjets, mm0jet_Yields.Wjets_StatErr
            ,ee1jet_Yields.Wjets, ee1jet_Yields.Wjets_StatErr
            ,mm1jet_Yields.Wjets, mm1jet_Yields.Wjets_StatErr
           );
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"total bkg."
            ,ee0jet_Yields.totBkg, ee0jet_Yields.totBkg_StatErr
            ,mm0jet_Yields.totBkg, mm0jet_Yields.totBkg_StatErr
            ,ee1jet_Yields.totBkg, ee1jet_Yields.totBkg_StatErr
            ,mm1jet_Yields.totBkg, mm1jet_Yields.totBkg_StatErr
           );

    fprintf(pFile,"\\hline\n");
    if(isdataBlinded) {
        fprintf(pFile,"%40s & %s & %s & %s & %s \\\\\n","Data","-","-","-","-");
    } else {
        fprintf(pFile,"%40s & %.0f & %.0f & %.0f & %.0f \\\\\n","Data",ee0jet_Yields.Data,mm0jet_Yields.Data,ee1jet_Yields.Data,mm1jet_Yields.Data);
    }
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");

    fprintf(pFile,"\\end{tabular}\n");
    fprintf(pFile,"\\caption{Observed(blinded) yields, pre-fit background estimates and signal predictions.}\n");
    fprintf(pFile,"\\label{tab:seltable}\n");
    fprintf(pFile,"\\end{center}\n");
    fprintf(pFile,"\\end{table}\n");
    fprintf(pFile,"\n");
    fclose(pFile);


    // printout <= 0,1 jet
    pFile = fopen("Yields_lesq1jet.tex","w");
    fprintf(pFile,"\n");
    fprintf(pFile,"\\begin{table}[h!]\n");
    fprintf(pFile,"\\begin{center}\n");
    fprintf(pFile,"\\begin{tabular}{c|c|c}\n");
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"Process & $ee$ & $\\mu\\mu$ \\\\\n");
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f \\\\\n"
            ,"Signal"
            ,eelesq1jet_Yields.Sig, eelesq1jet_Yields.Sig_StatErr
            ,mmlesq1jet_Yields.Sig, mmlesq1jet_Yields.Sig_StatErr
           );

    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$Z/\\gamma^*\\rightarrow\\ell^+\\ell^-$"
            ,eelesq1jet_Yields.Zjets, eelesq1jet_Yields.Zjets_StatErr
            ,mmlesq1jet_Yields.Zjets, mmlesq1jet_Yields.Zjets_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$WZ\\rightarrow 3\\ell\\nu$"
            ,eelesq1jet_Yields.WZ, eelesq1jet_Yields.WZ_StatErr
            ,mmlesq1jet_Yields.WZ, mmlesq1jet_Yields.WZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$ZZ\\rightarrow 2\\ell2\\nu$"
            ,eelesq1jet_Yields.ZZ, eelesq1jet_Yields.ZZ_StatErr
            ,mmlesq1jet_Yields.ZZ, mmlesq1jet_Yields.ZZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"Top/WW/$Z\\to\\tau^+\\tau^-$"
            ,eelesq1jet_Yields.WWtop, eelesq1jet_Yields.WWtop_StatErr
            ,mmlesq1jet_Yields.WWtop, mmlesq1jet_Yields.WWtop_StatErr
           );


    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"W+jets"
            ,eelesq1jet_Yields.Wjets, eelesq1jet_Yields.Wjets_StatErr
            ,mmlesq1jet_Yields.Wjets, mmlesq1jet_Yields.Wjets_StatErr
           );
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"total bkg."
            ,eelesq1jet_Yields.totBkg, eelesq1jet_Yields.totBkg_StatErr
            ,mmlesq1jet_Yields.totBkg, mmlesq1jet_Yields.totBkg_StatErr
           );

    fprintf(pFile,"\\hline\n");
    if(isdataBlinded) {
        fprintf(pFile,"%40s & %s & %s \\\\\n","Data","-","-");
    } else {
        fprintf(pFile,"%40s & %.0f & %.0f \\\\\n","Data",eelesq1jet_Yields.Data,mmlesq1jet_Yields.Data);
    }
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");

    fprintf(pFile,"\\end{tabular}\n");
    fprintf(pFile,"\\caption{Observed(blinded) yields, pre-fit background estimates and signal predictions.}\n");
    fprintf(pFile,"\\label{tab:seltable}\n");
    fprintf(pFile,"\\end{center}\n");
    fprintf(pFile,"\\end{table}\n");
    fprintf(pFile,"\n");
    fclose(pFile);



    cout << "########################## getYieldsFromShape ##########################" << endl;
    cout << endl;
}




std::vector<TString>  buildDataCard(Int_t mass, TString histo, TString url, TString Json)
{
    std::vector<TString> dcUrls;

    //get the datacard inputs
    DataCardInputs dci = convertHistosForLimits(mass,histo,url,Json);

    TString eecard = "";
    TString mumucard = "";
    TString combinedcard = "";
    TString combinedcardLL = "";

    //build the datacard separately for each channel
    for( auto const & channel : dci.ch ) {
        TString dcName=dci.shapesFile;
        dcName.ReplaceAll(".root","_"+channel+".dat");
        FILE* pFile = fopen(dcName.Data(),"w");

        if(!channel.Contains("ll")) combinedcard += channel+"="+dcName+" ";
        if(channel.Contains("lleq"))  combinedcardLL += channel+"="+dcName+" ";
        if(channel.Contains("ee"))    eecard += channel+"="+dcName+" ";
        if(channel.Contains("mumu"))  mumucard += channel+"="+dcName+" ";

        //header
        fprintf(pFile, "imax 1\n");
        fprintf(pFile, "jmax *\n");
        fprintf(pFile, "kmax *\n");
        fprintf(pFile, "-------------------------------\n");
        if(shape) {
            fprintf(pFile, "shapes * * %s %s/$PROCESS %s/$PROCESS_$SYSTEMATIC\n",dci.shapesFile.Data(), channel.Data(), channel.Data());
            fprintf(pFile, "-------------------------------\n");
        }
        //observations
        fprintf(pFile, "bin 1\n");
        fprintf(pFile, "Observation %f\n",dci.obs[RateKey_t("obs",channel)]);
        fprintf(pFile, "-------------------------------\n");

        fprintf(pFile,"%55s ", "bin");
        for(auto const & proc : dci.procs ) {
            if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
            fprintf(pFile,"%6i ", 1);
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%55s ", "process");
        for(auto const & proc : dci.procs ) {
            if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
            fprintf(pFile,"%6s ", (postfix+proc).Data());
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%55s ", "process");
        int procCtr(1-dci.nsignalproc);
        for(auto const & proc : dci.procs ) {
            if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
            fprintf(pFile,"%6i ", procCtr );
            procCtr++;
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%55s ", "rate");
        for(auto const & proc : dci.procs ) {
            if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
            fprintf(pFile,"%6f ", dci.rates[RateKey_t(proc,channel)] );
        }
        fprintf(pFile,"\n");
        fprintf(pFile, "-------------------------------\n");


        //systematics
        char sFile[2048];
        bool isSyst;
        if(runSystematics) {
            if(systpostfix.Contains("13")) {
                fprintf(pFile,"%45s %10s ", "lumi_13TeV", "lnN");
                for(auto const & proc : dci.procs ) {
                    if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
                    if(!proc.Contains("WWTopZtautau") && !proc.Contains("Zjets") && !proc.Contains("Wjets")) {
                        fprintf(pFile,"%6.5f ",1.0+normSysts["lumi_13TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }

            //leptont efficiency
            if(channel.Contains("ee")) {
                fprintf(pFile,"%45s %10s ", "CMS_eff_e", "lnN");
                for(auto const & proc : dci.procs ) {
                    if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
                    if(!proc.Contains("WWTopZtautau") && !proc.Contains("Zjets") /*&& !proc.Contains("Wjets")*/) {
                        fprintf(pFile,"%6.5f ",1.0+normSysts["CMS_eff_e"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else {
                fprintf(pFile,"%45s %10s ", "CMS_eff_m", "lnN");
                for(auto const & proc : dci.procs ) {
                    if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
                    if(!proc.Contains("WWTopZtautau") && !proc.Contains("Zjets") /*&& !proc.Contains("Wjets")*/) {
                        fprintf(pFile,"%6.5f ",1.0+normSysts["CMS_eff_m"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }


            for(std::map<TString, std::map<RateKey_t,Double_t> >::iterator it=dci.systs.begin(); it!=dci.systs.end(); it++) {
                if(!runSystematics && string(it->first.Data()).find("stat")>0 )continue;

                isSyst=false;
                if(it->first.Contains("_sys_")) {
                    sprintf(sFile,"%45s %10s ", it->first.Data(), "lnN");
                    for(auto const & proc : dci.procs ) {
                        if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;
                        if(it->second.find(RateKey_t(proc,channel)) != it->second.end()) {
                            Double_t systUnc = it->second[RateKey_t(proc,channel)];
                            if(systUnc<=0) {
                                sprintf(sFile,"%s%6s ",sFile,"-");
                            } else {
                                sprintf(sFile,"%s%6.5f ",sFile,(1.0+ (systUnc / dci.rates[RateKey_t(proc,channel)]) ));
                                isSyst=true;
                            }
                        } else {
                            sprintf(sFile,"%s%6s ",sFile,"-");
                        }
                    }
                    if(isSyst)fprintf(pFile,"%s\n",sFile);

                } else {
                    if(shape) {
                        sprintf(sFile,"%45s %10s ", it->first.Data(), "shape");
                    } else {
                        sprintf(sFile,"%45s %10s ", it->first.Data(), "lnN");
                    }
                    for(auto const & proc : dci.procs ) {
                        if(dci.rates.find(RateKey_t(proc,channel))==dci.rates.end()) continue;

                        if(it->second.find(RateKey_t(proc,channel)) != it->second.end()) {
                            //cout << "proc: " << proc << " it->first: " << it->first  << endl;
                            if(!shape) sprintf(sFile,"%s%6.5f ",sFile,it->second[RateKey_t(proc,channel)]);
                            else       sprintf(sFile,"%s%6s ",sFile,"1.0");
                            isSyst=true;
                        } else {
                            sprintf(sFile,"%s%6s ",sFile,"-");
                        }
                    }
                    if(isSyst)fprintf(pFile,"%s\n",sFile);
                }
            }//end of shape uncertainty


        }

        fclose(pFile);
        cout << "Data card for " << dci.shapesFile << " and " << channel << " channel @ " << dcName << endl;
        dcUrls.push_back(dcName);
    }

    FILE* pFile = fopen("combineCards.sh","w");
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcard + " > " + "card_combined.dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + eecard       + " > " + "card_ee.dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + mumucard     + " > " + "card_mumu.dat").Data());
    fprintf(pFile,"%s;\n",(TString("#combine -M Asymptotic --cl 0.95 --rRelAcc 0.00000001 --rAbsAcc 0.000000001 -m 1 -n Zwimps01jets --run expected --expectSignal=1 -t -1 card_combined.dat > COMB_asympt.log")).Data());
    fclose(pFile);

    return dcUrls;
}

//
DataCardInputs convertHistosForLimits(Int_t mass,TString histo,TString url,TString Json)
{
    DataCardInputs datacard_inputs;

    //init the json wrapper
    JSONWrapper::Object Root(Json.Data(), true);

    //init globalVariables
    TString massStr("");
    if(mass>0)massStr += mass;
    std::vector<TString> allCh,allProcs;

    TString MVStr("");
    if(MV>0)MVStr += MV;

    TString MAStr("");
    if(MA>0)MAStr += MA;


    TString K1Str("");
    std::ostringstream K1out;
    K1out << std::setprecision(2) << K1;
    std::string K1str = K1out.str();
    if(K1>0) K1Str = K1str;

    TString K2Str("");
    std::ostringstream K2out;
    K2out << std::setprecision(2) << K2;
    std::string K2str = K2out.str();
    if(K2>0) K2Str = K2str;




    //open input file
    TFile* input_file = TFile::Open(url);
    if( !input_file || input_file->IsZombie() ) {
        cout << "Invalid file name : " << url << endl;
        return datacard_inputs;
    }

    //get the shapes for each channel
    map<TString, Shape_t> shapes_map;
    std::vector<TString> channels = {"mumu","ee","emu","ll"};
    std::vector<TString> shapes;
    shapes.push_back(histo);
    if(subNRB2011 || subNRB2012)shapes.push_back(histo+"_NRBctrl");
    if(subNRB2012)		shapes.push_back(histo+"BTagSB");
    if(subWZ)			shapes.push_back(histo+"_3rdLepton");
    shapes.push_back(histo+"_NRBsyst");
    shapes.push_back("pfmet_minus_shapes");
    const size_t n_shapes=shapes.size();
    for ( auto & i_channel : channels ) {
        for( auto & i_bin : AnalysisBins ) {
            int indexcut_ = indexcut;
            double cutMin=shapeMin;
            double cutMax=shapeMax;
            for( auto & i_shape : shapes ) {
                shapes_map[i_channel+i_bin+i_shape]=getShapeFromFile(input_file, i_channel+i_bin,i_shape,indexcut_,Root,cutMin, cutMax);
                cout << "Loading shapes: " << i_channel+i_bin+i_shape << endl;
             }
        }
    }

    //all done with input file
    input_file->Close();

    //printf("done loading all shapes\n");

    //define vector for search
    std::vector<TString>& channels_for_limits = Channels;

    //// Non-resonant backgrounds are estimated from data.
    // Perform MC closure test, assign systematic uncertainties based on outcome
    dodataDrivenWWtW(channels_for_limits,"emu",shapes_map,histo,true);
    // Create the real estimate
    dodataDrivenWWtW(channels_for_limits,"emu",shapes_map,histo,false);

    // MC Z+jets -> DYExtrapolation
    doDYextrapolation(channels_for_limits,shapes_map,histo,"pfmet_minus_shapes", DYMET_EXTRAPOL, true);

    //replace data by total MC background
    if(blindData) BlindData(channels_for_limits,shapes_map,histo, blindWithSignal);

    //print event yields from the mt shapes
    //if(!fast)getYieldsFromShape(channels_for_limits,shapes_map,histo,true); //blind data
    if(!fast)getYieldsFromShape(channels_for_limits,shapes_map,histo,false); //unblind data



    //prepare the output
    datacard_inputs.shapesFile="zllwimps_"+massStr+systpostfix+".root";
    TFile *fout=TFile::Open(datacard_inputs.shapesFile,"recreate");

    //loop on channel/proc/systematics
    for(size_t ich=0; ich<channels_for_limits.size(); ich++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            TString chbin = channels_for_limits[ich]+AnalysisBins[b];
            fout->mkdir(chbin);
            fout->cd(chbin);
            allCh.push_back(chbin);
            Shape_t shapeSt = shapes_map.find(chbin+histo)->second;

            //signals
            datacard_inputs.nsignalproc = 0;

            for( auto signal_histo : shapeSt.signal ) {
                TString proc(signal_histo->GetTitle());

                cout << "############## Signal: " << proc << "##############" << endl;
                std::vector<TString> systs;
                std::vector<TH1*>    hshapes;
                systs.push_back("");
                hshapes.push_back(signal_histo);
                cout << "\n" << proc << " has rate: " << signal_histo->Integral() << endl;

                for( auto & signal_histo_variation : shapeSt.signalVars[signal_histo->GetTitle()] ) {
                    TString variation_name = signal_histo_variation.first;
                    TH1* variation_histo = signal_histo_variation.second;
                    printf("SYSTEMATIC FOR SIGNAL %s : %s: %f\n",proc.Data(), variation_name.Data(), variation_histo->Integral());
                    systs.push_back(signal_histo_variation.first);
                    hshapes.push_back(signal_histo_variation.second);
                }

                convertHistosForLimits_core(datacard_inputs, proc, AnalysisBins[b], chbin, systs, hshapes);
                if(ich==0 && b==0)allProcs.push_back(proc);
                datacard_inputs.nsignalproc++;
            }

            //backgrounds
            size_t nbckg=shapes_map.find(chbin+histo)->second.bckg.size();
            size_t nNonNullBckg=0;
            for(size_t ibckg=0; ibckg<nbckg; ibckg++) {
                TH1* h=shapeSt.bckg[ibckg];
                if(h->Integral()<=1e-10) continue;
                cout << "\n" << h->GetTitle() << " has rate: " << h->Integral() << endl;

                //check if any bin content has a negative value
                double tot_integralval = h->Integral();
                for(int bin=0; bin<h->GetXaxis()->GetNbins()+1; bin++) {
                    //cout << "bin: " << bin << " val: " << h->GetBinContent(bin) << endl;
                    if(h->GetBinContent(bin)<0) h->SetBinContent(bin,0);
                }
                h->Scale(tot_integralval/h->Integral());

                std::vector<std::pair<TString, TH1*> > vars = shapeSt.bckgVars[h->GetTitle()];

                std::vector<TString> systs;
                std::vector<TH1*>    hshapes;
                systs.push_back("");
                hshapes.push_back(shapeSt.bckg[ibckg]);
                for(size_t v=0; v<vars.size(); v++) {
                    printf("SYSTEMATIC FOR BACKGROUND %s : %s : %f\n",h->GetTitle(), vars[v].first.Data(), vars[v].second->Integral());

                    TH1* h= vars[v].second;

                    //check if any bin content has a negative value
                    double tot_integralval = h->Integral();
                    for(int bin=0; bin<h->GetXaxis()->GetNbins()+1; bin++) {
                        if(h->GetBinContent(bin)<0) h->SetBinContent(bin,0);
                    }
                    h->Scale(tot_integralval/h->Integral());

                    systs.push_back(vars[v].first);
                    hshapes.push_back(h/*vars[v].second*/);
                }

                TString proc(h->GetTitle());
                convertHistosForLimits_core(datacard_inputs, proc, AnalysisBins[b], chbin, systs, hshapes);
                if(ich==0 && b==0)allProcs.push_back(proc);

                //remove backgrounds with rate=0 (but keep at least one background)
                std::map<RateKey_t, Double_t>::iterator it = datacard_inputs.rates.find(RateKey_t(proc,chbin));
                if(it==datacard_inputs.rates.end()) {
                    printf("proc=%s not found --> THIS SHOULD NEVER HAPPENS.  PLEASE CHECK THE CODE\n",proc.Data());
                } else {
                    if(it->second>0) {
                        nNonNullBckg++;
                    } else if(ibckg<nbckg-1 || nNonNullBckg>0) {
                        datacard_inputs.rates.erase(datacard_inputs.rates.find(RateKey_t(proc,chbin)));
                    }
                }

            }

            //data
            TH1* h=shapeSt.data;
            std::vector<TString> systs;
            std::vector<TH1*>    hshapes;
            systs.push_back("");
            hshapes.push_back(h);
            TString proc(h->GetTitle());
            convertHistosForLimits_core(datacard_inputs, proc, AnalysisBins[b], chbin, systs, hshapes);

            //return to parent dir
            fout->cd("..");
        }
    }

    datacard_inputs.ch.resize(allCh.size());
    std::copy(allCh.begin(), allCh.end(),datacard_inputs.ch.begin());
    datacard_inputs.procs.resize(allProcs.size());
    std::copy(allProcs.begin(), allProcs.end(),datacard_inputs.procs.begin());



    //all done
    fout->Close();

    return datacard_inputs;
}

void rename_process( TString & proc ) {
    proc.ReplaceAll("#bar{t}","tbar");
    proc.ReplaceAll("Z-#gamma^{*}+jets#rightarrow ll","dy");
    proc.ReplaceAll("#rightarrow","");
    proc.ReplaceAll("(","");
    proc.ReplaceAll(")","");
    proc.ReplaceAll("+","");
    proc.ReplaceAll(" ","");
    proc.ReplaceAll("/","");
    proc.ReplaceAll("#","");
    proc.ReplaceAll("=","");
    proc.ReplaceAll(".","");
    proc.ReplaceAll("^","");
    proc.ReplaceAll("}","");
    proc.ReplaceAll("{","");
    proc.ToLower();

    proc.ReplaceAll("zh1052lmet","ZH_hinv");
    proc.ReplaceAll("zh1152lmet","ZH_hinv");
    proc.ReplaceAll("zh1252lmet","ZH_hinv");
    proc.ReplaceAll("zh1352lmet","ZH_hinv");
    proc.ReplaceAll("zh1452lmet","ZH_hinv");
    proc.ReplaceAll("zh1752lmet","ZH_hinv");
    proc.ReplaceAll("zh2002lmet","ZH_hinv");
    proc.ReplaceAll("zh3002lmet","ZH_hinv");

    proc.ReplaceAll("unpart101","UnPart1p01");
    proc.ReplaceAll("unpart102","UnPart1p02");
    proc.ReplaceAll("unpart104","UnPart1p04");
    proc.ReplaceAll("unpart106","UnPart1p06");
    proc.ReplaceAll("unpart109","UnPart1p09");
    proc.ReplaceAll("unpart110","UnPart1p10");
    proc.ReplaceAll("unpart120","UnPart1p20");
    proc.ReplaceAll("unpart130","UnPart1p30");
    proc.ReplaceAll("unpart140","UnPart1p40");
    proc.ReplaceAll("unpart150","UnPart1p50");
    proc.ReplaceAll("unpart160","UnPart1p60");
    proc.ReplaceAll("unpart170","UnPart1p70");
    proc.ReplaceAll("unpart180","UnPart1p80");
    proc.ReplaceAll("unpart190","UnPart1p90");
    proc.ReplaceAll("unpart200","UnPart2p00");
    proc.ReplaceAll("unpart220","UnPart2p20");

    proc.ReplaceAll("zz2l2nu","ZZ");
    proc.ReplaceAll("wz3lnu","WZ");
    proc.ReplaceAll("vvv","VVV");
    proc.ReplaceAll("zjetsdata","Zjets");
    proc.ReplaceAll("topwwztautaudata","WWTopZtautau");
    proc.ReplaceAll("wjetsdata","Wjets");
    proc.ReplaceAll("wjets","Wjets");
}

void rename_systematic( TString & syst, const TString & proc, const TString & channel ) {
    syst.ReplaceAll("down","Down");
    syst.ReplaceAll("up","Up");

    if(syst=="") {
        syst="";
    } else if(syst.BeginsWith("_jes")) {
        syst.ReplaceAll("_jes","_CMS_scale_j");
    } else if(syst.BeginsWith("_jer")) {
        syst.ReplaceAll("_jer","_CMS_res_j");
    } else if(syst.BeginsWith("_btag")) {
        syst.ReplaceAll("_btag","_CMS_eff_b");
    } else if(syst.BeginsWith("_pu" )) {
        syst.ReplaceAll("_pu", "_CMS_zllwimps_pu");
    } else if(syst.BeginsWith("_les" )) {
        if(channel.Contains("ee")) syst.ReplaceAll("_les", "_CMS_scale_e");
        if(channel.Contains("mumu")) syst.ReplaceAll("_les", "_CMS_scale_m");
    } else if(syst.BeginsWith("_umet" )) {
        syst.ReplaceAll("_umet", "_CMS_scale_met");
    } else if(syst.BeginsWith("_qcdscale")) {
        if ( proc.Contains("ZZ") || proc.Contains("WZ") )  syst.ReplaceAll("_qcdscale", "_QCDscale_VV");
        else if ( proc.Contains("VVV") ) 		       syst.ReplaceAll("_qcdscale", "_QCDscale_VVV");
        else if ( proc.Contains("ewk_s_dm") )              syst.ReplaceAll("_qcdscale", "_QCDscale_EWKDM");
        else if ( proc.Contains("UnPart") )		       syst.ReplaceAll("_qcdscale", "_QCDscale_Unpart");
        else if ( proc.Contains("add") )                   syst.ReplaceAll("_qcdscale", "_QCDscale_ADD");
        else if ( proc.Contains("dm") && (proc.Contains("mv") || proc.Contains("ma")) )                   syst.ReplaceAll("_qcdscale", "_QCDscale_VDM");
        else syst.ReplaceAll("_qcdscale", "_QCDscale");
    } else if(syst.BeginsWith("_pdf")) {
        syst.ReplaceAll("_pdf", "_pdf_qqbar");
    } else if(syst.BeginsWith("_qqZZewk")) {
        syst.ReplaceAll("_qqZZewk", "_qqZZewkcorr");
    } else {
        syst="_CMS_zllwimps"+syst;
    }
}
void convertHistosForLimits_core(DataCardInputs& dci, TString& proc, TString& bin, TString& ch, std::vector<TString>& systs, std::vector<TH1*>& hshapes)
{
    // Remove special characters, apply formatting, etc.
    rename_process(proc);

    for(unsigned int i=0; i<systs.size(); i++) {
        TString syst   = systs[i];
        TH1*    hshape = hshapes[i];
        hshape->SetDirectory(0);

        // Skip uncertainties we do not want to consider
        if(syst.BeginsWith("_umet" )) continue;

        //Do Renaming and cleaning
        rename_systematic( syst, proc, ch );
        double systUncertainty = hshape->GetBinError(0);
        hshape->SetBinError(0,0);

        //If cut&count keep only 1 bin in the histo
        if(!shape) {
            hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins());
            //make sure to also count the underflow and overflow
            double bin  = hshape->GetBinContent(0) + hshape->GetBinContent(1) + hshape->GetBinContent(2);
            double bine = sqrt(pow(hshape->GetBinError(0),2) + pow(hshape->GetBinError(1),2) + pow(hshape->GetBinError(2),2));
            hshape->SetBinContent(0,0);
            hshape->SetBinError  (0,0);
            hshape->SetBinContent(1,bin);
            hshape->SetBinError  (1,bine);
            hshape->SetBinContent(2,0);
            hshape->SetBinError  (2,0);
        }

        if(syst=="") {
            //central shape (for data call it data_obs)
            hshape->SetName(proc);
            if(proc=="data") {
                hshape->Write("data_obs");
            } else {
                hshape->Write(proc+postfix);

;

                if(hshape->Integral()>0) {
                    hshape->SetName(proc+syst);
                    TH1* statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Up");
                    TH1* statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Down");
                    if(proc.Contains("WWTopZtautau")) {
                        statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Up");
                        statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Down");
                    }
                    if(proc.Contains("Zjets")) {
                        TString zjetsch="";
                        if(ch.Contains("eq0jets")) zjetsch="_eq0jets";
                        if(ch.Contains("eq1jets")) zjetsch="_eq1jets";
                        if(ch.Contains("lesq1jets")) zjetsch="_lesq1jets";
                        statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+zjetsch+proc+"Up");
                        statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+zjetsch+proc+"Down");
                    }

                    for(int ibin=1; ibin<=statup->GetXaxis()->GetNbins(); ibin++) {
                        statup  ->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin),
                                                              std::max(0.01*hshape->GetBinContent(ibin), statup  ->GetBinContent(ibin) + statup  ->GetBinError(ibin))));
                        statdown->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin),
                                                              std::max(0.01*hshape->GetBinContent(ibin), statdown->GetBinContent(ibin) - statdown->GetBinError(ibin))));
                    }

                    //############################################################################
                    //### bin-by-bin uncertainty
                    bool useBinbyBin(true);
                    //############################################################################
                    useBinbyBin &= (proc.Contains("ZH") || proc.Contains("ZZ") || proc.Contains("WZ") || proc.Contains("WWTopZtautau") ||
                                    proc.Contains("ewk_s_dm") || proc.Contains("VVV") || proc.Contains("Zjets") ||
                                    proc.Contains("dm") ||  proc.Contains("add") || proc.Contains("UnPart")
                                   );
                    useBinbyBin &= (shape);

                    TString full_syst_name("CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix);
                    if(!useBinbyBin) {
                        statup  ->Write(proc+postfix+"_"+full_syst_name+"Up");
                        statdown->Write(proc+postfix+"_"+full_syst_name+"Down");
                        if(shape) {
                            dci.systs[full_syst_name][RateKey_t(proc,ch)]=1.0;
                        } else {
                            dci.systs[full_syst_name][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                        }
                    }

                    // bin-by-bin uncertainty
                    if(useBinbyBin) {

                        hshape->SetName(proc+syst);
                        //vector<int> saveBins;
                        for(int iBin=1; iBin<=hshape->GetXaxis()->GetNbins(); iBin++) {

                            if(hshape->GetBinContent(iBin)<=0) continue;

                            stringstream Bin_t;
                            Bin_t << iBin;
                            string sbin = Bin_t.str();
                            TString Bin = sbin;

                            TH1* statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"_Bin"+Bin+"Up");
                            TH1* statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"_Bin"+Bin+"Down");


                            for(int jbin=1; jbin<=statup->GetXaxis()->GetNbins(); jbin++) {
                                if(jbin!=iBin) continue;
                                statup  ->SetBinContent(jbin,std::min(2*hshape->GetBinContent(jbin),
                                                                      std::max(0.01*hshape->GetBinContent(jbin), statup  ->GetBinContent(jbin) + statup  ->GetBinError(jbin))));
                                statdown->SetBinContent(jbin,std::min(2*hshape->GetBinContent(jbin),
                                                                      std::max(0.01*hshape->GetBinContent(jbin), statdown->GetBinContent(jbin) - statdown->GetBinError(jbin))));
                            }

                            //#########Debug###########
                            double variance = (statup->Integral()/hshapes[0]->Integral());
                            TString name_variance = proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Up";
                            if(proc.Contains("WWTopZtautau")) name_variance = proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin+"Up";
                            cout << name_variance << ": " << variance << endl;


                            if(proc.Contains("WWTopZtautau")) {
                                statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin+"Up");
                                statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin+"Down");
                            } else {
                                statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Up");
                                statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Down");
                            }

                            if(shape) {
                                if(proc.Contains("WWTopZtautau")) dci.systs["CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin][RateKey_t(proc,ch)]=1.0;
                                else dci.systs[full_syst_name+"_Bin"+Bin][RateKey_t(proc,ch)]=1.0;

                            } else {
                                if(proc.Contains("WWTopZtautau")) dci.systs["CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                                else dci.systs[full_syst_name+"_Bin"+Bin][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                            }

                        }
                    }



                    //data driven Systematics
                    if(systUncertainty>0) {
                        printf("SYST in %s - %s - %s = %f\n",bin.Data(), ch.Data(), proc.Data(), systUncertainty);
                        //makesure that syst+stat error is never bigger than 100%
                        //double valerr, val  = hshape->IntegralAndError(1,hshape->GetXaxis()->GetNbins(),valerr);
                        //if(sqrt(pow(valerr,2)+pow(systUncertainty,2))>val){systUncertainty = sqrt(std::max(0.0, pow(val,2) - pow(valerr,2)));}

                        //add syst uncertainty as bin dependent or not
                        //dci.systs["CMS_zllwimps_sys_"+bin+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;

                        TString wjetsch="";
                        if(ch.Contains("ee")) wjetsch="ee";
                        if(ch.Contains("mumu")) wjetsch="mumu";

                        if(proc.Contains("WWTopZtautau") || proc.Contains("Zjets")) dci.systs["CMS_zllwimps_sys_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                        else if(proc.Contains("Wjets")) dci.systs["CMS_zllwimps_sys_"+wjetsch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                        else dci.systs["CMS_zllwimps_sys_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                    }
                }
            }
        } else if(runSystematics && proc!="data" && (syst.Contains("Up") || syst.Contains("Down"))) {
            //if empty histogram --> no variation is applied
            if(hshape->Integral()<hshapes[0]->Integral()*0.01 || isnan((float)hshape->Integral())) {
                hshape->Reset();
                hshape->Add(hshapes[0],1);
            }

            //write variation to file
            hshape->SetName(proc+syst);
            hshape->Write(proc+postfix+syst);
        } else if(runSystematics) {
            //for one sided systematics the down variation mirrors the difference bin by bin
            hshape->SetName(proc+syst);
            hshape->Write(proc+postfix+syst+"Up");
            TH1 *hmirrorshape=(TH1 *)hshape->Clone(proc+syst+"Down");
            for(int ibin=1; ibin<=hmirrorshape->GetXaxis()->GetNbins(); ibin++) {
                //double bin = hmirrorshape->GetBinContent(ibin);
                double bin = 2*hshapes[0]->GetBinContent(ibin)-hmirrorshape->GetBinContent(ibin);
                if(bin<0)bin=0;
                hmirrorshape->SetBinContent(ibin,bin);
            }
            if(hmirrorshape->Integral()<=0)hmirrorshape->SetBinContent(1, 1E-10);
            hmirrorshape->Write(proc+postfix+syst+"Down");
        }

        if(runSystematics && syst!="") {
            TString systName(syst);
            systName.ReplaceAll("Up","");
            systName.ReplaceAll("Down","");//  systName.ReplaceAll("_","");
            if(systName.First("_")==0)systName.Remove(0,1);


            TH1 *temp=(TH1*) hshape->Clone();
            temp->Add(hshapes[0],-1);
            if(temp->Integral()!=0) {
                double Unc = 1 + fabs(temp->Integral()/hshapes[0]->Integral());
                if(dci.systs.find(systName)==dci.systs.end() || dci.systs[systName].find(RateKey_t(proc,ch))==dci.systs[systName].end() ) {
                    dci.systs[systName][RateKey_t(proc,ch)]=Unc;
                } else {
                    dci.systs[systName][RateKey_t(proc,ch)]=(dci.systs[systName][RateKey_t(proc,ch)] + Unc)/2.0;
                }
            }


            delete temp;
        } else if(proc!="data" && syst=="") {
            dci.rates[RateKey_t(proc,ch)]= hshape->Integral();
        } else if(proc=="data" && syst=="") {
            dci.obs[RateKey_t("obs",ch)]=hshape->Integral();
        }
    }

}//convertHistosForLimits_core END

void doDYextrapolation(std::vector<TString>& signal_channels,map<TString, Shape_t>& allShapes, TString mainHisto, TString DY_EXTRAPOL_SHAPES, float DY_EXTRAPOL, bool isleftCR)
{
    cout << "\n#################### doDYextrapolation #######################\n" << endl;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t i=0; i<signal_channels.size(); i++) {

            THStack *stack = new THStack("stack","stack");

            Shape_t& shapeChan = allShapes.find(signal_channels[i]+AnalysisBins[b]+DY_EXTRAPOL_SHAPES)->second;
            TH1* hChan_data=shapeChan.data;
            cout << "DY extrapolation shapes: " << signal_channels[i]+AnalysisBins[b]+"_"+DY_EXTRAPOL_SHAPES << endl;
            cout << "Bins: " << hChan_data->GetXaxis()->GetNbins() << endl;

            TH1* hChan_MCnonDY = (TH1*)shapeChan.totalBckg->Clone("hChan_totMCnonDY");
            TH1* hChan_MCDY=NULL;//= (TH1*)shapeChan.totalBckg->Clone("hChan_MCDY");
            hChan_MCnonDY->Reset();
            //hChan_MCDY->Reset();
            for(size_t ibckg=0; ibckg<shapeChan.bckg.size(); ibckg++) {
                TString proc(shapeChan.bckg[ibckg]->GetTitle());
                //cout << "proc: " << proc << endl;
                if(proc.Contains("Z+jets")) {
                    hChan_MCDY=(TH1*)shapeChan.bckg[ibckg]->Clone("hChan_MCDY");
                    hChan_MCDY->Reset();
                    hChan_MCDY->Add(shapeChan.bckg[ibckg], 1);
                } else {
                    hChan_MCnonDY->Add(shapeChan.bckg[ibckg], 1);
                    stack->Add(shapeChan.bckg[ibckg],"HIST");
                }

                //stack->Add(shapeChan.bckg[ibckg],"HIST");
            }

            int thre_bin(0);
            int total_bin=hChan_MCnonDY->GetXaxis()->GetNbins();
            for(int bin=0; bin<=hChan_MCnonDY->GetXaxis()->GetNbins()+1; bin++) {
                double lowedge = hChan_MCnonDY->GetBinLowEdge(bin);
                if(lowedge >= DY_EXTRAPOL) {
                    thre_bin = bin;
                    break;
                }
            }

            double data_CR(0.);
            double subt_CR(0.);
            double zjet_SR(0.),zjet_CR(0.);
            double data_err2_CR(0.);
            double subt_err2_CR(0.);
            double zjet_err2_SR(0.),zjet_err2_CR(0.);

            if(isleftCR) {
                data_CR = hChan_data	->IntegralAndError(0,thre_bin-1,data_err2_CR);
                subt_CR = hChan_MCnonDY	->IntegralAndError(0,thre_bin-1,subt_err2_CR);
                zjet_CR = hChan_MCDY	->IntegralAndError(0,thre_bin-1,zjet_err2_CR);
                zjet_SR = hChan_MCDY	->IntegralAndError(thre_bin,total_bin+1,zjet_err2_SR);
            } else {
                data_CR = hChan_data	->IntegralAndError(thre_bin,total_bin+1,data_err2_CR);
                subt_CR = hChan_MCnonDY	->IntegralAndError(thre_bin,total_bin+1,subt_err2_CR);
                zjet_CR = hChan_MCDY	->IntegralAndError(thre_bin,total_bin+1,zjet_err2_CR);
                zjet_SR = hChan_MCDY	->IntegralAndError(0,thre_bin-1,zjet_err2_SR);
            }

            double sf=(data_CR-subt_CR)/zjet_CR;
            if(sf<0) sf=1;
            double zjet_Est = sf*zjet_SR;
            cout << "sf: " << sf << "\t" << "Est: " << zjet_Est << "\t" << "MC: " << zjet_SR << endl;

            hChan_MCDY->Scale(sf);
            stack->Add(hChan_MCDY,"HIST");


            TCanvas *c = new TCanvas("c", "c", 700, 700);
            TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
            t1->Draw();
            t1->SetLogy(true);
            t1->cd();
            //t1->SetBottomMargin(0.3);
            t1->SetRightMargin(0.05);
            //c->Divide(1,2);

            stack->Draw("");
            hChan_data->Draw("E1 same");

            stack->SetMinimum(5e-3);
            stack->GetXaxis()->SetTitle(hChan_data->GetXaxis()->GetTitle());
            stack->GetYaxis()->SetTitle("Events");

            c->SaveAs(signal_channels[i]+AnalysisBins[b]+"_"+DY_EXTRAPOL_SHAPES+"_DYEXTRAPOL.png");
            delete c;


            //add data-driven backgrounds
            Shape_t& shapeChan_SI = allShapes.find(signal_channels[i]+AnalysisBins[b]+mainHisto)->second;
            TH1* hChan_DATADY = NULL;

            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("Z+jets")) {
                    hChan_DATADY=(TH1*)shapeChan_SI.bckg[ibckg]->Clone("hChan_DATADY");
                    hChan_DATADY->SetTitle("Z+jets (data)");
                    //hChan_DATADY->Reset();
                    //remove MC Z+jets
                    shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                    ibckg--;
                }
            }

            cout << "Scale Factor: " << sf << " DY MC: " << hChan_DATADY->Integral() << endl;
            hChan_DATADY->Scale(sf);

            // 100% systematics
            hChan_DATADY->SetBinError(0, 1.0*hChan_DATADY->Integral());

            cout << "DY Extrapolation: " << hChan_DATADY->Integral() << endl;
            if(hChan_DATADY->Integral() > 1E-6) shapeChan_SI.bckg.push_back(hChan_DATADY);

            cout << "\n";


        }
    }
    cout << "\n#################### doDYextrapolation #######################\n" << endl;
}


void doWjetsBackground(std::vector<TString>& signal_channels,map<TString, Shape_t>& allShapes, TString mainHisto)
{
    cout << "\n#################### doWjetsBackground #######################\n" << endl;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t i=0; i<signal_channels.size(); i++) {
            TString label_Syst = signal_channels[i]+AnalysisBins[b];

            Shape_t& shapeChan_SI = allShapes.find(signal_channels[i]+AnalysisBins[b]+mainHisto)->second;
            Shape_t& shapeChan_Wjet = allShapes.find(signal_channels[i]+AnalysisBins[b]+"FR_WjetCtrl_"+mainHisto)->second;
            TH1* hChan_data=shapeChan_Wjet.data;

            double valerr_hChan_data, val_hChan_data;
            val_hChan_data = hChan_data->IntegralAndError(1,hChan_data->GetXaxis()->GetNbins(),valerr_hChan_data);

            cout << "Wjets Shapes: " << signal_channels[i]+AnalysisBins[b]+"_FR_WjetCtrl_"+mainHisto << "\t"
                 << val_hChan_data << endl;

            TH1* h_Wjet = NULL;
            h_Wjet = (TH1*)hChan_data->Clone("Wjets (data)");
            h_Wjet->SetTitle("Wjets (data)");

            //remove mc Wjet estimation
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("W+jets")) {
                    shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                    ibckg--;
                }
            }

            //add the background estimate
            if(val_hChan_data <= 0) continue;

            double minVal_bin(999.);
            for(int b=1; b<=h_Wjet->GetXaxis()->GetNbins()+1; b++) {
                double val = h_Wjet->GetBinContent(b);
                double err = h_Wjet->GetBinError(b);

                cout << "bin: " << b << "\t" << val << " +/- " << err << endl;
                if(val<minVal_bin) minVal_bin = val;
            }

            if(minVal_bin<0) {
                for(int b=1; b<=h_Wjet->GetXaxis()->GetNbins()+1; b++) {
                    double val = h_Wjet->GetBinContent(b);
                    if(val!=0) val -= minVal_bin;
                    //cout << "bin: " << b << "\t" << val << endl;
                    h_Wjet->SetBinContent(b, val);
                }
            }

            double newsum = h_Wjet->Integral();
            h_Wjet->Scale(val_hChan_data/newsum);

            cout << ">>>>>>>>>> remove negative bins <<<<<<<<<<" << endl;
            for(int b=1; b<=h_Wjet->GetXaxis()->GetNbins()+1; b++) {
                double val = h_Wjet->GetBinContent(b);
                double err = h_Wjet->GetBinError(b);
                cout << "bin: " << b << "\t" << val << " +/- " << err << endl;
            }


            double sysErr(0.);
            if(label_Syst.Contains("ee"))   sysErr = val_hChan_data*WjetsSyst_ee;
            else if(label_Syst.Contains("mumu")) sysErr = val_hChan_data*WjetsSyst_mm;
            h_Wjet->SetBinError(0,sysErr);

            shapeChan_SI.bckg.push_back(h_Wjet);
        }
    }
    cout << "\n#################### doWjetsBackground #######################\n" << endl;
}



void doQCDBackground(std::vector<TString>& signal_channels,map<TString, Shape_t>& allShapes, TString mainHisto)
{
    cout << "\n#################### doQCDBackground #######################\n" << endl;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t i=0; i<signal_channels.size(); i++) {


            cout << "QCD Shapes: " << signal_channels[i]+AnalysisBins[b]+mainHisto << endl;

            Shape_t& shapeChan_SI = allShapes.find(signal_channels[i]+AnalysisBins[b]+mainHisto)->second;
            TH1* hChan_data=shapeChan_SI.data;

            double valerr_hChan_data, val_hChan_data;
            val_hChan_data = hChan_data->IntegralAndError(1,hChan_data->GetXaxis()->GetNbins(),valerr_hChan_data);

            cout << "QCD Shapes: " << signal_channels[i]+AnalysisBins[b]+mainHisto << "\t"
                 << val_hChan_data << endl;

        }
    }
    cout << "\n#################### doQCDBackground #######################\n" << endl;
}








void dodataDrivenWWtW(std::vector<TString>& signal_channels,TString control_channel,map<TString, Shape_t>& allShapes, TString mainHisto, bool isMCclosureTest)
{

    if(isMCclosureTest) cout << "\n############ dodataDrivenWWtW (MCclosureTest) ###############\n" << endl;
    else 		cout << "\n#################### dodataDrivenWWtW #######################\n" << endl;
    FILE* pFile = NULL;
    if(isMCclosureTest) pFile = fopen("DataDrivenWWtW_MCclosureTest.tex","w");
    else 		pFile = fopen("DataDrivenWWtW.tex","w");
    fprintf(pFile,"\n\n\n");
    fprintf(pFile,"\\setlength{\\tabcolsep}{2pt}\n");
    fprintf(pFile,"\\begin{table}[ht!]\n");
    //fprintf(pFile,"\\begin{sidewaystable}[htp]\n");
    fprintf(pFile,"\\begin{center}\n");
    fprintf(pFile,"\\centering\n");
    fprintf(pFile,"\\scriptsize\\begin{tabular}{c|ccccc|cc|cc}\n");
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");
    if(isMCclosureTest) {
        fprintf(pFile,"Bin & $N_{ee}^{mc}$ & $N_{\\mu\\mu}^{mc}$ & $N_{e\\mu}^{mc}$ & $N_{e\\mu}^{mc,subtr}$ & $N_{e\\mu}^{corr}$ \n");
        fprintf(pFile,"& $N_{bkg,ee}^{closure}$ & $N_{bkg,ee}^{mc}$ & $N_{bkg,\\mu\\mu}^{closure}$ & $N_{bkg,\\mu\\mu}^{mc}$ \\\\\n");
    } else {
        fprintf(pFile,"Bin & $N_{ee}^{data}$ & $N_{\\mu\\mu}^{data}$ & $N_{e\\mu}^{data}$ & $N_{e\\mu}^{mc,subtr}$ & $N_{e\\mu}^{corr}$ \n");
        fprintf(pFile,"& $N_{bkg,ee}^{est}$ & $N_{bkg,ee}^{mc}$ & $N_{bkg,\\mu\\mu}^{est}$ & $N_{bkg,\\mu\\mu}^{mc}$ \\\\\n");
    }
    fprintf(pFile,"\\hline\n");

    for( auto const & bin : AnalysisBins ) {


        double N_data_ee(0.), N_data_mm(0.), N_data_ll(0.), N_data_em(0.);
        double N_mc_ee(0.), N_mc_mm(0.),N_mc_ll(0.);
        double N_mcsubtr_em(0.);

        double Err_data_ee(0.), Err_data_mm(0.),Err_data_ll(0.), Err_data_em(0.);
        double Err_mc_ee(0.), Err_mc_mm(0.),Err_mc_ll(0.);
        double Err_mcsubtr_em(0.);

        for(auto const & channel : signal_channels) {

            cout << control_channel+bin+mainHisto << " : " << channel+bin+mainHisto << endl;

            Shape_t& control_region_shapes = allShapes.find(control_channel+bin+mainHisto)->second;
            Shape_t& signal_region_shapes = allShapes.find(channel+bin+mainHisto)->second;

            TH1 *hCtrl_data(0), *hChan_data(0);
            if(isMCclosureTest) {
                hCtrl_data=control_region_shapes.totalBckg;
                hChan_data=signal_region_shapes.totalBckg;
            } else {
                hCtrl_data=control_region_shapes.data;
                hChan_data=signal_region_shapes.data;
            }


            TH1* hCtrl_MCnonNRB = (TH1*)control_region_shapes.totalBckg->Clone("hCtrl_MCnonNRB");
            TH1* hCtrl_MCNRB = (TH1*)control_region_shapes.totalBckg->Clone("hCtrl_MCNRB");
            TH1* hChan_MCnonNRB = (TH1*)signal_region_shapes.totalBckg->Clone("hChan_MCnonNRB");
            TH1* hChan_MCNRB = (TH1*)signal_region_shapes.totalBckg->Clone("hChan_MCNRB");
            hCtrl_MCnonNRB->Reset();
            hCtrl_MCNRB->Reset();
            hChan_MCnonNRB->Reset();
            hChan_MCNRB->Reset();


            //emu channel
            for(TH1* background_shape : control_region_shapes.bckg) {
                TString proc(background_shape->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    hCtrl_MCnonNRB->Add(background_shape, 1);
                }
                if(proc.Contains("t#bar{t}")
                        || proc.Contains("Single top")
                        || proc.Contains("WW#rightarrow l#nul#nu")
                        || proc.Contains("Z#rightarrow #tau#tau")
                  ) {
                    hCtrl_MCNRB->Add(background_shape, 1);
                }
            }

            //ee, mumu channel
            for(TH1* background_shape : signal_region_shapes.bckg) {
                TString proc(background_shape->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    hChan_MCnonNRB->Add(background_shape, 1);
                }
                if(proc.Contains("t#bar{t}")
                        || proc.Contains("Single top")
                        || proc.Contains("WW#rightarrow l#nul#nu")
                        || proc.Contains("Z#rightarrow #tau#tau")
                  ) {
                    hChan_MCNRB->Add(background_shape, 1);
                }
            }


            double valerr_hCtrl_data, val_hCtrl_data;
            val_hCtrl_data = hCtrl_data->IntegralAndError(1,hCtrl_data->GetXaxis()->GetNbins(),valerr_hCtrl_data);

            double valerr_hChan_data, val_hChan_data;
            val_hChan_data = hChan_data->IntegralAndError(1,hChan_data->GetXaxis()->GetNbins(),valerr_hChan_data);


            double valerr_hCtrl_MCnonNRB, val_hCtrl_MCnonNRB;
            double valerr_hCtrl_MCNRB, val_hCtrl_MCNRB;
            val_hCtrl_MCnonNRB = hCtrl_MCnonNRB->IntegralAndError(1,hCtrl_MCnonNRB->GetXaxis()->GetNbins(),valerr_hCtrl_MCnonNRB);
            val_hCtrl_MCNRB = hCtrl_MCNRB->IntegralAndError(1,hCtrl_MCNRB->GetXaxis()->GetNbins(),valerr_hCtrl_MCNRB);

            double valerr_hChan_MCnonNRB, val_hChan_MCnonNRB;
            double valerr_hChan_MCNRB, val_hChan_MCNRB;
            val_hChan_MCnonNRB = hChan_MCnonNRB->IntegralAndError(1,hChan_MCnonNRB->GetXaxis()->GetNbins(),valerr_hChan_MCnonNRB);
            val_hChan_MCNRB = hChan_MCNRB->IntegralAndError(1,hChan_MCNRB->GetXaxis()->GetNbins(),valerr_hChan_MCNRB);


            cout << "val_hCtrl_data: " << val_hCtrl_data << "\t"
                 << "val_hCtrl_MCnonNRB: " << val_hCtrl_MCnonNRB << "\t"
                 << "val_hCtrl_MCNRB:	" << val_hCtrl_MCNRB << endl;

            cout << "val_hChan_data: " << val_hChan_data << "\t"
                 << "val_hChan_MCnonNRB: " << val_hChan_MCnonNRB << "\t"
                 << "val_hChan_MCNRB:	" << val_hChan_MCNRB << endl;


            if(channel.Contains("ee")) {
                N_data_ee = val_hChan_data;
                N_mc_ee = val_hChan_MCNRB;
                Err_data_ee = valerr_hChan_data;
                Err_mc_ee = valerr_hChan_MCNRB;
            }
            if(channel.Contains("mumu")) {
                N_data_mm = val_hChan_data;
                N_mc_mm = val_hChan_MCNRB;
                Err_data_mm = valerr_hChan_data;
                Err_mc_mm = valerr_hChan_MCNRB;
            }
            if(channel.Contains("ll")) {
                N_data_ll = val_hChan_data;
                N_mc_ll = val_hChan_MCNRB;
                Err_data_ll = valerr_hChan_data;
                Err_mc_ll = valerr_hChan_MCNRB;
            }
            if(control_channel.Contains("emu")) {
                N_data_em = val_hCtrl_data;
                N_mcsubtr_em = val_hCtrl_MCnonNRB;
                Err_data_em = valerr_hCtrl_data;
                Err_mcsubtr_em = valerr_hCtrl_MCnonNRB;
            }


        } //ee, mumu, emu END


        cout << endl;
        cout << "N_data_ee: " << N_data_ee << "\t"
             << "N_data_mm: " << N_data_mm << "\t"
             << "N_data_ll: " << N_data_ll << "\t"
             << "N_data_em: " << N_data_em << "\t"
             << "N_mcsubtr_em: " << N_mcsubtr_em << endl;

        double k_ee = 0.5*sqrt(N_data_ee/N_data_mm);
        double k_mm = 0.5*sqrt(N_data_mm/N_data_ee);
        double N_data_corr = N_data_em-N_mcsubtr_em;
        double Err_data_corr = sqrt(pow(Err_data_em,2)+pow(Err_mcsubtr_em,2));
        double N_est_ee = N_data_corr*k_ee;
        double N_est_mm = N_data_corr*k_mm;
        double N_est_ll = N_est_ee + N_est_mm;

        double Err_est_ee = sqrt( pow(Err_data_corr,2)*N_data_ee/(4*N_data_mm)
                                  + pow(N_data_corr,2)*(pow(Err_data_ee/N_data_ee,2)+pow(Err_data_mm/N_data_mm,2))*(N_data_ee/(16*N_data_mm)) );

        double Err_est_mm = sqrt( pow(Err_data_corr,2)*N_data_mm/(4*N_data_ee)
                                  + pow(N_data_corr,2)*(pow(Err_data_ee/N_data_ee,2)+pow(Err_data_mm/N_data_mm,2))*(N_data_mm/(16*N_data_ee)) );

        double Err_est_ll = sqrt( pow(Err_data_corr,2)*N_data_mm/(4*N_data_ee)
                                  + pow(N_data_corr,2)*(pow(Err_data_ee/N_data_ee,2)+pow(Err_data_mm/N_data_mm,2))*(N_data_mm/(16*N_data_ee) + N_data_ee/(16*N_data_mm)) );

        cout << "N_est_ee: " << N_est_ee << "\t"
             << "N_est_mm: " << N_est_mm << "\t"
             << "N_est_ll: " << N_est_ll << endl;

        if(isMCclosureTest) {
            double sysEE = fabs(N_est_ee-N_mc_ee)/N_mc_ee;
            double sysMM = fabs(N_est_mm-N_mc_mm)/N_mc_mm;
            double sysLL = fabs(N_est_ll-N_mc_ll)/N_mc_ll;

            cout << "AnalysisBins: " << bin << " sysEE: " << sysEE << "\t" << "sysMM: " << sysMM << "\t" << "sysLL: " << sysLL << endl;

            if(bin == "eq0jets") {
                WWtopSyst_ee0jet = sysEE;
                WWtopSyst_mm0jet = sysMM;
                WWtopSyst_ll0jet = sysLL;
            }
            if(bin == "eq1jets") {
                WWtopSyst_ee1jet = sysEE;
                WWtopSyst_mm1jet = sysMM;
                WWtopSyst_ll1jet = sysLL;
            }
            if(bin == "lesq1jets") {
                WWtopSyst_eelesq1jet = sysEE;
                WWtopSyst_mmlesq1jet = sysMM;
                WWtopSyst_lllesq1jet = sysLL;
            }
        }





        if(isMCclosureTest) {
            fprintf(pFile,"%s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f  \\\\ \n",
                    bin.Data(),
                    N_data_ee,Err_data_ee,
                    N_data_mm,Err_data_mm,
                    N_data_em,Err_data_em,
                    N_mcsubtr_em,Err_mcsubtr_em,
                    N_data_corr,Err_data_corr,
                    N_est_ee,Err_est_ee, N_mc_ee,Err_mc_ee,
                    N_est_mm,Err_est_mm, N_mc_mm,Err_mc_mm
                   );
        } else {
            fprintf(pFile,"%s & %.0f $\\pm$ %.2f & %.0f $\\pm$ %.2f & %.0f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f  \\\\ \n",
                    bin.Data(),
                    N_data_ee,Err_data_ee,
                    N_data_mm,Err_data_mm,
                    N_data_em,Err_data_em,
                    N_mcsubtr_em,Err_mcsubtr_em,
                    N_data_corr,Err_data_corr,
                    N_est_ee,Err_est_ee, N_mc_ee,Err_mc_ee,
                    N_est_mm,Err_est_mm, N_mc_mm,Err_mc_mm
                   );
        }





        if(!isMCclosureTest) {
            //add data-driven backgrounds
            for(auto const & channel : signal_channels) {
                TString labelChan = channel+bin;
                Shape_t& signal_region_shapes = allShapes.find(channel+bin+mainHisto)->second;
                TH1* hChan_MCNRB = (TH1*)signal_region_shapes.totalBckg->Clone("hChan_MCNRB");
                hChan_MCNRB->Reset();

                //ee, mumu channel
                for(size_t ibckg=0; ibckg<signal_region_shapes.bckg.size(); ibckg++) {
                    TString proc(signal_region_shapes.bckg[ibckg]->GetTitle());
                    if(proc.Contains("t#bar{t}")
                            || proc.Contains("Single top")
                            || proc.Contains("WW#rightarrow l#nul#nu")
                            || proc.Contains("Z#rightarrow #tau#tau")
                      ) {
                        hChan_MCNRB->Add(signal_region_shapes.bckg[ibckg], 1);
                        //remove the separate parts
                        signal_region_shapes.bckg.erase(signal_region_shapes.bckg.begin()+ibckg);
                        ibckg--;

                    }
                }

                TH1* NonResonant = NULL;
                NonResonant = (TH1*)hChan_MCNRB->Clone("Top/WW/Ztautau (data)");
                NonResonant->SetTitle("Top/WW/Ztautau (data)");
                //add the background estimate

                double dataDrivenScale(1.0);
                if(channel.Contains("ee")) {
                    dataDrivenScale = N_est_ee/N_mc_ee;
                }
                if(channel.Contains("mumu")) {
                    dataDrivenScale = N_est_mm/N_mc_mm;
                }
                if(channel.Contains("ll")) {
                    dataDrivenScale = N_est_ll/N_mc_ll;
                }

                cout << "DATA/MC scale factors: " << dataDrivenScale << endl;

                for(int b=1; b<=NonResonant->GetXaxis()->GetNbins()+1; b++) {
                    double val = NonResonant->GetBinContent(b);
                    double newval = val*dataDrivenScale;
                    double err = NonResonant->GetBinError(b);
                    // Daniele's way (error on dataDrivenScale goes into the nirmalization error, later on)
                    double newerr = err*dataDrivenScale;

                    NonResonant->SetBinContent(b, newval );
                    if(newval!=0)NonResonant->SetBinError(b, newerr );
                }

                // Daniele's way: systematic uncertainty is from control region sample (possibly summed in quadrature with MC closure test)
                double Syst_WWTop(0.);
                if(labelChan.Contains("eeeq0jets")) 	Syst_WWTop = sqrt( pow(Err_est_ee/N_est_ee, 2) + pow(WWtopSyst_ee0jet, 2) );
                if(labelChan.Contains("mumueq0jets")) 	Syst_WWTop = sqrt( pow(Err_est_mm/N_est_mm, 2) + pow(WWtopSyst_mm0jet, 2) );
                if(labelChan.Contains("lleq0jets")) 	Syst_WWTop = sqrt( pow(Err_est_ll/N_est_ll, 2) + pow(WWtopSyst_ll0jet, 2) );
                if(labelChan.Contains("eeeq1jets")) 	Syst_WWTop = sqrt( pow(Err_est_ee/N_est_ee, 2) + pow(WWtopSyst_ee1jet, 2) );
                if(labelChan.Contains("mumueq1jets")) 	Syst_WWTop = sqrt( pow(Err_est_mm/N_est_mm, 2) + pow(WWtopSyst_mm1jet, 2) );
                if(labelChan.Contains("lleq1jets")) 	Syst_WWTop = sqrt( pow(Err_est_ll/N_est_ll, 2) + pow(WWtopSyst_ll1jet, 2) );
                if(labelChan.Contains("eelesq1jets")) 	Syst_WWTop = sqrt( pow(Err_est_ee/N_est_ee, 2) + pow(WWtopSyst_eelesq1jet, 2) );
                if(labelChan.Contains("mumulesq1jets")) Syst_WWTop = sqrt( pow(Err_est_mm/N_est_mm, 2) + pow(WWtopSyst_mmlesq1jet, 2) );
                if(labelChan.Contains("lllesq1jets"))   Syst_WWTop = sqrt( pow(Err_est_ll/N_est_ll, 2) + pow(WWtopSyst_lllesq1jet, 2) );


                cout << "labelchan: " << labelChan << " >>> Adding syst: " << Syst_WWTop << endl;
                NonResonant->SetBinError(0, Syst_WWTop*NonResonant->Integral());


                signal_region_shapes.bckg.push_back(NonResonant);
            }
        }

        cout << "\n\n";


    } //0jet, 1jet END


    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\end{tabular}\n");
    fprintf(pFile,"\\end{center}\n");
    //fprintf(pFile,"\\end{sidewaystable}\n");
    fprintf(pFile,"\\end{table}\n");
    fprintf(pFile,"\n\n\n");
    fclose(pFile);
    if(isMCclosureTest) cout << "############ dodataDrivenWWtW (MCclosureTest) ###############\n" << endl;
    else 		cout << "#################### dodataDrivenWWtW #######################\n" << endl;
}

//
void doBackgroundSubtraction(std::vector<TString>& signal_channels,TString control_channel,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto, TString url, JSONWrapper::Object &Root, bool isMCclosureTest)
{

    cout << "########################## doBackgroundSubtraction ##########################" << endl;
    string Lcol   = "\\begin{tabular}{|c";
    string Lchan  = "channel";
    string Lalph1 = "$\\alpha$ measured";
    string Lalph2 = "$\\alpha$ used";
    string Lyield   = "yield data";
    string LyieldMC = "yield mc";

    string Ccol   = "\\begin{tabular}{|c|c|c|c|c|c|";
    string Cname  = "channel & $\\alpha$ measured & $\\alpha$ used & yields $e\\mu$ & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("NonResonnant.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    for(size_t i=0; i<signal_channels.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Lcol += " |c";
            Lchan += string(" &")+signal_channels[i]+string(" - ")+AnalysisBins[b];
            Cval   = signal_channels[i]+string(" - ")+AnalysisBins[b];


            Shape_t& shapeCtrl_SB = allShapes.find(control_channel+AnalysisBins[b]+sideBandHisto)->second;
            TH1* hCtrl_SB=shapeCtrl_SB.data;

            Shape_t& shapeCtrl_SI = allShapes.find(control_channel+AnalysisBins[b]+mainHisto)->second;
            TH1* hCtrl_SI=shapeCtrl_SI.data;

            Shape_t& shapeChan_SB = allShapes.find(signal_channels[i]+AnalysisBins[b]+sideBandHisto)->second;
            TH1* hChan_SB=shapeChan_SB.data;

            Shape_t& shapeChan_SI = allShapes.find(signal_channels[i]+AnalysisBins[b]+mainHisto)->second;


            if(isMCclosureTest) {
                hCtrl_SB = shapeCtrl_SB.totalBckg;
                hCtrl_SI = shapeCtrl_SI.totalBckg;
                hChan_SB = shapeChan_SB.totalBckg;
            }

            cout << ">>>>>>>>>>>>>>>>> 2" << endl;

            //to subtract non WW/Top background from data
            TH1* tosubtrChan_SB = (TH1*)shapeChan_SB.totalBckg->Clone("tosubtrChan_SB");
            TH1* tosubtrCtrl_SB = (TH1*)shapeCtrl_SB.totalBckg->Clone("tosubtrCtrl_SB");
            tosubtrChan_SB->Reset();
            tosubtrCtrl_SB->Reset();


            for(size_t ibckg=0; ibckg<shapeChan_SB.bckg.size(); ibckg++) {
                TString proc(shapeChan_SB.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    tosubtrChan_SB->Add(shapeChan_SB.bckg[ibckg], 1);
                }
            }
            for(size_t ibckg=0; ibckg<shapeCtrl_SB.bckg.size(); ibckg++) {
                TString proc(shapeCtrl_SB.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    tosubtrCtrl_SB->Add(shapeCtrl_SB.bckg[ibckg], 1);
                }
            }

            double alpha=0 ,alpha_err=0;


            // >= 1 b-jets
            if(hCtrl_SB->GetBinContent(5)>0) {
                //new for monoZ analysis
                alpha = (hChan_SB->GetBinContent(5)-tosubtrChan_SB->GetBinContent(5))/(hCtrl_SB->GetBinContent(5)-tosubtrCtrl_SB->GetBinContent(5));
                double err1 = sqrt(pow(hChan_SB->GetBinError(5),2)+pow(tosubtrChan_SB->GetBinError(5),2));
                double err2 = sqrt(pow(hCtrl_SB->GetBinError(5),2)+pow(tosubtrCtrl_SB->GetBinError(5),2));
                alpha_err = alpha * sqrt(pow(err1/(hChan_SB->GetBinContent(5)-tosubtrChan_SB->GetBinContent(5)),2)
                                         + pow(err2/(hCtrl_SB->GetBinContent(5)-tosubtrCtrl_SB->GetBinContent(5)),2));
                cout << "alpha: " << alpha << " +/- " << alpha_err << endl;
            }

            // =0 b-jets
            int BinPosition = 2;
            double alpha_bveto=0 ,alpha_bveto_err=0;
            if(hCtrl_SB->GetBinContent(BinPosition)>0) {
                alpha_bveto = (hChan_SB->GetBinContent(BinPosition)-tosubtrChan_SB->GetBinContent(BinPosition))
                              /(hCtrl_SB->GetBinContent(BinPosition)-tosubtrCtrl_SB->GetBinContent(BinPosition));
                double err1 = sqrt(pow(hChan_SB->GetBinError(BinPosition),2)+pow(tosubtrChan_SB->GetBinError(BinPosition),2));
                double err2 = sqrt(pow(hCtrl_SB->GetBinError(BinPosition),2)+pow(tosubtrCtrl_SB->GetBinError(BinPosition),2));
                alpha_bveto_err = alpha_bveto * sqrt(pow(err1/(hChan_SB->GetBinContent(BinPosition)-tosubtrChan_SB->GetBinContent(BinPosition)),2)
                                                     + pow(err2/(hCtrl_SB->GetBinContent(BinPosition)-tosubtrCtrl_SB->GetBinContent(BinPosition)),2));
                cout << "alpha_bveto: " << alpha_bveto << " +/- " << alpha_bveto_err << endl;
            }

            cout << "<--------- Alpha systematics --------->" << endl;
            double syst_alpha = fabs(alpha_bveto-alpha)/alpha;
            double syst_alpha_err = pow(alpha_bveto_err/alpha_bveto,2);
            syst_alpha_err += pow(alpha_err/alpha,2);
            syst_alpha_err = sqrt(syst_alpha_err);
            cout << "syst_alpha: " << syst_alpha*100. << " +/- " << syst_alpha*100.*syst_alpha_err << " (%)"<< endl;
            cout << "<--------- Alpha systematics --------->" << endl;



            Lalph1 += string(" &") + toLatexRounded(alpha,alpha_err);
            Cval   += string(" &") + toLatexRounded(alpha,alpha_err);

            Lalph2 += string(" &") + toLatexRounded(alpha,alpha_err);
            Cval   += string(" &") + toLatexRounded(alpha,alpha_err);


            TH1* NonResonant = NULL;
            if(subNRB2011) {
                NonResonant = (TH1*)hCtrl_SI->Clone("Top/WW/W+Jets (data)");
                NonResonant->SetTitle("Top/WW/W+Jets (data)");
            } else if(subNRB2012) {
                Shape_t& shapeChan_BTag = allShapes.find(signal_channels[i]+mainHisto+"BTagSB")->second;
                NonResonant = (TH1*)shapeChan_BTag.data->Clone("Top (data)");
                NonResonant->SetTitle("Top (data)");
            } else {
                return;
            }



            TH1* emu_mcnonNRB = (TH1*)shapeCtrl_SI.totalBckg->Clone("emu_mcnonNRB");
            emu_mcnonNRB->Reset();

            for(size_t ibckg=0; ibckg<shapeCtrl_SI.bckg.size(); ibckg++) {
                TString proc(shapeCtrl_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu") ||
                        proc.Contains("WZ#rightarrow 3l#nu") ||
                        proc.Contains("Z+jets") ||
                        proc.Contains("W+jets")
                  ) {
                    emu_mcnonNRB->Add(shapeCtrl_SI.bckg[ibckg], 1);
                }
            }


            double valerr_MCnonNRB_emu;
            double val_MCnonNRB_emu = emu_mcnonNRB->IntegralAndError(1,emu_mcnonNRB->GetXaxis()->GetNbins(),valerr_MCnonNRB_emu);
            cout << "Total MCnonNRB in EMU: " << val_MCnonNRB_emu << endl;


            double valvalerr, valval;
            valval = NonResonant->IntegralAndError(1,NonResonant->GetXaxis()->GetNbins(),valvalerr);

            Cval   += string(" &") + toLatexRounded(valval,valvalerr);

            cout << control_channel+AnalysisBins[b]+mainHisto << ": " << valval << endl;
            for(int b=1; b<=NonResonant->GetXaxis()->GetNbins()+1; b++) {
                double val = NonResonant->GetBinContent(b);
                double subtr = emu_mcnonNRB->GetBinContent(b);
                double err = NonResonant->GetBinError(b);
                double newval = (val-subtr)*alpha;
                double newerr = sqrt(pow(err*alpha,2) + pow(val*alpha_err,2));
                NonResonant->SetBinContent(b, newval );
                NonResonant->SetBinError  (b, newerr );
            }
            NonResonant->Scale(DDRescale);

            Double_t valerr;
            Double_t val = NonResonant->IntegralAndError(1,NonResonant->GetXaxis()->GetNbins(),valerr);
            Double_t systError = val*NonResonnantSyst;
            NonResonant->SetBinError(0,systError);//save syst error in underflow bin that is always empty
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
                systError=-1;
            }
            Lyield += string(" &") + toLatexRounded(val,valerr,systError);
            Cval   += string(" &") + toLatexRounded(val,valerr,systError);

            //Clean background collection
            TH1* MCNRB = (TH1*)shapeChan_SI.totalBckg->Clone("MCNRB");
            MCNRB->Reset();

            TH1* MCnonNRB = (TH1*)shapeChan_SI.totalBckg->Clone("MCnonNRB");
            MCnonNRB->Reset();

            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                cout << "proc: " << proc << endl;
                if(proc.Contains("t#bar{t}")
                        || proc.Contains("Single top")
                        || proc.Contains("WW#rightarrow l#nul#nu")
                        || proc.Contains("Z#rightarrow #tau#tau")
                  ) {
                    MCNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                    NonResonant->SetFillColor(shapeChan_SI.bckg[ibckg]->GetFillColor());
                    NonResonant->SetLineColor(shapeChan_SI.bckg[ibckg]->GetLineColor());
                    shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                    ibckg--;
                } else {
                    MCnonNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                }
            }


            NonResonant->SetFillStyle(1001);
            NonResonant->SetFillColor(393);//592);
            NonResonant->SetLineColor(393);//592);


            //add the background estimate
            shapeChan_SI.bckg.push_back(NonResonant);


            //recompute total background
            shapeChan_SI.totalBckg->Reset();
            for(size_t i=0; i<shapeChan_SI.bckg.size(); i++) {
                shapeChan_SI.totalBckg->Add(shapeChan_SI.bckg[i]);
            }

            val = MCNRB->IntegralAndError(1,MCNRB->GetXaxis()->GetNbins(),valerr);
            cout << "MCNRB: " << val << endl;
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
            }
            LyieldMC += string(" &") + toLatexRounded(val,valerr);
            Cval     += string(" &") + toLatexRounded(val,valerr);

            if(pFile) {
                fprintf(pFile,"%s\\\\\n", Cval.c_str());
            }
        }
    }

    if(pFile) {
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fprintf(pFile,"\n\n\n\n");

        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s|}\\hline\n", Lcol.c_str());
        fprintf(pFile,"%s\\\\\n", Lchan.c_str());
        fprintf(pFile,"%s\\\\\n", Lalph1.c_str());
        fprintf(pFile,"%s\\\\\n", Lalph2.c_str());
        fprintf(pFile,"%s\\\\\n", Lyield.c_str());
        fprintf(pFile,"%s\\\\\n", LyieldMC.c_str());
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fclose(pFile);
    }


    cout << "########################## doBackgroundSubtraction ##########################" << endl;
}




void doWZSubtraction(std::vector<TString>& signal_channels,TString control_channel,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto)
{
    string Ccol   = "\\begin{tabular}{|l|c|c|c|c|";
    string Cname  = "channel & $\\alpha$ measured & $\\alpha$ used & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("WZ.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    for(size_t i=0; i<signal_channels.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Cval   = signal_channels[i]+string(" - ")+AnalysisBins[b];

            Shape_t& shapeCtrl_SB = allShapes.find(control_channel+AnalysisBins[b]+sideBandHisto)->second;
            Shape_t& shapeCtrl_SI = allShapes.find(control_channel+AnalysisBins[b]+mainHisto)->second;
            Shape_t& shapeChan_SB = allShapes.find(signal_channels[i]+AnalysisBins[b]+sideBandHisto)->second;
            Shape_t& shapeChan_SI = allShapes.find(signal_channels[i]+AnalysisBins[b]+mainHisto)->second;

            fprintf(pFile,"#############%s:\n",(string(" &")+signal_channels[i]+string(" - ")+AnalysisBins[b]).Data());
            fprintf(pFile,"MC: em 3leptons=%6.2E  em 2leptons=%6.2E  ll 3leptons=%6.2E  ll 2leptons=%6.2E\n",shapeCtrl_SB.totalBckg->Integral(), shapeCtrl_SI.totalBckg->Integral(), shapeChan_SB.totalBckg->Integral(), shapeChan_SI.totalBckg->Integral());

            TH1* histo1=NULL, *histo2=NULL, *histo3=NULL, *histo4=NULL;
            for(size_t ibckg=0; ibckg<shapeCtrl_SB.bckg.size(); ibckg++) {
                if(TString(shapeCtrl_SB.bckg[ibckg]->GetTitle()).Contains("WZ"))histo1=shapeCtrl_SB.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeCtrl_SI.bckg.size(); ibckg++) {
                if(TString(shapeCtrl_SI.bckg[ibckg]->GetTitle()).Contains("WZ"))histo2=shapeCtrl_SI.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeChan_SB.bckg.size(); ibckg++) {
                if(TString(shapeChan_SB.bckg[ibckg]->GetTitle()).Contains("WZ"))histo3=shapeChan_SB.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                if(TString(shapeChan_SI.bckg[ibckg]->GetTitle()).Contains("WZ"))histo4=shapeChan_SI.bckg[ibckg];
            }
            fprintf(pFile,"WZ: em 3leptons=%6.2E  em 2leptons=%6.2E  ll 3leptons=%6.2E  ll 2leptons=%6.2E\n",histo1->Integral(), histo2->Integral(), histo3->Integral(), histo4->Integral());

            double Num, Denom, NumError, DenomError;
            Num = histo4->IntegralAndError(1,histo4->GetXaxis()->GetNbins(),NumError);
            Denom = histo3->IntegralAndError(1,histo3->GetXaxis()->GetNbins(),DenomError);
            double ratio = Num/Denom;
            double ratio_err  = sqrt(pow(Num*DenomError,2) + pow(Denom*NumError,2))/ pow(Denom,2);
            double ratio_syst = fabs(histo3->Integral() - shapeChan_SB.totalBckg->Integral())/shapeChan_SB.totalBckg->Integral();
            fprintf(pFile,"Ratio = %s\n",toLatexRounded(ratio,ratio_err,ratio_syst).c_str());

            Double_t valerr;
            Double_t val = histo4->IntegralAndError(1,histo4->GetXaxis()->GetNbins(),valerr);

            Double_t valerr2, valsyst2;
            Double_t val2 = shapeChan_SB.totalBckg->IntegralAndError(1,shapeChan_SB.data->GetXaxis()->GetNbins(),valerr2);
            valerr2= sqrt(pow(valerr2*ratio,2) + pow(val2*ratio_err,2) );
            valsyst2 = val2*ratio_syst;
            val2=val2*ratio;
            fprintf(pFile,"WZ (MC closure test): %s --> %s\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val2,valerr2,valsyst2).c_str());

            bool noDataObserved=false;
            Double_t valerr3, valsyst3;
            Double_t val3 = shapeChan_SB.data->IntegralAndError(1,shapeChan_SB.data->GetXaxis()->GetNbins(),valerr3);
            if(val3<=0) {
                noDataObserved=true;
                val3=1.0;
                valerr3=1.0;
            }
            valerr3= sqrt(pow(valerr3*ratio,2) + pow(val3*ratio_err,2) );
            valsyst3 = val3*ratio_syst;
            val3=val3*ratio;
            if(!noDataObserved) {
                fprintf(pFile,"WZ (from data)      : %s --> %s\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val3,valerr3,valsyst3).c_str());
            } else {
                fprintf(pFile,"WZ (from data)      : %s --> smaller than %s (because no data was observed in 3dlepton SideBand--> assume 1+-1 observed data for rescale)\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val3,valerr3,valsyst3).c_str());
            }


        }
    }

    if(pFile) {
        fclose(pFile);
    }
}

void BlindData(std::vector<TString>& signal_channels, map<TString, Shape_t>& allShapes, TString mainHisto, bool addSignal)
{
    for(auto const & channel : signal_channels) {
        for(auto const & bin : AnalysisBins) {
            Shape_t& shapeChan_SI = allShapes.find(channel+bin+mainHisto)->second;
            shapeChan_SI.data->Reset();
            shapeChan_SI.data->Add(shapeChan_SI.totalBckg,1);
            if(addSignal) {
                for(auto const & signal_shape : shapeChan_SI.signal) {
                    shapeChan_SI.data->Add(signal_shape, 1);
                }
            }
        }
    }
}




