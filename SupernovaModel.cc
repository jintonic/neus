#include <TH1D.h>
#include <TH2D.h>

class SupernovaModel : public TObject
{
   protected:
      TH1D *fH1Empty;
      TH2D *fH2Empty;
      TH1D *fH1N1, *fH1N2, *fH1N3, *fH1E1, *fH1E2, *fH1E3;
      TH2D *fH2N1, *fH2N2, *fH2N3, *fH2E1, *fH2E2, *fH2E3;

      Float_t fInitialMass; ///< progenitor mass in Solar mass unit
      Float_t fMetallicity;
      Float_t fReviveTime; ///< shock revival time in ms

   public:
      SupernovaModel(
            Float_t initialMass=13, /* Solar mass */
            Float_t metallicity=0.02,
            Float_t reviveTime=100 /* ms */) : 
         TObject(), fH1Empty(0), fH2Empty(0),
         fH1N1(0), fH1N2(0), fH1N3(0), fH1E1(0), fH1E2(0), fH1E3(0), 
         fH2N1(0), fH2N2(0), fH2N3(0), fH2E1(0), fH2E2(0), fH2E3(0), 
         fInitialMass(initialMass), 
         fMetallicity(metallicity),
         fReviveTime(reviveTime) {};

      virtual ~SupernovaModel();

      //void LoadFullData(const char *databaseDir="./intpdata");
      void LoadIntegratedData(const char *databaseDir="./integdata");

      //TH2D* NumberSpectrum(const char *neutrino="#nu_{e}");
      TH1D* IntegratedNumberSpectrum(const char *neutrino="#nu_{e}");

      //TH2D* EnergySpectrum(const char *neutrino="#nu_{e}");
      TH1D* IntegratedEnergySpectrum(const char *neutrino="#nu_{e}");

      ClassDef(SupernovaModel,1);
};

//______________________________________________________________________________
//

SupernovaModel::~SupernovaModel()
{
   if (fH1Empty) delete fH1Empty;
   if (fH2Empty) delete fH2Empty;

   if (fH1N1) delete fH1N1;
   if (fH1N2) delete fH1N2;
   if (fH1N3) delete fH1N3;
   if (fH1E1) delete fH1E1;
   if (fH1E2) delete fH1E2;
   if (fH1E3) delete fH1E3;

   if (fH2N1) delete fH2N1;
   if (fH2N2) delete fH2N2;
   if (fH2N3) delete fH2N3;
   if (fH2E1) delete fH2E1;
   if (fH2E2) delete fH2E2;
   if (fH2E3) delete fH2E3;
}

//______________________________________________________________________________
//

#include <fstream>
using namespace std;

void SupernovaModel::LoadIntegratedData(const char *databaseDir)
{
   char *name = Form("%s/integ%.0f0%.0f.data",
         databaseDir, fInitialMass, fReviveTime/100);
   if (fMetallicity==0.004) 
      name = Form("%s/integ%.0f1%.0f.data",
            databaseDir, fInitialMass, fReviveTime/100);

   // check database
   ifstream file(name); 
   if (!(file.is_open())) {
      Warning("LoadData", "%s cannot be read!", name); 
      fH1Empty = new TH1D("h1Empty","Empty spectrum",0,0.,0.);
      fH2Empty = new TH2D("h2Empty","Empty spectrum",0,0.,0.,0,0.,0.);
      return;
   }

   // skip the first line
   char line[150];
   file.getline(line,150);

   // load data
   const UShort_t nbins = 20;
   Double_t binEdges[nbins+1]={0};
   Double_t number1[nbins];
   Double_t number2[nbins];
   Double_t number3[nbins];
   Double_t energy1[nbins];
   Double_t energy2[nbins];
   Double_t energy3[nbins];

   Double_t energy, n1, n2, n3, e1, e2, e3;

   UShort_t i=0;
   while(file>>energy>>energy>>n1>>n2>>n3>>e1>>e2>>e3) {
      number1[i]=n1; number2[i]=n2; number3[i]=n3; 
      energy1[i]=e1; energy2[i]=e2; energy3[i]=e3;
      i++;
      binEdges[i]=energy; 
   }

   file.close();

   // fill spectra
   fH1N1 = new TH1D("h1N1",
         ";Energy [MeV];Number luminosity/(1 MeV)",nbins,binEdges);
   fH1N2 = new TH1D("h1N2",
         ";Energy [MeV];Number luminosity/(1 MeV)",nbins,binEdges);
   fH1N3 = new TH1D("h1N3",
         ";Energy [MeV];Number luminosity/(1 MeV)",nbins,binEdges);

   fH1E1 = new TH1D("h1E1",
         ";Energy [MeV];Energy luminosity [erg]/(1 MeV)",nbins,binEdges);
   fH1E2 = new TH1D("h1E2",
         ";Energy [MeV];Energy luminosity [erg]/(1 MeV)",nbins,binEdges);
   fH1E3 = new TH1D("h1E3",
         ";Energy [MeV];Energy luminosity [erg]/(1 MeV)",nbins,binEdges);

   for (i=1; i<=nbins; i++) {
      fH1N1->SetBinContent(i,number1[i-1]);
      fH1N2->SetBinContent(i,number2[i-1]);
      fH1N3->SetBinContent(i,number3[i-1]);

      fH1E1->SetBinContent(i,energy1[i-1]);
      fH1E2->SetBinContent(i,energy2[i-1]);
      fH1E3->SetBinContent(i,energy3[i-1]);
   }

   fH1N1->SetStats(0);
   fH1N2->SetStats(0);
   fH1N3->SetStats(0);

   fH1E1->SetStats(0);
   fH1E2->SetStats(0);
   fH1E3->SetStats(0);

   fH1N1->SetLineColor(kBlack);
   fH1N2->SetLineColor(kRed);
   fH1N3->SetLineColor(kBlue);

   fH1E1->SetLineColor(kBlack);
   fH1E2->SetLineColor(kRed);
   fH1E3->SetLineColor(kBlue);

   fH1N1->SetMarkerColor(kBlack);
   fH1N2->SetMarkerColor(kRed);
   fH1N3->SetMarkerColor(kBlue);

   fH1E1->SetMarkerColor(kBlack);
   fH1E2->SetMarkerColor(kRed);
   fH1E3->SetMarkerColor(kBlue);

   fH1N1->SetMarkerStyle(20);
   fH1N2->SetMarkerStyle(20);
   fH1N3->SetMarkerStyle(20);

   fH1E1->SetMarkerStyle(20);
   fH1E2->SetMarkerStyle(20);
   fH1E3->SetMarkerStyle(20);

   fH1N1->SetMarkerSize(0.6);
   fH1N2->SetMarkerSize(0.6);
   fH1N3->SetMarkerSize(0.6);

   fH1E1->SetMarkerSize(0.6);
   fH1E2->SetMarkerSize(0.6);
   fH1E3->SetMarkerSize(0.6);
}

//______________________________________________________________________________
//

TH1D* SupernovaModel::IntegratedNumberSpectrum(const char *neutrino)
{
   if (!fH1N1) return fH1Empty;

   TString species(neutrino);
   if (species=="#nu_{e}") {
      fH1N1->SetTitle(Form("%s",neutrino));
      return fH1N1;
   } else if (species=="#bar{#nu}_{e}") {
      fH1N2->SetTitle(Form("%s",neutrino));
      return fH1N2;
   } else if (species=="#nu_{x}") {
      fH1N3->SetTitle(Form("%s",neutrino));
      return fH1N3;
   } else {
      Warning("IntegratedNumberSpectrum",
            "Neutrino species %s is not defined!",neutrino);
      Warning("IntegratedNumberSpectrum",
            "Please select one from #nu_{e}, #bar{#nu}_{e} or #nu_{x}.");
      Warning("IntegratedNumberSpectrum",
            "The spectrum of #nu_{e} is returned.");
      fH1N1->SetTitle("#nu_{e}");
      return fH1N1;
   }
}

//______________________________________________________________________________
//

TH1D* SupernovaModel::IntegratedEnergySpectrum(const char *neutrino)
{
   if (!fH1E1) return fH1Empty;

   TString species(neutrino);
   if (species=="#nu_{e}") {
      fH1E1->SetTitle(Form("%s",neutrino));
      return fH1E1;
   } else if (species=="#bar{#nu}_{e}") {
      fH1E2->SetTitle(Form("%s",neutrino));
      return fH1E2;
   } else if (species=="#nu_{x}") {
      fH1E3->SetTitle(Form("%s",neutrino));
      return fH1E3;
   } else {
      Warning("IntegratedEnergySpectrum",
            "Neutrino species %s is not defined!",neutrino);
      Warning("IntegratedEnergySpectrum",
            "Please select one from #nu_{e}, #bar{#nu}_{e} or #nu_{x}.");
      Warning("IntegratedEnergySpectrum",
            "The spectrum of #nu_{e} is returned.");
      fH1E1->SetTitle("#nu_{e}");
      return fH1E1;
   }
}
