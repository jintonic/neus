#include <TGraph.h>

class SupernovaModel : public TObject
{
   protected:
      TGraph *fGN1, *fGN2, *fGN3, *fGE1, *fGE2, *fGE3;

      Float_t fInitialMass; ///< progenitor mass in Solar mass unit
      Float_t fMetallicity;
      Float_t fReviveTime; ///< shock revival time in ms

   public:
      SupernovaModel(
            Float_t initialMass=13, /* Solar mass */
            Float_t metallicity=0.02,
            Float_t reviveTime=100, /* ms */
            const char *databaseDir="./integdata");
      virtual ~SupernovaModel();

      void LoadData(
            Float_t initialMass=13, /* Solar mass */
            Float_t metallicity=0.02,
            Float_t reviveTime=100, /* ms */
            const char *databaseDir="./integdata");

      TGraph* NumberSpectrum(const char *neutrino="#bar{#nu}_{e}");

      ClassDef(SupernovaModel,0);
};

//______________________________________________________________________________
//

SupernovaModel::SupernovaModel(
      Float_t initialMass, Float_t metallicity, Float_t reviveTime,
      const char *databaseDir) : TObject(), 
   fGN1(0), fGN2(0), fGN3(0), fGE1(0), fGE2(0), fGE3(0),
   fInitialMass(initialMass), fMetallicity(metallicity), fReviveTime(reviveTime)
{ 
   LoadData(fInitialMass,fMetallicity,fReviveTime,databaseDir); 
}

//______________________________________________________________________________
//

SupernovaModel::~SupernovaModel()
{
   if (fGN1) delete fGN1;
   if (fGN2) delete fGN2;
   if (fGN3) delete fGN3;
   if (fGE1) delete fGE1;
   if (fGE2) delete fGE2;
   if (fGE3) delete fGE3;
}

//______________________________________________________________________________
//

#include <fstream>
#include <string>
using namespace std;

void SupernovaModel::LoadData(
      Float_t initialMass, Float_t metallicity, Float_t reviveTime,
      const char *databaseDir)
{
   const char *name = Form("%s/integ%.0f%.0f%.0f.data",databaseDir,
         initialMass, metallicity*100/2, reviveTime/100);

   ifstream file(name); 
   if (!(file.is_open()))
      Fatal("LoadData", "%s cannot be read!", name); 

   char line[150];
   file.getline(line,150);

   const UShort_t fNp = 20;
   Double_t fEn[fNp];
   Double_t fN1[fNp];
   Double_t fN2[fNp];
   Double_t fN3[fNp];
   Double_t fE1[fNp];
   Double_t fE2[fNp];
   Double_t fE3[fNp];

   // load data
   Double_t energy, n1, n2, n3, e1, e2, e3;

   UShort_t i=0;
   while(file>>energy>>energy>>n1>>n2>>n3>>e1>>e2>>e3) {
      Printf("%e \t %e \t %e \t %e \t %e \t %e \t %e",
            energy, n1, n2, n3, e1, e2, e3);
      fEn[i]=energy; 
      fN1[i]=n1; fN2[i]=n2; fN3[i]=n3; 
      fE1[i]=e1; fE2[i]=e2; fE3[i]=e3;
      i++;
   }

   file.close();

   fGN1 = new TGraph(fNp,fEn,fN1);
   fGN2 = new TGraph(fNp,fEn,fN2);
   fGN3 = new TGraph(fNp,fEn,fN3);

   fGE1 = new TGraph(fNp,fEn,fE1);
   fGE2 = new TGraph(fNp,fEn,fE2);
   fGE3 = new TGraph(fNp,fEn,fE3);

   fGN1->SetLineColor(kBlack);
   fGN2->SetLineColor(kRed);
   fGN3->SetLineColor(kBlue);

   fGE1->SetLineColor(kBlack);
   fGE2->SetLineColor(kRed);
   fGE3->SetLineColor(kBlue);

   fGN1->SetMarkerColor(kBlack);
   fGN2->SetMarkerColor(kRed);
   fGN3->SetMarkerColor(kBlue);

   fGE1->SetMarkerColor(kBlack);
   fGE2->SetMarkerColor(kRed);
   fGE3->SetMarkerColor(kBlue);

   fGN1->SetMarkerStyle(20);
   fGN2->SetMarkerStyle(20);
   fGN3->SetMarkerStyle(20);

   fGE1->SetMarkerStyle(20);
   fGE2->SetMarkerStyle(20);
   fGE3->SetMarkerStyle(20);

   fGN1->SetMarkerSize(0.6);
   fGN2->SetMarkerSize(0.6);
   fGN3->SetMarkerSize(0.6);

   fGE1->SetMarkerSize(0.6);
   fGE2->SetMarkerSize(0.6);
   fGE3->SetMarkerSize(0.6);
}

//______________________________________________________________________________
//

TGraph* SupernovaModel::NumberSpectrum(const char *neutrino)
{
   TString species(neutrino);
   if (species=="#nu_{e}") {
      fGN1->SetTitle(
            Form(";Energy of %s [MeV];number of %s / MeV",neutrino,neutrino));
      return fGN1;
   } else if (species=="#bar{#nu}_{e}") {
      fGN2->SetTitle(
            Form(";Energy of %s [MeV];number of %s / MeV",neutrino,neutrino));
      return fGN2;
   } else if (species=="#nu_{x}") {
      fGN3->SetTitle(
            Form(";Energy of %s [MeV];number of %s / MeV",neutrino,neutrino));
      return fGN3;
   } else {
      Warning("NumberSpectrum","%s is not defined! Select one from",neutrino);
      Warning("NumberSpectrum", "#nu_{e}, #bar{#nu}_{e}, #nu_{x}");
      Warning("NumberSpectrum", "number spectrum of #nu_{e} is returned.");
      fGN1->SetTitle(";Energy of #nu_{e} [MeV];number of #nu_{e} / MeV");
      return fGN1;
   }
}
