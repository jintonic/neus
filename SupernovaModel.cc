#include "SupernovaModel.h"

SupernovaModel::SupernovaModel(
      Float_t initialMass, Float_t metallicity, Float_t reviveTime) :
   TNamed(), fH1Empty(0), fH2Empty(0),
   fH1N1(0), fH1N2(0), fH1N3(0), fH1E1(0), fH1E2(0), fH1E3(0),
   fH2N1(0), fH2N2(0), fH2N3(0), fH2E1(0), fH2E2(0), fH2E3(0),
   fInitialMass(initialMass), fMetallicity(metallicity), fReviveTime(reviveTime)
{
   if (metallicity>0.02-0.01) {
      fName = Form("model%.0f0%.0f",fInitialMass,fReviveTime/100);
      fTitle = Form("%.0f Solar mass, 0.02, %.0f ms",fInitialMass,fReviveTime);
   } else {
      fName = Form("model%.0f1%.0f",fInitialMass,fReviveTime/100);
      fTitle = Form("%.0f Solar mass, 0.004, %.0f ms",fInitialMass,fReviveTime);
   }
}

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
      Warning("LoadIntegratedData", "%s cannot be read!", name);
      fH1Empty = new TH1D(Form("h1Empty%.0f%.0f%.0f",
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            "Empty spectrum",0,0.,0.);
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
   fH1N1 = new TH1D(Form("h1N1%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Number luminosity/(1 MeV)",nbins,binEdges);
   fH1N2 = new TH1D(Form("h1N2%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Number luminosity/(1 MeV)",nbins,binEdges);
   fH1N3 = new TH1D(Form("h1N3%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Number luminosity/(1 MeV)",nbins,binEdges);

   fH1E1 = new TH1D(Form("h1E1%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Energy luminosity [erg]/(1 MeV)",nbins,binEdges);
   fH1E2 = new TH1D(Form("h1E2%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Energy luminosity [erg]/(1 MeV)",nbins,binEdges);
   fH1E3 = new TH1D(Form("h1E3%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Energy luminosity [erg]/(1 MeV)",nbins,binEdges);

   for (i=1; i<=nbins; i++) {
      fH1N1->SetBinContent(i,number1[i-1]);
      fH1N2->SetBinContent(i,number2[i-1]);
      fH1N3->SetBinContent(i,number3[i-1]);

      fH1E1->SetBinContent(i,energy1[i-1]);
      fH1E2->SetBinContent(i,energy2[i-1]);
      fH1E3->SetBinContent(i,energy3[i-1]);
   }

   // set properties
   fH1N1->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH1N2->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH1N3->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fH1E1->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH1E2->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH1E3->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

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
   if (species=="v_e") return fH1N1;
   else if (species=="anti-v_e") return fH1N2;
   else if (species=="v_x") return fH1N3;
   else {
      Warning("IntegratedNumberSpectrum",
            "Neutrino species %s is not defined!",neutrino);
      Warning("IntegratedNumberSpectrum",
            "Please select one from v_e, anti-v_e or v_x.");
      Warning("IntegratedNumberSpectrum",
            "The spectrum of #nu_{e} is returned.");
      return fH1N1;
   }
}

//______________________________________________________________________________
//

TH1D* SupernovaModel::IntegratedEnergySpectrum(const char *neutrino)
{
   if (!fH1E1) return fH1Empty;

   TString species(neutrino);
   if (species=="v_e") return fH1E1;
   else if (species=="anti-v_e") return fH1E2;
   else if (species=="v_x") return fH1E3;
   else {
      Warning("IntegratedEnergySpectrum",
            "Neutrino species %s is not defined!",neutrino);
      Warning("IntegratedEnergySpectrum",
            "Please select one from v_e, anti-v_e or v_x.");
      Warning("IntegratedEnergySpectrum",
            "The spectrum of #nu_{e} is returned.");
      return fH1E1;
   }
}

//______________________________________________________________________________
//

void SupernovaModel::LoadFullData(const char *databaseDir)
{
   char *name = Form("%s/intp%.0f0%.0f.data",
         databaseDir, fInitialMass, fReviveTime/100);
   if (fMetallicity==0.004)
      name = Form("%s/intp%.0f1%.0f.data",
            databaseDir, fInitialMass, fReviveTime/100);

   // check database
   ifstream file(name);
   if (!(file.is_open())) {
      Warning("LoadFullData", "%s cannot be read!", name);
      fH2Empty = new TH2D(Form("h2Empty%.0f%.0f%.0f",
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            "Empty spectrum",0,0.,0.,0,0.,0.);
      return;
   }

   // load data
   const UShort_t nbinx = 20, nbiny = 391;
   Double_t binEdgesx[nbinx+1]={0};
   Double_t binEdgesy[nbiny+1]={0};
   Double_t number1[nbinx][nbiny];
   Double_t number2[nbinx][nbiny];
   Double_t number3[nbinx][nbiny];
   Double_t energy1[nbinx][nbiny];
   Double_t energy2[nbinx][nbiny];
   Double_t energy3[nbinx][nbiny];

   Double_t time, energy, n1, n2, n3, e1, e2, e3;

   UShort_t ix=0, iy=0;
   while(file>>time) {
      ix=0;
      while(file>>energy>>energy>>n1>>n2>>n3>>e1>>e2>>e3) {
         number1[ix][iy]=n1; number2[ix][iy]=n2; number3[ix][iy]=n3;
         energy1[ix][iy]=e1; energy2[ix][iy]=e2; energy3[ix][iy]=e3;
         ix++;
         binEdgesx[ix]=energy;
         if (ix>=nbinx) break;
      }
      iy++;
      binEdgesy[iy]=time;
   }

   file.close();

   // The time axis in the database is not binned. In order to fill the data
   // into a 2D histogram, a time value in the database is regarded as the
   // center of a bin in the time axis of the histogram. Edges of bins are
   // set to the middle of two nearby time values.
   binEdgesy[0] = binEdgesy[1]-(binEdgesy[2] - binEdgesy[1])/2.;
   for (UShort_t i=1; i<nbiny; i++)
      binEdgesy[i] = (binEdgesy[i]+binEdgesy[i+1])/2.;
   binEdgesy[nbiny] = binEdgesy[nbiny]
      + (binEdgesy[nbiny] - binEdgesy[nbiny-1])/2.;

   // fill spectra
   fH2N1 = new TH2D(Form("h2N1%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Time [second];Number luminosity/(1 MeV)",
         nbinx,binEdgesx,nbiny,binEdgesy);
   fH2N2 = new TH2D(Form("h2N2%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Time [second];Number luminosity/(1 MeV)",
         nbinx,binEdgesx,nbiny,binEdgesy);
   fH2N3 = new TH2D(Form("h2N3%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Time [second];Number luminosity/(1 MeV)",
         nbinx,binEdgesx,nbiny,binEdgesy);

   fH2E1 = new TH2D(Form("h2E1%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Time [second];Energy luminosity [erg]/(1 MeV)",
         nbinx,binEdgesx,nbiny,binEdgesy);
   fH2E2 = new TH2D(Form("h2E2%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Time [second];Energy luminosity [erg]/(1 MeV)",
         nbinx,binEdgesx,nbiny,binEdgesy);
   fH2E3 = new TH2D(Form("h2E3%.0f%.0f%0.f",
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         ";Energy [MeV];Time [second];Energy luminosity [erg]/(1 MeV)",
         nbinx,binEdgesx,nbiny,binEdgesy);

   for (ix=0; ix<nbinx; ix++) {
      for (iy=0; iy<nbiny; iy++) {
         fH2N1->SetBinContent(ix+1,iy+1,number1[ix][iy]);
         fH2N2->SetBinContent(ix+1,iy+1,number2[ix][iy]);
         fH2N3->SetBinContent(ix+1,iy+1,number3[ix][iy]);

         fH2E1->SetBinContent(ix+1,iy+1,energy1[ix][iy]);
         fH2E2->SetBinContent(ix+1,iy+1,energy2[ix][iy]);
         fH2E3->SetBinContent(ix+1,iy+1,energy3[ix][iy]);
      }
   }

   // set properties
   fH2N1->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2N2->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2N3->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fH2E1->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2E2->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2E3->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fH2N1->SetStats(0);
   fH2N2->SetStats(0);
   fH2N3->SetStats(0);

   fH2E1->SetStats(0);
   fH2E2->SetStats(0);
   fH2E3->SetStats(0);

   fH2N1->SetLineColor(kBlack);
   fH2N2->SetLineColor(kRed);
   fH2N3->SetLineColor(kBlue);

   fH2E1->SetLineColor(kBlack);
   fH2E2->SetLineColor(kRed);
   fH2E3->SetLineColor(kBlue);
}

//______________________________________________________________________________
//

TH2D* SupernovaModel::NumberSpectrum(const char *neutrino)
{
   if (!fH2N1) return fH2Empty;

   TString species(neutrino);
   if (species=="v_e") return fH2N1;
   else if (species=="anti-v_e") return fH2N2;
   else if (species=="v_x") return fH2N3;
   else {
      Warning("NumberSpectrum","Neutrino species %s is not defined!",neutrino);
      Warning("NumberSpectrum","Please select one from v_e, anti-v_e or v_x.");
      Warning("NumberSpectrum","The spectrum of #nu_{e} is returned.");
      return fH2N1;
   }
}

//______________________________________________________________________________
//

TH2D* SupernovaModel::EnergySpectrum(const char *neutrino)
{
   if (!fH2E1) return fH2Empty;

   TString species(neutrino);
   if (species=="v_e") return fH2E1;
   else if (species=="anti-v_e") return fH2E2;
   else if (species=="v_x") return fH2E3;
   else {
      Warning("EnergySpectrum","Neutrino species %s is not defined!",neutrino);
      Warning("EnergySpectrum","Please select one from v_e, anti-v_e or v_x.");
      Warning("EnergySpectrum","The spectrum of #nu_{e} is returned.");
      return fH2E1;
   }
}

//______________________________________________________________________________
//

