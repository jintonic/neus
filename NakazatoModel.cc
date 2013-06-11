#include "NakazatoModel.h"

#include <TH2D.h>
#include <TDirectory.h>

#include <cmath>
#include <fstream>
using namespace std;

//______________________________________________________________________________
//

NEUS::NakazatoModel::NakazatoModel(
      Float_t initialMass, Float_t metallicity, Float_t reviveTime) : 
   SupernovaModel(), fInitialMass(initialMass), 
   fMetallicity(metallicity), fReviveTime(reviveTime)
{
   if (initialMass==30 && metallicity<0.01) fReviveTime=0;

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

void NEUS::NakazatoModel::LoadData(const char *dir)
{
   SupernovaModel::LoadData(dir); // set fDataLocation
   LoadIntegratedData();
   // no full data for the black hole
   if (fInitialMass!=30 || fMetallicity>0.01) LoadFullData();
}

//______________________________________________________________________________
//

void NEUS::NakazatoModel::LoadIntegratedData()
{
   char *name = Form("%s/integdata/integ%.0f0%.0f.data",
         fDataLocation.Data(), fInitialMass, fReviveTime/100);
   if (fMetallicity<0.01)
      name = Form("%s/integdata/integ%.0f1%.0f.data",
            fDataLocation.Data(), fInitialMass, fReviveTime/100);

   // check database
   ifstream file(name);
   if (!(file.is_open())) {
      Warning("LoadIntegratedData", "%s cannot be read!", name);
      return;
   }

   // skip the first line
   char line[150];
   file.getline(line,150);

   // load data
   Double_t binEdges[fNbinsE+1]={0};
   Double_t number1[fNbinsE];
   Double_t number2[fNbinsE];
   Double_t numberx[fNbinsE];
   Double_t energy1[fNbinsE];
   Double_t energy2[fNbinsE];
   Double_t energyx[fNbinsE];

   Double_t energy, n1, n2, nx, e1, e2, ex;

   UShort_t i=0;
   while(file>>energy>>energy>>n1>>n2>>nx>>e1>>e2>>ex) {
      number1[i]=n1; number2[i]=n2; numberx[i]=nx; 
      energy1[i]=e1; energy2[i]=e2; energyx[i]=ex;

      i++;
      binEdges[i]=energy;
   }

   file.close();

   fMinE = binEdges[0];
   fMaxE = binEdges[fNbinsE];
   fMaxT = 20.125; // used to set title of histograms for identification

   // fill spectra
   for (UShort_t i=1; i<=3; i++) {
      fHNe[i] = new TH1D(Form("hNe-%s-%d-%.4f", GetName(), i, fMaxT),
            ";energy [MeV];number of neutrinos [10^{50}/MeV]",fNbinsE,binEdges);

      fHLe[i] = new TH1D(Form("hLe-%s-%d-%.4f", GetName(), i, fMaxT),
            ";energy [MeV];luminosity [10^{50} erg/MeV]",fNbinsE,binEdges);
   }
   for (UShort_t i=4; i<=6; i++) {
      fHNe[i]=fHNe[3];
      fHLe[i]=fHLe[3];
   }

   for (i=1; i<=fNbinsE; i++) {
      fHNe[1]->SetBinContent(i,number1[i-1]/1e50);
      fHNe[2]->SetBinContent(i,number2[i-1]/1e50);
      fHNe[3]->SetBinContent(i,numberx[i-1]/1e50);

      fHLe[1]->SetBinContent(i,energy1[i-1]/1e50);
      fHLe[2]->SetBinContent(i,energy2[i-1]/1e50);
      fHLe[3]->SetBinContent(i,energyx[i-1]/1e50);
   }

   // set properties
   fHNe[1]->SetTitle(GetTitle());
   fHNe[2]->SetTitle(GetTitle());
   fHNe[3]->SetTitle(GetTitle());

   fHLe[1]->SetTitle(GetTitle());
   fHLe[2]->SetTitle(GetTitle());
   fHLe[3]->SetTitle(GetTitle());

   fHNe[1]->SetStats(0);
   fHNe[2]->SetStats(0);
   fHNe[3]->SetStats(0);

   fHLe[1]->SetStats(0);
   fHLe[2]->SetStats(0);
   fHLe[3]->SetStats(0);

   fHNe[1]->SetLineColor(kBlack);
   fHNe[2]->SetLineColor(kRed);
   fHNe[3]->SetLineColor(kBlue);

   fHLe[1]->SetLineColor(kBlack);
   fHLe[2]->SetLineColor(kRed);
   fHLe[3]->SetLineColor(kBlue);

   fHNe[1]->SetMarkerColor(kBlack);
   fHNe[2]->SetMarkerColor(kRed);
   fHNe[3]->SetMarkerColor(kBlue);

   fHLe[1]->SetMarkerColor(kBlack);
   fHLe[2]->SetMarkerColor(kRed);
   fHLe[3]->SetMarkerColor(kBlue);

   fHNe[1]->SetMarkerStyle(20);
   fHNe[2]->SetMarkerStyle(20);
   fHNe[3]->SetMarkerStyle(20);

   fHLe[1]->SetMarkerStyle(20);
   fHLe[2]->SetMarkerStyle(20);
   fHLe[3]->SetMarkerStyle(20);

   fHNe[1]->SetMarkerSize(0.6);
   fHNe[2]->SetMarkerSize(0.6);
   fHNe[3]->SetMarkerSize(0.6);

   fHLe[1]->SetMarkerSize(0.6);
   fHLe[2]->SetMarkerSize(0.6);
   fHLe[3]->SetMarkerSize(0.6);
}

//______________________________________________________________________________
//

void NEUS::NakazatoModel::LoadFullData()
{
   if (fInitialMass==30. && fMetallicity<0.01) {
      Warning("LoadFullData", "No full data for black hole");
      Warning("LoadFullData", "with inital mass = 30 Solar mass,");
      Warning("LoadFullData", "and metallicity of 0.004.");
      return;
   }

   char *name = Form("%s/intpdata/intp%.0f0%.0f.data",
         fDataLocation.Data(), fInitialMass, fReviveTime/100);
   if (fMetallicity<0.01)
      name = Form("%s/intpdata/intp%.0f1%.0f.data",
            fDataLocation.Data(), fInitialMass, fReviveTime/100);

   // check database
   ifstream file(name);
   if (!(file.is_open())) {
      Warning("LoadFullData", "%s cannot be read!", name);
      return;
   }

   // load data
   Double_t binEdgesx[fNbinsT+1]={0};
   Double_t binEdgesy[fNbinsE+1]={0};
   Double_t number1[fNbinsT][fNbinsE];
   Double_t number2[fNbinsT][fNbinsE];
   Double_t numberx[fNbinsT][fNbinsE];
   Double_t energy1[fNbinsT][fNbinsE];
   Double_t energy2[fNbinsT][fNbinsE];
   Double_t energyx[fNbinsT][fNbinsE];

   Double_t time, energy, n1, n2, nx, e1, e2, ex;

   UShort_t ix=0, iy=0;
   while(file>>time) {
      iy=0;
      while(file>>energy>>energy>>n1>>n2>>nx>>e1>>e2>>ex) {
         number1[ix][iy]=n1; number2[ix][iy]=n2; numberx[ix][iy]=nx;
         energy1[ix][iy]=e1; energy2[ix][iy]=e2; energyx[ix][iy]=ex;

         iy++;
         binEdgesy[iy]=energy;
         if (iy>=fNbinsE) break;
      }
      ix++;
      binEdgesx[ix]=time;
   }

   file.close();

   // The time axis in the database is not binned. In order to fill the data
   // into a 2D histogram, a time value in the database is regarded as the
   // center of a bin in the time axis of the histogram. Edges of bins are
   // set to the middle of two nearby time values.
   binEdgesx[0] = binEdgesx[1]-(binEdgesx[2] - binEdgesx[1])/2.;
   for (UShort_t i=1; i<fNbinsT; i++)
      binEdgesx[i] = (binEdgesx[i]+binEdgesx[i+1])/2.;
   binEdgesx[fNbinsT] = binEdgesx[fNbinsT]
      + (binEdgesx[fNbinsT] - binEdgesx[fNbinsT-1])/2.;

   fMinT = binEdgesx[0];
   fMaxT = binEdgesx[fNbinsT];
   fMinE = binEdgesy[0];
   fMaxE = binEdgesy[fNbinsE];

   // create histograms
   for (UShort_t i=1; i<=3; i++) {
      fHN2[i] = new TH2D(Form("hN2%d%.0f%.0f%.0f", i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";time [second];energy [MeV];",
            fNbinsT,binEdgesx,fNbinsE,binEdgesy);

      fHL2[i] = new TH2D(Form("hL2%d%.0f%.0f%.0f", i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";time [second];energy [MeV];",
            fNbinsT,binEdgesx,fNbinsE,binEdgesy);
   }
   for (UShort_t i=4; i<=6; i++) {
      fHN2[i]=fHN2[3];
      fHL2[i]=fHL2[3];
   }

   // fill spectra
   for (ix=0; ix<fNbinsT; ix++) {
      for (iy=0; iy<fNbinsE; iy++) {
         fHN2[1]->SetBinContent(ix+1,iy+1,number1[ix][iy]/1e50);
         fHN2[2]->SetBinContent(ix+1,iy+1,number2[ix][iy]/1e50);
         fHN2[3]->SetBinContent(ix+1,iy+1,numberx[ix][iy]/1e50);

         fHL2[1]->SetBinContent(ix+1,iy+1,energy1[ix][iy]/1e50);
         fHL2[2]->SetBinContent(ix+1,iy+1,energy2[ix][iy]/1e50);
         fHL2[3]->SetBinContent(ix+1,iy+1,energyx[ix][iy]/1e50);
      }
   }

   // set properties
   fHN2[1]->GetZaxis()->SetTitle("number of #nu_{e} [10^{50}/s/MeV]");
   fHN2[2]->GetZaxis()->SetTitle("number of #bar{#nu}_{e} [10^{50}/s/MeV]");
   fHN2[3]->GetZaxis()->SetTitle("number of #nu_{x} [10^{50}/s/MeV]");

   fHL2[1]->GetZaxis()->SetTitle("luminosity of #nu_{e} [10^{50} erg/s/MeV]");
   fHL2[2]->GetZaxis()->SetTitle("luminosity of #bar{#nu}_{e} [10^{50} erg/s/MeV]");
   fHL2[3]->GetZaxis()->SetTitle("luminosity of #nu_{x} [10^{50} erg/s/MeV]");

   fHN2[1]->GetZaxis()->SetTitleOffset(-0.5);
   fHN2[2]->GetZaxis()->SetTitleOffset(-0.5);
   fHN2[3]->GetZaxis()->SetTitleOffset(-0.5);

   fHL2[1]->GetZaxis()->SetTitleOffset(-0.5);
   fHL2[2]->GetZaxis()->SetTitleOffset(-0.5);
   fHL2[3]->GetZaxis()->SetTitleOffset(-0.5);

   fHN2[1]->GetZaxis()->CenterTitle();
   fHN2[2]->GetZaxis()->CenterTitle();
   fHN2[3]->GetZaxis()->CenterTitle();

   fHL2[1]->GetZaxis()->CenterTitle();
   fHL2[2]->GetZaxis()->CenterTitle();
   fHL2[3]->GetZaxis()->CenterTitle();

   fHN2[1]->SetTitle(GetTitle());
   fHN2[2]->SetTitle(GetTitle());
   fHN2[3]->SetTitle(GetTitle());

   fHL2[1]->SetTitle(GetTitle());
   fHL2[2]->SetTitle(GetTitle());
   fHL2[3]->SetTitle(GetTitle());

   fHN2[1]->SetStats(0);
   fHN2[2]->SetStats(0);
   fHN2[3]->SetStats(0);

   fHL2[1]->SetStats(0);
   fHL2[2]->SetStats(0);
   fHL2[3]->SetStats(0);

   fHN2[1]->SetLineColor(kBlack);
   fHN2[2]->SetLineColor(kRed);
   fHN2[3]->SetLineColor(kBlue);

   fHL2[1]->SetLineColor(kBlack);
   fHL2[2]->SetLineColor(kRed);
   fHL2[3]->SetLineColor(kBlue);
}

//______________________________________________________________________________
//

void NEUS::NakazatoModel::Print()
{
   Printf("%2.0f Solar mass, %.3f, %3.0f ms: N1=%1.2e, N2=%1.2e, Nx=%1.2e, N=%1.2e, L=%1.2e ergs",
         fInitialMass, fMetallicity, fReviveTime,
         Nall(1)*1e50, Nall(2)*1e50, Nall(3)*1e50, 
         Nall(1)*1e50 + Nall(2)*1e50 + Nall(3)*1e50*4,
         Lall(1)*1e50 + Lall(2)*1e50 + Lall(3)*1e50*4);
}
