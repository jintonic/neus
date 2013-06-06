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

   for (UShort_t i=0; i<fgNtype; i++) {
      fHN2[i]=NULL;
      fHL2[i]=NULL;
      fHNe[i]=NULL;
      fHNt[i]=NULL;
      fHLe[i]=NULL;
      fHLt[i]=NULL;
      fHEt[i]=NULL;
   }
}

//______________________________________________________________________________
//

NEUS::NakazatoModel::~NakazatoModel()
{
   for (UShort_t i=0; i<fgNtype; i++) {
      if (fHN2[i]) delete fHN2[i];
      if (fHL2[i]) delete fHL2[i];
      if (fHNe[i]) delete fHNe[i];
      if (fHNt[i]) delete fHNt[i];
      if (fHLe[i]) delete fHLe[i];
      if (fHLt[i]) delete fHLt[i];
      if (fHEt[i]) delete fHEt[i];
   }
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

   // fill spectra
   for (UShort_t i=1; i<=3; i++) {
      fHNe[i] = new TH1D(Form("hNe%d%.0f%.0f%.0f",i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";neutrino energy [MeV];",fNbinsE,binEdges);

      fHLe[i] = new TH1D(Form("hLe%d%.0f%.0f%.0f",i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";neutrino energy [MeV];",fNbinsE,binEdges);
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
   fHNe[1]->GetYaxis()->SetTitle("number of #nu_{e} [10^{50}/MeV]");
   fHNe[2]->GetYaxis()->SetTitle("number of #bar{#nu}_{e} [10^{50}/MeV]");
   fHNe[3]->GetYaxis()->SetTitle("number of #nu_{x} [10^{50}/MeV]");

   fHLe[1]->GetYaxis()->SetTitle("luminosity of #nu_{e} [10^{50} erg/MeV]");
   fHLe[2]->GetYaxis()->SetTitle("luminosity of #bar{#nu}_{e} [10^{50} erg/MeV]");
   fHLe[3]->GetYaxis()->SetTitle("luminosity of #nu_{x} [10^{50} erg/MeV]");

   fHNe[1]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHNe[2]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHNe[3]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fHLe[1]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHLe[2]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHLe[3]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

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
            ";time [second];neutrino energy [MeV];",
            fNbinsT,binEdgesx,fNbinsE,binEdgesy);

      fHL2[i] = new TH2D(Form("hL2%d%.0f%.0f%.0f", i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";time [second];neutrino energy [MeV];",
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

   fHN2[1]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHN2[2]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHN2[3]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fHL2[1]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHL2[2]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHL2[3]->SetTitle(Form(
            "progenitor: %.0f Solar mass, metallicity: %.3f, revive time: %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

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

Double_t NEUS::NakazatoModel::N2(UShort_t type, Double_t time, Double_t energy)
{
   return HN2(type)->Interpolate(time, energy);
}

//______________________________________________________________________________
//

Double_t NEUS::NakazatoModel::Ne(UShort_t type, Double_t energy)
{
   return HNe(type)->Interpolate(energy);
}

//______________________________________________________________________________
//

Double_t NEUS::NakazatoModel::Nt(UShort_t type, Double_t time)
{
   return HNt(type)->Interpolate(time);
}

//______________________________________________________________________________
//

Double_t NEUS::NakazatoModel::Nall(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("Nall","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Nall","Return 0!");
      return 0;
   }
   if (fTotalN[type]==0) fTotalN[type] = HNe(type)->Integral("width");
   return fTotalN[type];
}

//______________________________________________________________________________
//

Double_t NEUS::NakazatoModel::Lall(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("Lall","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Lall","Return 0!");
      return 0;
   }
   if (fTotalL[type]==0) fTotalL[type] = HLe(type)->Integral("width");
   return fTotalL[type];
}

//______________________________________________________________________________
//

Double_t NEUS::NakazatoModel::Eave(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("Eave","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Eave","Return 0!");
      return 0;
   }
   if (abs(fAverageE[type]-1.0)<0.01)
      fAverageE[type] = HLe(type)->Integral("width")/1.60217646e-6/Nall(type);
   return fAverageE[type];
}

//______________________________________________________________________________
//

TH2D* NEUS::NakazatoModel::HN2(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("HN2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HN2","NULL pointer is returned!");
      return NULL;
   }
   if (!fHN2[type]) {
      Warning("HN2","Spectrum does not exist!");
      Warning("HN2","Is the database correctly loaded?");
      Warning("HN2","NULL pointer is returned!");
   }
   return fHN2[type];
}

//______________________________________________________________________________
//

TH2D* NEUS::NakazatoModel::HL2(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("HL2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HL2","NULL pointer is returned!");
      return NULL;
   }
   if (!fHL2[type]) {
      Warning("HL2","Spectrum does not exist!");
      Warning("HL2","Is the database correctly loaded?");
      Warning("HL2","NULL pointer is returned!");
   }
   return fHL2[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::NakazatoModel::HNe(UShort_t type, Double_t tmax)
{
   if (type<1 || type>6) {
      Warning("HNe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HNe","NULL pointer is returned!");
      return NULL;
   }
   if (tmax==0) return fHNe[type]; // from database
   if (tmax>fMaxT) tmax=fMaxT;

   if (!fHNe[type]) {
      Warning("HNe","Spectrum does not exist!");
      Warning("HNe","Is the database correctly loaded?");
      Warning("HNe","NULL pointer is returned!");
      return NULL;
   }

   TH1D *h = (TH1D*) gDirectory->Get(Form("hNe%d%.0f%.0f%.0f%f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100,tmax));
   if (h) return h; // from memory

   // calculate integral in [0, tmax]
   fHNe[type]->Reset();
   for (UShort_t iy=1; iy<=fHN2[type]->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=fHN2[type]->GetNbinsX(); ix++) {
         if (tmax<fHN2[type]->GetXaxis()->GetBinCenter(ix)) break;
         content += fHN2[type]->GetBinContent(ix,iy) *
            fHN2[type]->GetXaxis()->GetBinWidth(ix);
      }
      fHNe[type]->SetBinContent(iy,content);
   }
   fHNe[type]->SetName(Form("%s%f",fHNe[type]->GetName(),tmax));
   return fHNe[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::NakazatoModel::HLe(UShort_t type, Double_t tmax)
{
   if (type<1 || type>6) {
      Warning("HLe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HLe","NULL pointer is returned!");
      return NULL;
   }
   if (tmax==0) return fHLe[type]; // from database
   if (tmax>fMaxT) tmax=fMaxT;

   if (!fHLe[type]) {
      Warning("HLe","Spectrum does not exist!");
      Warning("HLe","Is the database correctly loaded?");
      Warning("HLe","NULL pointer is returned!");
      return NULL;
   }

   TH1D *h = (TH1D*) gDirectory->Get(Form("hLe%d%.0f%.0f%.0f%f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100,tmax));
   if (h) return h; // from memory

   // calculate integral in [0, tmax]
   fHLe[type]->Reset();
   for (UShort_t iy=1; iy<=fHL2[type]->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=fHL2[type]->GetNbinsX(); ix++) {
         if (tmax<fHL2[type]->GetXaxis()->GetBinCenter(ix)) break;
         content += fHL2[type]->GetBinContent(ix,iy) *
            fHL2[type]->GetXaxis()->GetBinWidth(ix);
      }
      fHLe[type]->SetBinContent(iy,content);
   }
   fHLe[type]->SetName(Form("%s%f",fHLe[type]->GetName(),tmax));
   return fHLe[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::NakazatoModel::HNt(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("HNt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HNt","NULL pointer is returned!");
      return NULL;
   }
   if (!fHN2[type]) {
      Warning("HNt","Spectrum does not exist!");
      Warning("HNt","Is the database correctly loaded?");
      Warning("HNt","NULL pointer is returned!");
      return NULL;
   }
   if (fHNt[type]) return fHNt[type];

   fHNt[type] = fHN2[type]->ProjectionX(Form("hNt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHNt[type]->Reset();
   fHNt[type]->SetStats(0);

   // calculate integral
   for (UShort_t ix=1; ix<=fHN2[type]->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=fHN2[type]->GetNbinsY(); iy++)
         content += fHN2[type]->GetBinContent(ix,iy) *
            fHN2[type]->GetYaxis()->GetBinWidth(iy);
      fHNt[type]->SetBinContent(ix,content);
   }
   return fHNt[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::NakazatoModel::HLt(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("HLt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HLt","NULL pointer is returned!");
      return NULL;
   }
   if (!fHL2[type]) {
      Warning("HLt","Spectrum does not exist!");
      Warning("HLt","Is the database correctly loaded?");
      Warning("HLt","NULL pointer is returned!");
      return NULL;
   }
   if (fHLt[type]) return fHLt[type];

   fHLt[type] = fHL2[type]->ProjectionX(Form("hLt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHLt[type]->Reset();
   fHLt[type]->SetStats(0);
   fHLt[type]->GetYaxis()->SetTitle("Luminosity [10^{50} erg/second]");

   // calculate integral
   for (UShort_t ix=1; ix<=fHL2[type]->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=fHL2[type]->GetNbinsY(); iy++)
         content += fHL2[type]->GetBinContent(ix,iy) *
            fHL2[type]->GetYaxis()->GetBinWidth(iy);
      fHLt[type]->SetBinContent(ix,content);
   }
   return fHLt[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::NakazatoModel::HEt(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("HEt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HEt","NULL pointer is returned!");
      return NULL;
   }
   if (!fHN2[type]) {
      Warning("HEt","Spectrum does not exist!");
      Warning("HEt","Is the database correctly loaded?");
      Warning("HEt","NULL pointer is returned!");
      return NULL;
   }
   if (fHEt[type]) return fHEt[type];

   fHEt[type] = fHN2[type]->ProjectionX(Form("hEt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHEt[type]->Reset();
   fHEt[type]->SetStats(0);
   fHEt[type]->GetYaxis()->SetTitle("average neutrino energy [MeV/second]");

   // calculate average
   for (UShort_t ix=1; ix<=fHN2[type]->GetNbinsX(); ix++) {
      Double_t totalE=0, totalN=0;
      for (UShort_t iy=1; iy<=fHN2[type]->GetNbinsY(); iy++) {
         totalN += fHN2[type]->GetBinContent(ix,iy) *
            fHN2[type]->GetYaxis()->GetBinWidth(iy);
         totalE += fHN2[type]->GetBinContent(ix,iy) * 
            fHN2[type]->GetYaxis()->GetBinWidth(iy) *
            fHN2[type]->GetYaxis()->GetBinCenter(iy);
      }
      fHEt[type]->SetBinContent(ix,totalE/totalN);
   }
   return fHEt[type];
}

//______________________________________________________________________________
//

void NEUS::NakazatoModel::Print()
{
   Printf("%2.0f Solar mass, %.3f, %3.0f ms: N1=%1.2e, N2=%1.2e, Nx=%1.2e, N=%1.2e, L=%1.2e erg",
         fInitialMass, fMetallicity, fReviveTime,
         Nall(1)*1e50, Nall(2)*1e50, Nall(3)*1e50, 
         Nall(1)*1e50 + Nall(2)*1e50 + Nall(3)*1e50*4,
         Lall(1)*1e50 + Lall(2)*1e50 + Lall(3)*1e50*4);
}
