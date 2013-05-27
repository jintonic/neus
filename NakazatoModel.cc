#include "NakazatoModel.h"

#include <TH2D.h>
#include <TDirectory.h>

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

   for (UShort_t i=0; i<7; i++) {
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
   for (UShort_t i=0; i<7; i++) {
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
   const UShort_t nbins = 20;
   Double_t binEdges[nbins+1]={0};
   Double_t number1[nbins];
   Double_t number2[nbins];
   Double_t numberx[nbins];
   Double_t energy1[nbins];
   Double_t energy2[nbins];
   Double_t energyx[nbins];
   Double_t maxN1=0,maxN2=0, maxNx=0, maxL1=0, maxL2=0, maxLx=0;
   Double_t minN1=1e100, minN2=1e100, minNx=1e100;
   Double_t minL1=1e100, minL2=1e100, minLx=1e100;

   Double_t energy, n1, n2, nx, e1, e2, ex;

   UShort_t i=0;
   while(file>>energy>>energy>>n1>>n2>>nx>>e1>>e2>>ex) {
      number1[i]=n1; number2[i]=n2; numberx[i]=nx; 
      energy1[i]=e1; energy2[i]=e2; energyx[i]=ex;

      if (n1>maxN1) maxN1=n1;
      if (n2>maxN2) maxN2=n2;
      if (nx>maxNx) maxNx=nx;

      if (e1>maxL1) maxL1=e1;
      if (e2>maxL2) maxL2=e2;
      if (ex>maxLx) maxLx=ex;

      if (n1<minN1) minN1=n1;
      if (n2<minN2) minN2=n2;
      if (nx<minNx) minNx=nx;

      if (e1<minL1) minL1=e1;
      if (e2<minL2) minL2=e2;
      if (ex<minLx) minLx=ex;

      i++;
      binEdges[i]=energy;
   }

   file.close();

   // fill spectra
   for (UShort_t i=1; i<=3; i++) {
      fHNe[i] = new TH1D(Form("hNe%d%.0f%.0f%.0f",i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";Energy [MeV];Number luminosity/(1 MeV)",nbins,binEdges);

      fHLe[i] = new TH1D(Form("hLe%d%.0f%.0f%.0f",i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";Energy [MeV];Energy luminosity [erg]/(1 MeV)",nbins,binEdges);
   }
   for (UShort_t i=4; i<=6; i++) {
      fHNe[i]=fHNe[3];
      fHLe[i]=fHLe[3];
   }

   for (i=1; i<=nbins; i++) {
      fHNe[1]->SetBinContent(i,number1[i-1]);
      fHNe[2]->SetBinContent(i,number2[i-1]);
      fHNe[3]->SetBinContent(i,numberx[i-1]);

      fHLe[1]->SetBinContent(i,energy1[i-1]);
      fHLe[2]->SetBinContent(i,energy2[i-1]);
      fHLe[3]->SetBinContent(i,energyx[i-1]);
   }

   // set properties
   fHNe[1]->SetMaximum(maxN1);
   fHNe[1]->SetMinimum(minN1);
   fHNe[2]->SetMaximum(maxN2);
   fHNe[2]->SetMinimum(minN2);
   fHNe[3]->SetMaximum(maxNx);
   fHNe[3]->SetMinimum(minNx);

   fHLe[1]->SetMaximum(maxL1);
   fHLe[1]->SetMinimum(minL1);
   fHLe[2]->SetMaximum(maxL2);
   fHLe[2]->SetMinimum(minL2);
   fHLe[3]->SetMaximum(maxLx);
   fHLe[3]->SetMinimum(minLx);

   fHNe[1]->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHNe[2]->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHNe[3]->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fHLe[1]->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHLe[2]->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHLe[3]->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
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
   const UShort_t nbinx = 20, nbiny = 391;
   Double_t binEdgesx[nbinx+1]={0};
   Double_t binEdgesy[nbiny+1]={0};
   Double_t number1[nbinx][nbiny];
   Double_t number2[nbinx][nbiny];
   Double_t numberx[nbinx][nbiny];
   Double_t energy1[nbinx][nbiny];
   Double_t energy2[nbinx][nbiny];
   Double_t energyx[nbinx][nbiny];
   Double_t maxN1=0,maxN2=0, maxNx=0, maxL1=0, maxL2=0, maxLx=0;
   Double_t minN1=1e100, minN2=1e100, minNx=1e100;
   Double_t minL1=1e100, minL2=1e100, minLx=1e100;

   Double_t time, energy, n1, n2, nx, e1, e2, ex;

   UShort_t ix=0, iy=0;
   while(file>>time) {
      ix=0;
      while(file>>energy>>energy>>n1>>n2>>nx>>e1>>e2>>ex) {
         number1[ix][iy]=n1; number2[ix][iy]=n2; numberx[ix][iy]=nx;
         energy1[ix][iy]=e1; energy2[ix][iy]=e2; energyx[ix][iy]=ex;

         if (n1>maxN1) maxN1=n1;
         if (n2>maxN2) maxN2=n2;
         if (nx>maxNx) maxNx=nx;

         if (e1>maxL1) maxL1=e1;
         if (e2>maxL2) maxL2=e2;
         if (ex>maxLx) maxLx=ex;

         if (n1<minN1) minN1=n1;
         if (n2<minN2) minN2=n2;
         if (nx<minNx) minNx=nx;

         if (e1<minL1) minL1=e1;
         if (e2<minL2) minL2=e2;
         if (ex<minLx) minLx=ex;

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

   // create histograms
   for (UShort_t i=1; i<=3; i++) {
      fHN2[i] = new TH2D(Form("hN2%d%.0f%.0f%.0f", i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";Energy [MeV];Time [second];Number luminosity/(1 MeV)",
            nbinx,binEdgesx,nbiny,binEdgesy);

      fHL2[i] = new TH2D(Form("hL2%d%.0f%.0f%.0f", i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";Energy [MeV];Time [second];Energy luminosity [erg]/(1 MeV)",
            nbinx,binEdgesx,nbiny,binEdgesy);
   }
   for (UShort_t i=4; i<=6; i++) {
      fHN2[i]=fHN2[3];
      fHL2[i]=fHL2[3];
   }

   // fill spectra
   for (ix=0; ix<nbinx; ix++) {
      for (iy=0; iy<nbiny; iy++) {
         fHN2[1]->SetBinContent(ix+1,iy+1,number1[ix][iy]);
         fHN2[2]->SetBinContent(ix+1,iy+1,number2[ix][iy]);
         fHN2[3]->SetBinContent(ix+1,iy+1,numberx[ix][iy]);

         fHL2[1]->SetBinContent(ix+1,iy+1,energy1[ix][iy]);
         fHL2[2]->SetBinContent(ix+1,iy+1,energy2[ix][iy]);
         fHL2[3]->SetBinContent(ix+1,iy+1,energyx[ix][iy]);
      }
   }

   // set properties
   fHN2[1]->SetMaximum(maxN1);
   fHN2[1]->SetMinimum(minN1);
   fHN2[2]->SetMaximum(maxN2);
   fHN2[2]->SetMinimum(minN2);
   fHN2[3]->SetMaximum(maxNx);
   fHN2[3]->SetMinimum(minNx);

   fHL2[1]->SetMaximum(maxL1);
   fHL2[1]->SetMinimum(minL1);
   fHL2[2]->SetMaximum(maxL2);
   fHL2[2]->SetMinimum(minL2);
   fHL2[3]->SetMaximum(maxLx);
   fHL2[3]->SetMinimum(minLx);

   fHN2[1]->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHN2[2]->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHN2[3]->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fHL2[1]->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHL2[2]->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fHL2[3]->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
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

Double_t NEUS::NakazatoModel::N2(UShort_t type, Double_t energy, Double_t time)
{
   return HN2(type)->Interpolate(energy, time);
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
   if (time>20) time=20;

   // look for the right bin in time axis
   Double_t t=-5e-2;
   Int_t ibin=0;
   while (t<time) {
      ibin++;
      t += HNt(type)->GetBinContent(ibin);
   }

   return HNt(type)->Integral("width");
}

//______________________________________________________________________________
//

Double_t NEUS::NakazatoModel::Nall(UShort_t type)
{
   return HNe(type)->Integral("width");
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
   if (!fHNe[type]) {
      Warning("HNe","Spectrum does not exist!");
      Warning("HNe","Is the database correctly loaded?");
      Warning("HNe","NULL pointer is returned!");
      return NULL;
   }

   if (tmax==0) return fHNe[type]; // from database

   TH1D *h = (TH1D*) gDirectory->Get(Form("hNe%d%.0f%.0f%.0f%f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100,tmax));
   if (h) return h; // from memory

   // calculate integral in [0, tmax]
   fHNe[type]->Reset();
   for (UShort_t ix=1; ix<=fHN2[type]->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=fHN2[type]->GetNbinsY(); iy++) {
         if (tmax<fHN2[type]->GetYaxis()->GetBinCenter(iy)) break;
         content += fHN2[type]->GetBinContent(ix,iy) *
            fHN2[type]->GetYaxis()->GetBinWidth(iy);
      }
      fHNe[type]->SetBinContent(ix,content);
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
   if (!fHLe[type]) {
      Warning("HLe","Spectrum does not exist!");
      Warning("HLe","Is the database correctly loaded?");
      Warning("HLe","NULL pointer is returned!");
      return NULL;
   }

   if (tmax==0) return fHLe[type]; // from database

   TH1D *h = (TH1D*) gDirectory->Get(Form("hLe%d%.0f%.0f%.0f%f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100,tmax));
   if (h) return h; // from memory

   // calculate integral in [0, tmax]
   fHLe[type]->Reset();
   for (UShort_t ix=1; ix<=fHL2[type]->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=fHL2[type]->GetNbinsY(); iy++) {
         if (tmax<fHL2[type]->GetYaxis()->GetBinCenter(iy)) break;
         content += fHL2[type]->GetBinContent(ix,iy) *
            fHL2[type]->GetYaxis()->GetBinWidth(iy);
      }
      fHLe[type]->SetBinContent(ix,content);
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

   fHNt[type] = fHN2[type]->ProjectionY(Form("hNt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHNt[type]->Reset();
   fHNt[type]->SetStats(0);

   // calculate integral
   for (UShort_t iy=1; iy<=fHN2[type]->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=fHN2[type]->GetNbinsX(); ix++)
         content += fHN2[type]->GetBinContent(ix,iy) *
            fHN2[type]->GetXaxis()->GetBinWidth(ix);
      fHNt[type]->SetBinContent(iy,content);
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

   fHLt[type] = fHL2[type]->ProjectionY(Form("hLt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHLt[type]->Reset();
   fHLt[type]->SetStats(0);
   fHLt[type]->GetYaxis()->SetTitle("Luminosity [erg] / second");

   // calculate integral
   for (UShort_t iy=1; iy<=fHL2[type]->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=fHL2[type]->GetNbinsX(); ix++)
         content += fHL2[type]->GetBinContent(ix,iy) *
            fHL2[type]->GetXaxis()->GetBinWidth(ix);
      fHLt[type]->SetBinContent(iy,content);
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

   fHEt[type] = fHN2[type]->ProjectionY(Form("hEt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHEt[type]->Reset();
   fHEt[type]->SetStats(0);
   fHEt[type]->GetYaxis()->SetTitle("Average energy [MeV] / second");

   // calculate average
   for (UShort_t iy=1; iy<=fHN2[type]->GetNbinsY(); iy++) {
      Double_t totalE=0, totalN=0;
      for (UShort_t ix=1; ix<=fHN2[type]->GetNbinsX(); ix++) {
         totalN += fHN2[type]->GetBinContent(ix,iy) *
            fHN2[type]->GetXaxis()->GetBinWidth(ix);
         totalE += fHN2[type]->GetBinContent(ix,iy) * 
            fHN2[type]->GetXaxis()->GetBinWidth(ix) *
            fHN2[type]->GetXaxis()->GetBinCenter(ix);
      }
      fHEt[type]->SetBinContent(iy,totalE/totalN);
   }
   return fHEt[type];
}

//______________________________________________________________________________
//

void NEUS::NakazatoModel::Print()
{
   Printf("%2.0f Solar mass, %.3f, %3.0f ms: n1=%1.2e, n2=%1.2e, nx=%1.2e, total=%1.2e",
         fInitialMass, fMetallicity, fReviveTime, Nall(1), Nall(2), Nall(3),
         Nall(1) + Nall(2) + Nall(3)*4);
}
