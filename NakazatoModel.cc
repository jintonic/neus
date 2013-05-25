#include "NakazatoModel.h"

#include <TH2D.h>
#include <TDirectory.h>

#include <fstream>
using namespace std;

//______________________________________________________________________________
//

NEUS::NakazatoModel::NakazatoModel(
      Float_t initialMass, Float_t metallicity, Float_t reviveTime) : TNamed(),
   fInitialMass(initialMass), fMetallicity(metallicity), fReviveTime(reviveTime)
{
   if (metallicity>0.02-0.01) {
      fName = Form("model%.0f0%.0f",fInitialMass,fReviveTime/100);
      fTitle = Form("%.0f Solar mass, 0.02, %.0f ms",fInitialMass,fReviveTime);
   } else {
      fName = Form("model%.0f1%.0f",fInitialMass,fReviveTime/100);
      fTitle = Form("%.0f Solar mass, 0.004, %.0f ms",fInitialMass,fReviveTime);
   }

   for (UShort_t i=0; i<7; i++) {
      fH2N[i]=NULL;
      fH2L[i]=NULL;
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
      if (fH2N[i]) delete fH2N[i];
      if (fH2L[i]) delete fH2L[i];
      if (fHNe[i]) delete fHNe[i];
      if (fHNt[i]) delete fHNt[i];
      if (fHLe[i]) delete fHLe[i];
      if (fHLt[i]) delete fHLt[i];
      if (fHEt[i]) delete fHEt[i];
   }
}

//______________________________________________________________________________
//

void NEUS::NakazatoModel::LoadIntegratedData(const char *databaseDir)
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

void NEUS::NakazatoModel::LoadFullData(const char *databaseDir)
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
      fH2N[i] = new TH2D(Form("h2N%d%.0f%.0f%.0f", i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";Energy [MeV];Time [second];Number luminosity/(1 MeV)",
            nbinx,binEdgesx,nbiny,binEdgesy);

      fH2L[i] = new TH2D(Form("h2L%d%.0f%.0f%.0f", i,
               fInitialMass,fMetallicity*1000,fReviveTime/100),
            ";Energy [MeV];Time [second];Energy luminosity [erg]/(1 MeV)",
            nbinx,binEdgesx,nbiny,binEdgesy);
   }
   for (UShort_t i=4; i<=6; i++) {
      fH2N[i]=fH2N[3];
      fH2L[i]=fH2L[3];
   }

   // fill spectra
   for (ix=0; ix<nbinx; ix++) {
      for (iy=0; iy<nbiny; iy++) {
         fH2N[1]->SetBinContent(ix+1,iy+1,number1[ix][iy]);
         fH2N[2]->SetBinContent(ix+1,iy+1,number2[ix][iy]);
         fH2N[3]->SetBinContent(ix+1,iy+1,numberx[ix][iy]);

         fH2L[1]->SetBinContent(ix+1,iy+1,energy1[ix][iy]);
         fH2L[2]->SetBinContent(ix+1,iy+1,energy2[ix][iy]);
         fH2L[3]->SetBinContent(ix+1,iy+1,energyx[ix][iy]);
      }
   }

   // set properties
   fH2N[1]->SetMaximum(maxN1);
   fH2N[1]->SetMinimum(minN1);
   fH2N[2]->SetMaximum(maxN2);
   fH2N[2]->SetMinimum(minN2);
   fH2N[3]->SetMaximum(maxNx);
   fH2N[3]->SetMinimum(minNx);

   fH2L[1]->SetMaximum(maxL1);
   fH2L[1]->SetMinimum(minL1);
   fH2L[2]->SetMaximum(maxL2);
   fH2L[2]->SetMinimum(minL2);
   fH2L[3]->SetMaximum(maxLx);
   fH2L[3]->SetMinimum(minLx);

   fH2N[1]->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2N[2]->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2N[3]->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fH2L[1]->SetTitle(Form("#nu_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2L[2]->SetTitle(Form("#bar{#nu}_{e}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));
   fH2L[3]->SetTitle(Form("#nu_{x}, model: %.0f M_{#odot}, %.3f, %.0f ms",
            fInitialMass, fMetallicity, fReviveTime));

   fH2N[1]->SetStats(0);
   fH2N[2]->SetStats(0);
   fH2N[3]->SetStats(0);

   fH2L[1]->SetStats(0);
   fH2L[2]->SetStats(0);
   fH2L[3]->SetStats(0);

   fH2N[1]->SetLineColor(kBlack);
   fH2N[2]->SetLineColor(kRed);
   fH2N[3]->SetLineColor(kBlue);

   fH2L[1]->SetLineColor(kBlack);
   fH2L[2]->SetLineColor(kRed);
   fH2L[3]->SetLineColor(kBlue);
}

//______________________________________________________________________________
//

TH2D* NEUS::NakazatoModel::H2N(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("H2N","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("H2N","NULL pointer is returned!");
      return NULL;
   }
   if (!fH2N[type]) {
      Warning("H2N","Spectrum does not exist!");
      Warning("H2N","Is the database correctly loaded?");
      Warning("H2N","NULL pointer is returned!");
   }
   return fH2N[type];
}

//______________________________________________________________________________
//

TH2D* NEUS::NakazatoModel::H2L(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("H2L","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("H2L","NULL pointer is returned!");
      return NULL;
   }
   if (!fH2L[type]) {
      Warning("H2L","Spectrum does not exist!");
      Warning("H2L","Is the database correctly loaded?");
      Warning("H2L","NULL pointer is returned!");
   }
   return fH2L[type];
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
   for (UShort_t ix=1; ix<=fH2N[type]->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=fH2N[type]->GetNbinsY(); iy++) {
         if (tmax<fH2N[type]->GetYaxis()->GetBinCenter(iy)) break;
         content += fH2N[type]->GetBinContent(ix,iy) *
            fH2N[type]->GetYaxis()->GetBinWidth(iy);
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
   for (UShort_t ix=1; ix<=fH2L[type]->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=fH2L[type]->GetNbinsY(); iy++) {
         if (tmax<fH2L[type]->GetYaxis()->GetBinCenter(iy)) break;
         content += fH2L[type]->GetBinContent(ix,iy) *
            fH2L[type]->GetYaxis()->GetBinWidth(iy);
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
   if (!fH2N[type]) {
      Warning("HNt","Spectrum does not exist!");
      Warning("HNt","Is the database correctly loaded?");
      Warning("HNt","NULL pointer is returned!");
      return NULL;
   }
   if (fHNt[type]) return fHNt[type];

   fHNt[type] = fH2N[type]->ProjectionY(Form("hNt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHNt[type]->Reset();
   fHNt[type]->SetStats(0);

   // calculate integral
   for (UShort_t iy=1; iy<=fH2N[type]->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=fH2N[type]->GetNbinsX(); ix++)
         content += fH2N[type]->GetBinContent(ix,iy) *
            fH2N[type]->GetXaxis()->GetBinWidth(ix);
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
   if (!fH2L[type]) {
      Warning("HLt","Spectrum does not exist!");
      Warning("HLt","Is the database correctly loaded?");
      Warning("HLt","NULL pointer is returned!");
      return NULL;
   }
   if (fHLt[type]) return fHLt[type];

   fHLt[type] = fH2L[type]->ProjectionY(Form("hLt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHLt[type]->Reset();
   fHLt[type]->SetStats(0);
   fHLt[type]->GetYaxis()->SetTitle("Luminosity [erg] / second");

   // calculate integral
   for (UShort_t iy=1; iy<=fH2L[type]->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=fH2L[type]->GetNbinsX(); ix++)
         content += fH2L[type]->GetBinContent(ix,iy) *
            fH2L[type]->GetXaxis()->GetBinWidth(ix);
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
   if (!fH2N[type]) {
      Warning("HEt","Spectrum does not exist!");
      Warning("HEt","Is the database correctly loaded?");
      Warning("HEt","NULL pointer is returned!");
      return NULL;
   }
   if (fHEt[type]) return fHEt[type];

   fHEt[type] = fH2N[type]->ProjectionY(Form("hEt%d%.0f%.0f%.0f",type,
            fInitialMass,fMetallicity*1000,fReviveTime/100),
         1,1,"e");
   fHEt[type]->Reset();
   fHEt[type]->SetStats(0);
   fHEt[type]->GetYaxis()->SetTitle("Average energy [MeV] / second");

   // calculate average
   for (UShort_t iy=1; iy<=fH2N[type]->GetNbinsY(); iy++) {
      Double_t totalE=0, totalN=0;
      for (UShort_t ix=1; ix<=fH2N[type]->GetNbinsX(); ix++) {
         totalN += fH2N[type]->GetBinContent(ix,iy) *
            fH2N[type]->GetXaxis()->GetBinWidth(ix);
         totalE += fH2N[type]->GetBinContent(ix,iy) * 
            fH2N[type]->GetXaxis()->GetBinWidth(ix) *
            fH2N[type]->GetXaxis()->GetBinCenter(ix);
      }
      fHEt[type]->SetBinContent(iy,totalE/totalN);
   }
   return fHEt[type];
}

//______________________________________________________________________________
//

