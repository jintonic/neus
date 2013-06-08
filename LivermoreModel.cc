#include "LivermoreModel.h"

extern "C" {
   void wilson_nl_(Double_t*, Double_t*, Double_t*, Double_t*, Double_t*);
}

#include <TH2D.h>
#include <TSystem.h>

#include <cmath>
using namespace std;

//______________________________________________________________________________
//

NEUS::LivermoreModel::LivermoreModel(const char *name, const char *title) :
   SupernovaModel(name, title)
{
   fMinE= 2.5; // determined by wilson_NL_
   fMaxE=82.5; // no need to go higher
   fMinT= 0.0012; // determined by wilson_NL_
   fMaxT=17.9012; // determined by wilson_NL_
}

//______________________________________________________________________________
//

void NEUS::LivermoreModel::LoadData(const char *dir)
{
   SupernovaModel::LoadData(dir);
   gSystem->Setenv("TOTAL_DATA_DIR", dir);

   Double_t binEdgesx[400]={0};
   Double_t binEdgesy[400]={0};

   Double_t t=fMinT, dt;
   UShort_t nbinsx=0;
   binEdgesx[nbinsx] = fMinT;
   while (t<fMaxT) {
      if (t<0.1) dt = 1e-3;
      else if (t<0.5) dt = 5e-3;
      else if (t<1.0) dt = 1e-2;
      else if (t<4.0) dt = 5e-2;
      else if (t<10.) dt = 1e-1;
      else dt = 0.5;

      nbinsx++;
      binEdgesx[nbinsx] = binEdgesx[nbinsx-1] + dt;
      t+=dt;
   }
   binEdgesx[nbinsx] = fMaxT;

   Double_t e=fMinE, de;
   UShort_t nbinsy=0;
   binEdgesy[nbinsy] = fMinE;
   while (e<fMaxE) {
      if (e<10.) de = 0.5;
      else if (e<40.) de = 1.0;
      else de = 2.0;

      nbinsy++;
      binEdgesy[nbinsy] = binEdgesy[nbinsy-1] + de;
      e+=de;
   }
   binEdgesy[nbinsy] = fMaxE;

   // create histograms
   for (UShort_t i=1; i<=3; i++) {
      fHN2[i] = new TH2D(Form("hN2%s%d", GetName(), i),
            ";time [second];energy [MeV];",
            nbinsx,binEdgesx,nbinsy,binEdgesy);

      fHL2[i] = new TH2D(Form("hL2%s%d", GetName(), i),
            ";time [second];energy [MeV];",
            nbinsx,binEdgesx,nbinsy,binEdgesy);
   }
   for (UShort_t i=4; i<=6; i++) {
      fHN2[i]=fHN2[3];
      fHL2[i]=fHL2[3];
   }

   // fill spectra
   for (UShort_t ix=1; ix<=nbinsx; ix++) {
      for (UShort_t iy=1; iy<=nbinsy; iy++) {
         t = fHN2[1]->GetXaxis()->GetBinLowEdge(ix);
         e = fHN2[1]->GetYaxis()->GetBinLowEdge(iy);

         fHN2[1]->SetBinContent(ix,iy,N2(1,t,e));
         fHN2[2]->SetBinContent(ix,iy,N2(2,t,e));
         fHN2[3]->SetBinContent(ix,iy,N2(3,t,e));

         fHL2[1]->SetBinContent(ix,iy,L2(1,t,e));
         fHL2[2]->SetBinContent(ix,iy,L2(2,t,e));
         fHL2[3]->SetBinContent(ix,iy,L2(3,t,e));
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

void NEUS::LivermoreModel::UseDivariData()
{
   fTotalN[1] = 3.0e7; // * 1e50
   fTotalN[2] = 2.1e7; // * 1e50
   fTotalN[3] = 1.85e7; // * 1e50

   fAverageE[1] = 3.5*3; // MeV
   fAverageE[2] = 5.0*3; // MeV
   fAverageE[3] = 8.0*3; // MeV

   SetName("DivariApproximation");
   SetTitle("Divari approximation");
}

//______________________________________________________________________________
//

void NEUS::LivermoreModel::Clear(Option_t *option)
{
   fTotalN[1] = 0.;
   fTotalN[2] = 0.;
   fTotalN[3] = 0.;

   fAverageE[1] = 1.;
   fAverageE[2] = 1.;
   fAverageE[3] = 1.;

   SetName("LivermoreModel");
   SetTitle("Livermore model");
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::N2(UShort_t type, Double_t time, Double_t energy)
{
   if (type<1 || type>6) {
      Warning("N2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("N2","Return 0!");
      return 0;
   }
   if (time<fMinT || time>fMaxT) {
      Warning("N2","Time is out of range!");
      Warning("N2","Return 0!");
      return 0;
   }
   if (energy<fMinE || energy>fMaxE) {
      Warning("N2","Energy is out of range!");
      Warning("N2","Return 0!");
      return 0;
   }

   Double_t dNL1, dNL2, dNL3;
   wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);

   if (type==1) return dNL1/1e50;
   else if (type==2) return dNL2/1e50;
   else return dNL3/1e50;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::L2(UShort_t type, Double_t time, Double_t energy)
{
   if (type<1 || type>6) {
      Warning("L2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("L2","Return 0!");
      return 0;
   }
   if (time<fMinT || time>fMaxT) {
      Warning("L2","Time is out of range!");
      Warning("L2","Return 0!");
      return 0;
   }
   if (energy<fMinE || energy>fMaxE) {
      Warning("L2","Energy is out of range!");
      Warning("L2","Return 0!");
      return 0;
   }

   Double_t dNL1, dNL2, dNL3;
   wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);

   if (type==1) return dNL1/1e50*energy*1.60217646e-6;
   else if (type==2) return dNL2/1e50*energy*1.60217646e-6;
   else return dNL3/1e50*energy*1.60217646e-6;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Ne(UShort_t type, Double_t energy)
{
   if (type<1 || type>6) {
      Warning("Ne","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Ne","Return 0!");
      return 0;
   }
   if (energy<fMinE || energy>fMaxE) {
      Warning("Ne","Energy is out of range!");
      Warning("Ne","Return 0!");
      return 0;
   }

   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t time = fMinT; // start time [second]
   Double_t dt = 1.e-3; // time interval [second]

   while (time<=fMaxT) {
      if (time<0.1) dt = 1e-3;
      else if (time<0.5) dt = 5e-3;
      else if (time<1.0) dt = 1e-2;
      else if (time<4.0) dt = 5e-2;
      else if (time<10.) dt = 1e-1;
      else dt = 0.5;

      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1   += dNL1*dt;
      n2   += dNL2*dt;
      n3   += dNL3*dt;
      time += dt;
   }

   if (type==1) return n1/1e50;
   else if (type==2) return n2/1e50;
   else return n3/1e50;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Nt(UShort_t type, Double_t time)
{
   if (time<fMinT || time>fMaxT) {
      Warning("Nt","Time is out of range!");
      Warning("Nt","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("Nt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Nt","Return 0!");
      return 0;
   }

   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t energy = fMinE; // start energy
   Double_t dE = 2.; // energy inteval

   while (energy<fMaxE) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1 += dNL1*dE;
      n2 += dNL2*dE;
      n3 += dNL3*dE;
      energy += dE;
   }

   if (type==1) return n1/1e50;
   else if (type==2) return n2/1e50;
   else return n3/1e50;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Nall(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("Nall","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Nall","Return 0!");
      return 0;
   }
   if (fTotalN[type]!=0) return fTotalN[type];

   Double_t dNL1, dNL2, dNL3;

   Double_t energy = fMinE; // start energy [MeV]
   Double_t de = 2; // energy step [MeV]

   Double_t time = fMinT; // start time [second]
   Double_t dt = 1.e-3; // time interval [second]

   while (time<=fMaxT) {
      if (time<0.1) dt = 1e-3;
      else if (time<0.5) dt = 5e-3;
      else if (time<1.0) dt = 1e-2;
      else if (time<4.0) dt = 5e-2;
      else if (time<10.) dt = 1e-1;
      else dt = 0.5;

      energy = fMinE;
      while (energy<fMaxE) {
         wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
         fTotalN[1] += dNL1*dt*de;
         fTotalN[2] += dNL2*dt*de;
         fTotalN[3] += dNL3*dt*de;
         energy+=de;
      }
      time += dt;
   }

   fTotalN[1]/=1e50;
   fTotalN[2]/=1e50;
   fTotalN[3]/=1e50;

   if (type==1) return fTotalN[1];
   else if (type==2) return fTotalN[2];
   else return fTotalN[3];
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Lall(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("Lall","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Lall","Return 0!");
      return 0;
   }
   if (fTotalL[type]!=0) return fTotalL[type];

   Double_t dNL1, dNL2, dNL3;

   Double_t energy = fMinE; // start energy [MeV]
   Double_t de = 2; // energy step [MeV]

   Double_t time = fMinT; // start time [second]
   Double_t dt = 1.e-3; // time interval [second]

   while (time<=fMaxT) {
      if (time<0.1) dt = 1e-3;
      else if (time<0.5) dt = 5e-3;
      else if (time<1.0) dt = 1e-2;
      else if (time<4.0) dt = 5e-2;
      else if (time<10.) dt = 1e-1;
      else dt = 0.5;

      energy = fMinE;
      while (energy<fMaxE) {
         wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
         fTotalL[1] += dNL1*dt*de*energy*1.60217646e-6;
         fTotalL[2] += dNL2*dt*de*energy*1.60217646e-6;
         fTotalL[3] += dNL3*dt*de*energy*1.60217646e-6;
         energy+=de;
      }
      time += dt;
   }

   fTotalL[1]/=1e50;
   fTotalL[2]/=1e50;
   fTotalL[3]/=1e50;

   if (type==1) return fTotalL[1];
   else if (type==2) return fTotalL[2];
   else return fTotalL[3];
}

//______________________________________________________________________________
//

void NEUS::LivermoreModel::Print()
{
   Printf("20 Solar mass, Livermore:     N1=%1.2e, N2=%1.2e, Nx=%1.2e, N=%1.2e, L=%1.2e erg",
         Nall(1)*1e50, Nall(2)*1e50, Nall(3)*1e50, 
         Nall(1)*1e50 + Nall(2)*1e50 + Nall(3)*1e50*4,
         Lall(1)*1e50 + Lall(2)*1e50 + Lall(3)*1e50*4);
}

//______________________________________________________________________________
//

