#include "LivermoreModel.h"

extern "C" {
   void wilson_nl_(Double_t*, Double_t*, Double_t*, Double_t*, Double_t*);
}

#include <TF2.h>
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

   for (UShort_t i=0; i<fgNtype; i++) {
      fN2[i]=NULL;
      fL2[i]=NULL;
      fNe[i]=NULL;
      fNt[i]=NULL;
      fLe[i]=NULL;
      fLt[i]=NULL;
      fEt[i]=NULL;
   }
}

//______________________________________________________________________________
//

NEUS::LivermoreModel::~LivermoreModel()
{
   for (UShort_t i=0; i<fgNtype; i++) {
      if (fN2[i]) delete fN2[i];
      if (fL2[i]) delete fL2[i];
      if (fNe[i]) delete fNe[i];
      if (fNt[i]) delete fNt[i];
      if (fLe[i]) delete fLe[i];
      if (fLt[i]) delete fLt[i];
      if (fEt[i]) delete fEt[i];
   }
}

//______________________________________________________________________________
//

void NEUS::LivermoreModel::SetDataLocation(const char *dir)
{
   SupernovaModel::SetDataLocation(dir);
   gSystem->Setenv("TOTAL_DATA_DIR", dir);
}


//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::N2(UShort_t type, Double_t time, Double_t energy)
{
   Double_t x[2] = {time, energy};
   Double_t p[1] = {type};
   return WilsonN2(x,p);
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Ne(UShort_t type, Double_t energy)
{
   Double_t x[1] = {energy};
   Double_t p[2] = {type, fMaxT};
   return WilsonNe(x,p);
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Nt(UShort_t type, Double_t time)
{
   Double_t x[1] = {time};
   Double_t p[2] = {type};
   return WilsonNt(x,p);
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

   while (time<0.1) { // fine step before 0.1 second
      dt = 1.e-3;
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

   while (time>=0.1 && time<1) {
      dt = 1.e-2;
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

   while (time>=1 && time<10) {
      dt = 1.e-1;
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

   while (time>=10 && time<=fMaxT) { // coarse step after 10 second
      dt = 1; // time interval
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
   fTotalN[4] = fTotalN[3];
   fTotalN[5] = fTotalN[3];
   fTotalN[6] = fTotalN[3];

   if (type==1) return fTotalN[1];
   else if (type==2) return fTotalN[2];
   else return fTotalN[3];
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Eave(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("Eave","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("Eave","Return 0!");
      return 0;
   }
   if (abs(fAverageE[type]-1.0)>0.01) return fAverageE[type];

   fAverageE[type] = HLe(type)->Integral("width")/1.60217646e-6/Nall(type);

   return fAverageE[type];
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonN2(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   Double_t energy = x[1]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   if (time<fMinT || time>fMaxT) {
      Warning("WilsonN2","Time is out of range!");
      Warning("WilsonN2","Return 0!");
      return 0;
   }
   if (energy<fMinE || energy>fMaxE) {
      Warning("WilsonN2","Energy is out of range!");
      Warning("WilsonN2","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("WilsonN2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("WilsonN2","Return 0!");
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

Double_t NEUS::LivermoreModel::WilsonNe(Double_t *x, Double_t *parameter)
{
   Double_t energy = x[0]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino
   Double_t tmax = parameter[1]; // upper limit of time [second]
   if (energy<fMinE || energy>fMaxE) {
      Warning("WilsonNe","Energy is out of range!");
      Warning("WilsonNe","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("WilsonNe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("WilsonNe","Return 0!");
      return 0;
   }
   if (tmax<10) tmax = 10; // do not handle very small upper limit
   if (tmax>fMaxT) tmax=fMaxT; // can not handle very big time input

   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t time = fMinT; // start time [second]

   Double_t dt = 1.e-3; // time interval [second]
   while (time<0.1) { // fine step before 0.1 second
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1   += dNL1*dt;
      n2   += dNL2*dt;
      n3   += dNL3*dt;
      time += dt;
   }

   dt = 1.e-2;
   while (time>=0.1 && time<1) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1   += dNL1*dt;
      n2   += dNL2*dt;
      n3   += dNL3*dt;
      time += dt;
   }

   dt = 1.e-1;
   while (time>=1 && time<10) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1   += dNL1*dt;
      n2   += dNL2*dt;
      n3   += dNL3*dt;
      time += dt;
   }

   dt = 1.;
   while (time>=10 && time<=tmax) { // coarse step after 10 second
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

Double_t NEUS::LivermoreModel::WilsonNt(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   if (time<fMinT || time>fMaxT) {
      Warning("WilsonNt","Time is out of range!");
      Warning("WilsonNt","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("WilsonNt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("WilsonNt","Return 0!");
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

Double_t NEUS::LivermoreModel::WilsonL2(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   Double_t energy = x[1]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   if (time<fMinT || time>fMaxT) {
      Warning("WilsonL2","Time is out of range!");
      Warning("WilsonL2","Return 0!");
      return 0;
   }
   if (energy<fMinE || energy>fMaxE) {
      Warning("WilsonL2","Energy is out of range!");
      Warning("WilsonL2","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("WilsonL2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("WilsonL2","Return 0!");
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

Double_t NEUS::LivermoreModel::WilsonLe(Double_t *x, Double_t *parameter)
{
   Double_t energy = x[0]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino
   Double_t tmax = parameter[1]; // upper limit of time [second]
   if (energy<fMinE || energy>fMaxE) {
      Warning("WilsonLe","Energy is out of range!");
      Warning("WilsonLe","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("WilsonLe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("WilsonLe","Return 0!");
      return 0;
   }
   if (tmax<10) tmax = 10; // do not handle very small upper limit
   if (tmax>fMaxT) tmax=fMaxT; // can not handle very big time input

   Double_t l1=0, l2=0, l3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t time = 1.2e-3; // start time

   Double_t dt = 1.e-3; // time interval [second]
   while (time<0.1) { // fine step before 0.1 second
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      l1   += dNL1*dt*energy*1.60217646e-6;
      l2   += dNL2*dt*energy*1.60217646e-6;
      l3   += dNL3*dt*energy*1.60217646e-6;
      time += dt;
   }

   dt = 1.e-2;
   while (time>=0.1 && time<1) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      l1   += dNL1*dt*energy*1.60217646e-6;
      l2   += dNL2*dt*energy*1.60217646e-6;
      l3   += dNL3*dt*energy*1.60217646e-6;
      time += dt;
   }

   dt = 1.e-1;
   while (time>=1 && time<10) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      l1   += dNL1*dt*energy*1.60217646e-6;
      l2   += dNL2*dt*energy*1.60217646e-6;
      l3   += dNL3*dt*energy*1.60217646e-6;
      time += dt;
   }

   dt = 1.;
   while (time>=10 && time<=tmax) { // coarse step after 10 second
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      l1   += dNL1*dt*energy*1.60217646e-6;
      l2   += dNL2*dt*energy*1.60217646e-6;
      l3   += dNL3*dt*energy*1.60217646e-6;
      time += dt;
   }

   if (type==1) return l1/1e50;
   else if (type==2) return l2/1e50;
   else return l3/1e50;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonLt(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   if (time<fMinT || time>fMaxT) {
      Warning("WilsonLt","Time is out of range!");
      Warning("WilsonLt","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("WilsonLt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("WilsonLt","Return 0!");
      return 0;
   }

   Double_t l1=0, l2=0, l3=0;
   Double_t dNL1, dNL2, dNL3;
   Double_t energy = fMinE; // start energy
   Double_t dE = 2.; // energy inteval
   while (energy<fMaxE) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      l1 += dNL1*dE*energy*1.60217646e-6;
      l2 += dNL2*dE*energy*1.60217646e-6;
      l3 += dNL3*dE*energy*1.60217646e-6;
      energy += dE;
   }

   if (type==1) return l1/1e50;
   else if (type==2) return l2/1e50;
   else return l3/1e50;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonEt(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   if (time<fMinT || time>fMaxT) {
      Warning("WilsonEt","Time is out of range!");
      Warning("WilsonEt","Return 0!");
      return 0;
   }
   if (type<1 || type>6) {
      Warning("WilsonEt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("WilsonEt","Return 0!");
      return 0;
   }

   Double_t n1=0, n2=0, n3=0;
   Double_t e1=0, e2=0, e3=0;
   Double_t dNL1, dNL2, dNL3;
   Double_t energy = fMinE; // start energy
   Double_t dE = 2.; // energy inteval
   while (energy<fMaxE) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1 += dNL1*dE;
      n2 += dNL2*dE;
      n3 += dNL3*dE;
      e1 += dNL1*dE*energy;
      e2 += dNL2*dE*energy;
      e3 += dNL3*dE*energy;
      energy += dE;
   }

   if (type==1) return e1/n1;
   else if (type==2) return e2/n2;
   else return e3/n3;
}

//______________________________________________________________________________
//

TF2* NEUS::LivermoreModel::FN2(UShort_t type)
{
   if (fN2[type]) return fN2[type];

   fN2[type] = new TF2(Form("fN2%d",type), this, &LivermoreModel::WilsonN2,
         fMinT, fMaxT, fMinE, fMaxE, 1, "LivermoreModel", "FN2");
   fN2[type]->SetParameter(0,type);
   fN2[type]->SetNpx(179);
   fN2[type]->SetNpy(80);
   fN2[type]->GetXaxis()->SetTitle("time [second]");
   fN2[type]->GetYaxis()->SetTitle("neutrino energy [MeV]");
   if (type==1) {
      fN2[type]->SetTitle("number of #nu_{e} [10^{50}/second/MeV]");
      fN2[type]->SetLineColor(kBlack);
   } else if (type==2) {
      fN2[type]->SetTitle("number of #bar{#nu}_{e} [10^{50}/second/MeV]");
      fN2[type]->SetLineColor(kRed);
   } else {
      fN2[type]->SetTitle("number of #nu_{x} [10^{50}/second/MeV]");
      fN2[type]->SetLineColor(kBlue);
   }
   return fN2[type];
}

//______________________________________________________________________________
//

TF2* NEUS::LivermoreModel::FL2(UShort_t type)
{
   if (fL2[type]) return fL2[type];

   fL2[type] = new TF2(Form("fL2%d",type), this, &LivermoreModel::WilsonL2,
         fMinT, fMaxT, fMinE, fMaxE, 1, "LivermoreModel", "FL2");
   fL2[type]->SetParameter(0,type);
   fL2[type]->SetNpx(179);
   fL2[type]->SetNpy(80);
   fL2[type]->GetXaxis()->SetTitle("time [second]");
   fL2[type]->GetYaxis()->SetTitle("neutrino energy [MeV]");
   if (type==1) {
      fL2[type]->SetTitle("luminosity of #nu_{e} [10^{50}/second/erg]");
      fL2[type]->SetLineColor(kBlack);
   } else if (type==2) {
      fL2[type]->SetTitle("luminosity of #bar{#nu}_{e} [10^{50}/second/erg]");
      fL2[type]->SetLineColor(kRed);
   } else {
      fN2[type]->SetTitle("luminosity of #nu_{x} [10^{50}/second/erg]");
      fL2[type]->SetLineColor(kBlue);
   }
   return fL2[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FNe(UShort_t type, Double_t tmax)
{
   if (fNe[type]) return fNe[type];

   fNe[type] = new TF1(Form("fNe%d",type), this, &LivermoreModel::WilsonNe,
         fMinE, fMaxE, 2, "LivermoreModel", "FNe");
   fNe[type]->SetParameter(0,type);
   fNe[type]->SetParameter(1,tmax);

   fNe[type]->SetNpx(80);
   fNe[type]->GetXaxis()->SetTitle("neutrino energy [MeV]");
   if (type==1) {
      fNe[type]->GetYaxis()->SetTitle("number of #nu_{e} [10^{50}/MeV]");
      fNe[type]->SetLineColor(kBlack);
   } else if (type==2) {
      fNe[type]->GetYaxis()
         ->SetTitle("number of #bar{#nu}_{e} [10^{50}/MeV]");
      fNe[type]->SetLineColor(kRed);
   }else {
      fNe[type]->GetYaxis()->SetTitle("number of #nu_{x} [10^{50}/MeV]");
      fNe[type]->SetLineColor(kBlue);
   }
   return fNe[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FLe(UShort_t type, Double_t tmax)
{
   if (fLe[type]) return fLe[type];

   fLe[type] = new TF1(Form("fLe%d",type), this, &LivermoreModel::WilsonLe,
         fMinE, fMaxE, 2, "LivermoreModel", "FLe");
   fLe[type]->SetParameter(0,type);
   fLe[type]->SetParameter(1,tmax);

   fLe[type]->SetNpx(80);
   fLe[type]->GetXaxis()->SetTitle("neutrino energy [MeV]");
   if (type==1) {
      fLe[type]->GetYaxis()
      ->SetTitle("luminosity of #nu_{e} [10^{50} erg/MeV]");
      fLe[type]->SetLineColor(kBlack);
   } else if (type==2) {
      fLe[type]->GetYaxis()
      ->SetTitle("luminosity of #bar{#nu}_{e} [10^{50} erg/MeV]");
      fLe[type]->SetLineColor(kRed);
   } else {
      fLe[type]->GetYaxis()
      ->SetTitle("luminosity of #nu_{x} [10^{50} erg/MeV]");
      fLe[type]->SetLineColor(kBlue);
   }
   return fLe[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FNt(UShort_t type)
{
   if (fNt[type]) return fNt[type];

   fNt[type] = new TF1(Form("fNt%d",type), this, &LivermoreModel::WilsonNt,
         fMinT, fMaxT, 1, "LivermoreModel", "FNt");
   fNt[type]->SetParameter(0,type);

   fNt[type]->SetNpx(179);
   fNt[type]->GetXaxis()->SetTitle("time [second]");
   if (type==1) {
      fNt[type]->GetYaxis()->SetTitle("number of #nu_{e} [10^{50}/second]");
      fNt[type]->SetLineColor(kBlack);
   } else if (type==2) {
      fNt[type]->GetYaxis()
      ->SetTitle("number of #bar{#nu}_{e} [10^{50}/second]");
      fNt[type]->SetLineColor(kRed);
   } else {
      fNt[type]->GetYaxis()->SetTitle("number of #nu_{x} [10^{50}/second]");
      fNt[type]->SetLineColor(kBlue);
   }
   return fNt[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FLt(UShort_t type)
{
   if (fLt[type]) return fLt[type];

   fLt[type] = new TF1(Form("fLt%d",type), this, &LivermoreModel::WilsonLt,
         fMinT, fMaxT, 1, "LivermoreModel", "FLt");
   fLt[type]->SetParameter(0,type);

   fLt[type]->SetNpx(179);
   fLt[type]->GetXaxis()->SetTitle("time [second]");
   if (type==1) {
      fLt[type]->GetYaxis()
      ->SetTitle("luminosity of #nu_{e} [10^{50} erg/second]");
      fLt[type]->SetLineColor(kBlack);
   } else if (type==2) {
      fLt[type]->GetYaxis()
      ->SetTitle("luminosity of #bar{#nu}_{e} [10^{50} erg/second]");
      fLt[type]->SetLineColor(kRed);
   } else {
      fLt[type]->GetYaxis()
      ->SetTitle("luminosity of #nu_{x} [10^{50} erg/second]");
      fLt[type]->SetLineColor(kBlue);
   }
   return fLt[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FEt(UShort_t type)
{
   if (fEt[type]) return fEt[type];

   fEt[type] = new TF1(Form("fEt%d",type), this, &LivermoreModel::WilsonEt,
         fMinT, fMaxT, 1, "LivermoreModel", "FEt");
   fEt[type]->SetParameter(0,type);

   fEt[type]->SetNpx(179);
   fEt[type]->GetXaxis()->SetTitle("time [second]");
   if (type==1) {
      fEt[type]->GetYaxis()->SetTitle("<E> of #nu_{e} [MeV/second]");
      fEt[type]->SetLineColor(kBlack);
   } else if (type==2) {
      fEt[type]->GetYaxis()->SetTitle("<E> of #bar{#nu}_{e} [MeV/second]");
      fEt[type]->SetLineColor(kRed);
   } else {
      fEt[type]->GetYaxis()->SetTitle("<E> of #nu_{x} [MeV/second]");
      fEt[type]->SetLineColor(kBlue);
   }
   return fEt[type];
}

//______________________________________________________________________________
//

TH2D* NEUS::LivermoreModel::HN2(UShort_t type)
{
   return (TH2D*) FN2(type)->GetHistogram();
}

//______________________________________________________________________________
//

TH1D* NEUS::LivermoreModel::HNe(UShort_t type, Double_t tmax)
{
   return (TH1D*) FNe(type, tmax)->GetHistogram();
}

//______________________________________________________________________________
//

TH1D* NEUS::LivermoreModel::HNt(UShort_t type)
{
   return (TH1D*) FNt(type)->GetHistogram();
}

//______________________________________________________________________________
//

TH2D* NEUS::LivermoreModel::HL2(UShort_t type)
{
   return (TH2D*) FL2(type)->GetHistogram();
}

//______________________________________________________________________________
//

TH1D* NEUS::LivermoreModel::HLe(UShort_t type, Double_t tmax)
{
   return (TH1D*) FLe(type, tmax)->GetHistogram();
}

//______________________________________________________________________________
//

TH1D* NEUS::LivermoreModel::HLt(UShort_t type)
{
   return (TH1D*) FLt(type)->GetHistogram();
}

//______________________________________________________________________________
//

TH1D* NEUS::LivermoreModel::HEt(UShort_t type)
{
   return (TH1D*) FEt(type)->GetHistogram();
}

//______________________________________________________________________________
//

void NEUS::LivermoreModel::Print()
{
   Printf("20 Solar mass, Livermore:     n1=%1.2e, n2=%1.2e, nx=%1.2e, total=%1.2e",
         Nall(1)*1e50, Nall(2)*1e50, Nall(3)*1e50, 
         Nall(1)*1e50 + Nall(2)*1e50 + Nall(3)*1e50*4);
}
