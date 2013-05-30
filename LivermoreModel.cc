#include "LivermoreModel.h"

extern "C" {
   void wilson_nl_(Double_t*, Double_t*, Double_t*, Double_t*, Double_t*);
}

#include <TF2.h>
#include <TH2D.h>
#include <TSystem.h>

//______________________________________________________________________________
//

NEUS::LivermoreModel::LivermoreModel(const char *name, const char *title) :
   SupernovaModel(name, title)
{
   for (UShort_t i=0; i<7; i++) {
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
   for (UShort_t i=0; i<7; i++) {
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

Double_t NEUS::LivermoreModel::N2(UShort_t type, Double_t energy, Double_t time)
{
   Double_t dNL1, dNL2, dNL3;
   wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);

   if (type==1) return dNL1;
   else if (type==2) return dNL2;
   else return dNL3;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Ne(UShort_t type, Double_t energy)
{
   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t time = 1.2e-3; // start time [second]
   Double_t tmax=17.9012; // end time [second]

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
   while (time>=10 && time<tmax) { // coarse step after 10 second
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1   += dNL1*dt;
      n2   += dNL2*dt;
      n3   += dNL3*dt;
      time += dt;
   }

   if (type==1) return n1;
   else if (type==2) return n2;
   else return n3;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Nt(UShort_t type, Double_t time)
{
   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t energy = 2.5; // start energy [MeV]
   Double_t de = 2; // energy step [MeV]
   Double_t emax=80; // max energy [MeV]

   while (energy<emax) {
      wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
      n1 += dNL1*de;
      n2 += dNL2*de;
      n3 += dNL3*de;
      energy+=de;
   }
   if (type==1) return n1;
   else if (type==2) return n2;
   else return n3;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::Nall(UShort_t type)
{
   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t energy = 2.5; // start energy [MeV]
   Double_t de = 2; // energy step [MeV]
   Double_t emax=80; // max energy [MeV]

   Double_t time = 1.2e-3; // start time [second]
   Double_t dt = 1.e-3; // time interval [second]
   Double_t tmax=17.9012; // end time [second]

   while (time<0.1) { // fine step before 0.1 second
      dt = 1.e-3;
      energy = 2.5;
      while (energy<emax) {
         wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
         n1 += dNL1*dt*de;
         n2 += dNL2*dt*de;
         n3 += dNL3*dt*de;
         energy+=de;
      }
      time += dt;
   }

   while (time>=0.1 && time<1) {
      dt = 1.e-2;
      energy = 2.5;
      while (energy<emax) {
         wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
         n1 += dNL1*dt*de;
         n2 += dNL2*dt*de;
         n3 += dNL3*dt*de;
         energy+=de;
      }
      time += dt;
   }

   while (time>=1 && time<10) {
      dt = 1.e-1;
      energy = 2.5;
      while (energy<emax) {
         wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
         n1 += dNL1*dt*de;
         n2 += dNL2*dt*de;
         n3 += dNL3*dt*de;
         energy+=de;
      }
      time += dt;
   }

   while (time>=10 && time<tmax) { // coarse step after 10 second
      dt = 1; // time interval
      energy = 2.5;
      while (energy<emax) {
         wilson_nl_(&time, &energy, &dNL1, &dNL2, &dNL3);
         n1 += dNL1*dt*de;
         n2 += dNL2*dt*de;
         n3 += dNL3*dt*de;
         energy+=de;
      }
      time += dt;
   }

   if (type==1) return n1;
   else if (type==2) return n2;
   else return n3;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonN2(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   Double_t Ev = x[1]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   Double_t dNL1, dNL2, dNL3;
   wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);

   if (type==1) return dNL1/1e50;
   else if (type==2) return dNL2/1e50;
   else return dNL3/1e50;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonNe(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino
   Double_t tmax = parameter[1]; // upper limit of time [second]

   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;

   Double_t time = 1.2e-3; // start time [second]
   Double_t dt = 1.e-3; // time inteval
   if (tmax<10) {
      while (time<=tmax) { // fine step before 10 second
         wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);
         n1   += dNL1*dt;
         n2   += dNL2*dt;
         n3   += dNL3*dt;
         time += dt;
      }
   } else {
      if (tmax>17.9012) tmax=17.9012; // end time [second]
      while (time<=10) { // fine step before 10 second
         wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);
         n1   += dNL1*dt;
         n2   += dNL2*dt;
         n3   += dNL3*dt;
         time += dt;
      }
      dt = 1; // time inteval
      while (time>=10 && time<tmax) { // coarse step after 10 second
         wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);
         n1   += dNL1*dt;
         n2   += dNL2*dt;
         n3   += dNL3*dt;
         time += dt;
      }
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
   Double_t maxEv= parameter[1]; // upper limit of neutrino energy [MeV]
   if (maxEv>100.) maxEv=100.;

   Double_t n1=0, n2=0, n3=0;
   Double_t dNL1, dNL2, dNL3;
   Double_t Ev = 2.5; // start energy
   Double_t dE = 1.; // energy inteval
   while (Ev<maxEv) {
      wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);
      n1 += dNL1*dE;
      n2 += dNL2*dE;
      n3 += dNL3*dE;
      Ev += dE;
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
   Double_t Ev = x[1]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino

   Double_t dNL1, dNL2, dNL3;
   wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);

   if (type==1) return dNL1*Ev*1.60217646e-6;
   else if (type==2) return dNL2*Ev*1.60217646e-6;
   else return dNL3*Ev*1.60217646e-6;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonLe(Double_t *x, Double_t *parameter)
{
   Double_t Ev = x[0]; // neutrino energy [MeV]
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino
   Double_t tmax = parameter[1]; // upper limit of time [second]
   if (tmax>17.9012) tmax=17.9012;

   Double_t l1=0, l2=0, l3=0;
   Double_t dNL1, dNL2, dNL3;
   Double_t time = 1.2e-3; // start time
   Double_t dt = 1.e-3; // time inteval
   while (time<tmax) {
      wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);
      l1   += dNL1*dt*Ev*1.60217646e-6;
      l2   += dNL2*dt*Ev*1.60217646e-6;
      l3   += dNL3*dt*Ev*1.60217646e-6;
      time += dt;
   }

   if (type==1) return l1;
   else if (type==2) return l2;
   else return l3;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonLt(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino
   Double_t maxEv= parameter[1]; // upper limit of neutrino energy [MeV]
   if (maxEv>100.) maxEv=100.;

   Double_t l1=0, l2=0, l3=0;
   Double_t dNL1, dNL2, dNL3;
   Double_t Ev = 2.5; // start energy
   Double_t dE = 1.; // energy inteval
   while (Ev<maxEv) {
      wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);
      l1 += dNL1*dE*Ev*1.60217646e-6;
      l2 += dNL2*dE*Ev*1.60217646e-6;
      l3 += dNL3*dE*Ev*1.60217646e-6;
      Ev += dE;
   }

   if (type==1) return l1;
   else if (type==2) return l2;
   else return l3;
}

//______________________________________________________________________________
//

Double_t NEUS::LivermoreModel::WilsonEt(Double_t *x, Double_t *parameter)
{
   Double_t time = x[0]; // second
   UShort_t type = static_cast<UShort_t>(parameter[0]); // type of neutrino
   Double_t maxEv= parameter[1]; // upper limit of neutrino energy [MeV]
   if (maxEv>100.) maxEv=100.;

   Double_t n1=0, n2=0, n3=0;
   Double_t e1=0, e2=0, e3=0;
   Double_t dNL1, dNL2, dNL3;
   Double_t Ev = 2.5; // start energy
   Double_t dE = 1.; // energy inteval
   while (Ev<maxEv) {
      wilson_nl_(&time, &Ev, &dNL1, &dNL2, &dNL3);
      n1 += dNL1;
      n2 += dNL2;
      n3 += dNL3;
      e1 += dNL1*dE*Ev;
      e2 += dNL2*dE*Ev;
      e3 += dNL3*dE*Ev;
      Ev += dE;
   }

   if (type==1) return e1/n1;
   else if (type==2) return e2/n2;
   else return e3/n3;
}

//______________________________________________________________________________
//

TF2* NEUS::LivermoreModel::FN2(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("FN2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("FN2","NULL pointer is returned!");
      return NULL;
   }
   if (fN2[type]) return fN2[type];

   fN2[type] = new TF2(Form("fN2%d",type), this, &LivermoreModel::WilsonN2,
         1.2e-3, 17.9012, 2.5, 62.5, 1, "LivermoreModel", "FN2");
   fN2[type]->SetParameter(0,type);
   fN2[type]->SetNpx(179);
   fN2[type]->SetNpy(60);
   fN2[type]->SetMaximum(1e10);
   fN2[type]->SetMinimum(1e-10);
   fN2[type]->GetXaxis()->SetTitle("time [second]");
   fN2[type]->GetYaxis()->SetTitle("neutrino energy [MeV]");
   if (type==1)
      fN2[type]->SetTitle("number of #nu_{e}/10^{50} [/second/MeV]");
   return fN2[type];
}

//______________________________________________________________________________
//

TF2* NEUS::LivermoreModel::FL2(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("FL2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("FL2","NULL pointer is returned!");
      return NULL;
   }
   if (fL2[type]) return fL2[type];

   fL2[type] = new TF2(Form("fL2%d",type), this, &LivermoreModel::WilsonL2,
         1.2e-3, 18, 2.5, 100, 1, "LivermoreModel", "FL2");
   fL2[type]->SetParameter(0,type);
   return fL2[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FNe(UShort_t type, Double_t tmax)
{
   if (type<1 || type>6) {
      Warning("FNe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("FNe","NULL pointer is returned!");
      return NULL;
   }
   if (fNe[type]) return fNe[type];

   fNe[type] = new TF1(Form("fNe%d",type), this, &LivermoreModel::WilsonNe,
         2.5, 100, 2, "LivermoreModel", "FNe");
   fNe[type]->SetParameter(0,type);
   fNe[type]->SetParameter(1,tmax);

   return fNe[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FLe(UShort_t type, Double_t tmax)
{
   if (type<1 || type>6) {
      Warning("FLe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("FLe","NULL pointer is returned!");
      return NULL;
   }
   if (fLe[type]) return fLe[type];

   fLe[type] = new TF1(Form("fLe%d",type), this, &LivermoreModel::WilsonLe,
         2.5, 100, 2, "LivermoreModel", "FLe");
   fLe[type]->SetParameter(0,type);
   fLe[type]->SetParameter(1,tmax);

   return fLe[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FNt(UShort_t type, Double_t maxEv)
{
   if (type<1 || type>6) {
      Warning("FNt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("FNt","NULL pointer is returned!");
      return NULL;
   }
   if (fNt[type]) return fNt[type];

   fNt[type] = new TF1(Form("fNt%d",type), this, &LivermoreModel::WilsonNt,
         1.2e-3, 18., 2, "LivermoreModel", "FNt");
   fNt[type]->SetParameter(0,type);
   fNt[type]->SetParameter(1,maxEv);

   return fNt[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FLt(UShort_t type, Double_t maxEv)
{
   if (type<1 || type>6) {
      Warning("FLt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("FLt","NULL pointer is returned!");
      return NULL;
   }
   if (fLt[type]) return fLt[type];

   fLt[type] = new TF1(Form("fLt%d",type), this, &LivermoreModel::WilsonLt,
         1.2e-3, 18., 2, "LivermoreModel", "FLt");
   fLt[type]->SetParameter(0,type);
   fLt[type]->SetParameter(1,maxEv);

   return fLt[type];
}

//______________________________________________________________________________
//

TF1* NEUS::LivermoreModel::FEt(UShort_t type, Double_t maxEv)
{
   if (type<1 || type>6) {
      Warning("HEt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HEt","NULL pointer is returned!");
      return NULL;
   }
   if (fEt[type]) return fEt[type];

   fEt[type] = new TF1(Form("fEt%d",type), this, &LivermoreModel::WilsonEt,
         1.2e-3, 18., 2, "LivermoreModel", "FEt");
   fEt[type]->SetParameter(0,type);
   fEt[type]->SetParameter(1,maxEv);

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

TH1D* NEUS::LivermoreModel::HNt(UShort_t type, Double_t maxEv)
{
   return (TH1D*) FNt(type,maxEv)->GetHistogram();
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

TH1D* NEUS::LivermoreModel::HLt(UShort_t type, Double_t maxEv)
{
   return (TH1D*) FLt(type, maxEv)->GetHistogram();
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
         Nall(1), Nall(2), Nall(3), Nall(1) + Nall(2) + Nall(3)*4);
}
