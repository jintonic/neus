#include "SupernovaModel.h"

#include <TF1.h>
#include <TH2D.h>
#include <TAxis.h>

#include <cmath>
using namespace std;

//______________________________________________________________________________
//

NEUS::SupernovaModel::SupernovaModel() : TNamed(), fDataLocation(),
   fMinE(0), fMaxE(0), fMinT(0), fMaxT(0)
{
   for (UShort_t i=0; i<fgNtype; i++) {
      fTotalN[i] = 0;
      fTotalL[i] = 0;
      fAverageE[i] = 1.0; // 0 cannot be in the denominator
      fHN2[i] = 0;
      fHL2[i] = 0;
      fHNe[i] = 0;
      fHNt[i] = 0;
      fHLe[i] = 0;
      fHLt[i] = 0;
      fHEt[i] = 0;
      fNeFD[i]= 0;
   }
}

//______________________________________________________________________________
//

NEUS::SupernovaModel::SupernovaModel(const char *name, const char *title) : 
   TNamed(name, title), fDataLocation(), fMinE(0), fMaxE(0), fMinT(0), fMaxT(0)
{
   for (UShort_t i=0; i<fgNtype; i++) {
      fTotalN[i] = 0;
      fTotalL[i] = 0;
      fAverageE[i] = 1.0;
      fHN2[i] = 0;
      fHL2[i] = 0;
      fHNe[i] = 0;
      fHNt[i] = 0;
      fHLe[i] = 0;
      fHLt[i] = 0;
      fHEt[i] = 0;
      fNeFD[i]= 0;
   }
}

//______________________________________________________________________________
//

void NEUS::SupernovaModel::Clear(Option_t *option)
{
   for (UShort_t i=0; i<fgNtype; i++) {
      fTotalN[i] = 0;
      fTotalL[i] = 0;
      fAverageE[i] = 1.; // 0 cannot be in the denominator
      if (fNeFD[i]) delete fNeFD[i];
      if (fHN2[i]) delete fHN2[i];
      if (fHL2[i]) delete fHL2[i];
      if (fHNe[i]) delete fHNe[i];
      if (fHNt[i]) delete fHNt[i];
      if (fHLe[i]) delete fHLe[i];
      if (fHLt[i]) delete fHLt[i];
      if (fHEt[i]) delete fHEt[i];
      fHN2[i] = 0;
      fHL2[i] = 0;
      fHNe[i] = 0;
      fHNt[i] = 0;
      fHLe[i] = 0;
      fHLt[i] = 0;
      fHEt[i] = 0;
      fNeFD[i] = 0;
   }
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::NeFD(UShort_t type, Double_t energy)
{
   Double_t co = 0.55 * Nall(type); // coefficient
   Double_t kT = Eave(type)*2/6.; // <E> = NDF*kT/2 => kT = <E>*2/NDF
   return co/kT/kT/kT * energy*energy/(1.+exp(energy/kT));
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::NeFermiDirac(Double_t *x, Double_t *parameter)
{
   Double_t energy = x[0];
   UShort_t type = static_cast<UShort_t>(parameter[0]);
   return NeFD(type, energy);
}

//______________________________________________________________________________
//

TF1* NEUS::SupernovaModel::FNeFD(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("FNeFD","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("FNeFD","NULL pointer is returned!");
      return 0;
   }
   if (fNeFD[type]) return fNeFD[type];

   fNeFD[type] = new TF1(Form("fNeFD%d",type), this,
         &SupernovaModel::NeFermiDirac, fMinE, fMaxE, 1);
   fNeFD[type]->SetParameter(0,type);

   // fHistogram in TF1 has not yet been created at this moment,
   // preperties should be saved in TF1 instead of fHistogram.
   // That is why fNeFD[type]->SetTitle() is called instead of
   // fNeFD[type]->GetXaxis()->SetTitle(), 
   // which returns the x axis of fHistogram.
   if (type==1) {
      fNeFD[type]->SetTitle(
            ";neutrino energy [MeV];number of #nu_{e} [10^{50}/second];");
      fNeFD[type]->SetLineColor(kBlack);
      fNeFD[type]->SetLineStyle(kDashed);
   } else if (type==2) {
      fNeFD[type]->SetTitle(
            ";neutrino energy [MeV];number of #bar{#nu}_{e} [10^{50}/second];");
      fNeFD[type]->SetLineColor(kRed);
      fNeFD[type]->SetLineStyle(kDashed);
   } else {
      fNeFD[type]->SetTitle(
            ";neutrino energy [MeV];number of #nu_{x} [10^{50}/second];");
      fNeFD[type]->SetLineColor(kBlue);
      fNeFD[type]->SetLineStyle(kDashed);
   }
   return fNeFD[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::SupernovaModel::HNeFD(UShort_t type)
{
   return (TH1D*) FNeFD(type)->GetHistogram();
}

//______________________________________________________________________________
//

TH2D* NEUS::SupernovaModel::HN2(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("HN2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HN2","NULL pointer is returned!");
      return 0;
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

TH2D* NEUS::SupernovaModel::HL2(UShort_t type)
{
   if (type<1 || type>6) {
      Warning("HL2","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HL2","NULL pointer is returned!");
      return 0;
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

TH1D* NEUS::SupernovaModel::HNe(UShort_t type, Double_t tmax)
{
   if (type<1 || type>6) {
      Warning("HNe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HNe","NULL pointer is returned!");
      return 0;
   }
   if (tmax>fMaxT) tmax=fMaxT;

   TString name = Form("hNe-%s-%d-%.4f", GetName(), type, tmax);
   if (fHNe[type]) {
      if (name.CompareTo(fHNe[type]->GetName())==0) return fHNe[type];
      else fHNe[type]->Reset();
   } else {
      fHNe[type] = new TH1D(name.Data(),
            ";energy [MeV];number of neutrinos [10^{50}/MeV]",
            HN2(type)->GetNbinsY(),
            HN2(type)->GetYaxis()->GetXbins()->GetArray());
      fHNe[type]->SetStats(0);
      fHNe[type]->SetLineColor(HN2(type)->GetLineColor());
      fHNe[type]->SetTitle(GetTitle());
   }

   // calculate integral in [0, tmax]
   for (UShort_t iy=1; iy<=HN2(type)->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=HN2(type)->GetNbinsX(); ix++) {
         if (tmax<HN2(type)->GetXaxis()->GetBinCenter(ix)) break;
         content += HN2(type)->GetBinContent(ix,iy) *
            HN2(type)->GetXaxis()->GetBinWidth(ix);
      }
      fHNe[type]->SetBinContent(iy,content);
   }
   return fHNe[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::SupernovaModel::HLe(UShort_t type, Double_t tmax)
{
   if (type<1 || type>6) {
      Warning("HLe","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HLe","NULL pointer is returned!");
      return 0;
   }
   if (tmax>fMaxT) tmax=fMaxT;

   TString name = Form("hLe-%s-%d-%.4f", GetName(), type, tmax);
   if (fHLe[type]) {
      if (name.CompareTo(fHLe[type]->GetName())==0) return fHLe[type];
      else fHLe[type]->Reset();
   } else {
      fHLe[type] = new TH1D(name.Data(),
            ";energy [MeV];luminosity [10^{50} erg/MeV]",
            HL2(type)->GetNbinsY(),
            HL2(type)->GetYaxis()->GetXbins()->GetArray());
      fHLe[type]->SetStats(0);
      fHLe[type]->SetLineColor(HL2(type)->GetLineColor());
      fHLe[type]->SetTitle(GetTitle());
   }

   // calculate integral in [0, tmax]
   for (UShort_t iy=1; iy<=HL2(type)->GetNbinsY(); iy++) {
      Double_t content=0;
      for (UShort_t ix=1; ix<=HL2(type)->GetNbinsX(); ix++) {
         if (tmax<HL2(type)->GetXaxis()->GetBinCenter(ix)) break;
         content += HL2(type)->GetBinContent(ix,iy) *
            HL2(type)->GetXaxis()->GetBinWidth(ix);
      }
      fHLe[type]->SetBinContent(iy,content);
   }
   return fHLe[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::SupernovaModel::HNt(UShort_t type, Double_t emax)
{
   if (type<1 || type>6) {
      Warning("HNt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HNt","NULL pointer is returned!");
      return 0;
   }
   if (emax>fMaxE) emax=fMaxE;

   TString name = Form("hNt-%s-%d-%.1f", GetName(), type, emax);
   if (fHNt[type]) {
     if (name.CompareTo(fHNt[type]->GetName())==0) return fHNt[type];
     else fHNt[type]->Reset();
   } else {
      Info("HNt", "Create %s",name.Data());
      fHNt[type] = new TH1D(name.Data(),
            ";time [second];number of neutrinos [10^{50}/second]",
            HN2(type)->GetNbinsX(),
            HN2(type)->GetXaxis()->GetXbins()->GetArray());
      fHNt[type]->SetStats(0);
      fHNt[type]->SetLineColor(HN2(type)->GetLineColor());
      fHNt[type]->SetTitle(GetTitle());
   }

   // calculate integral
   for (UShort_t ix=1; ix<=HN2(type)->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=HN2(type)->GetNbinsY(); iy++) {
         if (emax<HN2(type)->GetYaxis()->GetBinLowEdge(iy)) break;
         content += HN2(type)->GetBinContent(ix,iy) *
            HN2(type)->GetYaxis()->GetBinWidth(iy);
      }
      fHNt[type]->SetBinContent(ix,content);
   }
   return fHNt[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::SupernovaModel::HLt(UShort_t type, Double_t emax)
{
   if (type<1 || type>6) {
      Warning("HLt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HLt","NULL pointer is returned!");
      return 0;
   }
   if (emax>fMaxE) emax=fMaxE;

   TString name = Form("hLt-%s-%d-%.1f", GetName(), type, emax);
   if (fHLt[type]) {
      if (name.CompareTo(fHLt[type]->GetName())==0) return fHLt[type];
      else fHLt[type]->Reset();
   } else {
      Info("HLt", "Create %s",name.Data());
      fHLt[type] = new TH1D(name.Data(),
            ";time [second];luminosity [10^{50} erg/second]",
            HL2(type)->GetNbinsX(),
            HL2(type)->GetXaxis()->GetXbins()->GetArray());
      fHLt[type]->SetStats(0);
      fHLt[type]->SetLineColor(HL2(type)->GetLineColor());
      fHLt[type]->SetTitle(GetTitle());
   }

   // calculate integral
   for (UShort_t ix=1; ix<=HL2(type)->GetNbinsX(); ix++) {
      Double_t content=0;
      for (UShort_t iy=1; iy<=HL2(type)->GetNbinsY(); iy++) {
         if (emax<HL2(type)->GetYaxis()->GetBinLowEdge(iy)) break;
         content += HL2(type)->GetBinContent(ix,iy) *
            HL2(type)->GetYaxis()->GetBinWidth(iy);
      }
      fHLt[type]->SetBinContent(ix,content);
   }
   return fHLt[type];
}

//______________________________________________________________________________
//

TH1D* NEUS::SupernovaModel::HEt(UShort_t type, Double_t emax)
{
   if (type<1 || type>6) {
      Warning("HEt","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HEt","NULL pointer is returned!");
      return 0;
   }
   if (emax>fMaxE) emax=fMaxE;

   TString name = Form("hEt-%s-%d-%.1f", GetName(), type, emax);
   if (fHEt[type]) {
      if (name.CompareTo(fHEt[type]->GetName())==0) return fHEt[type];
      else fHEt[type]->Reset();
   } else {
      Info("HEt", "Create %s",name.Data());
      fHEt[type] = new TH1D(name.Data(),
            ";time [second];average energy [MeV/second]",
            HN2(type)->GetNbinsX(),
            HN2(type)->GetXaxis()->GetXbins()->GetArray());
      fHEt[type]->SetStats(0);
      fHEt[type]->SetLineColor(HN2(type)->GetLineColor());
      fHEt[type]->SetTitle(GetTitle());
   }

   // calculate average
   for (UShort_t ix=1; ix<=HN2(type)->GetNbinsX(); ix++) {
      Double_t totalE=0, totalN=0;
      for (UShort_t iy=1; iy<=HN2(type)->GetNbinsY(); iy++) {
         if (emax<HN2(type)->GetYaxis()->GetBinLowEdge(iy)) break;
         totalN += HN2(type)->GetBinContent(ix,iy) *
            HN2(type)->GetYaxis()->GetBinWidth(iy);
         totalE += HN2(type)->GetBinContent(ix,iy) * 
            HN2(type)->GetYaxis()->GetBinWidth(iy) *
            HN2(type)->GetYaxis()->GetBinCenter(iy);
      }
      fHEt[type]->SetBinContent(ix,totalE/totalN);
   }
   return fHEt[type];
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::Eave(UShort_t type)
{
   if (abs(fAverageE[type]-1.0)<0.01)
      fAverageE[type] = Lall(type)/Nall(type)/1.60217646e-6;
   return fAverageE[type];
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::N2(UShort_t type, Double_t time, Double_t energy)
{
   return HN2(type)->Interpolate(time, energy);
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::L2(UShort_t type, Double_t time, Double_t energy)
{
   return HL2(type)->Interpolate(time, energy);
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::Ne(UShort_t type, Double_t energy)
{
   return HNe(type)->Interpolate(energy);
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::Nt(UShort_t type, Double_t time)
{
   return HNt(type)->Interpolate(time);
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::Nall(UShort_t type)
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

Double_t NEUS::SupernovaModel::Lall(UShort_t type)
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

