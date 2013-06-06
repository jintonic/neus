#include "SupernovaModel.h"

#include <TF1.h>
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
      fAverageE[i] = 1.0; // 0 cannot be in the denominator
      fNeFD[i] = NULL;
   }
}

//______________________________________________________________________________
//

NEUS::SupernovaModel::SupernovaModel(const char *name, const char *title) : 
   TNamed(name, title), fDataLocation(), fMinE(0), fMaxE(0), fMinT(0), fMaxT(0)
{
   for (UShort_t i=0; i<fgNtype; i++) {
      fTotalN[i] = 0;
      fAverageE[i] = 1;
      fNeFD[i] = NULL;
   }
}

//______________________________________________________________________________
//

NEUS::SupernovaModel::~SupernovaModel()
{
   for (UShort_t i=0; i<fgNtype; i++) if (fNeFD[i]) delete fNeFD[i];
}

//______________________________________________________________________________
//

Double_t NEUS::SupernovaModel::NeFD(UShort_t type, Double_t energy)
{
   if (type<1 || type>6) {
      Warning("NeFD","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("NeFD","0 is returned!");
      return 0;
   }
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
      Warning("HNeFD","Type of neutrino must be one of 1, 2, 3, 4, 5, 6!");
      Warning("HNeFD","NULL pointer is returned!");
      return NULL;
   }
   if (fNeFD[type]) return fNeFD[type];

   fNeFD[type] = new TF1(Form("fNeFD%d",type), this, 
         &SupernovaModel::NeFermiDirac,
         fMinE, fMaxE, 1, "SupernovaModel", "HNeFD");
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
