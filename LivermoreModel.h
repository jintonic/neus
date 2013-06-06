#ifndef LIVERMOREMODEL_H
#define LIVERMOREMODEL_H

#include "SupernovaModel.h"

class TF1;
class TF2;

namespace NEUS { class LivermoreModel; }

class NEUS::LivermoreModel : public SupernovaModel
{
   protected:
      TF2 *fN2[fgNtype], *fL2[fgNtype];
      TF1 *fNe[fgNtype], *fNt[fgNtype];
      TF1 *fLe[fgNtype], *fLt[fgNtype];
      TF1 *fEt[fgNtype];

      Double_t WilsonN2(Double_t *x, Double_t *parameter);
      Double_t WilsonNe(Double_t *x, Double_t *parameter);
      Double_t WilsonNt(Double_t *x, Double_t *parameter);

      Double_t WilsonL2(Double_t *x, Double_t *parameter);
      Double_t WilsonLe(Double_t *x, Double_t *parameter);
      Double_t WilsonLt(Double_t *x, Double_t *parameter);

      Double_t WilsonEt(Double_t *x, Double_t *parameter);

   public:
      LivermoreModel(const char *name="LivermoreModel",
            const char *title="Livermore model");
      ~LivermoreModel();

      void SetDataLocation(const char *dir);
      void UseDivariData();
      void Clear(Option_t *option="");

      Double_t N2(UShort_t type, Double_t time, Double_t energy); // N(t, E) [1e50/s/MeV]
      Double_t Ne(UShort_t type, Double_t energy); // N(t) [1e50/s]
      Double_t Nt(UShort_t type, Double_t time); // N(E) [1e50/MeV]
      Double_t Nall(UShort_t type); // N [1e50]
      Double_t Eave(UShort_t type);

      TH2D* HN2(UShort_t type=1); // N(t, E) [1e50/s/MeV]
      TH1D* HNt(UShort_t type=1); // N(t) [1e50/s]
      TH1D* HNe(UShort_t type=1, Double_t tmax=17.9012/*second*/); // N(E) [1e50/MeV]

      TH2D* HL2(UShort_t type=1); // L(t, E) [1e50 erg/s/MeV]
      TH1D* HLt(UShort_t type=1); // L(t) [1e50 erg/s]
      TH1D* HLe(UShort_t type=1, Double_t tmax=17.9012/*second*/); // L(E) [1e50 erg/MeV]

      TH1D* HEt(UShort_t type=1); // <E>(t)/s

      TF2* FN2(UShort_t type=1); // N(t, E) [1e50/s/MeV]
      TF1* FNt(UShort_t type=1); // N(t) [1e50/s]
      TF1* FNe(UShort_t type=1, Double_t tmax=17.9012/*second*/); // N(E) [1e50/MeV]

      TF2* FL2(UShort_t type=1); // L(t, E) [1e50 erg/s/MeV]
      TF1* FLt(UShort_t type=1); // L(t) [1e50 erg/s]
      TF1* FLe(UShort_t type=1, Double_t tmax=17.9012/*second*/); // L(E) [1e50 erg/MeV]

      TF1* FEt(UShort_t type=1); // <E>(t)/s

      void Print();

      ClassDef(LivermoreModel,1);
};

#endif
