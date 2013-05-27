#ifndef LIVERMOREMODEL_H
#define LIVERMOREMODEL_H

#include "SupernovaModel.h"

class TF1;
class TF2;

namespace NEUS { class LivermoreModel; }

class NEUS::LivermoreModel : public SupernovaModel
{
   protected:
      TF2 *fN2[7], *fL2[7];
      TF1 *fNe[7], *fNt[7], *fLe[7], *fLt[7], *fEt[7];

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

      Double_t N2(UShort_t type, Double_t energy, Double_t time);
      Double_t Ne(UShort_t type, Double_t energy);
      Double_t Nt(UShort_t type, Double_t time);
      Double_t Nall(UShort_t type);

      TH2D* HN2(UShort_t type=1); // N(E, t)
      TH1D* HNt(UShort_t type=1, Double_t maxEv=100/*MeV*/); // N(t)
      TH1D* HNe(UShort_t type=1, Double_t tmax=18/*second*/); // N(E)

      TH2D* HL2(UShort_t type=1); // L(E, t)
      TH1D* HLt(UShort_t type=1, Double_t maxEv=100/*MeV*/); // L(t)
      TH1D* HLe(UShort_t type=1, Double_t tmax=18/*second*/); // L(E)

      TH1D* HEt(UShort_t type=1); // <E>(t)

      TF2* FN2(UShort_t type=1); // N(E, t)
      TF1* FNt(UShort_t type=1, Double_t maxEv=100/*MeV*/); // N(t)
      TF1* FNe(UShort_t type=1, Double_t tmax=18/*second*/); // N(E)

      TF2* FL2(UShort_t type=1); // L(E, t)
      TF1* FLt(UShort_t type=1, Double_t maxEv=100/*MeV*/); // L(t)
      TF1* FLe(UShort_t type=1, Double_t tmax=18/*second*/); // L(E)

      TF1* FEt(UShort_t type=1, Double_t maxEv=100/*MeV*/); // <E>(t)

      void Print();

      ClassDef(LivermoreModel,1);
};

#endif
