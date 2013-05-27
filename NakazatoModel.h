#ifndef NAKAZATOMODEL_H
#define NAKAZATOMODEL_H

#include "SupernovaModel.h"

namespace NEUS { class NakazatoModel; }

class NEUS::NakazatoModel : public SupernovaModel
{
   private:
      TH2D *fHN2[7], *fHL2[7];
      TH1D *fHNe[7], *fHNt[7], *fHLe[7], *fHLt[7], *fHEt[7];

      Float_t fInitialMass; // progenitor mass in Solar mass unit
      Float_t fMetallicity;
      Float_t fReviveTime; // shock revival time in ms

   public:
      NakazatoModel(
            Float_t initialMass=13, /* Solar mass */
            Float_t metallicity=0.02,
            Float_t reviveTime=100 /* ms */);
      ~NakazatoModel();

      Double_t InitialMass() { return fInitialMass; }
      Double_t Metallicity() { return fMetallicity; }
      Double_t ReviveTime() { return fReviveTime; }

      void LoadFullData();
      void LoadIntegratedData();

      Double_t N2(UShort_t type, Double_t energy, Double_t time);
      Double_t Ne(UShort_t type, Double_t energy);
      Double_t Nall(UShort_t type);

      TH2D* HN2(UShort_t type=1); // N(E, t)
      TH1D* HNt(UShort_t type=1); // N(t)
      TH1D* HNe(UShort_t type=1, Double_t tmax=0); // N(E)

      TH2D* HL2(UShort_t type=1); // L(E, t)
      TH1D* HLt(UShort_t type=1); // L(t)
      TH1D* HLe(UShort_t type=1, Double_t tmax=0); // L(E)

      TH1D* HEt(UShort_t type=1); // <E>(t)

      ClassDef(NakazatoModel,1);
};

#endif
