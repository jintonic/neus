#ifndef NAKAZATOMODEL_H
#define NAKAZATOMODEL_H

#include "SupernovaModel.h"

namespace NEUS { class NakazatoModel; }

class NEUS::NakazatoModel : public SupernovaModel
{
   private:
      TH2D *fHN2[fgNtype], *fHL2[fgNtype];
      TH1D *fHNe[fgNtype], *fHNt[fgNtype];
      TH1D *fHLe[fgNtype], *fHLt[fgNtype];
      TH1D *fHEt[fgNtype];

      Float_t fInitialMass; // progenitor mass in Solar mass unit
      Float_t fMetallicity;
      Float_t fReviveTime; // shock revival time in ms

      static const UShort_t fNbinsT = 391;
      static const UShort_t fNbinsE = 20;

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

      Double_t N2(UShort_t type, Double_t time, Double_t energy); // N(t, E) [1e50/s/MeV]
      Double_t Nt(UShort_t type, Double_t time); // N(t) [1e50/s]
      Double_t Ne(UShort_t type, Double_t energy); // N(E) [1e50/MeV]
      Double_t Nall(UShort_t type); // N [1e50]
      Double_t Eave(UShort_t type);

      TH2D* HN2(UShort_t type=1); // N(t, E) [1e50/s/MeV]
      TH1D* HNt(UShort_t type=1); // N(t) [1e50/s]
      TH1D* HNe(UShort_t type=1, Double_t tmax=0); // N(E) [1e50/MeV]

      TH2D* HL2(UShort_t type=1); // L(t, E) [1e50 erg/s/MeV]
      TH1D* HLt(UShort_t type=1); // L(t) [1e50 erg/s]
      TH1D* HLe(UShort_t type=1, Double_t tmax=0); // L(E) [1e50 erg/MeV]

      TH1D* HEt(UShort_t type=1); // <E>(t)/s

      void Print();

      ClassDef(NakazatoModel,1);
};

#endif
