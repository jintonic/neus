#ifndef NAKAZATOMODEL_HH
#define NAKAZATOMODEL_HH

#include <TH2D.h>

class NakazatoModel : public TNamed
{
   protected:
      TH2D *fH2N[7], *fH2L[7];
      TH1D *fHNe[7], *fHNt[7], *fHLe[7], *fHLt[7], *fHEt[7];

      Float_t fInitialMass; // progenitor mass in Solar mass unit
      Float_t fMetallicity;
      Float_t fReviveTime; // shock revival time in ms

   public:
      NakazatoModel(
            Float_t initialMass=13, /* Solar mass */
            Float_t metallicity=0.02,
            Float_t reviveTime=100 /* ms */);
      virtual ~NakazatoModel();

      void LoadFullData(const char *databaseDir="./intpdata");
      void LoadIntegratedData(const char *databaseDir="./integdata");

      TH2D* H2N(UShort_t type=1);
      TH1D* HNt(UShort_t type=1);
      TH1D* HNe(UShort_t type=1, Double_t tmax=0);

      Double_t N2(UShort_t type, Double_t energy, Double_t time);
      Double_t Nt(UShort_t type, Double_t time);
      Double_t Ne(UShort_t type, Double_t energy, Double_t tmax);

      TH2D* H2L(UShort_t type=1);
      TH1D* HLt(UShort_t type=1);
      TH1D* HLe(UShort_t type=1, Double_t tmax=0);

      Double_t E2(UShort_t type, Double_t energy, Double_t time);
      Double_t Lt(UShort_t type, Double_t time);
      Double_t Le(UShort_t type, Double_t energy, Double_t tmax);

      TH1D* HEt(UShort_t type=1);

      Double_t InitialMass() { return fInitialMass; }
      Double_t Metallicity() { return fMetallicity; }
      Double_t ReviveTime() { return fReviveTime; }

      ClassDef(NakazatoModel,1);
};

#endif // NAKAZATOMODEL_HH
