#include <TH1D.h>
#include <TH2D.h>

class SupernovaModel : public TNamed
{
   protected:
      TH1D *fH1Empty;
      TH2D *fH2Empty;
      TH1D *fH1N1, *fH1N2, *fH1N3, *fH1E1, *fH1E2, *fH1E3;
      TH2D *fH2N1, *fH2N2, *fH2N3, *fH2E1, *fH2E2, *fH2E3;

      Float_t fInitialMass; // progenitor mass in Solar mass unit
      Float_t fMetallicity;
      Float_t fReviveTime; // shock revival time in ms

   public:
      SupernovaModel(
            Float_t initialMass=13, /* Solar mass */
            Float_t metallicity=0.02,
            Float_t reviveTime=100 /* ms */);
      virtual ~SupernovaModel();

      void LoadFullData(const char *databaseDir="./intpdata");
      void LoadIntegratedData(const char *databaseDir="./integdata");

      TH2D* NumberSpectrum(const char *neutrino="v_e");
      TH1D* IntegratedNumberSpectrum(const char *neutrino="anti-v_e");

      TH2D* EnergySpectrum(const char *neutrino="v_e");
      TH1D* IntegratedEnergySpectrum(const char *neutrino="anti-v_e");

      Double_t InitialMass() { return fInitialMass; }
      Double_t Metallicity() { return fMetallicity; }
      Double_t ReviveTime() { return fReviveTime; }

      ClassDef(SupernovaModel,1);
};
