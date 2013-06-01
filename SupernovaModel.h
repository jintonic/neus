#ifndef SUPERNOVAMODEL_H
#define SUPERNOVAMODEL_H

#include <TNamed.h>

class TF1;
class TH1D;
class TH2D;

namespace NEUS { class SupernovaModel; }

class NEUS::SupernovaModel : public TNamed
{
   protected:
      TString fDataLocation;

      Double_t fMinE, fMaxE;
      Double_t fMinT, fMaxT;

      static const UShort_t fgNtype = 7;
      Double_t fTotalN[fgNtype];
      Double_t fAverageE[fgNtype];

      TF1 *fNeFD[fgNtype];

      virtual Double_t NeFermiDirac(Double_t *x, Double_t *parameter);

   public:
      SupernovaModel();
      SupernovaModel(const char *name, const char *title);
      virtual ~SupernovaModel();

      virtual void SetDataLocation(const char *dir) { fDataLocation=dir; }
      const char* DataLocation() { return fDataLocation; }

      virtual Double_t N2(UShort_t type/*1:v_e, 2: anti-v_e, 3, 4, 5, 6: v_x*/,
            Double_t time/*second*/, Double_t energy/*MeV*/) { return 0; }
      virtual Double_t Ne(UShort_t type, Double_t energy/*MeV*/) { return 0; }
      virtual Double_t Nt(UShort_t type, Double_t time/*second*/) { return 0; }
      virtual Double_t Nall(UShort_t type) { return fTotalN[type]; }
      virtual Double_t Eave(UShort_t type) { return fAverageE[type]; }

      // Fermi-Dirac approximation of N(E)
      virtual Double_t NeFD(UShort_t type, Double_t energy/*MeV*/);

      virtual TH2D* HN2(UShort_t type) { return 0; } // N(E, t)
      virtual TH1D* HNt(UShort_t type) { return 0; } // N(t)
      virtual TH1D* HNe(UShort_t type, Double_t tmax/*second*/) // N(E)
      { return 0; }

      virtual TH2D* HL2(UShort_t type) { return 0; } // L(E, t)
      virtual TH1D* HLt(UShort_t type) { return 0; } // L(t)
      virtual TH1D* HLe(UShort_t type, Double_t tmax/*second*/) // L(E)
      { return 0; }

      virtual TH1D* HEt(UShort_t type) { return 0; } // <E>(t)

      TF1* FNeFD(UShort_t type); // Fermi-Dirac approximation of N(E)
      TH1D* HNeFD(UShort_t type); // Fermi-Dirac approximation of N(E)

      ClassDef(SupernovaModel,1);
};

#endif
