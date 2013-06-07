#ifndef SUPERNOVAMODEL_H
#define SUPERNOVAMODEL_H

#include <TNamed.h>

class TF1;
class TH1D;
class TH2D;

namespace NEUS { class SupernovaModel; }

/**
 * Base class of all models.
 * The Fermi-Dirac approximation of energy spectra is implemenated here.
 */
class NEUS::SupernovaModel : public TNamed
{
   protected:
      TString fDataLocation;

      Double_t fMinE, fMaxE;
      Double_t fMinT, fMaxT;

      /**
       * Numbers of types of neutrinos.
       * There are 6 types of neutrinos:
       * type 1: v_e
       * type 2: anti-v_e
       * type 3: v_mu
       * type 4: anti-v_mu
       * type 5: v_tau
       * type 6: anti-v_tau
       * type 0: reserved for special uses
       */
      static const UShort_t fgNtype = 7;
      Double_t fTotalN[fgNtype], fTotalL[fgNtype];
      Double_t fAverageE[fgNtype];

      TF1 *fNeFD[fgNtype];

      virtual Double_t NeFermiDirac(Double_t *x, Double_t *parameter);

   public:
      SupernovaModel();
      SupernovaModel(const char *name, const char *title);
      virtual ~SupernovaModel();

      virtual void SetDataLocation(const char *dir) { fDataLocation=dir; }
      const char* DataLocation() { return fDataLocation; }

      Double_t TMax() { return fMaxT; }
      Double_t TMin() { return fMinT; }
      Double_t EMax() { return fMaxT; }
      Double_t EMin() { return fMinE; }

      /**
       * Number of neutrinos as a function of type, time and energy.
       * It is in unit of 1e50/MeV/second.
       * type 1: v_e
       * type 2: anti-v_e
       * type 3, 4, 5, 6: v_x, or types other than type 1 and 2
       * time: second after core collapse
       * energy: neutrino energy in unit of MeV/c2
       */
      virtual Double_t N2(UShort_t type, Double_t time, Double_t energy)
      { return 0; }
      /**
       * Number of neutrinos integrated over time.
       * It is in unit of 1e50/MeV.
       */
      virtual Double_t Ne(UShort_t type, Double_t energy) { return 0; }
      /**
       * Number of neutrinos integrated over energy.
       * It is in unit of 1e50/second.
       */
      virtual Double_t Nt(UShort_t type, Double_t time) { return 0; }
      /**
       * Fermi-Dirac approximation of N(E)
       * It is in unit of 1e50/MeV/second.
       */
      virtual Double_t NeFD(UShort_t type, Double_t energy);

      virtual Double_t Nall(UShort_t type) { return fTotalN[type]; }
      virtual Double_t Lall(UShort_t type) { return fTotalL[type]; }
      virtual Double_t Eave(UShort_t type) { return fAverageE[type]; }

      /**
       * Number of neutrinos in TH2D format.
       * x axis: second after the core collapse, in [TMin(), TMax()].
       * y axis: neutrino energy in unit of MeV/c2, in [EMin(), EMax()].
       * z axis: number of neutrinos in unit of 1e50/MeV/second.
       */
      virtual TH2D* HN2(UShort_t type) { return 0; }
      /**
       * Number of neutrinos integrated from [EMin(), emax]
       * If emax>EMax(), emax is set to be EMax().
       */
      virtual TH1D* HNt(UShort_t type, Double_t emax) { return 0; } // N(t)
      /**
       * Number of neutrinos integrated from [TMin(), tmax]
       * If tmax>TMax(), tmax is set to be TMax().
       */
      virtual TH1D* HNe(UShort_t type, Double_t tmax) { return 0; }

      TF1* FNeFD(UShort_t type); // Fermi-Dirac approximation of N(E)
      TH1D* HNeFD(UShort_t type); // Fermi-Dirac approximation of N(E)

      /**
       * Luminosity of neutrinos in TH2D format.
       * x axis: second after the core collapse, in [TMin(), TMax()].
       * y axis: neutrino energy in unit of MeV/c2, in [EMin(), EMax()].
       * z axis: luminosity of neutrinos in unit of 1e50*erg/MeV/second.
       */
      virtual TH2D* HL2(UShort_t type) { return 0; }
      virtual TH1D* HLt(UShort_t type, Double_t emax) { return 0; } // L(t)
      virtual TH1D* HLe(UShort_t type, Double_t tmax) { return 0; } // L(E)

      /**
       * Average energy of neutrinos in TH1D format.
       * The average is calculated in [Emin(), emax].
       * x axis: second after the core collapse, in [TMin(), TMax()].
       * y axis: average energy in unit of MeV/second.
       */
      virtual TH1D* HEt(UShort_t type, Double_t emax) { return 0; } // <E>(t)

      ClassDef(SupernovaModel,1);
};

#endif
