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
   public:
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

   protected:
      TString fDataLocation;

      Double_t fMinE, fMaxE;
      Double_t fMinT, fMaxT;

      Double_t fTotalN[fgNtype], fTotalL[fgNtype];
      Double_t fAverageE[fgNtype];

      /**
       * Histograms holding data.
       * They are used for both visualization and fast interpolation.
       */
      TH2D *fHN2[fgNtype], *fHL2[fgNtype];
      TH1D *fHNt[fgNtype], *fHNe[fgNtype];
      TH1D *fHLt[fgNtype], *fHLe[fgNtype];
      TH1D *fHEt[fgNtype];

      TF1 *fNeFD[fgNtype];

      Double_t NeFermiDirac(Double_t *x, Double_t *parameter);

   public:
      SupernovaModel();
      SupernovaModel(const char *name, const char *title);
      virtual ~SupernovaModel() { Clear(); }

      virtual void Clear(Option_t *option="");

      /**
       * Load data to histograms.
       */
      virtual void LoadData(const char *dir) { fDataLocation=dir; }
      const char* DataLocation() { return fDataLocation; }

      Double_t TMax() { return fMaxT; }
      Double_t TMin() { return fMinT; }
      Double_t EMax() { return fMaxE; }
      Double_t EMin() { return fMinE; }
      void SetEMin(double E) { fMinE=E; }
      void SetEMax(double E) { fMaxE=E; }

      /**
       * Number of neutrinos as a function of time and energy, N(t, E).
       * It is in unit of 1e50/MeV/second.
       * type 1: v_e
       * type 2: anti-v_e
       * type 3, 4, 5, 6: v_x, or types other than type 1 and 2
       * time: second after core collapse
       * energy: neutrino energy in unit of MeV/c2
       */
      Double_t N2(UShort_t type, Double_t time, Double_t energy);
      Double_t L2(UShort_t type, Double_t time, Double_t energy);
      /**
       * Number of neutrinos integrated over energy, N(t).
       * It is in unit of 1e50/second.
       */
      Double_t Nt(UShort_t type, Double_t time);
      /**
       * Number of neutrinos integrated over time, N(E).
       * It is in unit of 1e50/MeV.
       */
      Double_t Ne(UShort_t type, Double_t energy);
      /**
       * Fermi-Dirac approximation of N(E)
       * It is in unit of 1e50/MeV/second.
       */
      Double_t NeFD(UShort_t type, Double_t energy);

      virtual Double_t Nall(UShort_t type);
      virtual Double_t Lall(UShort_t type);
      Double_t Eave(UShort_t type);

      /**
       * N(t, E) in TH2D format.
       * x axis: second after the core collapse, in [TMin(), TMax()].
       * y axis: neutrino energy in unit of MeV/c2, in [EMin(), EMax()].
       * z axis: number of neutrinos in unit of 1e50/MeV/second.
       */
      TH2D* HN2(UShort_t type=1); 
      /**
       * N(t) in TH1D format.
       * It is the integration of N(t, E) over [EMin(), emax].
       * If emax>EMax(), emax is set to be EMax().
       * x axis: second after the core collapse, in [TMin(), TMax()].
       * y axis: number of neutrinos in unit of 1e50/second.
       */
      TH1D* HNt(UShort_t type=1, Double_t emax=999.);
      /**
       * N(E) in TH1D format.
       * It is the integration of N(t, E) over [TMin(), tmax].
       * If tmax>TMax(), tmax is set to be TMax().
       * x axis: neutrino energy, in [EMin(), EMax()].
       * y axis: number of neutrinos in unit of 1e50/MeV.
       */
      TH1D* HNe(UShort_t type=1, Double_t tmax=999.);

      TF1* FNeFD(UShort_t type=1); // Fermi-Dirac approximation of N(E)
      TH1D* HNeFD(UShort_t type=1); // Fermi-Dirac approximation of N(E)

      /**
       * Luminosity of neutrinos, L(t, E), in TH2D format.
       * x axis: second after the core collapse, in [TMin(), TMax()].
       * y axis: neutrino energy in unit of MeV, in [EMin(), EMax()].
       * z axis: luminosity of neutrinos in unit of 1e50*erg/MeV/second.
       */
      TH2D* HL2(UShort_t type=1);
      TH1D* HLt(UShort_t type=1, Double_t emax=999.); // L(t)
      TH1D* HLe(UShort_t type=1, Double_t tmax=999.); // L(E)

      /**
       * Average energy of neutrinos, <E>(t) in TH1D format.
       * The average is calculated in [Emin(), emax].
       * If emax>EMax(), emax is set to be EMax().
       * x axis: second after the core collapse, in [TMin(), TMax()].
       * y axis: average energy in unit of MeV/second.
       */
      TH1D* HEt(UShort_t type=1, Double_t emax=999.);

      ClassDef(SupernovaModel,1);
};

#endif
