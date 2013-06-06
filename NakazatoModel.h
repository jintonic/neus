#ifndef NAKAZATOMODEL_H
#define NAKAZATOMODEL_H

#include "SupernovaModel.h"

namespace NEUS { class NakazatoModel; }

/**
 * Nakazato's supernova model.
 * The data are saved in TH2D and TH1D format so that they can be easily
 * visualized and interpolated. The data are not equally binned so do the
 * histograms. They can be nicely displayed both in linear and in logrithmic
 * scale.
 */
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

      /**
       * Load data as function of energy and time.
       * The data are divided by 1e50 and then loaded into TH2D objects.
       */
      void LoadFullData();
      /**
       * Load data as function of E.
       * They are integrated in [TMin(), TMax()].
       * The data are divided by 1e50 and then loaded into TH2D objects.
       */
      void LoadIntegratedData();

      Double_t N2(UShort_t type, Double_t time, Double_t energy);
      Double_t Nt(UShort_t type, Double_t time);
      Double_t Ne(UShort_t type, Double_t energy);
      Double_t Nall(UShort_t type);
      Double_t Lall(UShort_t type);
      Double_t Eave(UShort_t type);

      TH2D* HN2(UShort_t type=1);
      TH1D* HNt(UShort_t type=1);
      /**
       * Number of neutrinos integrated from [TMin(), tmax]
       * If tmax>TMax(), tmax is set to be TMax().
       * If tmax=0, the integrated database is used instead of the full one.
       */
      TH1D* HNe(UShort_t type=1, Double_t tmax=0);

      TH2D* HL2(UShort_t type=1);
      TH1D* HLt(UShort_t type=1);
      TH1D* HLe(UShort_t type=1, Double_t tmax=0);

      TH1D* HEt(UShort_t type=1);

      void Print();

      ClassDef(NakazatoModel,1);
};

#endif
