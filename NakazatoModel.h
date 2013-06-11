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
      Float_t fInitialMass; // progenitor mass in Solar mass unit
      Float_t fMetallicity;
      Float_t fReviveTime; // shock revival time in ms

      static const UShort_t fNbinsT = 391;
      static const UShort_t fNbinsE = 20;

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

   public:
      NakazatoModel(
            Float_t initialMass=13, /* Solar mass */
            Float_t metallicity=0.02,
            Float_t reviveTime=100 /* ms */);
      ~NakazatoModel() {};

      Double_t InitialMass() { return fInitialMass; }
      Double_t Metallicity() { return fMetallicity; }
      Double_t ReviveTime() { return fReviveTime; }

      void LoadData(const char *dir);

      void Print();

      ClassDef(NakazatoModel,1);
};

#endif
