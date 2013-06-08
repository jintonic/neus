#ifndef LIVERMOREMODEL_H
#define LIVERMOREMODEL_H

#include "SupernovaModel.h"

namespace NEUS { class LivermoreModel; }

/**
 * Livermore model.
 * Totani's interpolator written in Fortran is used to fill histograms, which
 * are used for visualization and fast interpolation.
 */
class NEUS::LivermoreModel : public SupernovaModel
{
   public:
      LivermoreModel(const char *name="LivermoreModel",
            const char *title="Livermore model");
      ~LivermoreModel() {};

      void LoadData(const char *dir);
      /**
       * Use <E> and N given in Divari 2012.
       * Divari et al. claim that they use the Livermore model for their
       * calculation in http://stacks.iop.org/0954-3899/39/i=9/a=095204. But
       * the average energy and total number of neutrinos used in their paper
       * are different from Totani's Fortran code. This function is used to set
       * <E> and N to the values in their paper.
       */
      void UseDivariData();
      /**
       * Clear internal data.
       */
      void Clear(Option_t *option="");

      Double_t N2(UShort_t type, Double_t time, Double_t energy);
      Double_t L2(UShort_t type, Double_t time, Double_t energy);
      Double_t Ne(UShort_t type, Double_t energy);
      Double_t Nt(UShort_t type, Double_t time);
      Double_t Nall(UShort_t type);
      Double_t Lall(UShort_t type);

      void Print();

      ClassDef(LivermoreModel,1);
};

#endif
