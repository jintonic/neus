#ifndef LIVERMOREMODEL_H
#define LIVERMOREMODEL_H

#include "SupernovaModel.h"

class TF1;
class TF2;

namespace NEUS { class LivermoreModel; }

/**
 * C++ wrapper of Totani's interpolator Written in Fortran.
 * TF1 and TF2 are used to visualize the results. TH1D and TH2D are created
 * inside TF1 and TF2 for drawing. They are recreated if the scale of x axis is
 * switched between linear and logarithmic for better visualization. They can
 * also be used to get a value without going through the interpolation again.
 */
class NEUS::LivermoreModel : public SupernovaModel
{
   protected:
      TF2 *fN2[fgNtype], *fL2[fgNtype];
      TF1 *fNe[fgNtype], *fNt[fgNtype];
      TF1 *fLe[fgNtype], *fLt[fgNtype];
      TF1 *fEt[fgNtype];

      Double_t WilsonN2(Double_t *x, Double_t *parameter);
      Double_t WilsonNe(Double_t *x, Double_t *parameter);
      Double_t WilsonNt(Double_t *x, Double_t *parameter);

      Double_t WilsonL2(Double_t *x, Double_t *parameter);
      Double_t WilsonLe(Double_t *x, Double_t *parameter);
      Double_t WilsonLt(Double_t *x, Double_t *parameter);

      Double_t WilsonEt(Double_t *x, Double_t *parameter);

   public:
      LivermoreModel(const char *name="LivermoreModel",
            const char *title="Livermore model");
      ~LivermoreModel();

      /**
       * Set path to directory holding wilson-*.dat
       */
      void SetDataLocation(const char *dir);
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
      Double_t Ne(UShort_t type, Double_t energy);
      Double_t Nt(UShort_t type, Double_t time);
      Double_t Nall(UShort_t type);
      Double_t Lall(UShort_t type);
      Double_t Eave(UShort_t type);

      TH2D* HN2(UShort_t type=1);
      TH1D* HNt(UShort_t type=1, Double_t emax=99);
      /**
       * Number of neutrinos integrated from [TMin(), tmax]
       * If tmax>TMax(), tmax is set to be TMax().
       * If tmax<10*second, tmax is set to be 10*second.
       */
      TH1D* HNe(UShort_t type=1, Double_t tmax=20);

      TH2D* HL2(UShort_t type=1);
      TH1D* HLt(UShort_t type=1, Double_t emax=99);
      TH1D* HLe(UShort_t type=1, Double_t tmax=20);

      TH1D* HEt(UShort_t type=1, Double_t emax=99);

      TF2* FN2(UShort_t type=1);
      TF1* FNt(UShort_t type=1, Double_t emax=99);
      TF1* FNe(UShort_t type=1, Double_t tmax=20);

      TF2* FL2(UShort_t type=1);
      TF1* FLt(UShort_t type=1, Double_t emax=99);
      TF1* FLe(UShort_t type=1, Double_t tmax=20);

      TF1* FEt(UShort_t type=1, Double_t emax=99);

      /**
       * Set number of bins in energy axis.
       * It sets precision of histograms created inside TF1 or TF2.
       */
      void SetNbinsE(Int_t n=1000);
      void SetNbinsT(Int_t n=1000);

      void Print();

      ClassDef(LivermoreModel,1);
};

#endif
