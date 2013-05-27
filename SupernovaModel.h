#ifndef SUPERNOVAMODEL_H
#define SUPERNOVAMODEL_H

#include <TNamed.h>

class TH1D;
class TH2D;

namespace NEUS { class SupernovaModel; }

class NEUS::SupernovaModel : public TNamed
{
   protected:
      TString fDataLocation;

   public:
      SupernovaModel() : TNamed(), fDataLocation() {};
      SupernovaModel(const char *name, const char *title);
      virtual ~SupernovaModel() {};

      virtual void SetDataLocation(const char *dir) { fDataLocation=dir; }
      const char* DataLocation() { return fDataLocation; }

      virtual Double_t N2(UShort_t type/*1:v_e, 2: anti-v_e, 3, 4, 5, 6: v_x*/,
            Double_t energy/*MeV*/, Double_t time/*second*/);
      virtual Double_t Ne(UShort_t type, Double_t energy/*MeV*/) { return 0; }
      virtual Double_t Nt(UShort_t type, Double_t time/*second*/) { return 0; }
      virtual Double_t Nall(UShort_t type) { return 0; }

      virtual TH2D* HN2(UShort_t type=1); // N(E, t)
      virtual TH1D* HNt(UShort_t type=1, Double_t maxEv=100/*MeV*/); // N(t)
      virtual TH1D* HNe(UShort_t type=1, Double_t tmax=18/*second*/); // N(E)

      virtual TH2D* HL2(UShort_t type=1); // L(E, t)
      virtual TH1D* HLt(UShort_t type=1, Double_t maxEv=100/*MeV*/); // L(t)
      virtual TH1D* HLe(UShort_t type=1, Double_t tmax=18/*second*/); // L(E)

      virtual TH1D* HEt(UShort_t type=1); // <E>(t)

      ClassDef(SupernovaModel,1);
};

#endif
