#include "SupernovaModel.h"

NEUS::SupernovaModel::SupernovaModel(const char *name, const char *title)
   : TNamed(name, title), fDataLocation() {};

Double_t NEUS::SupernovaModel::N2(UShort_t type, Double_t energy, Double_t time)
{ return 0.; }

TH2D* NEUS::SupernovaModel::HN2(UShort_t type) { return 0; }
TH1D* NEUS::SupernovaModel::HNt(UShort_t type, Double_t maxEv) { return 0; }
TH1D* NEUS::SupernovaModel::HNe(UShort_t type, Double_t tmax) { return 0; }

TH2D* NEUS::SupernovaModel::HL2(UShort_t type) { return 0; }
TH1D* NEUS::SupernovaModel::HLt(UShort_t type, Double_t maxEv) { return 0; }
TH1D* NEUS::SupernovaModel::HLe(UShort_t type, Double_t tmax) { return 0; }

TH1D* NEUS::SupernovaModel::HEt(UShort_t type) { return 0; }
