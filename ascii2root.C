#include "NakazatoModel.h"
#include "LivermoreModel.h"
using namespace NEUS;

#include <TFile.h>

int main()
{
   const UShort_t nm = 21;
   Float_t mass[nm] = {13,13,13,13,13,13, 20,20,20,20,20,20, 
      30,30,30, 50,50,50,50,50,50};
   Float_t meta[nm] = {0.02,0.02,0.02,0.004,0.004,0.004,
      0.02,0.02,0.02,0.004,0.004,0.004, 0.02,0.02,0.02,
      0.02,0.02,0.02,0.004,0.004,0.004};
   Float_t trev[nm] = {100,200,300,100,200,300, 100,200,300,100,200,300,
   100,200,300, 100,200,300,100,200,300};
   NakazatoModel *model[nm];

   for (UShort_t i=0; i<nm; i++) {
      model[i] = new NakazatoModel(mass[i],meta[i],trev[i]);
      model[i]->SetDataLocation(".");
      model[i]->LoadIntegratedData();
      model[i]->LoadFullData();
      model[i]->Print();
   }

   NakazatoModel *blackHole = new NakazatoModel(30,0.004);
   blackHole->SetDataLocation(".");
   blackHole->LoadIntegratedData();
   blackHole->Print();

   LivermoreModel *totani = new LivermoreModel;
   totani->SetDataLocation("../total");
   totani->Print();

   TFile *output = new TFile("models.root","recreate");
   for (UShort_t i=0; i<nm; i++) model[i]->Write();
   blackHole->Write();
   totani->Write();
   output->Close();
}
