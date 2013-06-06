#include "NakazatoModel.h"
using namespace NEUS;

#include <TF1.h>
#include <TH2D.h>
#include <TFile.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
using namespace std;

int main()
{
   // load database
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
   }

   NakazatoModel *blackHole = new NakazatoModel(30,0.004);
   blackHole->SetDataLocation(".");
   blackHole->LoadIntegratedData();

   // draw
   TCanvas *can = new TCanvas;
   can->Print("nakazato.ps[");

   // N(E) in 20 second
   TH1D *hNe1 = model[14]->HNe(1);
   TH1D *hNe2 = model[14]->HNe(2);
   TH1D *hNex = model[14]->HNe(3);

   can->SetLogy();
   hNe1->GetXaxis()->SetRangeUser(0,60);
   hNe1->GetYaxis()->SetRangeUser(1e3,1e7);
   hNe1->GetYaxis()->SetTitle("number of neutrino [10^{50}/MeV]");
   hNe1->Draw();
   hNe2->Draw("same");
   hNex->Draw("same");

   model[14]->FNeFD(1)->Draw("same");
   model[14]->FNeFD(2)->Draw("same");
   model[14]->FNeFD(3)->Draw("same");

   TLegend *leg = new TLegend(0.45,0.60,0.88,0.88);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetHeader("total number of #nu within 20 second:");
   leg->AddEntry(hNe1, Form("#nu_{e}: %.1e",model[14]->Nall(1)*1e50),"l");
   leg->AddEntry(hNe2,Form("#bar{#nu}_{e}: %.1e",model[14]->Nall(2)*1e50),"l");
   leg->AddEntry(hNex, Form("#nu_{x}: %.1e",model[14]->Nall(3)*1e50),"l");
   leg->Draw();

   can->Print("nakazato.ps");

   // N(E) in Livermore time window
   Double_t tmax=17.9012; // second
   hNe1 = model[14]->HNe(1, tmax);
   hNe2 = model[14]->HNe(2, tmax);
   hNex = model[14]->HNe(3, tmax);

   hNe1->GetXaxis()->SetRangeUser(0,60);
   hNe1->GetYaxis()->SetTitle("number of neutrino [10^{50}/MeV]");
   hNe1->Draw();
   hNe2->Draw("same");
   hNex->Draw("same");

   leg->Clear();
   leg->SetHeader(Form("total number of #nu within %.1f second:",tmax));
   leg->AddEntry(hNe1, Form("#nu_{e}: %.1e",hNe1->Integral("width")*1e50),"l");
   leg->AddEntry(hNe2,
         Form("#bar{#nu}_{e}: %.1e",hNe2->Integral("width")*1e50),"l");
   leg->AddEntry(hNex, Form("#nu_{x}: %.1e",hNex->Integral("width")*1e50),"l");
   leg->Draw();

   can->Print("nakazato.ps");

   // <E>(t)
   can->SetLogy(0);
   TH1D *hEt1 = model[14]->HEt(1);
   TH1D *hEt2 = model[14]->HEt(2);
   TH1D *hEt3 = model[14]->HEt(3);

   hEt1->GetXaxis()->SetRangeUser(-0.05,0.9);
   hEt1->GetYaxis()->SetRangeUser(0,25);
   hEt1->Draw();
   hEt2->Draw("same");
   hEt3->Draw("same");

   leg->Clear();
   leg->SetHeader("");
   leg->SetX1NDC(0.6);
   leg->AddEntry(hNe1, "#nu_{e}");
   leg->AddEntry(hNe2, "#bar{#nu}_{e}");
   leg->AddEntry(hNex, "#nu_{x}");
   leg->Draw();
   can->Print("nakazato.ps");

   hEt1->GetXaxis()->SetRangeUser(0,20);
   hEt1->GetYaxis()->SetRangeUser(5,20);
   hEt1->Draw();
   hEt2->Draw("same");
   hEt3->Draw("same");

   leg->Draw();
   can->Print("nakazato.ps");

   // L(t)
   can->SetLogy();
   TH1D *hLt1 = model[14]->HLt(1);
   TH1D *hLt2 = model[14]->HLt(2);
   TH1D *hLt3 = model[14]->HLt(3);

   hLt1->GetXaxis()->SetRangeUser(-0.5,0.9);
   hLt1->Draw();
   hLt2->Draw("same");
   hLt3->Draw("same");

   leg->Draw();
   can->Print("nakazato.ps");

   hLt1->GetXaxis()->SetRangeUser(0.,20.);
   hLt1->Draw();
   hLt2->Draw("same");
   hLt3->Draw("same");

   leg->Draw();
   can->Print("nakazato.ps");

   can->Print("nakazato.ps]");
}
