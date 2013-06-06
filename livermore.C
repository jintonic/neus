#include "LivermoreModel.h"
using namespace NEUS;

#include <TF2.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>

#include <iostream>
using namespace std;

int main()
{
   LivermoreModel *sn = new LivermoreModel;
   sn->SetDataLocation("../total");

   cout<<"average neutrino energy:"<<endl;
   cout<<"type 1: "<<sn->Eave(1)<<" MeV"<<endl;
   cout<<"type 2: "<<sn->Eave(2)<<" MeV"<<endl;
   cout<<"type x: "<<sn->Eave(3)<<" MeV"<<endl;

   // draw spectra
   TCanvas *can = new TCanvas;
   can->Print("livermore.ps[");

   can->SetLogz();
   sn->FN2(1)->Draw("colz");
   can->Print("livermore.ps");

   // N(E) in 18 second
   TF1 *fNe1 = sn->FNe(1);
   TF1 *fNe2 = sn->FNe(2);
   TF1 *fNex = sn->FNe(3);

   can->SetLogy();
   fNe1->GetXaxis()->SetRangeUser(0,60);
   fNe1->GetYaxis()->SetRangeUser(1e3,1e7);
   fNe1->GetYaxis()->SetTitle("number of neutrino [10^{50}/MeV]");
   fNe1->Draw();
   fNe2->Draw("same");
   fNex->Draw("same");

   sn->FNeFD(1)->Draw("same");
   sn->FNeFD(2)->Draw("same");
   sn->FNeFD(3)->Draw("same");

   TLegend *leg = new TLegend(0.5,0.60,0.88,0.88);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetHeader("total number of #nu within 18 second:");
   leg->AddEntry(fNe1, Form("#nu_{e}: %.1e",sn->Nall(1)*1e50),"l");
   leg->AddEntry(fNe2,Form("#bar{#nu}_{e}: %.1e",sn->Nall(2)*1e50),"l");
   leg->AddEntry(fNex, Form("#nu_{x}: %.1e",sn->Nall(3)*1e50),"l");
   leg->Draw();

   can->Print("livermore.ps");

   // L(t)
   TF1 *fLt1 = sn->FLt(1);
   TF1 *fLt2 = sn->FLt(2);
   TF1 *fLt3 = sn->FLt(3);

   can->SetLogx();
   fLt1->SetRange(0.02,17.9012); // works, but recreate histogram :(
   fLt1->GetXaxis()->SetTitle("time [second]");
   fLt1->GetYaxis()->SetRangeUser(3,1e4);
   fLt1->GetYaxis()->SetTitle("luminosity [10^{50} erg/second]");
   fLt1->SetTitle(sn->GetTitle());
   fLt1->Draw();
   fLt2->Draw("same");
   fLt3->Draw("same");

   leg->Clear();
   leg->SetHeader("");
   leg->SetX1NDC(0.6);
   leg->AddEntry(fNe1, "#nu_{e}");
   leg->AddEntry(fNe2, "#bar{#nu}_{e}");
   leg->AddEntry(fNex, "#nu_{x}");
   leg->Draw();
   can->Print("livermore.ps");

   // <E>(t)
   TF1 *fEt1 = sn->FEt(1);
   TF1 *fEt2 = sn->FEt(2);
   TF1 *fEt3 = sn->FEt(3);

   can->SetLogy(0);
   fEt1->SetRange(2.5e-2,17);
   fEt1->GetXaxis()->SetTitle("time [second]");
   fEt1->GetYaxis()->SetRangeUser(5,30);
   fEt1->GetYaxis()->SetTitle("average energy [MeV/second]");
   fEt1->Draw();
   fEt2->Draw("same");
   fEt3->Draw("same");

   leg->SetX1NDC(0.3);
   leg->SetX2NDC(0.5);
   leg->Draw();
   can->Print("livermore.ps");

   can->Print("livermore.ps]");
}
