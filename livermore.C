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

   sn->SetNbinsT(200);
   can->SetLogz();
   sn->FN2(1)->Draw("colz");
   can->Print("livermore.ps");

   // N(E) in 17.9012 second
   sn->SetNbinsE(200);
   can->SetLogy(); // set before creating histograms
   TH1D *hNe1 = sn->HNe(1);
   TH1D *hNe2 = sn->HNe(2);
   TH1D *hNex = sn->HNe(3);

   hNe1->GetXaxis()->SetRangeUser(0,60);
   hNe1->GetYaxis()->SetRangeUser(1e3,1e7);
   hNe1->GetYaxis()->SetTitle("number of neutrino [10^{50}/MeV]");
   hNe1->Draw();
   hNe2->Draw("same");
   hNex->Draw("same");

   sn->FNeFD(1)->Draw("same");
   sn->FNeFD(2)->Draw("same");
   sn->FNeFD(3)->Draw("same");

   TLegend *leg = new TLegend(0.5,0.60,0.88,0.88);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetHeader("total number of #nu within 17.9 second:");
   leg->AddEntry(hNe1, Form("#nu_{e}: %.1e",sn->Nall(1)*1e50),"l");
   leg->AddEntry(hNe2,Form("#bar{#nu}_{e}: %.1e",sn->Nall(2)*1e50),"l");
   leg->AddEntry(hNex, Form("#nu_{x}: %.1e",sn->Nall(3)*1e50),"l");
   leg->Draw();

   can->Print("livermore.ps");

   // L(t)
   sn->SetNbinsT(1000); // set nbins
   can->SetLogx(); // set before creating histograms
   can->SetLogy(); // set before creating histograms

   TH1D *hLt1 = sn->HLt(1); // create internal histogram
   TH1D *hLt2 = sn->HLt(2);
   TH1D *hLt3 = sn->HLt(3);

   hLt1->GetXaxis()->SetRangeUser(0.0212,17.9012);
   hLt1->GetYaxis()->SetRangeUser(3,1e4);
   hLt1->GetYaxis()->SetTitle("luminosity [10^{50} erg/second]");
   hLt1->SetTitle(sn->GetTitle());
   hLt1->Draw();
   hLt2->Draw("same");
   hLt3->Draw("same");

   leg->Clear();
   leg->SetHeader("");
   leg->SetX1NDC(0.6);
   leg->AddEntry(hNe1, "#nu_{e}");
   leg->AddEntry(hNe2, "#bar{#nu}_{e}");
   leg->AddEntry(hNex, "#nu_{x}");
   leg->Draw();
   can->Print("livermore.ps");

   // draw linear scale after log scale
   // the x axis is still binned in log scale,
   // it looks better also in linear scale
   can->SetLogx(0);
   can->SetLogy(0);

   hLt1->GetXaxis()->SetRangeUser(0.0012,1.5012);
   hLt1->GetYaxis()->SetRangeUser(0,1000);
   hLt1->Draw();
   hLt2->Draw("same");
   hLt3->Draw("same");

   leg->Draw();
   can->Print("livermore.ps");

   // <E>(t)
   sn->SetNbinsT(1000);
   can->SetLogx(1); // set before creating histograms
   can->SetLogy(0); // set before creating histograms

   TH1D *hEt1 = sn->HEt(1);
   TH1D *hEt2 = sn->HEt(2);
   TH1D *hEt3 = sn->HEt(3);

   hEt1->GetXaxis()->SetRangeUser(2.5e-2,17);
   hEt1->GetYaxis()->SetRangeUser(5,30);
   hEt1->GetYaxis()->SetTitle("average energy [MeV/second]");
   hEt1->Draw();
   hEt2->Draw("same");
   hEt3->Draw("same");

   leg->SetX1NDC(0.3);
   leg->SetX2NDC(0.5);
   leg->Draw();
   can->Print("livermore.ps");

   // draw linear scale after log scale
   // the x axis is still binned in log scale,
   // it looks better also in linear scale
   can->SetLogx(0); // set before creating histograms
   can->SetLogy(0); // set before creating histograms

   hEt1->GetXaxis()->SetRangeUser(0.0012,1.5012);
   hEt1->GetYaxis()->SetRangeUser(5,27);
   hEt1->Draw();
   hEt2->Draw("same");
   hEt3->Draw("same");

   leg->SetX1NDC(0.15);
   leg->SetX2NDC(0.35);
   leg->SetY1NDC(0.65);
   leg->SetY2NDC(0.88);
   leg->Draw();
   can->Print("livermore.ps");

   can->Print("livermore.ps]");
}
