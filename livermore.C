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
   sn->SetDataLocation("/home/jingliu/total");

   // draw spectra
   const Int_t nRGBs = 5;
   const Int_t nContours = 255;

   Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(nRGBs,stops,red,green,blue,nContours);
   gStyle->SetNumberContours(nContours);

   TCanvas *can = new TCanvas;
   can->Print("livermore.ps[");

   can->SetLogz();
   sn->FN2(1)->Draw("colz");
   can->Print("livermore.ps");

   //TH1D *hNe1 = sn->HNe(1);
   //TH1D *hNe2 = sn->HNe(2);
   //TH1D *hNex = sn->HNe(3);

   //hNe1->GetXaxis()->SetRangeUser(0,40);
   //hNe1->GetYaxis()->SetRangeUser(0,4.8e56);
   //hNe1->SetTitle(sn->GetTitle());
   //hNe1->Draw();
   //hNe2->Draw("same");
   //hNex->Scale(4);
   //hNex->Draw("same");

   //TLegend *leg = new TLegend(0.45,0.60,0.88,0.88);
   //leg->SetFillColor(0);
   //leg->SetBorderSize(0);
   //leg->SetHeader("total number of #nu within 18 second:");
   //leg->AddEntry(hNe1, Form("#nu_{e}: %.1e",hNe1->Integral()),"l");
   //leg->AddEntry(hNe2,Form("#bar{#nu}_{e}: %.1e",hNe2->Integral()),"l");
   //leg->AddEntry(hNex, Form("#nu_{x}: %.1e",hNex->Integral()*4),"l");
   //leg->Draw();

   //can->Print("livermore.ps");

   //Double_t tmax=18;/*second*/
   //hNe1 = sn->HNe(1, tmax);
   //hNe2 = sn->HNe(2, tmax);
   //hNex = sn->HNe(3, tmax);

   //hNe1->GetXaxis()->SetRangeUser(0,40);
   //hNe1->GetYaxis()->SetRangeUser(0,4.8e56);
   //hNe1->SetTitle(sn->GetTitle());
   //hNe1->Draw();
   //hNe2->Draw("same");
   //hNex->Scale(4);
   //hNex->Draw("same");

   //leg->Clear();
   //leg->SetHeader(Form("total number of #nu within %.0f second:",tmax));
   //leg->AddEntry(hNe1, Form("#nu_{e}: %.1e",hNe1->Integral()),"l");
   //leg->AddEntry(hNe2,Form("#bar{#nu}_{e}: %.1e",hNe2->Integral()),"l");
   //leg->AddEntry(hNex, Form("#nu_{x}: %.1e",hNex->Integral()*4),"l");
   //leg->Draw();

   //can->Print("livermore.ps");

   //TH1D *hEt1 = sn->HEt(1);
   //TH1D *hEt2 = sn->HEt(2);
   //TH1D *hEt3 = sn->HEt(3);

   //hEt1->GetXaxis()->SetRangeUser(-0.05,0.55);
   //hEt1->GetYaxis()->SetRangeUser(0,25);
   //hEt1->SetTitle(sn->GetTitle());
   //hEt1->Draw();
   //hEt2->Draw("same");
   //hEt3->Draw("same");

   //leg->Clear();
   //leg->SetHeader("");
   //leg->SetX1NDC(0.6);
   //leg->AddEntry(hNe1, "#nu_{e}");
   //leg->AddEntry(hNe2, "#bar{#nu}_{e}");
   //leg->AddEntry(hNex, "#nu_{x}");
   //leg->Draw();
   //can->Print("livermore.ps");

   //hEt1->GetXaxis()->SetRangeUser(0,20);
   //hEt1->GetYaxis()->SetRangeUser(5,20);
   //hEt1->Draw();
   //hEt2->Draw("same");
   //hEt3->Draw("same");

   //leg->Draw();
   //can->Print("livermore.ps");

   //can->SetLogy();
   //TH1D *hLt1 = sn->HLt(1);
   //TH1D *hLt2 = sn->HLt(2);
   //TH1D *hLt3 = sn->HLt(3);

   //hLt1->GetXaxis()->SetRangeUser(-0.5,0.55);
   //hLt1->GetYaxis()->SetRangeUser(1e50,3e53);
   //hLt1->SetTitle(sn->GetTitle());
   //hLt1->Draw();
   //hLt2->Draw("same");
   //hLt3->Draw("same");

   //leg->Draw();
   //can->Print("livermore.ps");

   //hLt1->GetXaxis()->SetRangeUser(0.,20.);
   //hLt1->GetYaxis()->SetRangeUser(1e50,1e52);
   //hLt1->Draw();
   //hLt2->Draw("same");
   //hLt3->Draw("same");

   //leg->Draw();
   //can->Print("livermore.ps");

   can->Print("livermore.ps]");
}
