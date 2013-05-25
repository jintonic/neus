#include "NakazatoModel.h"
using namespace NEUS;

#include <TH2D.h>
#include <TFile.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
using namespace std;

int main()
{
   // load data
   TFile *input = new TFile("models.root");
   if (!input->IsOpen()) {
      cout<<".x ascii2root.C"<<endl;
      gROOT->Macro("ascii2root.C");
      input = new TFile("models.root");
   }
   NakazatoModel *model1301 = (NakazatoModel*) input->Get("model1301");
   //NakazatoModel *model1302 = (NakazatoModel*) input->Get("model1302");
   //NakazatoModel *model1303 = (NakazatoModel*) input->Get("model1303");
   //NakazatoModel *model1311 = (NakazatoModel*) input->Get("model1311");
   //NakazatoModel *model1312 = (NakazatoModel*) input->Get("model1312");
   //NakazatoModel *model1313 = (NakazatoModel*) input->Get("model1313");
   //NakazatoModel *model2001 = (NakazatoModel*) input->Get("model2001");
   //NakazatoModel *model2002 = (NakazatoModel*) input->Get("model2002");
   //NakazatoModel *model2003 = (NakazatoModel*) input->Get("model2003");
   //NakazatoModel *model2011 = (NakazatoModel*) input->Get("model2011");
   //NakazatoModel *model2012 = (NakazatoModel*) input->Get("model2012");
   //NakazatoModel *model2013 = (NakazatoModel*) input->Get("model2013");
   //NakazatoModel *model3001 = (NakazatoModel*) input->Get("model3001");
   //NakazatoModel *model3002 = (NakazatoModel*) input->Get("model3002");
   //NakazatoModel *model3003 = (NakazatoModel*) input->Get("model3003");
   //NakazatoModel *model3011 = (NakazatoModel*) input->Get("model3011");
   //NakazatoModel *model3012 = (NakazatoModel*) input->Get("model3012");
   //NakazatoModel *model3013 = (NakazatoModel*) input->Get("model3013");
   //NakazatoModel *model5001 = (NakazatoModel*) input->Get("model5001");
   //NakazatoModel *model5002 = (NakazatoModel*) input->Get("model5002");
   //NakazatoModel *model5003 = (NakazatoModel*) input->Get("model5003");
   //NakazatoModel *model5011 = (NakazatoModel*) input->Get("model5011");
   //NakazatoModel *model5012 = (NakazatoModel*) input->Get("model5012");
   //NakazatoModel *model5013 = (NakazatoModel*) input->Get("model5013");

   // draw spectra
   TCanvas *can = new TCanvas;
   can->Print("spectra.ps[");

   TH1D *hNe1 = model1301->HNe(1);
   TH1D *hNe2 = model1301->HNe(2);
   TH1D *hNex = model1301->HNe(3);

   hNe1->GetXaxis()->SetRangeUser(0,40);
   hNe1->GetYaxis()->SetRangeUser(0,4.8e56);
   hNe1->SetTitle(model1301->GetTitle());
   hNe1->Draw();
   hNe2->Draw("same");
   hNex->Scale(4);
   hNex->Draw("same");

   TLegend *leg = new TLegend(0.45,0.60,0.88,0.88);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetHeader("total number of #nu within 20 second:");
   leg->AddEntry(hNe1, Form("#nu_{e}: %.1e",hNe1->Integral()),"l");
   leg->AddEntry(hNe2,Form("#bar{#nu}_{e}: %.1e",hNe2->Integral()),"l");
   leg->AddEntry(hNex, Form("#nu_{x}: %.1e",hNex->Integral()*4),"l");
   leg->Draw();

   can->Print("spectra.ps");

   Double_t tmax=18;/*second*/
   hNe1 = model1301->HNe(1, tmax);
   hNe2 = model1301->HNe(2, tmax);
   hNex = model1301->HNe(3, tmax);

   hNe1->GetXaxis()->SetRangeUser(0,40);
   hNe1->GetYaxis()->SetRangeUser(0,4.8e56);
   hNe1->SetTitle(model1301->GetTitle());
   hNe1->Draw();
   hNe2->Draw("same");
   hNex->Scale(4);
   hNex->Draw("same");

   leg->Clear();
   leg->SetHeader(Form("total number of #nu within %.0f second:",tmax));
   leg->AddEntry(hNe1, Form("#nu_{e}: %.1e",hNe1->Integral()),"l");
   leg->AddEntry(hNe2,Form("#bar{#nu}_{e}: %.1e",hNe2->Integral()),"l");
   leg->AddEntry(hNex, Form("#nu_{x}: %.1e",hNex->Integral()*4),"l");
   leg->Draw();

   can->Print("spectra.ps");

   TH1D *hEt1 = model1301->HEt(1);
   TH1D *hEt2 = model1301->HEt(2);
   TH1D *hEt3 = model1301->HEt(3);

   hEt1->GetXaxis()->SetRangeUser(-0.05,0.55);
   hEt1->GetYaxis()->SetRangeUser(0,25);
   hEt1->SetTitle(model1301->GetTitle());
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
   can->Print("spectra.ps");

   hEt1->GetXaxis()->SetRangeUser(0,20);
   hEt1->GetYaxis()->SetRangeUser(5,20);
   hEt1->Draw();
   hEt2->Draw("same");
   hEt3->Draw("same");

   leg->Draw();
   can->Print("spectra.ps");

   can->SetLogy();
   TH1D *hLt1 = model1301->HLt(1);
   TH1D *hLt2 = model1301->HLt(2);
   TH1D *hLt3 = model1301->HLt(3);

   hLt1->GetXaxis()->SetRangeUser(-0.5,0.55);
   hLt1->GetYaxis()->SetRangeUser(1e50,3e53);
   hLt1->SetTitle(model1301->GetTitle());
   hLt1->Draw();
   hLt2->Draw("same");
   hLt3->Draw("same");

   leg->Draw();
   can->Print("spectra.ps");

   hLt1->GetXaxis()->SetRangeUser(0.,20.);
   hLt1->GetYaxis()->SetRangeUser(1e50,1e52);
   hLt1->Draw();
   hLt2->Draw("same");
   hLt3->Draw("same");

   leg->Draw();
   can->Print("spectra.ps");

   can->Print("spectra.ps]");
}
