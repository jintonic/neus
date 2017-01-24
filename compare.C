#include "NakazatoModel.h"
#include "LivermoreModel.h"
using namespace NEUS;

#include <TH1D.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLegend.h>

int main()
{
   NakazatoModel *nakazato = new NakazatoModel(20,0.02,200);
   nakazato->LoadData("../neus");

   LivermoreModel *totani = new LivermoreModel;
   totani->LoadData("../total");

   gROOT->SetStyle("Plain");
   TCanvas *c = new TCanvas;
   c->SetRightMargin(0.01);
   c->SetLeftMargin(0.12);
   c->SetTopMargin(0.03);
   c->SetBottomMargin(0.11);
   c->SetLogy();

   // N(E) in 20 second
   TH1D *hNe1 = nakazato->HNe(1);
   TH1D *hNe2 = nakazato->HNe(2);
   TH1D *hNex = nakazato->HNe(3);
   TH1D *hNv1 = totani->HNe(1);
   TH1D *hNv2 = totani->HNe(2);
   TH1D *hNvx = totani->HNe(3);

   hNe1->SetTitle("");
   hNe1->GetXaxis()->SetRangeUser(0,68);
   hNe1->GetXaxis()->SetLabelFont(22);
   hNe1->GetXaxis()->SetLabelSize(0.05);
   hNe1->GetXaxis()->SetTitleFont(22);
   hNe1->GetXaxis()->SetTitleSize(0.05);
   hNe1->GetYaxis()->SetLabelFont(22);
   hNe1->GetYaxis()->SetLabelSize(0.05);
   hNe1->GetYaxis()->SetTitleFont(22);
   hNe1->GetYaxis()->SetTitleSize(0.05);
   hNe1->GetYaxis()->SetTitleOffset(1.2);
   hNe1->GetYaxis()->SetRangeUser(1e3,1e7);
   hNe1->GetYaxis()->SetTitle("total number of neutrino [10^{50}/MeV]");
   hNe1->GetXaxis()->SetTitle("neutrino energy [MeV]");
   hNe1->SetLineWidth(6);
   hNe1->SetLineStyle(kDashed);
   hNe1->Draw();
   hNe2->SetLineWidth(4);
   hNe2->SetLineStyle(kDashed);
   hNe2->Draw("same");
   hNex->SetLineWidth(2);
   hNex->SetLineStyle(kDashed);
   hNex->Draw("same");
   hNvx->SetLineWidth(4);
   hNvx->Draw("same");
   hNv2->SetLineWidth(2);
   hNv2->Draw("same");
   hNv1->Draw("same");

   TLegend *l = new TLegend(0.75,0.5,0.98,0.96);
   l->SetBorderSize(0);
   l->SetTextFont(22);
   l->AddEntry(hNe1, "#nu_{e}, Nakazato","l");
   l->AddEntry(hNe2, "#bar{#nu}_{e}, Nakazato","l");
   l->AddEntry(hNex, "#nu_{x}, Nakazato","l");
   l->AddEntry(hNv1, "#nu_{e}, Totani","l");
   l->AddEntry(hNv2, "#bar{#nu}_{e}, Totani","l");
   l->AddEntry(hNvx, "#nu_{x}, Totani","l");
   l->Draw();

   c->Print("snve.eps");
   c->Print("snve.pdf");
}
