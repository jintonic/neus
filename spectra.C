{
   if (!TClass::GetDict("NakazatoModel")) gSystem->Load("libSNVD.so");

   // load data
   TFile *input = new TFile("models.root");
   if (!input->IsOpen()) {
      cout<<".x ascii2root.C"<<endl;
      gROOT->Macro("ascii2root.C");
      input = new TFile("models.root");
   }
   NakazatoModel *model1301 = (NakazatoModel*) input->Get("model1301");
   NakazatoModel *model1302 = (NakazatoModel*) input->Get("model1302");
   NakazatoModel *model1303 = (NakazatoModel*) input->Get("model1303");
   NakazatoModel *model1311 = (NakazatoModel*) input->Get("model1311");
   NakazatoModel *model1312 = (NakazatoModel*) input->Get("model1312");
   NakazatoModel *model1313 = (NakazatoModel*) input->Get("model1313");
   NakazatoModel *model2001 = (NakazatoModel*) input->Get("model2001");
   NakazatoModel *model2002 = (NakazatoModel*) input->Get("model2002");
   NakazatoModel *model2003 = (NakazatoModel*) input->Get("model2003");
   NakazatoModel *model2011 = (NakazatoModel*) input->Get("model2011");
   NakazatoModel *model2012 = (NakazatoModel*) input->Get("model2012");
   NakazatoModel *model2013 = (NakazatoModel*) input->Get("model2013");
   NakazatoModel *model3001 = (NakazatoModel*) input->Get("model3001");
   NakazatoModel *model3002 = (NakazatoModel*) input->Get("model3002");
   NakazatoModel *model3003 = (NakazatoModel*) input->Get("model3003");
   NakazatoModel *model3011 = (NakazatoModel*) input->Get("model3011");
   NakazatoModel *model3012 = (NakazatoModel*) input->Get("model3012");
   NakazatoModel *model3013 = (NakazatoModel*) input->Get("model3013");
   NakazatoModel *model5001 = (NakazatoModel*) input->Get("model5001");
   NakazatoModel *model5002 = (NakazatoModel*) input->Get("model5002");
   NakazatoModel *model5003 = (NakazatoModel*) input->Get("model5003");
   NakazatoModel *model5011 = (NakazatoModel*) input->Get("model5011");
   NakazatoModel *model5012 = (NakazatoModel*) input->Get("model5012");
   NakazatoModel *model5013 = (NakazatoModel*) input->Get("model5013");

   // draw spectra
   TCanvas *can = new TCanvas;
   can->Print("spectra.ps[");

   TH1D *hVe = model1301->IntegratedNumberSpectrum("v_e");
   TH1D *haVe= model1301->IntegratedNumberSpectrum("anti-v_e");
   TH1D *hVx = model1301->IntegratedNumberSpectrum("v_x");

   hVe->GetXaxis()->SetRangeUser(0,50);
   hVe->GetYaxis()->SetRangeUser(0,480e54);
   hVe->SetTitle(model1301->GetTitle());
   hVe->Draw();
   haVe->Draw("same");
   hVx->Scale(4);
   hVx->Draw("same");

   TLegend *leg = new TLegend(0.5,0.30,0.8,0.58);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->SetHeader("total number of #nu:");
   leg->AddEntry(hVe, Form("#nu_{e}: %.1e",hVe->Integral()),"l");
   leg->AddEntry(haVe,Form("#bar{#nu}_{e}: %.1e",haVe->Integral()),"l");
   leg->AddEntry(hVx, Form("#nu_{x}: %.1e",hVx->Integral()*3),"l");
   leg->Draw();

   can->Print("spectra.ps");

   TH2D *hVe2 = model1301->NumberSpectrum("anti-v_e");
   haVe->Draw();
   can->Print("spectra.ps");
   hVe2->ProjectionX()->Draw();
   can->Print("spectra.ps");

   can->Print("spectra.ps]");
}
