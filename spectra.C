{
   if (!TClass::GetDict("SupernovaModel")) gSystem->Load("libSNVD.so");

   // load data
   TFile *input = new TFile("models.root");
   if (!input->IsOpen()) {
      cout<<".x ascii2root.C"<<endl;
      gROOT->Macro("ascii2root.C");
      input = new TFile("models.root");
   }
   SupernovaModel *model1301 = (SupernovaModel*) input->Get("model1301");
   SupernovaModel *model1302 = (SupernovaModel*) input->Get("model1302");
   SupernovaModel *model1303 = (SupernovaModel*) input->Get("model1303");
   SupernovaModel *model1311 = (SupernovaModel*) input->Get("model1311");
   SupernovaModel *model1312 = (SupernovaModel*) input->Get("model1312");
   SupernovaModel *model1313 = (SupernovaModel*) input->Get("model1313");
   SupernovaModel *model2001 = (SupernovaModel*) input->Get("model2001");
   SupernovaModel *model2002 = (SupernovaModel*) input->Get("model2002");
   SupernovaModel *model2003 = (SupernovaModel*) input->Get("model2003");
   SupernovaModel *model2011 = (SupernovaModel*) input->Get("model2011");
   SupernovaModel *model2012 = (SupernovaModel*) input->Get("model2012");
   SupernovaModel *model2013 = (SupernovaModel*) input->Get("model2013");
   SupernovaModel *model3001 = (SupernovaModel*) input->Get("model3001");
   SupernovaModel *model3002 = (SupernovaModel*) input->Get("model3002");
   SupernovaModel *model3003 = (SupernovaModel*) input->Get("model3003");
   SupernovaModel *model3011 = (SupernovaModel*) input->Get("model3011");
   SupernovaModel *model3012 = (SupernovaModel*) input->Get("model3012");
   SupernovaModel *model3013 = (SupernovaModel*) input->Get("model3013");
   SupernovaModel *model5001 = (SupernovaModel*) input->Get("model5001");
   SupernovaModel *model5002 = (SupernovaModel*) input->Get("model5002");
   SupernovaModel *model5003 = (SupernovaModel*) input->Get("model5003");
   SupernovaModel *model5011 = (SupernovaModel*) input->Get("model5011");
   SupernovaModel *model5012 = (SupernovaModel*) input->Get("model5012");
   SupernovaModel *model5013 = (SupernovaModel*) input->Get("model5013");

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

   can->Print("spectra.ps]");
}
