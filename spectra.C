{
   if (!TClass::GetDict("Supernova")) gSystem->Load("libSNVD.so");

   // load data
   TFile *input = new TFile("models.root");
   if (!input->IsOpen()) {
      cout<<".x ascii2root.C"<<endl;
      gROOT->Macro("ascii2root.C");
      input = new TFile("models.root");
   }
   Supernova *model1301 = (Supernova*) input->Get("model1301");
   Supernova *model1302 = (Supernova*) input->Get("model1302");
   Supernova *model1303 = (Supernova*) input->Get("model1303");
   Supernova *model1311 = (Supernova*) input->Get("model1311");
   Supernova *model1312 = (Supernova*) input->Get("model1312");
   Supernova *model1313 = (Supernova*) input->Get("model1313");
   Supernova *model2001 = (Supernova*) input->Get("model2001");
   Supernova *model2002 = (Supernova*) input->Get("model2002");
   Supernova *model2003 = (Supernova*) input->Get("model2003");
   Supernova *model2011 = (Supernova*) input->Get("model2011");
   Supernova *model2012 = (Supernova*) input->Get("model2012");
   Supernova *model2013 = (Supernova*) input->Get("model2013");
   Supernova *model3001 = (Supernova*) input->Get("model3001");
   Supernova *model3002 = (Supernova*) input->Get("model3002");
   Supernova *model3003 = (Supernova*) input->Get("model3003");
   Supernova *model3011 = (Supernova*) input->Get("model3011");
   Supernova *model3012 = (Supernova*) input->Get("model3012");
   Supernova *model3013 = (Supernova*) input->Get("model3013");
   Supernova *model5001 = (Supernova*) input->Get("model5001");
   Supernova *model5002 = (Supernova*) input->Get("model5002");
   Supernova *model5003 = (Supernova*) input->Get("model5003");
   Supernova *model5011 = (Supernova*) input->Get("model5011");
   Supernova *model5012 = (Supernova*) input->Get("model5012");
   Supernova *model5013 = (Supernova*) input->Get("model5013");

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
