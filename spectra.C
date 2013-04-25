{
   if (!TClass::GetDict("SupernovaModel")) gSystem->Load("libSNVD.so");

   // load data
   TFile *input = new TFile("models.root");
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

   model1301->IntegratedEnergySpectrum("v_e")->Draw();
   model1301->IntegratedEnergySpectrum("anti-v_e")->Draw("same");
   model1301->IntegratedEnergySpectrum("v_x")->Draw("same");
   can->Print("spectra.ps");

   //model1301->NumberSpectrum("v_e")->ProfileY()->Draw();
   //model1301->NumberSpectrum("anti-v_e")->ProfileY()->Draw("same");
   //model1301->NumberSpectrum("v_x")->ProfileY()->Draw("same");
   //can->Print("spectra.ps");

   //model1301->IntegratedNumberSpectrum("v_e")->Draw();
   //model1301->IntegratedNumberSpectrum("anti-v_e")->Draw("same");
   //model1301->IntegratedNumberSpectrum("v_x")->Draw("same");
   //can->Print("spectra.ps");

   can->Print("spectra.ps]");
}
