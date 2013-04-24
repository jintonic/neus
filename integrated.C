{
   // check if class SupernovaModel is already loaded
   if (!TClass::GetDict("SupernovaModel"))
      gROOT->ProcessLine(".L SupernovaModel.cc++");

   // initialize the model
   Double_t initialMass, metallicity, reviveTime;
   SupernovaModel *model = new SupernovaModel(
         initialMass=13,/*Solar mass*/
         metallicity=0.02,
         reviveTime=100/*ms*/);
   model->LoadIntegratedData("./integdata");

   // draw spectra
   model->IntegratedNumberSpectrum("#nu_{e}")->Draw();
   model->IntegratedNumberSpectrum("#bar{#nu}_{e}")->Draw("same");
   model->IntegratedNumberSpectrum("#nu_{x}")->Draw("same");

   model->IntegratedEnergySpectrum("#nu_{e}")->Draw();
   model->IntegratedEnergySpectrum("#bar{#nu}_{e}")->Draw("same");
   model->IntegratedEnergySpectrum("#nu_{x}")->Draw("same");
}
