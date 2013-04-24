{
   // Check if SupernovaModel is already loaded
   if (!TClass::GetDict("SupernovaModel")) {
      gROOT->ProcessLine(".L SupernovaModel.cc++");
   }

   // Use the Class
   SupernovaModel *data = new SupernovaModel;
   data->NumberSpectrum("#nu_{e}")->Draw("apl");
   data->NumberSpectrum("#bar{#nu}_{e}")->Draw("pl");
   data->NumberSpectrum("#nu_{x}")->Draw("pl");
}
