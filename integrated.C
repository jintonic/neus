{
   // Check if TimeIntegrated is already loaded
   if (!TClass::GetDict("TimeIntegrated")) {
      gROOT->ProcessLine(".L TimeIntegrated.cc++");
   }

   // Use the Class
   TimeIntegrated *data = new TimeIntegrated;
   data->NumberSpectrum("#nu_{e}")->Draw("apl");
   data->NumberSpectrum("#bar{#nu}_{e}")->Draw("pl");
   data->NumberSpectrum("#nu_{x}")->Draw("pl");
}
