{
   if (!TClass::GetDict("NakazatoModel")) gSystem->Load("libSNVD.so");

   NakazatoModel *model1301 = new NakazatoModel(13,0.02,100);
   model1301->LoadIntegratedData("./integdata");
   model1301->LoadFullData("./intpdata");

   NakazatoModel *model1302 = new NakazatoModel(13,0.02,200);
   model1302->LoadIntegratedData("./integdata");
   model1302->LoadFullData("./intpdata");

   NakazatoModel *model1303 = new NakazatoModel(13,0.02,300);
   model1303->LoadIntegratedData("./integdata");
   model1303->LoadFullData("./intpdata");

   NakazatoModel *model1311 = new NakazatoModel(13,0.004,100);
   model1311->LoadIntegratedData("./integdata");
   model1311->LoadFullData("./intpdata");

   NakazatoModel *model1312 = new NakazatoModel(13,0.004,200);
   model1312->LoadIntegratedData("./integdata");
   model1312->LoadFullData("./intpdata");

   NakazatoModel *model1313 = new NakazatoModel(13,0.004,300);
   model1313->LoadIntegratedData("./integdata");
   model1313->LoadFullData("./intpdata");

   NakazatoModel *model2001 = new NakazatoModel(20,0.02,100);
   model2001->LoadIntegratedData("./integdata");
   model2001->LoadFullData("./intpdata");

   NakazatoModel *model2002 = new NakazatoModel(20,0.02,200);
   model2002->LoadIntegratedData("./integdata");
   model2002->LoadFullData("./intpdata");

   NakazatoModel *model2003 = new NakazatoModel(20,0.02,300);
   model2003->LoadIntegratedData("./integdata");
   model2003->LoadFullData("./intpdata");

   NakazatoModel *model2011 = new NakazatoModel(20,0.004,100);
   model2011->LoadIntegratedData("./integdata");
   model2011->LoadFullData("./intpdata");

   NakazatoModel *model2012 = new NakazatoModel(20,0.004,200);
   model2012->LoadIntegratedData("./integdata");
   model2012->LoadFullData("./intpdata");

   NakazatoModel *model2013 = new NakazatoModel(20,0.004,300);
   model2013->LoadIntegratedData("./integdata");
   model2013->LoadFullData("./intpdata");

   NakazatoModel *model3001 = new NakazatoModel(30,0.02,100);
   model3001->LoadIntegratedData("./integdata");
   model3001->LoadFullData("./intpdata");

   NakazatoModel *model3002 = new NakazatoModel(30,0.02,200);
   model3002->LoadIntegratedData("./integdata");
   model3002->LoadFullData("./intpdata");

   NakazatoModel *model3003 = new NakazatoModel(30,0.02,300);
   model3003->LoadIntegratedData("./integdata");
   model3003->LoadFullData("./intpdata");

   NakazatoModel *model3011 = new NakazatoModel(30,0.004,100);
   model3011->LoadIntegratedData("./integdata");
   model3011->LoadFullData("./intpdata");

   NakazatoModel *model3012 = new NakazatoModel(30,0.004,200);
   model3012->LoadIntegratedData("./integdata");
   model3012->LoadFullData("./intpdata");

   NakazatoModel *model3013 = new NakazatoModel(30,0.004,300);
   model3013->LoadIntegratedData("./integdata");
   model3013->LoadFullData("./intpdata");

   NakazatoModel *model5001 = new NakazatoModel(50,0.02,100);
   model5001->LoadIntegratedData("./integdata");
   model5001->LoadFullData("./intpdata");

   NakazatoModel *model5002 = new NakazatoModel(50,0.02,200);
   model5002->LoadIntegratedData("./integdata");
   model5002->LoadFullData("./intpdata");

   NakazatoModel *model5003 = new NakazatoModel(50,0.02,300);
   model5003->LoadIntegratedData("./integdata");
   model5003->LoadFullData("./intpdata");

   NakazatoModel *model5011 = new NakazatoModel(50,0.004,100);
   model5011->LoadIntegratedData("./integdata");
   model5011->LoadFullData("./intpdata");

   NakazatoModel *model5012 = new NakazatoModel(50,0.004,200);
   model5012->LoadIntegratedData("./integdata");
   model5012->LoadFullData("./intpdata");

   NakazatoModel *model5013 = new NakazatoModel(50,0.004,300);
   model5013->LoadIntegratedData("./integdata");
   model5013->LoadFullData("./intpdata");

   TFile *output = new TFile("models.root","recreate");
   model1301->Write();
   model1302->Write();
   model1303->Write();
   model1311->Write();
   model1312->Write();
   model1313->Write();
   model2001->Write();
   model2002->Write();
   model2003->Write();
   model2011->Write();
   model2012->Write();
   model2013->Write();
   model3001->Write();
   model3002->Write();
   model3003->Write();
   model3011->Write();
   model3012->Write();
   model3013->Write();
   model5001->Write();
   model5002->Write();
   model5003->Write();
   model5011->Write();
   model5012->Write();
   model5013->Write();
   output->Close();
}
