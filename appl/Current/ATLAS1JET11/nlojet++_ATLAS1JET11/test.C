


{

  TH1D* hnew = (TH1D*)_file0->Get("grid/reference");
  TH1D* hold = (TH1D*)_file1->Get("grid/reference");
  TH1D* hh = (TH1D*)hold->Clone("hh");
    
  hh->Divide(hnew);
  hh->Draw();
    
}


