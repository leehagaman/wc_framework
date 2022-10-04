void test_fid_boundary(){
  TFile *file = new TFile("processed_files/checkout_run1_intrinsic_nue_POT1.2E23.root");
  // TFile *file = new TFile("processed_files/checkout_run1_intrinsic_nue_LowE_POT6.1E23.root");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");

  Float_t truth_vtxX, truth_vtxY, truth_vtxZ;
  T_eval->SetBranchStatus("*",0);
  T_eval->SetBranchStatus("truth_vtxX",1);
  T_eval->SetBranchStatus("truth_vtxY",1);
  T_eval->SetBranchStatus("truth_vtxZ",1);
  T_eval->SetBranchAddress("truth_vtxX",&truth_vtxX);
  T_eval->SetBranchAddress("truth_vtxY",&truth_vtxY);
  T_eval->SetBranchAddress("truth_vtxZ",&truth_vtxZ);

  float bx[2]={100,100};
  float by[2]={0,0};
  float bz[2]={500,500};
  for (Int_t i=0;i!=T_eval->GetEntries();i++){
    T_eval->GetEntry(i);
    if (truth_vtxX < bx[0]) bx[0] = truth_vtxX;
    if (truth_vtxX > bx[1]) bx[1] = truth_vtxX;
    if (truth_vtxY < by[0]) by[0] = truth_vtxY;
    if (truth_vtxY > by[1]) by[1] = truth_vtxY;
    if (truth_vtxZ < bz[0]) bz[0] = truth_vtxZ;
    if (truth_vtxZ > bz[1]) bz[1] = truth_vtxZ;
  }

  std::cout << bx[0] << " " << bx[1] << std::endl;
  std::cout << by[0] << " " << by[1] << std::endl;
  std::cout << bz[0] << " " << bz[1] << std::endl;
  
  
}
