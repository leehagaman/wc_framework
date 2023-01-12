#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>
#include<vector>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

//////////////////////////////////////////////////////////////////////////////////

//
// How to run: 
//              root -l cal_CL.cc+
//

void cal_CL_Asimov_Wilks_simplified()
{
  TString str_file_data = "./fc_files/sub_fit_Asimov.root ";

  TFile *file_data = new TFile(str_file_data, "read");

  /////////////////////////////////

  // Declaration of leaf types
  Int_t           grid_Np;
  Int_t           grid_0p;
  Double_t        true_Np;
  Double_t        true_0p;
  vector<int>     *vec_min_status;
  vector<double>  *vec_chi2_var;
  vector<double>  *vec_min_chi2;
  vector<double>  *vec_dchi2;
  vector<double>  *vec_min_fNp_val;
  vector<double>  *vec_min_fNp_err;
  vector<double>  *vec_min_f0p_val;
  vector<double>  *vec_min_f0p_err;

  // List of branches
  TBranch        *b_grid_Np;   //!
  TBranch        *b_grid_0p;   //!
  TBranch        *b_true_Np;   //!
  TBranch        *b_true_0p;   //!
  TBranch        *b_vec_min_status;   //!
  TBranch        *b_vec_chi2_var;   //!
  TBranch        *b_vec_min_chi2;   //!
  TBranch        *b_vec_dchi2;   //!
  TBranch        *b_vec_min_fNp_val;   //!
  TBranch        *b_vec_min_fNp_err;   //!
  TBranch        *b_vec_min_f0p_val;   //!
  TBranch        *b_vec_min_f0p_err;   //!
  
  // Set object pointer
  vec_min_status = 0;
  vec_chi2_var = 0;
  vec_min_chi2 = 0;
  vec_dchi2 = 0;
  vec_min_fNp_val = 0;
  vec_min_fNp_err = 0;
  vec_min_f0p_val = 0;
  vec_min_f0p_err = 0;

  ///////

  TTree *tree_data = (TTree*)file_data->Get("tree");
  tree_data->SetBranchAddress("grid_Np", &grid_Np, &b_grid_Np);
  tree_data->SetBranchAddress("grid_0p", &grid_0p, &b_grid_0p);
  tree_data->SetBranchAddress("true_Np", &true_Np, &b_true_Np);
  tree_data->SetBranchAddress("true_0p", &true_0p, &b_true_0p);
  tree_data->SetBranchAddress("vec_min_status", &vec_min_status, &b_vec_min_status);
  tree_data->SetBranchAddress("vec_chi2_var", &vec_chi2_var, &b_vec_chi2_var);
  tree_data->SetBranchAddress("vec_min_chi2", &vec_min_chi2, &b_vec_min_chi2);
  tree_data->SetBranchAddress("vec_dchi2", &vec_dchi2, &b_vec_dchi2);
  tree_data->SetBranchAddress("vec_min_fNp_val", &vec_min_fNp_val, &b_vec_min_fNp_val);
  tree_data->SetBranchAddress("vec_min_fNp_err", &vec_min_fNp_err, &b_vec_min_fNp_err);
  tree_data->SetBranchAddress("vec_min_f0p_val", &vec_min_f0p_val, &b_vec_min_f0p_val);
  tree_data->SetBranchAddress("vec_min_f0p_err", &vec_min_f0p_err, &b_vec_min_f0p_err);

  ///////////////////////////////// data
  ///////////////////////////////// (1) find the real minimum chi2
  ///////////////////////////////// (2) map chi2 at each grid

  int entries_data = tree_data->GetEntries();
  cout<<endl<<" ---> entries_data: "<<entries_data<<endl<<endl;
  
  double data_min_chi2 = 1e8;
  double data_min_fNp_val = 0;
  double data_min_f0p_val = 0;
  double data_min_fNp_err = 0;
  double data_min_f0p_err = 0;
  
  map<int, map<int, double> >map_data_chi2_var;
  map<int, map<int, double> >map_data_dchi2;
  
  for(int ientry=0; ientry<entries_data; ientry++) {
    tree_data->GetEntry(ientry);

    int vec_size = vec_min_status->size();
    for(int isize=0; isize<vec_size; isize++) {
      if( vec_min_status->at(isize)==0 && (vec_min_chi2->at(isize) < data_min_chi2) ) {
	data_min_chi2 = vec_min_chi2->at(isize);
	if (abs(data_min_chi2) < 1e-9) data_min_chi2 = 0.;  // added to fix a suspected numerical issue
	data_min_fNp_val = vec_min_fNp_val->at(isize);
	data_min_f0p_val = vec_min_f0p_val->at(isize);
	data_min_fNp_err = vec_min_fNp_err->at(isize);
	data_min_f0p_err = vec_min_f0p_err->at(isize);	
      }

      map_data_chi2_var[grid_Np][grid_0p] = vec_chi2_var->at(isize);
      
    }// for(int isize=0; isize<vec_size; isize++)        
  }// for(int ientry=0; ientry<entries_data; ientry++)
  
  if( data_min_chi2==1e8 ) {
    cout<<endl<<" ERROR: No valid data min_chi2"<<endl<<endl;
  }
  
  int number_grid_Np = map_data_chi2_var.size();
  int number_grid_0p = map_data_chi2_var.begin()->second.size();
  
  TH2D *h2d_space = new TH2D("h2d_space", "x for index of Np and y for index of 0p", number_grid_Np, 0, number_grid_Np, number_grid_0p, 0, number_grid_0p);
 
  cout<<number_grid_Np<<"\t"<<number_grid_0p<<endl;
  
  for(auto it_Np=map_data_chi2_var.begin(); it_Np!=map_data_chi2_var.end(); it_Np++) {
    for(auto it_0p=it_Np->second.begin(); it_0p!=it_Np->second.end(); it_0p++) {
      grid_Np = it_Np->first;
      grid_0p = it_0p->first;

      map_data_dchi2[grid_Np][grid_0p] = map_data_chi2_var[grid_Np][grid_0p] - data_min_chi2;

      double dchi2_val = map_data_dchi2[grid_Np][grid_0p];
      if (dchi2_val < 0.) dchi2_val = 0.; // similar line is in read_TLee_v20.cxx, this should only matter for Asimov on the exact grid point
     
      //cout << "lhagaman debug, grids, dchi2, prob: " << grid_Np << ", " << grid_0p << ", " << map_data_dchi2[grid_Np][grid_0p] << ", " << TMath::Prob(map_data_dchi2[grid_Np][grid_0p],2) << "\n";

      double CL_value = 1. - TMath::Prob(dchi2_val,2);

      //cout<<"lhagaman debug, "<<grid_Np<<"\t"<<grid_0p<<"\t"<<map_data_chi2_var[grid_Np][grid_0p]<<"\t"<<data_min_chi2<<"\t"<<CL_value<<endl;


      h2d_space->SetBinContent(grid_Np, grid_0p, CL_value); 
    }
  }

  cout<<TString::Format(" ---> The parameters space is divided into %dx%d grids", number_grid_Np, number_grid_0p)<<endl;
  
  //cout<<endl;
  cout<<" ---> fitting: data_min_chi2 "<<data_min_chi2<<endl;
  cout<<" ---> fitting: data_min_fNp  "<<data_min_fNp_val<<" +/- "<<data_min_fNp_err<<" (error by default Minuit2)"<<endl;
  cout<<" ---> fitting: data_min_f0p  "<<data_min_f0p_val<<" +/- "<<data_min_f0p_err<<" (error by default Minuit2)"<<endl;  
  cout<<endl;

  //TH2D *h2d_space = new TH2D("h2d_space", "x for index of Np and y for index of 0p", number_grid_Np, 0, number_grid_Np, number_grid_0p, 0, number_grid_0p);
  

  h2d_space->SetStats(0);
  h2d_space->Draw("colz");
  h2d_space->SaveAs("file_map_CL_Asimov_Wilks.root");
  
}


