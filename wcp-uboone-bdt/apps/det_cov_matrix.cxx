#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" 
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TVectorD.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  
  if (argc < 2){
    std::cout << "./det_cov_matrix -r[#sys 1-10] -g[#gp 0,1]" << std::endl;
  }
  int run = 1; // run 1 ...
  bool flag_osc = false;
  int flag_gp = 0; // gaussian process smoothing
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run = atoi(&argv[i][2]); // which run period
      break;
    case 'o':
      flag_osc = atoi(&argv[i][2]); // which run period
      break;
    case 'g':
      flag_gp = atoi(&argv[i][2]); // 1: on, 0: off.  2 = off but do GPSmoothing debugging, 3 = smoothing and debugging
      break;
    }
  }
  
  CovMatrix cov("./configurations/cov_input.txt", "./configurations/det_input.txt", "./configurations/det_file_ch.txt", "./configurations/rw_cv_input.txt");
  
  // appling the detvar background approximation to the signal channels:
  cov.add_disabled_ch_name("wc_1gNp_overlay");
  cov.add_disabled_ch_name("wc_1g0p_overlay");
  cov.add_disabled_ch_name("glee_1g1p_overlay");
  cov.add_disabled_ch_name("glee_1g0p_overlay");
  
  cov.add_disabled_ch_name("wc_1gNp_cc_pi0_overlay");
  cov.add_disabled_ch_name("wc_1g0p_cc_pi0_overlay");
  cov.add_disabled_ch_name("glee_1g1p_cc_pi0_overlay");
  cov.add_disabled_ch_name("glee_1g0p_cc_pi0_overlay");
  
  cov.add_disabled_ch_name("wc_1gNp_nue_overlay");
  cov.add_disabled_ch_name("wc_1g0p_nue_overlay");
  cov.add_disabled_ch_name("glee_1g1p_nue_overlay");
  cov.add_disabled_ch_name("glee_1g0p_nue_overlay");
  
  
  //cov.add_disabled_ch_name("nc_delta_Np_01_overlay");
  //cov.add_disabled_ch_name("nc_delta_0p_01_overlay");

  //cov.add_disabled_ch_name("joint_BG_nueCC_Enu_FC_nu_overlay");
  //cov.add_disabled_ch_name("joint_BG_nueCC_Enu_PC_nu_overlay");

  
  //cov.add_disabled_ch_name("joint2_BG_nueCC_Enu_FC_nu_overlay");
  //cov.add_disabled_ch_name("joint2_BG_nueCC_Enu_PC_nu_overlay");
  
  
  //cov.add_disabled_ch_name("BG_nueCC_FC_overlay");
  //cov.add_disabled_ch_name("BG_nueCC_PC_overlay");
  // cov.add_disabled_ch_name("nueCC_FC_nueoverlay");
  // cov.add_disabled_ch_name("nueCC_PC_nueoverlay");
  // cov.add_disabled_ch_name("LEE_FC_nueoverlay");
  // cov.add_disabled_ch_name("LEE_PC_nueoverlay");
  
  //cov.add_disabled_ch_name("BG_nueCC2_FC_overlay");
  //cov.add_disabled_ch_name("BG_nueCC2_PC_overlay");
  //cov.add_disabled_ch_name("BG_nueCC3_FC_overlay");
  //cov.add_disabled_ch_name("BG_nueCC3_PC_overlay");
  //cov.add_disabled_ch_name("BG_nueCC_extra_FC_overlay");
  //cov.add_disabled_ch_name("BG_nueCC_extra_PC_overlay");
  
  
  if (flag_osc) cov.add_osc_config();

  cov.print_rw(cov.get_rw_info());

  // Get the file based on runno ...
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();
  // Construct the histogram ...

  // outfilename ...
  TString outfile_name;

  TH1F *htemp;
  std::map<TString, TH1F*> map_histoname_hist;
  std::map<int, TH1F*> map_covch_hist;
  
  for (auto it = map_inputfile_info.begin(); it!=map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);

    
    if (period == run){
      outfile_name = out_filename;
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename, 0);
      
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	TString histoname = std::get<0>(*it1);
	Int_t nbin = std::get<1>(*it1);
	float llimit = std::get<2>(*it1);
	float hlimit = std::get<3>(*it1);
	TString var_name = std::get<4>(*it1);
	TString ch_name = std::get<5>(*it1);
	TString add_cut = std::get<6>(*it1);
	TString weight = std::get<7>(*it1);
	
	//	std::cout << std::get<0>(*it1)  << " " << std::get<1>(*it1) << " " << std::get<4>(*it1) << " " << std::get<5>(*it1) << " " << std::get<6>(*it1) << " " << std::get<7>(*it1) << std::endl;
	htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
	map_histoname_hist[histoname] = htemp;

	int covch = cov.get_covch_name(ch_name);
	if (map_covch_hist.find(covch) == map_covch_hist.end()){
	  TH1F *htemp1 = (TH1F*)htemp->Clone(Form("pred_covch_%d",covch));
	  map_covch_hist[covch] = htemp1;
	}
      }
      //  std::cout << input_filename << " " << filetype << " " << out_filename << std::endl; 
    }
  }
  std::cout << outfile_name << std::endl;

  TMatrixD* cov_add_mat = cov.get_add_cov_matrix();
  // create a covariance matrix for bootstrapping ...
  TMatrixD* cov_mat_bootstrapping = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  // create a covariance matrix for det systematics ...
  TMatrixD* cov_det_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  TVectorD* vec_mean_diff = new TVectorD(cov_add_mat->GetNrows());
  TVectorD* vec_mean = new TVectorD(cov_add_mat->GetNrows());
  
  cov.gen_det_cov_matrix(run, map_covch_hist, map_histoname_hist, vec_mean, vec_mean_diff, cov_mat_bootstrapping, cov_det_mat, flag_gp);

  TMatrixD* frac_cov_det_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  for (size_t i=0; i!= frac_cov_det_mat->GetNrows(); i++){
    double val_1 = (*vec_mean)(i);
    for (size_t j=0; j!=frac_cov_det_mat->GetNrows();j++){
      double val_2 = (*vec_mean)(j);
      double val = (*cov_det_mat)(i,j);
      if (val_1 ==0 && val_2 == 0){
	(*frac_cov_det_mat)(i,j) = 0;
      }else if (val_1 ==0 || val_2 ==0){
	if (val !=0){
	  if (i==j){
	    (*frac_cov_det_mat)(i,j) = 1./16.; // 25% uncertainties ...
	  }else{
	    (*frac_cov_det_mat)(i,j) = 0;
	  }
	}else{
	  (*frac_cov_det_mat)(i,j) = 0;
	}
      }else{
	(*frac_cov_det_mat)(i,j) = val/val_1/val_2;
      }
    }
  }
  
  
  TFile *file = new TFile(outfile_name,"RECREATE");
  vec_mean->Write(Form("vec_mean_%d",run));
  vec_mean_diff->Write(Form("vec_mean_diff_%d",run));
  
  cov_mat_bootstrapping->Write(Form("cov_mat_boostrapping_%d",run));
  cov_det_mat->Write(Form("cov_det_mat_%d",run));
  frac_cov_det_mat->Write(Form("frac_cov_det_mat_%d",run));
  

  // save central ... results ...
  // for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
  //  ((TH1F*)it->second)->SetDirectory(file);
  // }
  for (auto it = map_covch_hist.begin(); it != map_covch_hist.end(); it++){
    ((TH1F*)it->second)->SetDirectory(file);
  }
  
  file->Write();
  file->Close();
}
