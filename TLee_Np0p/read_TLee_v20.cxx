#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "WCPLEEANA/TLee.h"

#include "WCPLEEANA/Configure_Lee.h"

#include "TApplication.h"

#include <chrono>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*  
usage:

make clean
make  
./read_TLee_v20 -f 1 -p 1

---> README:
---> Makefile: comment the line "ROOTSYS=/home/xji/data0/software/root_build", if you have your own "ROOTSYS"
---> minuit2 is in the ROOT
*/

int main(int argc, char** argv)
{
	TString roostr = "";

	double scaleF_POT = 1;
	int ifile = 1;

	for(int i=1; i<argc; i++) {    
		if( strcmp(argv[i],"-p")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>scaleF_POT ) ) { cerr<<" ---> Error scaleF_POT !"<<endl; exit(1); }
		}

		if( strcmp(argv[i],"-f")==0 ) {
			stringstream convert( argv[i+1] );
			if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
		}
	}

	cout<<endl<<" ---> check, scaleF_POT "<<scaleF_POT<<", ifile "<<ifile<<endl<<endl;

	//////////////////////////////////////////////////////////////////////////////////////// Draw style

	gStyle->SetOptStat(0);
	//gStyle->SetPalette(kBird);

	double snWidth = 2;

	// use medium bold lines and thick markers
	gStyle->SetLineWidth(snWidth);
	gStyle->SetFrameLineWidth(snWidth);
	gStyle->SetHistLineWidth(snWidth);
	gStyle->SetFuncWidth(snWidth);
	gStyle->SetGridWidth(snWidth);
	gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetEndErrorSize(4);
	gStyle->SetEndErrorSize(0);

	if( !config_Lee::flag_display_graphics ) {
		gROOT->SetBatch( 1 );
	}

	////////////////////////////////////////////////////////////////////////////////////////

	TApplication theApp("theApp",&argc,argv);

	TLee *Lee_test = new TLee();

	////////// just do it one time in the whole procedure

	Lee_test->channels_observation   = config_Lee::channels_observation;
	Lee_test->syst_cov_flux_Xs_begin = config_Lee::syst_cov_flux_Xs_begin;
	Lee_test->syst_cov_flux_Xs_end   = config_Lee::syst_cov_flux_Xs_end;
	Lee_test->syst_cov_mc_stat_begin = config_Lee::syst_cov_mc_stat_begin;
	Lee_test->syst_cov_mc_stat_end   = config_Lee::syst_cov_mc_stat_end;  

	Lee_test->flag_lookelsewhere     = config_Lee::flag_lookelsewhere;

	Lee_test->Set_config_file_directory(config_Lee::spectra_file, config_Lee::flux_Xs_directory,
			config_Lee::detector_directory, config_Lee::mc_directory);

	int size_array_LEE_ch = sizeof(config_Lee::array_LEE_ch)/sizeof(config_Lee::array_LEE_ch[0]);
	for(int idx=0; idx<size_array_LEE_ch; idx++) {
		if( config_Lee::array_LEE_ch[idx]!=0 ) Lee_test->map_Lee_ch[config_Lee::array_LEE_ch[idx]] = 1;
	}

	int size_array_LEE_Np_ch = sizeof(config_Lee::array_LEE_Np_ch)/sizeof(config_Lee::array_LEE_Np_ch[0]);
	for(int idx=0; idx<size_array_LEE_Np_ch; idx++) {
		if( config_Lee::array_LEE_Np_ch[idx]!=0 ) Lee_test->map_Lee_Np_ch[config_Lee::array_LEE_Np_ch[idx]] = 1;
	}

	int size_array_LEE_0p_ch = sizeof(config_Lee::array_LEE_0p_ch)/sizeof(config_Lee::array_LEE_0p_ch[0]);
	for(int idx=0; idx<size_array_LEE_0p_ch; idx++) {
		if( config_Lee::array_LEE_0p_ch[idx]!=0 ) Lee_test->map_Lee_0p_ch[config_Lee::array_LEE_0p_ch[idx]] = 1;
	}

	Lee_test->scaleF_POT = scaleF_POT;

	Lee_test->flag_syst_flux_Xs    = config_Lee::flag_syst_flux_Xs;
	Lee_test->flag_syst_detector   = config_Lee::flag_syst_detector;
	Lee_test->flag_syst_additional = config_Lee::flag_syst_additional;
	Lee_test->flag_syst_mc_stat    = config_Lee::flag_syst_mc_stat;
	Lee_test->flag_syst_mc_data_stat_cor    = config_Lee::flag_syst_mc_data_stat_cor;


	// declaring Lee_test as an empty int[4]
	Lee_test->array_no_stat_bins = new int[4];	
	for (int i = 0; i < 4; i++) {
		Lee_test->array_no_stat_bins[i] = config_Lee::array_no_stat_bins[i];
	}

	Lee_test->num_no_stat_bins = config_Lee::num_no_stat_bins;

	Lee_test->Set_Spectra_MatrixCov();
	Lee_test->Set_POT_implement();
	Lee_test->Set_TransformMatrix();

	////////// can do any times

	Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
	Lee_test->scaleF_Lee_Np = config_Lee::Lee_Np_strength_for_outputfile_covariance_matrix;
	Lee_test->scaleF_Lee_0p = config_Lee::Lee_0p_strength_for_outputfile_covariance_matrix;

	// weird that this was just overwritten, commented it out

	//Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
	//Lee_test->scaleF_Lee_Np = config_Lee::Lee_Np_strength_for_GoF;
	//Lee_test->scaleF_Lee_0p = config_Lee::Lee_0p_strength_for_GoF;

	Lee_test->Set_Collapse();

	//////////

	if( 0 ) {// shape only covariance

		int nbins = Lee_test->bins_newworld;
		TMatrixD matrix_pred = Lee_test->matrix_pred_newworld;
		TMatrixD matrix_syst = Lee_test->matrix_absolute_flux_cov_newworld;    
		TMatrixD matrix_shape(nbins, nbins);
		TMatrixD matrix_mixed(nbins, nbins);
		TMatrixD matrix_norm(nbins, nbins);

		///
		double N_T = 0;
		for(int idx=0; idx<nbins; idx++) N_T += matrix_pred(0, idx);

		///
		double M_kl = 0;
		for(int i=0; i<nbins; i++) {
			for(int j=0; j<nbins; j++) {
				M_kl += matrix_syst(i,j);
			}
		}

		///
		for(int i=0; i<nbins; i++) {
			for(int j=0; j<nbins; j++) {      
				double N_i = matrix_pred(0, i);
				double N_j = matrix_pred(0, j);

				double M_ij = matrix_syst(i,j);      
				double M_ik = 0; for(int k=0; k<nbins; k++) M_ik += matrix_syst(i,k);
				double M_kj = 0; for(int k=0; k<nbins; k++) M_kj += matrix_syst(k,j);

				matrix_shape(i,j) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
				matrix_mixed(i,j) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2*N_i*N_j*M_kl/N_T/N_T;       
				matrix_norm(i,j) = N_i*N_j*M_kl/N_T/N_T;
			}
		}

		Lee_test->matrix_absolute_flux_cov_newworld = matrix_shape;

	}// shape only covariance

	/////////////////////////////
	/////////////////////////////

	TFile *file_collapsed_covariance_matrix = new TFile("file_collapsed_covariance_matrix.root", "recreate");

	TTree *tree_config = new TTree("tree", "configure information");
	int flag_syst_flux_Xs = config_Lee::flag_syst_flux_Xs;
	int flag_syst_detector = config_Lee::flag_syst_detector;
	int flag_syst_additional = config_Lee::flag_syst_additional;
	int flag_syst_mc_stat = config_Lee::flag_syst_mc_stat;
	int flag_syst_mc_data_stat_cor = config_Lee::flag_syst_mc_data_stat_cor;
	double user_Lee_strength_for_output_covariance_matrix = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
	double user_scaleF_POT = scaleF_POT;
	vector<double>vc_val_GOF;
	vector<int>vc_val_GOF_NDF;
	tree_config->Branch("flag_syst_flux_Xs", &flag_syst_flux_Xs, "flag_syst_flux_Xs/I" );
	tree_config->Branch("flag_syst_detector", &flag_syst_detector, "flag_syst_detector/I" );
	tree_config->Branch("flag_syst_additional", &flag_syst_additional, "flag_syst_additional/I" );
	tree_config->Branch("flag_syst_mc_stat", &flag_syst_mc_stat, "flag_syst_mc_stat/I" );
	tree_config->Branch("flag_syst_mc_data_stat_cor", &flag_syst_mc_data_stat_cor, "flag_syst_mc_data_stat_cor/I" );
	tree_config->Branch("user_Lee_strength_for_output_covariance_matrix", &user_Lee_strength_for_output_covariance_matrix,
			"user_Lee_strength_for_output_covariance_matrix/D" );
	tree_config->Branch("user_scaleF_POT", &user_scaleF_POT, "user_scaleF_POT/D" );
	tree_config->Branch("vc_val_GOF", &vc_val_GOF);
	tree_config->Branch("vc_val_GOF_NDF", &vc_val_GOF_NDF);
	file_collapsed_covariance_matrix->cd();

	Lee_test->matrix_absolute_cov_newworld.Write("matrix_absolute_cov_newworld");// (bins, bins)
	Lee_test->matrix_absolute_flux_cov_newworld.Write("matrix_absolute_flux_cov_newworld");
	Lee_test->matrix_absolute_Xs_cov_newworld.Write("matrix_absolute_Xs_cov_newworld");
	Lee_test->matrix_absolute_detector_cov_newworld.Write("matrix_absolute_detector_cov_newworld");
	Lee_test->matrix_absolute_mc_stat_cov_newworld.Write("matrix_absolute_mc_stat_cov_newworld");
	Lee_test->matrix_absolute_additional_cov_newworld.Write("matrix_absolute_additional_cov_newworld");

	for(auto it=Lee_test->matrix_input_cov_detector_sub.begin(); it!=Lee_test->matrix_input_cov_detector_sub.end(); it++) {
		int idx = it->first;
		roostr = TString::Format("matrix_absolute_detector_sub_cov_newworld_%02d", idx);
		Lee_test->matrix_absolute_detector_sub_cov_newworld[idx].Write(roostr);
	}

	Lee_test->matrix_pred_newworld.Write("matrix_pred_newworld");// (1, bins)
	Lee_test->matrix_data_newworld.Write("matrix_data_newworld");// (1, bins)

	for(auto it_sub=Lee_test->matrix_sub_flux_geant4_Xs_newworld.begin(); it_sub!=Lee_test->matrix_sub_flux_geant4_Xs_newworld.end(); it_sub++) {
		int index = it_sub->first;
		roostr = TString::Format("matrix_sub_flux_geant4_Xs_newworld_%d", index);
		Lee_test->matrix_sub_flux_geant4_Xs_newworld[index].Write(roostr);
	}

	//file_collapsed_covariance_matrix->Close();

	//////////

	if( config_Lee::flag_plotting_systematics ) Lee_test->Plotting_systematics();

	//////////////////////////////////////////////////////////////////////////////////////// Goodness of fit

	//Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
	//Lee_test->Set_Collapse();

	if( config_Lee::flag_GoF_output2file_default_0 ) {
		file_collapsed_covariance_matrix->cd();

		for(auto it=Lee_test->map_data_spectrum_ch_bin.begin(); it!=Lee_test->map_data_spectrum_ch_bin.end(); it++) {
			int val_ch = it->first;
			int size_map = it->second.size();
			int size_before = 0;
			for(int idx=1; idx<val_ch; idx++) {
				int size_current = Lee_test->map_data_spectrum_ch_bin[idx].size();
				size_before += size_current;
			}

			vector<int>vc_target_chs;
			for(int ibin=1; ibin<size_map; ibin++) {
				vc_target_chs.push_back( size_before + ibin -1 );
			}

			vector<int>vc_support_chs;

			Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 100 + val_ch );

			vc_val_GOF.push_back( Lee_test->val_GOF_noConstrain );
			vc_val_GOF_NDF.push_back( Lee_test->val_GOF_NDF );
			//cout<<" ---> "<<Lee_test->val_GOF_NDF<<"\t"<<Lee_test->val_GOF_noConstrain<<endl;
		}

		tree_config->Fill();
		tree_config->Write();
		file_collapsed_covariance_matrix->Close();
	}

	// bool flag_both_numuCC            = config_Lee::flag_both_numuCC;// 1
	// bool flag_CCpi0_FC_by_numuCC     = config_Lee::flag_CCpi0_FC_by_numuCC;// 2
	// bool flag_CCpi0_PC_by_numuCC     = config_Lee::flag_CCpi0_PC_by_numuCC;// 3
	// bool flag_NCpi0_by_numuCC        = config_Lee::flag_NCpi0_by_numuCC;// 4
	// bool flag_nueCC_PC_by_numuCC_pi0 = config_Lee::flag_nueCC_PC_by_numuCC_pi0;// 5
	// bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = config_Lee::flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC;// 6, HghE>800 MeV
	// bool flag_nueCC_LowE_FC_by_all   = config_Lee::flag_nueCC_LowE_FC_by_all;// 7
	// bool flag_nueCC_FC_by_all        = config_Lee::flag_nueCC_FC_by_all;// 8

	//////////////////////////////////////////////////////////////////////////////////////// user's goodness-of-fit test

	// Lee_test->scaleF_Lee = ?;
	// Lee_test->scaleF_Lee_Np = ?;
	// Lee_test->scaleF_Lee_0p = ?;

	// Lee_test->Set_Collapse();

	///////// validated by 1d = 5, and 2d = [5,5], [1,5], [5,1]

	if( 0 ) {
		Lee_test->scaleF_Lee = 5;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back( 0 );
		vc_target_chs.push_back( 2 );

		vector<int>vc_support_chs;
		for(int ibin=1; ibin<=16*4; ibin++) vc_support_chs.push_back( 4+ibin -1 );

		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 12301 );    
	}

	if( 0 ) {
		Lee_test->scaleF_Lee_Np = 5;
		Lee_test->scaleF_Lee_0p = 5;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back( 0 );
		vc_target_chs.push_back( 2 );

		vector<int>vc_support_chs;
		for(int ibin=1; ibin<=16*4; ibin++) vc_support_chs.push_back( 4+ibin -1 );

		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 12302 );  
	}   

	//////////////////////////////////////////////////////////////////////////////////////// Asimov/Data fitting

	if( 0 ) {
		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 15;
		Lee_test->Set_Collapse();

		///////// Four options:
		//Lee_test->Set_measured_data();// (1), really data
		// Lee_test->Set_toy_Asimov();// (2), Asiomv sample
		// Lee_test->Set_Variations(int num_toy);// (3), generate many toy-MC
		// Lee_test->Set_toy_Variation(int itoy);// which toy-MC to be used   
		// Lee_test->Set_fakedata(TMatrixD matrix_fakedata);// (4), user's defined fakedata

		Lee_test->Set_toy_Asimov();
		
		Lee_test->scaleF_Lee_Np = 1;
                Lee_test->scaleF_Lee_0p = 1;
                Lee_test->Set_Collapse();

		Lee_test->Minimization_Lee_Np_0p_strength_FullCov(4, 4, "");
		// Lee_test->Minimization_Lee_Np_0p_strength_FullCov(5, 4, "Np");// fix Np if the string containts "Np", e.g "NNp"
		// Lee_test->Minimization_Lee_Np_0p_strength_FullCov(4, 5, "0p");// fix 0p if the string containts "0p", e.g. "A_0p"
		// Lee_test->Minimization_Lee_Np_0p_strength_FullCov(5, 5, "Np_0p");// fit both if the string containts "Np" and "0p", e.g. "Np_abc_0p"

		///////// the fitting results/info are saved in
		// Lee_test->minimization_status;// (int type), should be 0 for a successful fitting
		// Lee_test->minimization_chi2;
		// Lee_test->minimization_Lee_Np_strength_val;
		// Lee_test->minimization_Lee_Np_strength_err;// error from the default Minuit2 setting, which is at the points delta_chi2 = chi2_var - chi2_min = 1
		// Lee_test->minimization_Lee_0p_strength_val;
		// Lee_test->minimization_Lee_0p_strength_err;// error from the default Minuit2 setting, which is at the points delta_chi2 = chi2_var - chi2_min = 1

		cout<<endl;
		cout<<" ---> Fitting results"<<endl;
		cout<<" Lee_test->minimization_status: "<<Lee_test->minimization_status<<endl;
		cout<<" Lee_test->minimization_chi2:   "<<Lee_test->minimization_chi2<<endl;
		cout<<" Lee_test->minimization_Lee_Np_strength_val: "<<Lee_test->minimization_Lee_Np_strength_val
			<<" +/- "<<Lee_test->minimization_Lee_Np_strength_err<<endl;
		cout<<" Lee_test->minimization_Lee_0p_strength_val: "<<Lee_test->minimization_Lee_0p_strength_val
			<<" +/- "<<Lee_test->minimization_Lee_0p_strength_err<<endl;        
	}

	//////////////////////////////////////////////////////////////////////////////////////// Feldman-Cousins approach

	/// The (Np, 0p) space are divided into many grid. e.g. 10x10: Np(0, 10), 0p(0, 10) ---> TH2D *h2d_space = new TH2D("h2d_space", "", 10, 0, 10, 10, 0, 10);
	///
	/// At each space point (f_Np, f_0p), we will cacluate the confidence level following Feldman-Cousins approach as following:
	///
	/// (1) Generate the delta_chi2 distribution "distribution_dchi2" by many pseudo experiments, which is corresponding the the pdf of dchi2 at one space point.
	///     * dchi2 = chi2_(f_Np, f_0p) - chi2_min
	///     * pseduo experiments are generated considering both statistical and systematic uncertainties: 
	///
	///     For each pseudo experiment, we will caculate the chi2_(f_Np, f_0p) and chi2_min
	///     ...
	///     By many pseudo experiments, we will have the distribution of dchi2
	///
	/// (2) Calculate the dchi2 of "data", "data" can be real data, or Asimov sample, or others wanted to study
	///     dchi2_data = chi2_(f_Np, f_0p)_data - chi2_min_data
	///
	/// (3) Calculte the pvalue or the confidence value at one space point (f_Np, f_0p)
	///     pvalue: N( distribution_dchi2 > dchi2_data ) / N( distribution_dchi2 )
	///     confidence level = 1 - pvalue
	///

	//////////////////////////////// Scripts:

	bool chi2_grid_checks = false;
	if (chi2_grid_checks) cout << "checking for negative chi2 values...\n";
	if (chi2_grid_checks) {
		for (float i = 0; i <= 20; i += 1){
			for (float j = 0; j <= 20; j += 1) {
				Lee_test->scaleF_Lee_Np = 1;
		                Lee_test->scaleF_Lee_0p = 1;
		                Lee_test->Set_Collapse();// apply the values
                		Lee_test->Set_toy_Asimov();
                		double pars_2d[2] = {i, j};
                		double chi2_var = Lee_test->FCN_Np_0p( pars_2d );// calcualte the chi2 value at the point
				if (chi2_var < 0) cout << "\ntesting chi2 value for asimov data at " << i << ", " << j << "\n";
                		if (chi2_var < 0) cout<<"   negative chi2! "<<chi2_var<<endl;
			}
		}
	}	
	if (chi2_grid_checks) {
		for (float i = 0; i <= 20; i += 1){
			for (float j = 0; j <= 20; j += 1) {
				Lee_test->scaleF_Lee_Np = 1;
		                Lee_test->scaleF_Lee_0p = 1;
		                Lee_test->Set_Collapse();// apply the values
                		Lee_test->Set_measured_data();
				double pars_2d[2] = {i, j};
                		double chi2_var = Lee_test->FCN_Np_0p( pars_2d );// calcualte the chi2 value at the point
				if (chi2_var < 0) cout << "\ntesting chi2 value for real data at " << i << ", " << j << "\n";
                		if (chi2_var < 0) cout<<"   negative chi2! "<<chi2_var<<endl;
			}
		}
	}	
	if (chi2_grid_checks) cout << "done checking for negative chi2 values\n";

	if (0) { 

		cout << "getting chi2 at CV point, all bins...\n";

		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		//Lee_test->Set_measured_data(); // measurement is ready, real data
		Lee_test->Set_toy_Asimov(); // measurement is ready, asimov data, hidden input is prediction

		Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
		// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method

		double pars_2d[2] = {5, 5};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_target_chs.push_back(i);
		}
		vector<int>vc_support_chs;
		vc_support_chs.push_back(1); // empty bin
				
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 111001 );

		cout << "done getting chi2 at CV point\n";

	}

	if (0) { // test GoF vals, coarser. When doing this test, realized that I didn't have Set_Collapse for WC_only

		cout << "printing chi2 at each grid point, WC only...\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}

		for (float i = 0; i <= 4; i += 2){
			for (float j = 0; j <= 4; j += 2) {
				// I think these aren't used, and are reset by pars_2d below
				Lee_test->scaleF_Lee_Np = 1;
		        Lee_test->scaleF_Lee_0p = 1;
				Lee_test->Set_Collapse();
				Lee_test->Set_toy_Asimov();
				Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
				// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method

				double pars_2d[2] = {i, j};
                double chi2_var = Lee_test->FCN_Np_0p( pars_2d );

				cout << "pars_2d: " << i << ", " << j << "\n";
				
				// get constrained chi2 value for all four signal channels, no overflow bins, all constraining bins
				Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );

			}
		}

		cout << "done printing chi2 at each grid point\n";

	}

	bool get_data_grid_point_chi2s = false;
	bool get_asimov_grid_point_chi2s = false;

	if (get_data_grid_point_chi2s or get_asimov_grid_point_chi2s) { // getting GoF at each grid point

		cout << "printing chi2 at each grid point, all bins...\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_target_chs.push_back(i);
		}
		vector<int>vc_support_chs;
		vc_support_chs.push_back(1); // empty bin

		for (float i = 0; i <= 15; i += 0.5){
			for (float j = 0; j <= 15; j += 0.5) {
				// I think these aren't used, and are reset by pars_2d below
				Lee_test->scaleF_Lee_Np = 1;
		        Lee_test->scaleF_Lee_0p = 1;
		        Lee_test->Set_Collapse();// apply the values

				if (get_data_grid_point_chi2s) {
                	Lee_test->Set_measured_data();
				} else {
					Lee_test->Set_toy_Asimov();
					Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
					// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
				}

				double pars_2d[2] = {i, j};
				cout << "i, j = " << i << ", " << j << "\n";
                double chi2_var = Lee_test->FCN_Np_0p( pars_2d );
				cout << "all bins: chi2 = " << chi2_var << "\n";
				
				// get constrained chi2 value for all four signal channels, no overflow bins, all constraining bins
				Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );

			}
		}

		cout << "done printing chi2 at each grid point\n";

	}

	if (get_data_grid_point_chi2s or get_asimov_grid_point_chi2s) { // getting GoF at each grid point

		cout << "printing chi2 at each grid point, WC+gLEE...\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}

		for (float i = 0; i <= 15; i += 0.5){
			for (float j = 0; j <= 15; j += 0.5) {
				// I think these aren't used, and are reset by pars_2d below
				Lee_test->scaleF_Lee_Np = 1;
		        Lee_test->scaleF_Lee_0p = 1;
		        Lee_test->Set_Collapse();// apply the values
				if (get_data_grid_point_chi2s) {
                	Lee_test->Set_measured_data();
				} else {
					Lee_test->Set_toy_Asimov();
					Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
					// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
				}
				double pars_2d[2] = {i, j};
				double chi2_var = Lee_test->FCN_Np_0p( pars_2d );
				
				// get constrained chi2 value for all four signal channels, no overflow bins, all constraining bins
				Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
			}
		}

		cout << "done printing chi2 at each grid point\n";

	}

	if (get_data_grid_point_chi2s or get_asimov_grid_point_chi2s) { // getting GoF at each grid point

		cout << "printing chi2 at each grid point, WC only...\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}

		for (float i = 0; i <= 15; i += 0.5){
			for (float j = 0; j <= 15; j += 0.5) {
				// I think these aren't used, and are reset by pars_2d below
				Lee_test->scaleF_Lee_Np = 1;
		        Lee_test->scaleF_Lee_0p = 1;
				Lee_test->Set_Collapse();// apply the values
				if (get_data_grid_point_chi2s) {
                	Lee_test->Set_measured_data();
				} else {
					Lee_test->Set_toy_Asimov();
					Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
					// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
				}
				double pars_2d[2] = {i, j};
                double chi2_var = Lee_test->FCN_Np_0p( pars_2d );
				
				// get constrained chi2 value for all four signal channels, no overflow bins, all constraining bins
				Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );

			}
		}

		cout << "done printing chi2 at each grid point\n";

	}

	if (get_data_grid_point_chi2s or get_asimov_grid_point_chi2s) { // getting GoF at each grid point

		cout << "printing chi2 at each grid point, gLEE only...\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}

		for (float i = 0; i <= 15; i += 0.5){
			for (float j = 0; j <= 15; j += 0.5) {
				// I think these aren't used, and are reset by pars_2d below
				Lee_test->scaleF_Lee_Np = 1;
		        Lee_test->scaleF_Lee_0p = 1;
		        Lee_test->Set_Collapse();// apply the values
				if (get_data_grid_point_chi2s) {
                	Lee_test->Set_measured_data();
				} else {
					Lee_test->Set_toy_Asimov();
					Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
					// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
				}
				double pars_2d[2] = {i, j};
                double chi2_var = Lee_test->FCN_Np_0p( pars_2d );
				
				// get constrained chi2 value for all four signal channels, no overflow bins, all constraining bins
				Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );

			}
		}

		cout << "done printing chi2 at each grid point\n";

	}

	bool make_sig_bkg_constr_v3_plot = false;
	if (make_sig_bkg_constr_v3_plot) {
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		for (int i=0; i < 2*6*4; i++){
			vc_target_chs.push_back(i);
		}
		vector<int>vc_support_chs;
		for (int i=2*6*4; i < 2*6*4 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 123123 );
	}


	bool make_one_bin_plots = true;

	// WC 1gNp
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333001 );
	}

	// WC 1g0p
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back(2);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333002 );
	}

	// gLEE 1g1p
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back(4);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333003 );
	}

	// gLEE 1g0p
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333004 );
	}

	// WC 1gNp and 1g0p
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333006 );
	}

	// gLEE 1g1p and 1g0p
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333007 );
	}

	// WC 1gNp and 1g0p and gLEE 1g1p and 1g0p
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 15;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333008 );
	}

	// Constraining channels (constrained by signal channels, but shouldn't actually use or care about the constrained result)
	if (make_one_bin_plots){

		Lee_test->scaleF_Lee = 1;
		Lee_test->Set_Collapse();

		vector<int>vc_target_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
		vc_target_chs.push_back(i);
		}
		vector<int>vc_support_chs;
		vc_support_chs.push_back(0);
		vc_support_chs.push_back(2);
		vc_support_chs.push_back(4);
		vc_support_chs.push_back(6);

		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333009 );
	}

	if (make_one_bin_plots) {
		cout << "getting constrained predictions at (1, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333101 );
	}
	if (make_one_bin_plots) {
		cout << "getting constrained predictions at (0, 0) point...\n";
		Lee_test->scaleF_Lee_Np = 0;
		Lee_test->scaleF_Lee_0p = 0;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {0, 0};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333102 );
	}
	if (make_one_bin_plots) {
		cout << "getting constrained predictions at (1, 0) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 0;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 0};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333103 );
	}
	if (make_one_bin_plots) {
		cout << "getting constrained predictions at (0, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 0;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {0, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 333104 );
	}


	bool make_small_set_of_four_bin_excess_plots = false;
	if (make_small_set_of_four_bin_excess_plots) {
		cout << "getting constrained predictions at (1, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 999001001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_excess_plots) {
		cout << "getting constrained predictions at (15, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {15, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 999015001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_excess_plots) {
		cout << "getting constrained predictions at (1, 15) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 15;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 15};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 999001015 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_excess_plots) {
		cout << "getting constrained predictions at (15, 15) point...\n";
		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 15;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {15, 15};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 999015015 );
		cout << "done\n";
	}


	bool make_small_set_of_four_bin_excess_plots_asimov = false;
	if (make_small_set_of_four_bin_excess_plots_asimov) {
		cout << "getting asimov constrained predictions at (1, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_toy_Asimov(); // measurement is ready, asimov data, hidden input is prediction
		Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
		// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
		double pars_2d[2] = {1, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 998001001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_excess_plots_asimov) {
		cout << "getting asimov constrained predictions at (15, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_toy_Asimov(); // measurement is ready, asimov data, hidden input is prediction
		Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
		// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
		double pars_2d[2] = {15, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 998015001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_excess_plots_asimov) {
		cout << "getting asimov constrained predictions at (1, 15) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 15;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_toy_Asimov(); // measurement is ready, asimov data, hidden input is prediction
		Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
		// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
		double pars_2d[2] = {1, 15};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 998001015 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_excess_plots_asimov) {
		cout << "getting asimov constrained predictions at (15, 15) point...\n";
		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 15;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_toy_Asimov(); // measurement is ready, asimov data, hidden input is prediction
		Lee_test->matrix_data_newworld = Lee_test->matrix_pred_newworld; // updating so that constrained chi2s use asimov data
		// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
		double pars_2d[2] = {15, 15};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 998015015 );
		cout << "done\n";
	}

	bool make_small_set_of_four_bin_super_excess_plots = false;
	if (make_small_set_of_four_bin_super_excess_plots) {
		cout << "getting constrained predictions at (1, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 997001001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_super_excess_plots) {
		cout << "getting constrained predictions at (1000, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 1000;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1000, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 997015001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_super_excess_plots) {
		cout << "getting constrained predictions at (1, 1000) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1000;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 1000};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 997001015 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_super_excess_plots) {
		cout << "getting constrained predictions at (1000, 1000) point...\n";
		Lee_test->scaleF_Lee_Np = 1000;
		Lee_test->scaleF_Lee_0p = 1000;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1000, 1000};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 997015015 );
		cout << "done\n";
	}


	bool make_small_set_of_four_bin_super_excess_plots_limited_constraints = false;
	if (make_small_set_of_four_bin_super_excess_plots_limited_constraints) {
		cout << "getting constrained predictions at (1, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8+16*2; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 996001001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_super_excess_plots_limited_constraints) {
		cout << "getting constrained predictions at (1000, 1) point...\n";
		Lee_test->scaleF_Lee_Np = 100;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {100, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8+16*2; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 996015001 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_super_excess_plots_limited_constraints) {
		cout << "getting constrained predictions at (1, 1000) point...\n";
		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 100;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {1, 100};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8+16*2; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 996001015 );
		cout << "done\n";
	}
	if (make_small_set_of_four_bin_super_excess_plots_limited_constraints) {
		cout << "getting constrained predictions at (1000, 1000) point...\n";
		Lee_test->scaleF_Lee_Np = 100;
		Lee_test->scaleF_Lee_0p = 100;
		Lee_test->Set_Collapse(); // prediction is ready
		Lee_test->Set_measured_data(); // measurement is ready, real data
		double pars_2d[2] = {100, 100};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8+16*2; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 996015015 );
		cout << "done\n";
	}


	/*

	if ( constraints_at_certain_bins ) {

		cout << "getting constrained predictions at (0, 1) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {0, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
				
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 112003 );

		cout << "done\n";
	}

	if ( constraints_at_certain_bins ) {

		cout << "getting constrained predictions at (1, 0) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1, 0};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
				
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 112004 );

		cout << "done\n";
	}
	*/


	bool calculate_chi2_15_5_data = false;
	if ( calculate_chi2_15_5_data ) {

		cout << "in calculate_chi2_15_5_data...\ngetting chi2 at (15, 5) point...\n";

		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 5;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {15, 5};

		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		cout << "data vs (15, 5) chi2 = " << chi2_var << "\n";
	}

	bool calculate_chi2_1_1_data = false;
	if ( calculate_chi2_1_1_data ) {

		cout << "in calculate_chi2_1_1_data...\ngetting chi2 at (1, 1) point...\n";

		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1, 1};

		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		cout << "data vs (1, 1) chi2 = " << chi2_var << "\n";
	}

	bool calculate_chi2_1248_0_data = false;
	if ( calculate_chi2_1248_0_data ) {

		cout << "in calculate_chi2_1248_0_data...\ngetting chi2 at (1248, 0) point...\n";

		Lee_test->scaleF_Lee_Np = 1248;
		Lee_test->scaleF_Lee_0p = 0;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1248, 0};

		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		cout << "data vs (1248, 0) chi2 = " << chi2_var << "\n";
	}

	bool calculate_chi2_200_0_data = false;
	if ( calculate_chi2_200_0_data ) {

		cout << "in calculate_chi2_200_0_data...\ngetting chi2 at (200, 0) point...\n";

		Lee_test->scaleF_Lee_Np = 200;
		Lee_test->scaleF_Lee_0p = 0;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {200, 0};

		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		cout << "data vs (200, 0) chi2 = " << chi2_var << "\n";
	}

	bool calculate_chi2_weird_fake_data_point = false;
	if (calculate_chi2_weird_fake_data_point) {
		/*
		From 1000 pseudo-experiments, example of weird point
		(15, 5) toy #283 fake data: 109 0 205 0 155 0 269 0 0 20 251 490 606 571 515 390 386 69 87 67 55 52 20 46 0 208 1028 945 802 480 383 161 51 82 39 36 2 21 5 32 0 9 126 768 1643 2591 3162 3168 2964 2397 1755 1505 1084 744 472 1096 15 28 675 1758 1998 2056 1818 1710 1529 1249 1018 894 518 253 383 614 
		(15, 5) toy #283 vs (15, 5) chi2 = 89.7803
		(15, 5) toy #283 minimization point = (1248.43, 1e-06)
		(15, 5) toy #283 vs (15, 5) chi2_min = 50.7839
		(15, 5) toy #283 vs (15, 5) dchi2 = 38.9964
		*/

		cout << "investigating weird fake data point...\n";

		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 5;
		Lee_test->Set_Collapse(); // prediction is ready
		cout << "set the prediction\n";

		Double_t weird_data[72] = {109, 0, 205, 0, 155, 0, 269, 0, 0, 20, 251, 490, 606, 571, 515, 390, 386, 69, 87, 67, 55, 52, 20, 46, 0, 208, 1028, 945, 802, 480, 383, 161, 51, 82, 39, 36, 2, 21, 5, 32, 0, 9, 126, 768, 1643, 2591, 3162, 3168, 2964, 2397, 1755, 1505, 1084, 744, 472, 1096, 15, 28, 675, 1758, 1998, 2056, 1818, 1710, 1529, 1249, 1018, 894, 518, 253, 383, 614};
		TMatrixD weird_data_TMatrixD = TMatrixD(1, 72, weird_data);
		Lee_test->matrix_data_newworld = weird_data_TMatrixD;
		Lee_test->Set_measured_data(); // measurement is ready, real data
		cout << "set the data\n";

		double pars_2d[2] = {15, 5}; // phase space point for the prediction
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse

		cout << "weird_data vs (15, 5) chi2 = " << chi2_var << "\n";

		// minimization
		/*double initial_Np = 1;
		double initial_0p = 1;
		Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");
		cout << "weird_data vs (15, 5) minimization point = (" << Lee_test->minimization_Lee_Np_strength_val << ", " << Lee_test->minimization_Lee_0p_strength_val << ")\n"; 
		double chi2_min = Lee_test->minimization_chi2;	
		cout << "weird_data vs (15, 5) chi2_min = " << chi2_min << "\n";
		double dchi2 = chi2_var - chi2_min;
		cout << "weird_data vs (15, 5) dchi2 = " << dchi2 << "\n";
		*/

	}


	bool create_chi2_map = false;
	if (create_chi2_map) {
		
		cout << "creating real data chi2 map...\n";

		ofstream data_chi2_map;
		data_chi2_map.open("data_chi2_map.txt");

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1, 1};
		double chi2_var = -1;
		float scale_Np = 0;
		float scale_0p = 0;

		int num_Np_bins = 100;
		int num_0p_bins = 100;
		for (int i = 0; i <= num_Np_bins; i++){
			cout << "done with row " << i << "\n";
			for (int j = 0; j <= num_0p_bins; j++) {
				//scale_Np = 10*i;
				//scale_0p = 10*j;
				scale_Np = i;
				scale_0p = j;
				Lee_test->scaleF_Lee_Np = scale_Np;
				Lee_test->scaleF_Lee_0p = scale_0p;
				Lee_test->Set_Collapse(); // prediction is ready
				pars_2d[0] = scale_Np; pars_2d[1] = scale_0p;
				chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
				data_chi2_map << "(" << scale_Np << ", " << scale_0p << ") : " << chi2_var << "\n";
			}
		}

		data_chi2_map.close();
	}

	bool create_weird_data_chi2_map = false;
	if (create_weird_data_chi2_map) {
		
		cout << "creating weird data chi2 map...\n";

		ofstream weird_data_chi2_map;
		weird_data_chi2_map.open("weird_data_chi2_map.txt");

		Double_t weird_data[72] = {109, 0, 205, 0, 155, 0, 269, 0, 0, 20, 251, 490, 606, 571, 515, 390, 386, 69, 87, 67, 55, 52, 20, 46, 0, 208, 1028, 945, 802, 480, 383, 161, 51, 82, 39, 36, 2, 21, 5, 32, 0, 9, 126, 768, 1643, 2591, 3162, 3168, 2964, 2397, 1755, 1505, 1084, 744, 472, 1096, 15, 28, 675, 1758, 1998, 2056, 1818, 1710, 1529, 1249, 1018, 894, 518, 253, 383, 614};
		TMatrixD weird_data_TMatrixD = TMatrixD(1, 72, weird_data);
		Lee_test->matrix_data_newworld = weird_data_TMatrixD;
		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1, 1};
		double chi2_var = -1;
		float scale_Np = 0;
		float scale_0p = 0;

		int num_Np_bins = 100;
		int num_0p_bins = 100;
		for (int i = 0; i <= num_Np_bins; i++){
			cout << "done with row " << i << "\n";
			for (int j = 0; j <= num_0p_bins; j++) {
				//scale_Np = 10*i;
				//scale_0p = 10*j;
				scale_Np = i;
				scale_0p = j;
				Lee_test->scaleF_Lee_Np = scale_Np;
				Lee_test->scaleF_Lee_0p = scale_0p;
				Lee_test->Set_Collapse(); // prediction is ready
				pars_2d[0] = scale_Np; pars_2d[1] = scale_0p;
				chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
				weird_data_chi2_map << "(" << scale_Np << ", " << scale_0p << ") : " << chi2_var << "\n";
			}
		}

		weird_data_chi2_map.close();
	}


	bool calculate_chi2_and_chi2min_15_5_data_and_toys = false;
	int num_toys = 1000;
	ofstream chi2_and_chi2min_15_5_data_and_toys;
	if (calculate_chi2_and_chi2min_15_5_data_and_toys) {
		chi2_and_chi2min_15_5_data_and_toys.open("chi2_and_chi2min_15_5_data_and_toys.txt");
	}
	if ( calculate_chi2_and_chi2min_15_5_data_and_toys ) { // doing data
		cout << "doing data in calculate_chi2_and_chi2min_15_5_data_and_toys...\n";

		Lee_test->Set_measured_data(); // measurement is ready, real data

		chi2_and_chi2min_15_5_data_and_toys << "data: ";
		for (int i = 0; i < Lee_test->matrix_data_newworld.GetNcols(); i++){
			chi2_and_chi2min_15_5_data_and_toys << Lee_test->matrix_data_newworld[0][i] << " ";
		}
		chi2_and_chi2min_15_5_data_and_toys << "\n";


		double pars_2d[2] = {15, 5}; // phase space point for the prediction

		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		chi2_and_chi2min_15_5_data_and_toys << "data vs (15, 5) chi2 = " << chi2_var << "\n";

		// minimization
		double initial_Np = 1;
		double initial_0p = 1;
		Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");
		chi2_and_chi2min_15_5_data_and_toys  << "data vs (15, 5) minimization point = (" << Lee_test->minimization_Lee_Np_strength_val << ", " << Lee_test->minimization_Lee_0p_strength_val << ")\n"; 
		double chi2_min = Lee_test->minimization_chi2;	
		chi2_and_chi2min_15_5_data_and_toys << "data vs (15, 5) chi2_min = " << chi2_min << "\n";
		double dchi2 = chi2_var - chi2_min;
		chi2_and_chi2min_15_5_data_and_toys << "data vs (15, 5) dchi2 = " << dchi2 << "\n";
	}
	if ( calculate_chi2_and_chi2min_15_5_data_and_toys ) { // doing (15, 5) toys
		cout << "doing (15, 5) toys in calculate_chi2_and_chi2min_15_5_data_and_toys...\n";

		Lee_test->scaleF_Lee_Np = 15;
		Lee_test->scaleF_Lee_0p = 5;
		Lee_test->Set_Collapse(); // prediction is ready

		//int num_toys = 10;

		Lee_test->Set_Variations( num_toys ); // fake data number of toys

		double pars_2d[2] = {15, 5};

		cout << "starting loop over toys...\n";
		for (int itoy = 1; itoy <= num_toys; itoy++){
			Lee_test->Set_toy_Variation( itoy ); //  look at ith toy

			chi2_and_chi2min_15_5_data_and_toys << "(15, 5) toy #" << itoy << " fake data: ";
			for (int i = 0; i < Lee_test->matrix_data_newworld.GetNcols(); i++){
				chi2_and_chi2min_15_5_data_and_toys << Lee_test->map_fake_data[i] << " ";
			}
			chi2_and_chi2min_15_5_data_and_toys << "\n";

			double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
			chi2_and_chi2min_15_5_data_and_toys << "(15, 5) toy #" << itoy << " vs (15, 5) chi2 = " << chi2_var << "\n";

			// minimization	
			double initial_Np = 15;
			double initial_0p = 5;
			Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");
			chi2_and_chi2min_15_5_data_and_toys  << "(15, 5) toy #" << itoy << " minimization point = (" << Lee_test->minimization_Lee_Np_strength_val << ", " << Lee_test->minimization_Lee_0p_strength_val << ")\n"; 
			double chi2_min = Lee_test->minimization_chi2;	
			chi2_and_chi2min_15_5_data_and_toys << "(15, 5) toy #" << itoy << " vs (15, 5) chi2_min = " << chi2_min << "\n";
			double dchi2 = chi2_var - chi2_min;
			chi2_and_chi2min_15_5_data_and_toys << "(15, 5) toy #" << itoy << " vs (15, 5) dchi2 = " << dchi2 << "\n";
		}
		cout << "done\n";
	}
	if ( calculate_chi2_and_chi2min_15_5_data_and_toys ) { // doing (1, 1) toys
		cout << "doing (1, 1) toys in calculate_chi2_and_chi2min_15_5_data_and_toys...\n";

		Lee_test->scaleF_Lee_Np = 1;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		//int num_toys = 10;

		Lee_test->Set_Variations( num_toys ); // fake data number of toys

		double pars_2d[2] = {15, 5};

		cout << "starting loop over toys...\n";
		for (int itoy = 1; itoy <= num_toys; itoy++){
			Lee_test->Set_toy_Variation( itoy ); //  look at ith toy

			chi2_and_chi2min_15_5_data_and_toys << "(1, 1) toy #" << itoy << " fake data: ";
			for (int i = 0; i < Lee_test->matrix_data_newworld.GetNcols(); i++){
				chi2_and_chi2min_15_5_data_and_toys << Lee_test->map_fake_data[i] << " ";
			}
			chi2_and_chi2min_15_5_data_and_toys << "\n";

			double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
			chi2_and_chi2min_15_5_data_and_toys << "(1, 1) toy #" << itoy << " vs (15, 5) chi2 = " << chi2_var << "\n";

			// minimization			
			double initial_Np = 1;
			double initial_0p = 1;		
			Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");
			chi2_and_chi2min_15_5_data_and_toys  << "(1, 1) toy #" << itoy << " minimization point = (" << Lee_test->minimization_Lee_Np_strength_val << ", " << Lee_test->minimization_Lee_0p_strength_val << ")\n"; 
			double chi2_min = Lee_test->minimization_chi2;	
			chi2_and_chi2min_15_5_data_and_toys << "(1, 1) toy #" << itoy << " vs (15, 5) chi2_min = " << chi2_min << "\n";
			double dchi2 = chi2_var - chi2_min;
			chi2_and_chi2min_15_5_data_and_toys << "(1, 1) toy #" << itoy << " vs (15, 5) dchi2 = " << dchi2 << "\n";
		}
		cout << "done\n";
	}
	chi2_and_chi2min_15_5_data_and_toys.close();


	bool test_15_5 = false;

	if ( test_15_5 ) {

		cout << "getting constrained predictions at (15, 5) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {15, 5};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 112011 );

		cout << "done\n";
	}

	if ( test_15_5 ) {

		cout << "(wc_only) getting constrained predictions at (15, 5) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {15, 5};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		//vc_target_chs.push_back(2);
		//vc_target_chs.push_back(4);
		//vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 112012 );

		cout << "done\n";
	}

	if ( test_15_5 ) {

		cout << "(glee_only) getting constrained predictions at (15, 5) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {15, 5};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		//vc_target_chs.push_back(0);
		//vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		//vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 112013 );

		cout << "done\n";
	}

	if ( test_15_5 ) {

		cout << "(wc + constraining bins) getting predictions at (15, 5) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {15, 5};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		//vc_target_chs.push_back(4);
		//vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_target_chs.push_back(i);
		}
		vc_support_chs.push_back(1); // empty bin, just for testing
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 112022 );

		cout << "done\n";
	}


	bool test_15_5_min = false;
	if ( test_15_5_min ) {

		cout << "(wc + constraining bins) getting predictions at (15, 10) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {15, 10};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(0);
		vc_target_chs.push_back(2);
		//vc_target_chs.push_back(4);
		//vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_target_chs.push_back(i);
		}
		vc_support_chs.push_back(1); // empty bin, just for testing
		//Lee_test->Get_Constrained_Chi2_detailed( vc_target_chs, vc_support_chs, 111001 );
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 112022 );

		cout << "done\n";
	}

	bool collapse_1_1000 = false;
	if (collapse_1_1000) { // collapsing to get the error breakdown at an excess point
		cout << "getting predictions at (1, 1000) point...\n";

		Lee_test->scaleF_Lee_Np = 1000;
		Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Plotting_systematics();

		cout << "done\n";
	}

	bool constrained_true_Np_0p_sig_channels = false; // to be used with the sig bkg cnstr setup
	if (constrained_true_Np_0p_sig_channels) {
		cout << "constraining true Np and 0p signal channels...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1000, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		for (int i=0; i < 40; i+=2){
			vc_target_chs.push_back(i);
		}

		vector<int>vc_support_chs;
		for (int i=40; i < 40 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 888002 );

		cout << "done\n";
	}

	bool constrained_sig_channels_to_compare_above = false; // to be used with the normal setup
	if (constrained_sig_channels_to_compare_above) {
		cout << "constraining signal channels...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		for (int i=0; i < 8; i+=2){
			vc_target_chs.push_back(i);
		}

		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 888003 );

		cout << "done\n";
	}

	bool constraining_glee1g0p_with_glee1g1p = false; // to be used with the normal setup
	if (constraining_glee1g0p_with_glee1g1p) {
		cout << "constraining gLEE 1g0p with gLEE 1g1p...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1000, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		vc_target_chs.push_back(6); // gLEE 1g0p
		vector<int>vc_support_chs;
		vc_support_chs.push_back(4); // gLEE 1g1p
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 888004 );

		cout << "done\n";
	}

	bool constraining_glee_only = false;
	if (constraining_glee_only) {
		cout << "constraining gLEE only...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1000, 1};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		//vc_target_chs.push_back(0);
		//vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 888005 );

		cout << "done\n";
	}

	bool constraining_glee_with_nc_pi0 = false;
	if (constraining_glee_with_nc_pi0) {
		cout << "constraining gLEE only using NC Pi0...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1, 1000};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		//vc_target_chs.push_back(0);
		//vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 4*16; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 888006 );

		cout << "done\n";
	}



	bool constrained_glee_15_5_excess = false;
	if (constrained_glee_15_5_excess) {
		cout << "(wc + constraining bins) getting predictions at (15, 10) point...\n";

		//Lee_test->scaleF_Lee_Np = 1;
		//Lee_test->scaleF_Lee_0p = 1;
		Lee_test->Set_Collapse(); // prediction is ready

		Lee_test->Set_measured_data(); // measurement is ready, real data

		double pars_2d[2] = {1, 1000};
		double chi2_var = Lee_test->FCN_Np_0p( pars_2d ); // this re-does Set_Collapse
		//cout << "all bins: chi2 = " << chi2_var << "\n";

		vector<int>vc_target_chs;
		//vc_target_chs.push_back(0);
		//vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}
		Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 777001 );

		cout << "done\n";
	}


	bool various_chi2_checks = false;

	if (various_chi2_checks) {
		int Np_scale_points[] = {0, 7, 15};
		int zero_p_scale_points[] = {0, 7, 15};
		int num_Np_scale_points = sizeof(Np_scale_points) / sizeof(Np_scale_points[0]);
		int num_zero_p_scale_points = sizeof(zero_p_scale_points) / sizeof(zero_p_scale_points[0]);

		double pars_2d[2] = {0, 0};

		double data_chi2 = -1;
		double data_chi2_min = -1;
		double data_dchi2 = -1;
		double data_min_Np_strength_val = -1;
		double data_min_0p_strength_val = -1;

		double asimov_chi2 = -1;
		double asimov_chi2_min = -1;
		double asimov_dchi2 = -1;
		double asimov_min_Np_strength_val = -1;
		double asimov_min_0p_strength_val = -1;

		double scaled_asimov_chi2 = -1;
		double scaled_asimov_chi2_min = -1;
		double scaled_asimov_dchi2 = -1;
		double scaled_asimov_min_Np_strength_val = -1;
		double scaled_asimov_min_0p_strength_val = -1;

		// gLEE only constrained
		vector<int>vc_target_chs;
		//vc_target_chs.push_back(0);
		//vc_target_chs.push_back(2);
		vc_target_chs.push_back(4);
		vc_target_chs.push_back(6);
		vector<int>vc_support_chs;
		for (int i=8; i < 8 + 16 * 4; i++){
			vc_support_chs.push_back(i);
		}

		for (int i = 0; i <num_Np_scale_points; i++) {
			for (int j = 0; j < num_zero_p_scale_points; j++) {
				cout << "chi2 value for (" << Np_scale_points[i] << ", " << zero_p_scale_points[j] << "):\n";

				Lee_test->scaleF_Lee_Np = 1;
				Lee_test->scaleF_Lee_0p = 1;
				Lee_test->Set_Collapse();
				Lee_test->Set_toy_Asimov();
				pars_2d[0] = Np_scale_points[i]; 
				pars_2d[1] = zero_p_scale_points[j];
				asimov_chi2 = Lee_test->FCN_Np_0p( pars_2d );
				Lee_test->Minimization_Lee_Np_0p_strength_FullCov(1, 1, "");
				asimov_chi2_min = Lee_test->minimization_chi2;
				asimov_dchi2 = asimov_chi2 - asimov_chi2_min;
				asimov_min_Np_strength_val = Lee_test->minimization_Lee_Np_strength_val;
				asimov_min_0p_strength_val = Lee_test->minimization_Lee_0p_strength_val;

				Lee_test->scaleF_Lee_Np = 0;
				Lee_test->scaleF_Lee_0p = 15;
				Lee_test->Set_Collapse();
				Lee_test->Set_toy_Asimov();
				pars_2d[0] = Np_scale_points[i]; 
				pars_2d[1] = zero_p_scale_points[j];
				scaled_asimov_chi2 = Lee_test->FCN_Np_0p( pars_2d );
				Lee_test->Minimization_Lee_Np_0p_strength_FullCov(0, 15, "");
				scaled_asimov_chi2_min = Lee_test->minimization_chi2;
				scaled_asimov_dchi2 = scaled_asimov_chi2 - scaled_asimov_chi2_min;
				scaled_asimov_min_Np_strength_val = Lee_test->minimization_Lee_Np_strength_val;
				scaled_asimov_min_0p_strength_val = Lee_test->minimization_Lee_0p_strength_val;

				Lee_test->scaleF_Lee_Np = 1;
				Lee_test->scaleF_Lee_0p = 1;
				Lee_test->Set_Collapse();
				Lee_test->Set_measured_data();
				pars_2d[0] = Np_scale_points[i]; 
				pars_2d[1] = zero_p_scale_points[j];
				data_chi2 = Lee_test->FCN_Np_0p( pars_2d );
				Lee_test->Minimization_Lee_Np_0p_strength_FullCov(1, 1, "");
				data_chi2_min = Lee_test->minimization_chi2;
				data_dchi2 = data_chi2 - data_chi2_min;
				data_min_Np_strength_val = Lee_test->minimization_Lee_Np_strength_val;
				data_min_0p_strength_val = Lee_test->minimization_Lee_0p_strength_val;


				cout << "\n                                     asimov: " << asimov_chi2 << "\n";
				cout << "                                 asimov min: " << asimov_chi2_min << "\n";
				cout << "                               asimov dchi2: " << asimov_dchi2 << "\n";
				cout << "                  asimov min (Np, 0p) scale: (" << asimov_min_Np_strength_val << ", " << asimov_min_0p_strength_val << ")\n";

				cout << "\n                             (0, 15) asimov: " << scaled_asimov_chi2 << "\n";
				cout << "                         (0, 15) asimov min: " << scaled_asimov_chi2_min << "\n";
				cout << "                       (0, 15) asimov dchi2: " << scaled_asimov_dchi2 << "\n";
				cout << "          (0, 15) asimov min (Np, 0p) scale: (" << scaled_asimov_min_Np_strength_val << ", " << scaled_asimov_min_0p_strength_val << ")\n";

				cout << "\n                                     data: " << data_chi2 << "\n";
				cout << "                                 data min: " << data_chi2_min << "\n";
				cout << "                               data dchi2: " << data_dchi2 << "\n";
				cout << "                  data min (Np, 0p) scale: (" << data_min_Np_strength_val << ", " << data_min_0p_strength_val << ")\n";

			}
		}
	}



	// small test before doing full grid
	if( 0 ) {

		///////
		
		int grid_Np = 0;
		int grid_0p = 0;
		double true_Np = 0;
		double true_0p = 0;
		vector<int>vec_min_status;
		vector<double>vec_chi2_var;
		vector<double>vec_min_chi2;
		vector<double>vec_dchi2;
		vector<double>vec_min_fNp_val;
		vector<double>vec_min_fNp_err;
		vector<double>vec_min_f0p_val;
		vector<double>vec_min_f0p_err;    
	
		std::cout << "i file: " << ifile << "\n";
		
		///////
		int Ntoys = 1; // number of toy-MC used to generate the distribution_dchi2
		// changed this from 100 to 1 since it shouldn't matter for Asimov
		// changed this from 100 to 10 to speed it up, same number of toys as for 1D initially


		if (ifile==0 or ifile==-1) { // run it with data or Asimov
			Ntoys = 1;
		}
		
		/////// 2d space of (Np, 0p)
		int bins_Np = 1;
		int bins_0p = 1;
		
		TH2D *h2d_space = new TH2D("h2d_space", "", bins_Np, -0.5, 15.5, bins_0p, -0.5, 15.5);

		double pars_2d[2] = {0};
	

		/////// scan the entrie space defined
		for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++) {
			for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++) {

				cout<<TString::Format(" processing Np/0p: %3d %3d", bin_Np, bin_0p)<<endl;
				
				/// give values of one space point
				double grid_Np_val = h2d_space->GetXaxis()->GetBinCenter( bin_Np );
				double grid_0p_val = h2d_space->GetYaxis()->GetBinCenter( bin_0p );

				grid_Np = bin_Np;
				grid_0p = bin_0p;
				true_Np = grid_Np_val;
				true_0p = grid_0p_val;

				// temp, testing 15,5 point
				grid_Np = 14; // bin 0 means a value of 1
				grid_0p = 4;
				true_Np = 15;
				true_0p = 5;
				grid_Np_val = 15;
				grid_0p_val = 5;

				/*
				grid_Np = bin_Np;
				grid_0p = bin_0p;
				true_Np = grid_Np_val;
				true_0p = grid_0p_val;
				*/
				
				pars_2d[0] = grid_Np_val;
				pars_2d[1] = grid_0p_val;
				cout << "    pars2d: " << pars_2d[0] << ", " << pars_2d[1] << "\n";

			
				/// generate pseudo experiments
				
				if (ifile==0) { // make the data file

					bool do_hybrid_asimov = false;
					if (do_hybrid_asimov) {

						// doing a hybrid data-asimov
						cout << "initial matrix_data_newworld, indices 0, 2, 4, 6, 8, 9, 10, 11, 12:\n";
						cout << Lee_test->matrix_data_newworld[0][0] << ", " << Lee_test->matrix_data_newworld[0][2] << ", " << Lee_test->matrix_data_newworld[0][4] << ", " << Lee_test->matrix_data_newworld[0][6] << ", " << Lee_test->matrix_data_newworld[0][8] << ", " << Lee_test->matrix_data_newworld[0][9] << ", " << Lee_test->matrix_data_newworld[0][10] << ", " << Lee_test->matrix_data_newworld[0][11] << ", " << Lee_test->matrix_data_newworld[0][12] << "\n";
						for (int bin_num_asimov = 8; bin_num_asimov < 2*4+16*4; bin_num_asimov += 1) {
							//cout << "current bin: " << bin_num_asimov << "\n";
							Lee_test->matrix_data_newworld[0][bin_num_asimov] = Lee_test->matrix_pred_newworld[0][bin_num_asimov];
						}
						cout << "updated matrix_data_newworld, indices 0, 2, 4, 6, 8, 9, 10, 11, 12:\n";
						cout << Lee_test->matrix_data_newworld[0][0] << ", " << Lee_test->matrix_data_newworld[0][2] << ", " << Lee_test->matrix_data_newworld[0][4] << ", " << Lee_test->matrix_data_newworld[0][6] << ", " << Lee_test->matrix_data_newworld[0][8] << ", " << Lee_test->matrix_data_newworld[0][9] << ", " << Lee_test->matrix_data_newworld[0][10] << ", " << Lee_test->matrix_data_newworld[0][11] << ", " << Lee_test->matrix_data_newworld[0][12] << "\n";
						// dangerous, could mess up things later in the program, does not save the real measurement, have to do one at a time with this method
					}

					Lee_test->Set_measured_data();
				} else if (ifile==-1) { // make the Asimov file at the Standard Model
					Lee_test->scaleF_Lee_Np = 1.;
					Lee_test->scaleF_Lee_0p = 1.;
					Lee_test->Set_Collapse();// apply the values
					Lee_test->Set_toy_Asimov();
				} else { // make the variations file
					Lee_test->scaleF_Lee_Np = grid_Np_val;
					Lee_test->scaleF_Lee_0p = grid_0p_val;
					Lee_test->Set_Collapse();// apply the values
					Lee_test->Set_Variations( Ntoys );
				}
				
				for(int itoy=1; itoy<=Ntoys; itoy++) {// scan each pseudo experiment
				
					if (ifile > 0) { // only for variation files	
						Lee_test->Set_toy_Variation(itoy);
					}

					std::cout << "getting chi2 for this toy...\n";
					double chi2_var = Lee_test->FCN_Np_0p( pars_2d );// calcualte the chi2 value at the point
					std::cout << "done\n";
								
					double initial_Np = grid_Np_val;
					double initial_0p = grid_0p_val;

					std::cout << "starting minimization...\n";
					
					Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");

					std::cout << "done with minimization\n";


					double chi2_min = Lee_test->minimization_chi2;
					
					double dchi2 = chi2_var - chi2_min;

					cout << "    chi2_var, chi2_min, dchi2: " << chi2_var << ", " << chi2_min << ", " << dchi2 << "\n";

				}// for(int itoy=1; itoy<=Ntoys; itoy++)
			}// for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++)
		}// for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++)
      
  	}




	// sparse grid for debugging
    if( 0 ) {

		///////
		
		int grid_Np = 0;
		int grid_0p = 0;
		double true_Np = 0;
		double true_0p = 0;
		vector<int>vec_min_status;
		vector<double>vec_chi2_var;
		vector<double>vec_min_chi2;
		vector<double>vec_dchi2;
		vector<double>vec_min_fNp_val;
		vector<double>vec_min_fNp_err;
		vector<double>vec_min_f0p_val;
		vector<double>vec_min_f0p_err;    
	
		std::cout << "i file: " << ifile << "\n";

		if (ifile==-2){
			std::cout << "creating mixed Asimov/data root file\n";
			roostr = TString("sub_fit_mixed_Asimov_data.root");
		} else if (ifile==-1) {
			std::cout << "creating asimov root file\n";
			roostr = TString("sub_fit_Asimov.root");
		} else if (ifile==0) {
			roostr = TString("sub_fit_data.root");
		} else {
			roostr = TString::Format("sub_fit_%06d.root", ifile);
		}
		
		TFile *subroofile = new TFile(roostr, "recreate");
		TTree *tree = new TTree("tree", "tree");

		tree->Branch( "grid_Np",          &grid_Np,          "grid_Np/I" );
		tree->Branch( "grid_0p",          &grid_0p,          "grid_0p/I" );
		tree->Branch( "true_Np",          &true_Np,          "true_Np/D" );
		tree->Branch( "true_0p",          &true_0p,          "true_0p/D" );    
		tree->Branch( "vec_min_status",   &vec_min_status );    
		tree->Branch( "vec_chi2_var",     &vec_chi2_var );
		tree->Branch( "vec_min_chi2",     &vec_min_chi2 );
		tree->Branch( "vec_dchi2",        &vec_dchi2 );
		tree->Branch( "vec_min_fNp_val",  &vec_min_fNp_val );
		tree->Branch( "vec_min_fNp_err",  &vec_min_fNp_err );
		tree->Branch( "vec_min_f0p_val",  &vec_min_f0p_val );
		tree->Branch( "vec_min_f0p_err",  &vec_min_f0p_err );

		
		///////
		int Ntoys = 100; // number of toy-MC used to generate the distribution_dchi2

		if (ifile==0 || ifile==-1 || ifile==-2) { // run it with data or Asimov
			Ntoys = 1;
		}
		
		/////// 2d space of (Np, 0p)
		int bins_Np = 3;
		int bins_0p = 3;

		double pars_2d[2] = {0};

		std::cout << "opening spectra text file at " << "spectra_" + std::to_string(ifile) + ".txt" << "\n";
		ofstream spectra_text_file;
		spectra_text_file.open("spectra_" + std::to_string(ifile) + ".txt");

		/////// scan the entrie space defined
		for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++) {
			for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++) {

				cout<<TString::Format(" processing Np/0p: %3d %3d", bin_Np, bin_0p)<<endl;

				double grid_Np_val = -1;
				double grid_0p_val = -1;

				if (bin_Np == 1) grid_Np_val = 1;
				if (bin_Np == 2) grid_Np_val = 7;
				if (bin_Np == 3) grid_Np_val = 15;
				if (bin_0p == 1) grid_0p_val = 1;
				if (bin_0p == 2) grid_0p_val = 7;
				if (bin_0p == 3) grid_0p_val = 15;

				/// 
				grid_Np = bin_Np;
				grid_0p = bin_0p;
				true_Np = grid_Np_val;
				true_0p = grid_0p_val;
				vec_min_status.clear();
				vec_chi2_var.clear();
				vec_min_chi2.clear();
				vec_dchi2.clear();
				vec_min_fNp_val.clear();
				vec_min_fNp_err.clear();
				vec_min_f0p_val.clear();
				vec_min_f0p_err.clear();
				
				pars_2d[0] = grid_Np_val;
				pars_2d[1] = grid_0p_val;
				
				/// generate pseudo experiments
				
				if (ifile==-2) {
					Lee_test->scaleF_Lee_Np = 1.;
					Lee_test->scaleF_Lee_0p = 1.;
					Lee_test->Set_Collapse();
					//Lee_test->Set_first_eight_bins_asimov_rest_measured();
					Lee_test->Set_first_eight_bins_constr_asimov_rest_measured();
				} else if (ifile==0) { // make the data file
					Lee_test->Set_measured_data();
				} else if (ifile==-1) { // make the Asimov file at the Standard Model
					Lee_test->scaleF_Lee_Np = 1.;
					Lee_test->scaleF_Lee_0p = 1.;
					Lee_test->Set_Collapse();// apply the values
					Lee_test->Set_toy_Asimov();
				} else { // make the variations file
					Lee_test->scaleF_Lee_Np = grid_Np_val;
					Lee_test->scaleF_Lee_0p = grid_0p_val;
					Lee_test->Set_Collapse();// apply the values
					Lee_test->Set_Variations( Ntoys );
				}
			
				for(int itoy=1; itoy<=Ntoys; itoy++) {// scan each pseudo experiment
				
					spectra_text_file << "new pseudo-experiment" << grid_Np_val << " " << grid_0p_val << "\n";
				
					if (ifile > 0) { // only for variation files	
						Lee_test->Set_toy_Variation(itoy);
					}

					for (int ibin=0; ibin<72; ibin++) {
						spectra_text_file << Lee_test->matrix_pred_newworld(0, ibin) << " ";
					}
					spectra_text_file << "\n";

					double chi2_var = Lee_test->FCN_Np_0p( pars_2d );// calcualte the chi2 value at the point
				
					//cout << "    pars2d: " << pars_2d[0] << ", " << pars_2d[1] << "\n";

					/////// do minimization: initial value is important. May find local minimum if the values are not suitable
					
					double initial_Np = grid_Np_val;
					double initial_0p = grid_0p_val;
					
					/// a simple way: the true values as the initial value
					//if( 1 ) {
					//  initial_Np = grid_Np_val;
					//  initial_0p = grid_0p_val;
					//}

					/// a more exact way: scan the space to find suitable initial values	  
					
					///////

					Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");

					double chi2_min = Lee_test->minimization_chi2;
					
					double dchi2 = chi2_var - chi2_min;
					
					if (dchi2 < 0.) dchi2 = 0.; // protects against small negative values, mostly relevant for Asimov at the exact grid point

					spectra_text_file << "chi2_var, chi2_min, dchi2: " << chi2_var << ", " << chi2_min << ", " << dchi2 << "\n";

					/////// 
					vec_min_status.push_back( Lee_test->minimization_status );
					vec_chi2_var.push_back( chi2_var );
					vec_min_chi2.push_back( chi2_min );
					vec_dchi2.push_back( dchi2 );
					vec_min_fNp_val.push_back( Lee_test->minimization_Lee_Np_strength_val );
					vec_min_fNp_err.push_back( Lee_test->minimization_Lee_Np_strength_err);
					vec_min_f0p_val.push_back( Lee_test->minimization_Lee_0p_strength_val );
					vec_min_f0p_err.push_back( Lee_test->minimization_Lee_0p_strength_err);
			

				}// for(int itoy=1; itoy<=Ntoys; itoy++)

				////// save the information into root-file for each space point
		
				tree->Fill();
		
			}// for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++)
		}// for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++)

		spectra_text_file.close();

		tree->Write();
		subroofile->Close();
      
  	}
	


	// full grid
    if( 0 ) {

		bool thirty_by_thirty = false;

		///////
		
		int grid_Np = 0;
		int grid_0p = 0;
		double true_Np = 0;
		double true_0p = 0;
		vector<int>vec_min_status;
		vector<double>vec_chi2_var;
		vector<double>vec_min_chi2;
		vector<double>vec_dchi2;
		vector<double>vec_min_fNp_val;
		vector<double>vec_min_fNp_err;
		vector<double>vec_min_f0p_val;
		vector<double>vec_min_f0p_err;    
	
		std::cout << "i file: " << ifile << "\n";

		if (ifile==-2){
			std::cout << "creating mixed Asimov/data root file\n";
			roostr = TString("sub_fit_mixed_Asimov_data_v2.root");
		} else if (ifile==-1) {
			std::cout << "creating asimov root file\n";
			roostr = TString("sub_fit_Asimov.root");
		} else if (ifile==0) {
			roostr = TString("sub_fit_data.root");
		} else {
			roostr = TString::Format("sub_fit_%06d.root", ifile);
		}
		
		TFile *subroofile = new TFile(roostr, "recreate");
		TTree *tree = new TTree("tree", "tree");

		tree->Branch( "grid_Np",          &grid_Np,          "grid_Np/I" );
		tree->Branch( "grid_0p",          &grid_0p,          "grid_0p/I" );
		tree->Branch( "true_Np",          &true_Np,          "true_Np/D" );
		tree->Branch( "true_0p",          &true_0p,          "true_0p/D" );    
		tree->Branch( "vec_min_status",   &vec_min_status );    
		tree->Branch( "vec_chi2_var",     &vec_chi2_var );
		tree->Branch( "vec_min_chi2",     &vec_min_chi2 );
		tree->Branch( "vec_dchi2",        &vec_dchi2 );
		tree->Branch( "vec_min_fNp_val",  &vec_min_fNp_val );
		tree->Branch( "vec_min_fNp_err",  &vec_min_fNp_err );
		tree->Branch( "vec_min_f0p_val",  &vec_min_f0p_val );
		tree->Branch( "vec_min_f0p_err",  &vec_min_f0p_err );

		
		///////
		int Ntoys = 50; // number of toy-MC used to generate the distribution_dchi2
		// changed this from 100 to 1 since it shouldn't matter for Asimov
		// changed this from 100 to 10 to speed it up, same number of toys as for 1D initially


		if (ifile==0 || ifile==-1 || ifile==-2) { // run it with data or Asimov
			Ntoys = 1;
		}
		
		/////// 2d space of (Np, 0p)
		int bins_Np = 16;
		int bins_0p = 16;

		if (thirty_by_thirty) { // actually keeping the same number of bins; 0, 2, 4, ..., 28, 30
			bins_Np = 16;
			bins_0p = 16;
		}
		
		TH2D *h2d_space = new TH2D("h2d_space", "", bins_Np, -0.5, 15.5, bins_0p, -0.5, 15.5);

		if (thirty_by_thirty) {
			h2d_space = new TH2D("h2d_space", "", bins_Np, -1, 31, bins_0p, -1, 31);
		}

		double pars_2d[2] = {0};

		
		/////// scan the entrie space defined
		for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++) {
			for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++) {

				cout<<TString::Format(" processing Np/0p: %3d %3d", bin_Np, bin_0p)<<endl;
				
				/// give values of one space point
				double grid_Np_val = h2d_space->GetXaxis()->GetBinCenter( bin_Np );
				double grid_0p_val = h2d_space->GetYaxis()->GetBinCenter( bin_0p );

				/// 
				grid_Np = bin_Np;
				grid_0p = bin_0p;
				true_Np = grid_Np_val;
				true_0p = grid_0p_val;
				vec_min_status.clear();
				vec_chi2_var.clear();
				vec_min_chi2.clear();
				vec_dchi2.clear();
				vec_min_fNp_val.clear();
				vec_min_fNp_err.clear();
				vec_min_f0p_val.clear();
				vec_min_f0p_err.clear();
				
				pars_2d[0] = grid_Np_val;
				pars_2d[1] = grid_0p_val;
				
				/// generate pseudo experiments
				
				if (ifile==-2) {
					Lee_test->scaleF_Lee_Np = 1.;
					Lee_test->scaleF_Lee_0p = 1.;
					Lee_test->Set_Collapse();
					//Lee_test->Set_first_eight_bins_asimov_rest_measured();
					Lee_test->Set_first_eight_bins_constr_asimov_rest_measured();
				} else if (ifile==0) { // make the data file
					Lee_test->Set_measured_data();
				} else if (ifile==-1) { // make the Asimov file at the Standard Model
					Lee_test->scaleF_Lee_Np = 1.;
					Lee_test->scaleF_Lee_0p = 1.;
					Lee_test->Set_Collapse();// apply the values
					Lee_test->Set_toy_Asimov();
				} else { // make the variations file
					Lee_test->scaleF_Lee_Np = grid_Np_val;
					Lee_test->scaleF_Lee_0p = grid_0p_val;
					Lee_test->Set_Collapse();// apply the values
					Lee_test->Set_Variations( Ntoys );
				}
			
				for(int itoy=1; itoy<=Ntoys; itoy++) {// scan each pseudo experiment

					if (ifile > 0) { // only for variation files	
						Lee_test->Set_toy_Variation(itoy);
					}

					double chi2_var = Lee_test->FCN_Np_0p( pars_2d );// calcualte the chi2 value at the point
				
					//cout << "    pars2d: " << pars_2d[0] << ", " << pars_2d[1] << "\n";

					/////// do minimization: initial value is important. May find local minimum if the values are not suitable
					
					double initial_Np = grid_Np_val;
					double initial_0p = grid_0p_val;
					
					/// a simple way: the true values as the initial value
					//if( 1 ) {
					//  initial_Np = grid_Np_val;
					//  initial_0p = grid_0p_val;
					//}

					/// a more exact way: scan the space to find suitable initial values	  
					
					///////

					Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");

					double chi2_min = Lee_test->minimization_chi2;
					
					double dchi2 = chi2_var - chi2_min;
					
					if (dchi2 < 0.) dchi2 = 0.; // protects against small negative values, mostly relevant for Asimov at the exact grid point

					/////// 
					vec_min_status.push_back( Lee_test->minimization_status );
					vec_chi2_var.push_back( chi2_var );
					vec_min_chi2.push_back( chi2_min );
					vec_dchi2.push_back( dchi2 );
					vec_min_fNp_val.push_back( Lee_test->minimization_Lee_Np_strength_val );
					vec_min_fNp_err.push_back( Lee_test->minimization_Lee_Np_strength_err);
					vec_min_f0p_val.push_back( Lee_test->minimization_Lee_0p_strength_val );
					vec_min_f0p_err.push_back( Lee_test->minimization_Lee_0p_strength_err);
			

				}// for(int itoy=1; itoy<=Ntoys; itoy++)

				////// save the information into root-file for each space point
		
				tree->Fill();
		
			}// for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++)
		}// for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++)

		tree->Write();
		subroofile->Close();
      
  	}
  
	////////////////////////////////////////////////////////////////////////////////////////

	cout<<endl;
	cout<<" -----------------------------------"<<endl;
	cout<<" Check the initialization at the end"<<endl;
	cout<<" -----------------------------------"<<endl;
	
	cout<<endl;
	cout<<" ---> flag_syst_flux_Xs    "<<Lee_test->flag_syst_flux_Xs<<endl;
	cout<<" ---> flag_syst_detector   "<<Lee_test->flag_syst_detector<<endl;
	cout<<" ---> flag_syst_additional "<<Lee_test->flag_syst_additional<<endl;
	cout<<" ---> flag_syst_mc_stat    "<<Lee_test->flag_syst_mc_stat<<endl;
	cout<<" ---> flag_syst_mc_data_stat_cor    "<<Lee_test->flag_syst_mc_data_stat_cor<<endl;
	
	cout<<endl;

	cout<<" ---> LEE channel size (set by array_LEE_ch in config): "<<Lee_test->map_Lee_ch.size()<<endl;
	if( (int)(Lee_test->map_Lee_ch.size()) ) {
		for(auto it_map_Lee=Lee_test->map_Lee_ch.begin(); it_map_Lee!=Lee_test->map_Lee_ch.end(); it_map_Lee++) {
		cout<<" ---> LEE channel: "<< it_map_Lee->first<<endl;
		}
	}
	cout<<endl;

	cout<<" ---> LEE_Np channel size (set by array_LEE_Np_ch in config): "<<Lee_test->map_Lee_Np_ch.size()<<endl;
	if( (int)(Lee_test->map_Lee_Np_ch.size()) ) {
		for(auto it_map_Lee_Np=Lee_test->map_Lee_Np_ch.begin(); it_map_Lee_Np!=Lee_test->map_Lee_Np_ch.end(); it_map_Lee_Np++) {
		cout<<" ---> LEE_Np channel: "<< it_map_Lee_Np->first<<endl;
		}
	}
	cout<<endl;

	cout<<" ---> LEE_0p channel size (set by array_LEE_0p_ch in config): "<<Lee_test->map_Lee_0p_ch.size()<<endl;
	if( (int)(Lee_test->map_Lee_0p_ch.size()) ) {
		for(auto it_map_Lee_0p=Lee_test->map_Lee_0p_ch.begin(); it_map_Lee_0p!=Lee_test->map_Lee_0p_ch.end(); it_map_Lee_0p++) {
		cout<<" ---> LEE_0p channel: "<< it_map_Lee_0p->first<<endl;
		}
	}
	cout<<endl;

	cout<<" ---> MC stat file: "<<Lee_test->syst_cov_mc_stat_begin<<".log - "<<Lee_test->syst_cov_mc_stat_end<<".log"<<endl;
	cout<<endl;
	
	cout<<" ---> check: scaleF_POT "<<scaleF_POT<<", file_flag "<<ifile<<endl<<endl;

	cout<<endl<<endl;
	cout<<" ---> Complete all the program"<<endl;
	cout<<endl<<endl;
	
	if( config_Lee::flag_display_graphics ) {
		cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
		cout<<" Enter Ctrl+c to end the program"<<endl<<endl;
		
		theApp.Run();
	}
	
	return 0;
}
