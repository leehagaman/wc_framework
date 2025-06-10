#include "TRandom.h"
#include "TRandom3.h"

void LEEana::CovMatrix::gen_xf_cov_matrix(int run, std::map<int, TH1F*>& map_covch_hist, std::map<TString, TH1F*>& map_histoname_hist, TVectorD* vec_mean,  TMatrixD* cov_xf_mat){
  // prepare the maps ... name --> no,  covch, lee
  std::map<TString, std::tuple<int, int, int, TString>> map_histoname_infos ; 
  std::map<int, TString> map_no_histoname; 
  
  int ncount = 0;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    
    if (period != run) continue;
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);
    
    for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
      int ch = map_name_ch[std::get<5>(*it1)];
      int obsch = get_obsch_name(std::get<5>(*it1));
      int covch = get_covch_name(std::get<5>(*it1));
      int flag_lee = std::get<7>(map_ch_hist[ch]);
      flag_lee = 0; // turning off all lee scaling for XsFlux/DetVar, since it uses eLEE weights and we're making uncollapsed cov matrices
      
      TString histoname = std::get<0>(*it1);
      TH1F *htemp = map_histoname_hist[histoname];
      //
      map_histoname_infos[histoname] = std::make_tuple(ncount, covch, flag_lee, input_filename);
      map_no_histoname[ncount] = histoname;
      ncount ++;

      // std::cout << histoname << " " << obsch << " " << covch << " " << flag_lee << std::endl;
    }
  }

  // now prepare the output ...
  // filename ...   # events               #weight  #leeweight #difference   #different types
  std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > > map_passed_events; // last one is variable name ...
  std::map<TString, double> map_filename_pot;
  std::vector<int> max_lengths;
  std::vector<int> max_sup_lengths;
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    //int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    if (period != run) continue;

    
    //map_all_events[input_filename];
    std::pair<std::vector<int>, std::vector<int>> lengths_pair = get_events_weights(input_filename, map_passed_events, map_filename_pot, map_histoname_infos);
    std::vector<int> lengths = lengths_pair.first;
    std::vector<int> sup_lengths = lengths_pair.second;

    if (lengths.size() > max_lengths.size()) max_lengths.resize(lengths.size());
    for (size_t i = 0; i != lengths.size();i++){
      if (lengths.at(i) > max_lengths.at(i)) max_lengths.at(i) = lengths.at(i);
    }
    
    if (sup_lengths.size() > max_sup_lengths.size()) max_sup_lengths.resize(sup_lengths.size());
    for (size_t i = 0; i != sup_lengths.size();i++){
      if (sup_lengths.at(i) > max_sup_lengths.at(i)) max_sup_lengths.at(i) = sup_lengths.at(i);
    }
    
    //std::cout << input_filename << " " << lengths.size() << std::endl;
  }
  
  double data_pot = 5e19;
  const int rows = cov_xf_mat->GetNcols();
  float x[rows];
  (*cov_xf_mat).Zero();
 
  int acc_no = 0;
  for (size_t j = 0; j!=max_lengths.size(); j++){

    //    if (j>0) continue;
    
    int nsize = max_lengths.at(j);
    int sup_nsize = max_sup_lengths.at(j);

    //    std::cout << j << " " << nsize << " " << sup_nsize << std::endl;
    
    TMatrixD temp_mat(rows, rows);
    temp_mat.Zero(); 
    
    
    for (int i=0;i!=nsize;i++){

      //      if (i>=10) continue;
      
      for (int k = 0; k!= rows;k++){
	x[k] = 0;
      }
      fill_xf_histograms(j, max_lengths.size(), acc_no, i, nsize,  map_passed_events, map_histoname_infos, map_no_histoname, map_histoname_hist);

      // merge histograms according to POTs ...
      for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
        std::cout << it->first << std::endl;
        int covch = it->first;
        //std::cout << "looking for " << covch << " in map_covch_hist...";
        TH1F *hpred = map_covch_hist[covch];
        //std::cout << " found it\n";
        hpred->Reset();
        
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
          TH1F *htemp = (TH1F*)hpred->Clone("htemp");
          htemp->Reset();
          std::map<int, double> temp_map_mc_acc_pot;
          
          for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
            TString histoname = (*it2).first;
            TString input_filename = map_histogram_inputfile[histoname];
            auto it3 = map_inputfile_info.find(input_filename);
            int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
            int norm_period = std::get<6>(it3->second);
            double mc_pot = map_filename_pot[input_filename];
            //std::cout << mc_pot << std::endl;
            if (temp_map_mc_acc_pot.find(norm_period) == temp_map_mc_acc_pot.end()){
              temp_map_mc_acc_pot[norm_period] = mc_pot;
            }else{
              temp_map_mc_acc_pot[norm_period] += mc_pot;
            }
          }
          
          for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
            TString histoname = (*it2).first;
            TString input_filename = map_histogram_inputfile[histoname];
            auto it3 = map_inputfile_info.find(input_filename);
            int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
            int norm_period = std::get<6>(it3->second);
            data_pot = std::get<5>(map_inputfile_info[input_filename]);
            double ratio = data_pot/temp_map_mc_acc_pot[norm_period];
            
            TH1F *hmc = map_histoname_hist[histoname];
            htemp->Add(hmc, ratio);
            //	std::cout << covch << " " << histoname << " " << ratio << std::endl;
          }
          
          hpred->Add(htemp);
          delete htemp;
        }

        int start_bin = map_covch_startbin[covch];
        for (int k=0;k!=hpred->GetNbinsX()+1;k++){
          x[start_bin+k] = hpred->GetBinContent(k+1) ;
          //	  std::cout << i << " " << x[start_bin+i] << std::endl;
        }
      }
      // add covariance matrix ...
      for (size_t n = 0;n!=rows; n++){
        for (size_t m =0; m!=rows;m++){
          temp_mat(n,m) += x[n] * x[m];
        }
      }
    } // i

    
    if (nsize==2){  // second check
      temp_mat *= 1./sup_nsize;
    }else{
      temp_mat *= 1./nsize;
    }

    //std::cout << j << " " << temp_mat(26+26+5,26+26+26) << std::endl;
    
    (*cov_xf_mat) +=  temp_mat;

    //std::cout << "lhagaman debug, j, temp_mat(0,0), temp_mat(2,2), temp_mat(4,4), temp_mat(6,6): " << j << ", " << temp_mat(0,0) << ", " << temp_mat(2,2) << ", " << temp_mat(4,4) << ", " << temp_mat(6,6) << "\n";
    
    
    acc_no += nsize;
  } //j
  
  // calculate the CV ...
  for (int i=0;i!=rows;i++){
    (*vec_mean)(i) = 0;
  }

  fill_xf_histograms(map_passed_events, map_histoname_infos, map_no_histoname, map_histoname_hist);

  // merge histograms according to POTs ...
  for (auto it = map_pred_covch_histos.begin(); it!=map_pred_covch_histos.end();it++){
    //std::cout << it->first << std::endl;
    int covch = it->first;
    TH1F *hpred = map_covch_hist[covch];
    hpred->Reset();
    
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      TH1F *htemp = (TH1F*)hpred->Clone("htemp");
      htemp->Reset();
      std::map<int, double> temp_map_mc_acc_pot;
      
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	TString histoname = (*it2).first;
	TString input_filename = map_histogram_inputfile[histoname];
	auto it3 = map_inputfile_info.find(input_filename);
	int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
	int norm_period = std::get<6>(it3->second);
	double mc_pot = map_filename_pot[input_filename];
	//std::cout << mc_pot << std::endl;
	if (temp_map_mc_acc_pot.find(norm_period) == temp_map_mc_acc_pot.end()){
	  temp_map_mc_acc_pot[norm_period] = mc_pot;
	}else{
	  temp_map_mc_acc_pot[norm_period] += mc_pot;
	}
	//std::cout << histoname << " " << input_filename << " " << mc_pot << " " << period << std::endl;
      }

      //      std::cout << "haha " << std::endl;
      
      for (auto it2 = it1->begin(); it2 != it1->end(); it2++){
	TString histoname = (*it2).first;
	TString input_filename = map_histogram_inputfile[histoname];
	auto it3 = map_inputfile_info.find(input_filename);
	int period = std::get<1>(it3->second);  if (period != run) continue; // skip ...
	int norm_period = std::get<6>(it3->second);
	data_pot = std::get<5>(map_inputfile_info[input_filename]);
	double ratio = data_pot/temp_map_mc_acc_pot[norm_period];
	
	TH1F *hmc = map_histoname_hist[histoname];
	htemp->Add(hmc, ratio);
	
	//	std::cout << covch << " " << histoname << " " << ratio << " " << data_pot << std::endl;
      }
      
      hpred->Add(htemp);
      delete htemp;
    }
    
    int start_bin = map_covch_startbin[covch];
    for (int k=0;k!=hpred->GetNbinsX()+1;k++){
      (*vec_mean)(start_bin+k) = hpred->GetBinContent(k+1) ;
      	  //std::cout << i << " " << x[start_bin+i] << std::endl;
    }
  }
  
  
  
}




void LEEana::CovMatrix::fill_xf_histograms(int num, int tot_num, int acc_no, int no, int tot_no, std::map<TString, 
    std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > >& map_passed_events, 
    std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  
    std::map<TString, TH1F*>& map_histoname_hist){
  for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
     it->second->Reset();
   }

  //  std::cout << acc_no << " " << no << std::endl;
  for (auto it = map_passed_events.begin(); it != map_passed_events.end(); it++){
    TString filename = it->first;
    // loop over events ...
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      float weight = std::get<0>(*it1);
      float weight_lee = std::get<1>(*it1);
      if (std::get<3>(*it1).size() != tot_num) std::cout << "Incorrect Match Sys No! " << std::endl;
      if (std::get<3>(*it1).at(num) != tot_no) std::cout << "Mismatch No of Universe! " << std::endl;
      float rel_weight_diff = std::get<2>(*it1).at(acc_no+no);
      for (auto it2 = std::get<4>(*it1).begin(); it2 != std::get<4>(*it1).end(); it2++){
	int no = (*it2).first;
	float val = (*it2).second;

	TString histoname = map_no_histoname[no];
	TH1F *htemp = map_histoname_hist[histoname];
	int flag_lee = std::get<2>(map_histoname_infos[histoname]);
  flag_lee = 0; // turning off all lee scaling for XsFlux/DetVar, since it uses eLEE weights and we're making uncollapsed cov matrices

	if (std::isnan(rel_weight_diff) || std::isinf(rel_weight_diff)) continue;
	// seems to have extremely small cv weight
	if (fabs(rel_weight_diff)>100) continue;
	

	if (flag_lee){
	  htemp->Fill(val, rel_weight_diff * weight * weight_lee);
	}else{
	  htemp->Fill(val, rel_weight_diff * weight);
	}
	
      }
    }
  }
  
}

void LEEana::CovMatrix::fill_xf_histograms(std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > >& map_passed_events, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos, std::map<int, TString>& map_no_histoname,  std::map<TString, TH1F*>& map_histoname_hist){
  for (auto it = map_histoname_hist.begin(); it != map_histoname_hist.end(); it++){
     it->second->Reset();
   }

  //std::cout << acc_no << " " << no << std::endl;
  for (auto it = map_passed_events.begin(); it != map_passed_events.end(); it++){
    TString filename = it->first;
    // loop over events ...
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      float weight = std::get<0>(*it1);
      float weight_lee = std::get<1>(*it1);
      for (auto it2 = std::get<4>(*it1).begin(); it2 != std::get<4>(*it1).end(); it2++){
	int no = (*it2).first;
	float val = (*it2).second;

	TString histoname = map_no_histoname[no];
	TH1F *htemp = map_histoname_hist[histoname];
	int flag_lee = std::get<2>(map_histoname_infos[histoname]);
  flag_lee = 0; // turning off all lee scaling for XsFlux/DetVar, since it uses eLEE weights and we're making uncollapsed cov matrices
	
	if (flag_lee){
	  htemp->Fill(val, weight * weight_lee);
	}else{
	  htemp->Fill(val,  weight);
	}
	
      }
    }
  }
  
}



std::pair<std::vector<int>, std::vector<int> > LEEana::CovMatrix::get_events_weights(TString input_filename, std::map<TString, std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > > >& map_passed_events, std::map<TString, double>& map_filename_pot, std::map<TString, std::tuple<int, int, int, TString>>& map_histoname_infos){
  TFile *file = new TFile(input_filename);

  TTree *T_BDTvars = (TTree*)file->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file->Get("wcpselection/T_KINEvars");

  EvalInfo eval;
  POTInfo pot;
  TaggerInfo tagger;
  PFevalInfo pfeval;
  KineInfo kine;

  kine.kine_energy_particle = new std::vector<float>;
  kine.kine_energy_info = new std::vector<int>;
  kine.kine_particle_type = new std::vector<int>;
  kine.kine_energy_included = new std::vector<int>;
  
  
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger.pio_2_v_dis2 = new std::vector<float>;
  tagger.pio_2_v_angle2 = new std::vector<float>;
  tagger.pio_2_v_acc_length = new std::vector<float>;
  tagger.pio_2_v_flag = new std::vector<float>;
  tagger.sig_1_v_angle = new std::vector<float>;
  tagger.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger.sig_1_v_energy = new std::vector<float>;
  tagger.sig_1_v_energy_1 = new std::vector<float>;
  tagger.sig_1_v_flag = new std::vector<float>;
  tagger.sig_2_v_energy = new std::vector<float>;
  tagger.sig_2_v_shower_angle = new std::vector<float>;
  tagger.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger.sig_2_v_flag = new std::vector<float>;
  tagger.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_energy = new std::vector<float>;
  tagger.stw_2_v_angle = new std::vector<float>;
  tagger.stw_2_v_dir_length = new std::vector<float>;
  tagger.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger.stw_2_v_flag = new std::vector<float>;
  tagger.stw_3_v_angle = new std::vector<float>;
  tagger.stw_3_v_dir_length = new std::vector<float>;
  tagger.stw_3_v_energy = new std::vector<float>;
  tagger.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger.stw_3_v_flag = new std::vector<float>;
  tagger.stw_4_v_angle = new std::vector<float>;
  tagger.stw_4_v_dis = new std::vector<float>;
  tagger.stw_4_v_energy = new std::vector<float>;
  tagger.stw_4_v_flag = new std::vector<float>;
  tagger.br3_3_v_energy = new std::vector<float>;
  tagger.br3_3_v_angle = new std::vector<float>;
  tagger.br3_3_v_dir_length = new std::vector<float>;
  tagger.br3_3_v_length = new std::vector<float>;
  tagger.br3_3_v_flag = new std::vector<float>;
  tagger.br3_5_v_dir_length = new std::vector<float>;
  tagger.br3_5_v_total_length = new std::vector<float>;
  tagger.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger.br3_5_v_n_seg = new std::vector<float>;
  tagger.br3_5_v_angle = new std::vector<float>;
  tagger.br3_5_v_sg_length = new std::vector<float>;
  tagger.br3_5_v_energy = new std::vector<float>;
  tagger.br3_5_v_n_main_segs = new std::vector<float>;
  tagger.br3_5_v_n_segs = new std::vector<float>;
  tagger.br3_5_v_shower_main_length = new std::vector<float>;
  tagger.br3_5_v_shower_total_length = new std::vector<float>;
  tagger.br3_5_v_flag = new std::vector<float>;
  tagger.br3_6_v_angle = new std::vector<float>;
  tagger.br3_6_v_angle1 = new std::vector<float>;
  tagger.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger.br3_6_v_direct_length = new std::vector<float>;
  tagger.br3_6_v_length = new std::vector<float>;
  tagger.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger.br3_6_v_energy = new std::vector<float>;
  tagger.br3_6_v_flag = new std::vector<float>;
  tagger.tro_1_v_particle_type = new std::vector<float>;
  tagger.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger.tro_1_v_min_dis = new std::vector<float>;
  tagger.tro_1_v_sg1_length = new std::vector<float>;
  tagger.tro_1_v_shower_main_length = new std::vector<float>;
  tagger.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger.tro_1_v_tmp_length = new std::vector<float>;
  tagger.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger.tro_1_v_flag = new std::vector<float>;
  tagger.tro_2_v_energy = new std::vector<float>;
  tagger.tro_2_v_stem_length = new std::vector<float>;
  tagger.tro_2_v_iso_angle = new std::vector<float>;
  tagger.tro_2_v_max_length = new std::vector<float>;
  tagger.tro_2_v_angle = new std::vector<float>;
  tagger.tro_2_v_flag = new std::vector<float>;
  tagger.tro_4_v_dir2_mag = new std::vector<float>;
  tagger.tro_4_v_angle = new std::vector<float>;
  tagger.tro_4_v_angle1 = new std::vector<float>;
  tagger.tro_4_v_angle2 = new std::vector<float>;
  tagger.tro_4_v_length = new std::vector<float>;
  tagger.tro_4_v_length1 = new std::vector<float>;
  tagger.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger.tro_4_v_energy = new std::vector<float>;
  tagger.tro_4_v_shower_main_length = new std::vector<float>;
  tagger.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger.tro_4_v_flag = new std::vector<float>;
  tagger.tro_5_v_max_angle = new std::vector<float>;
  tagger.tro_5_v_min_angle = new std::vector<float>;
  tagger.tro_5_v_max_length = new std::vector<float>;
  tagger.tro_5_v_iso_angle = new std::vector<float>;
  tagger.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger.tro_5_v_min_count = new std::vector<float>;
  tagger.tro_5_v_max_count = new std::vector<float>;
  tagger.tro_5_v_energy = new std::vector<float>;
  tagger.tro_5_v_flag = new std::vector<float>;
  tagger.lol_1_v_energy = new std::vector<float>;
  tagger.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_1_v_nseg = new std::vector<float>;
  tagger.lol_1_v_angle = new std::vector<float>;
  tagger.lol_1_v_flag = new std::vector<float>;
  tagger.lol_2_v_length = new std::vector<float>;
  tagger.lol_2_v_angle = new std::vector<float>;
  tagger.lol_2_v_type = new std::vector<float>;
  tagger.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger.lol_2_v_energy = new std::vector<float>;
  tagger.lol_2_v_shower_main_length = new std::vector<float>;
  tagger.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger.lol_2_v_flag = new std::vector<float>;
  tagger.cosmict_flag_10 = new std::vector<float>;
  tagger.cosmict_10_flag_inside = new std::vector<float>;
  tagger.cosmict_10_vtx_z = new std::vector<float>;
  tagger.cosmict_10_flag_shower = new std::vector<float>;
  tagger.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger.cosmict_10_angle_beam = new std::vector<float>;
  tagger.cosmict_10_length = new std::vector<float>;
  tagger.numu_cc_flag_1 = new std::vector<float>;
  tagger.numu_cc_1_particle_type = new std::vector<float>;
  tagger.numu_cc_1_length = new std::vector<float>;
  tagger.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger.numu_cc_1_direct_length = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger.numu_cc_flag_2 = new std::vector<float>;
  tagger.numu_cc_2_length = new std::vector<float>;
  tagger.numu_cc_2_total_length = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger.numu_cc_2_n_daughter_all = new std::vector<float>;

  set_tree_address(T_BDTvars, tagger,2 );
  set_tree_address(T_eval, eval);
  set_tree_address(T_PFeval, pfeval);
  set_tree_address(T_pot, pot);
  set_tree_address(T_KINEvars, kine);

  double total_pot = 0;
  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    total_pot += pot.pot_tor875;
  }
  // total POT calculations ...
  map_filename_pot[input_filename] = total_pot;

  // fill histogram ...
  T_BDTvars->SetBranchStatus("*",0);
  T_BDTvars->SetBranchStatus("numu_cc_flag",1);
  T_BDTvars->SetBranchStatus("numu_score",1);
  T_BDTvars->SetBranchStatus("nue_score",1);
  T_BDTvars->SetBranchStatus("cosmict_flag",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_0",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_1",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_2",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_3",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_4",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_5",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_6",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_7",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_8",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_9",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_10",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_11",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_12",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_13",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_14",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_15",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_16",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_17",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_18",1);
  T_BDTvars->SetBranchStatus("mip_vec_dQ_dx_19",1);
  T_BDTvars->SetBranchStatus("mip_energy",1);
  T_BDTvars->SetBranchStatus("mip_angle_beam",1);
  T_BDTvars->SetBranchStatus("spt_angle_vertical",1);
  T_BDTvars->SetBranchStatus("mip_quality_n_tracks",1);
  T_BDTvars->SetBranchStatus("mip_quality_n_showers",1);
  T_BDTvars->SetBranchStatus("gap_n_bad",1);
  T_BDTvars->SetBranchStatus("spt_angle_beam",1);
  T_BDTvars->SetBranchStatus("spt_angle_vertical",1);

   if (tagger.flag_nc_gamma_bdt){
    T_BDTvars->SetBranchStatus("nc_delta_score", 1);
    T_BDTvars->SetBranchStatus("nc_pio_score", 1);
  }
  
  T_eval->SetBranchStatus("*",0);
  T_eval->SetBranchStatus("run",1);
  T_eval->SetBranchStatus("subrun",1);
  T_eval->SetBranchStatus("event",1);
  T_eval->SetBranchStatus("match_energy",1);
  T_eval->SetBranchStatus("match_isFC",1);
  T_eval->SetBranchStatus("match_found",1);
  if (T_eval->GetBranch("match_found_asInt")) T_eval->SetBranchStatus("match_found_asInt",1); 
  T_eval->SetBranchStatus("stm_eventtype",1);
  T_eval->SetBranchStatus("stm_lowenergy",1);
  T_eval->SetBranchStatus("stm_LM",1);
  T_eval->SetBranchStatus("stm_TGM",1);
  T_eval->SetBranchStatus("stm_STM",1);
  T_eval->SetBranchStatus("stm_FullDead",1);
  T_eval->SetBranchStatus("stm_clusterlength",1);
  
  T_eval->SetBranchStatus("weight_spline",1);
  T_eval->SetBranchStatus("weight_cv",1);
  //T_eval->SetBranchStatus("weight_lee",1);
  T_eval->SetBranchStatus("weight_change",1);
  // MC enable truth information ...
  T_eval->SetBranchStatus("truth_isCC",1);
  T_eval->SetBranchStatus("truth_nuPdg",1);
  T_eval->SetBranchStatus("truth_vtxInside",1);
  T_eval->SetBranchStatus("truth_nuEnergy",1);
  T_eval->SetBranchStatus("truth_energyInside",1);
  T_eval->SetBranchStatus("truth_vtxX",1);
  T_eval->SetBranchStatus("truth_vtxY",1);
  T_eval->SetBranchStatus("truth_vtxZ",1);
  T_eval->SetBranchStatus("match_completeness_energy",1);

  bool flag_glee_merge = true;
  if (flag_glee_merge) {
      T_eval->SetBranchStatus("gl_sel_type",1);
      T_eval->SetBranchStatus("gl_sel_type",1);
      T_eval->SetBranchStatus("gl_true_Enu",1);
      T_eval->SetBranchStatus("gl_true_Elep",1);
      T_eval->SetBranchStatus("gl_reco_Eshower",1);
      T_eval->SetBranchStatus("gl_reco_shower_dirz",1);
      T_eval->SetBranchStatus("gl_reco_shower_implied_dirz",1);
      T_eval->SetBranchStatus("gl_simple_pot_weight",1);
      T_eval->SetBranchStatus("gl_rem_orig_wc_pot_weight",1);
      T_eval->SetBranchStatus("gl_new_pot_weight",1);
      T_eval->SetBranchStatus("gl_overlap_weight",1);
      T_eval->SetBranchStatus("gl_overlap_weight_tweaked",1);
      T_eval->SetBranchStatus("gl_wc_total_overlapped_weight",1);
  }
  T_eval->SetBranchStatus("run",1);
  T_eval->SetBranchStatus("subrun",1);
  T_eval->SetBranchStatus("event",1);

  
  T_KINEvars->SetBranchStatus("*",0);
  T_KINEvars->SetBranchStatus("kine_reco_Enu",1);
  T_KINEvars->SetBranchStatus("kine_energy_particle",1);
  T_KINEvars->SetBranchStatus("kine_particle_type",1);
  T_KINEvars->SetBranchStatus("kine_energy_info",1);
  T_KINEvars->SetBranchStatus("kine_energy_included",1);
  T_KINEvars->SetBranchStatus("kine_reco_add_energy",1);
  T_KINEvars->SetBranchStatus("kine_pio_mass",1);
  T_KINEvars->SetBranchStatus("kine_pio_flag",1);
  T_KINEvars->SetBranchStatus("kine_pio_vtx_dis",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_1",1);
  T_KINEvars->SetBranchStatus("kine_pio_energy_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_theta_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_phi_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_dis_2",1);
  T_KINEvars->SetBranchStatus("kine_pio_angle",1);
  if (T_KINEvars->GetBranch("vlne_v4_numu_full_primaryE")) {
    T_KINEvars->SetBranchStatus("vlne_v4_numu_full_primaryE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_full_totalE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_partial_primaryE",1);
    T_KINEvars->SetBranchStatus("vlne_v4_numu_partial_totalE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_full_primaryE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_full_totalE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_partial_primaryE",1);
    // T_KINEvars->SetBranchStatus("vlne_nue_partial_totalE",1);
  }

  T_PFeval->SetBranchStatus("*",0);
  T_PFeval->SetBranchStatus("reco_nuvtxX",1);
  T_PFeval->SetBranchStatus("reco_nuvtxY",1);
  T_PFeval->SetBranchStatus("reco_nuvtxZ",1);
  T_PFeval->SetBranchStatus("reco_showervtxX",1);
  T_PFeval->SetBranchStatus("reco_showervtxY",1);
  T_PFeval->SetBranchStatus("reco_showervtxZ",1);
  T_PFeval->SetBranchStatus("reco_muonMomentum",1);
  T_PFeval->SetBranchStatus("reco_showerKE",1);
  T_PFeval->SetBranchStatus("nuvtx_diff",1);
  T_PFeval->SetBranchStatus("showervtx_diff",1);
  T_PFeval->SetBranchStatus("muonvtx_diff",1);
  T_PFeval->SetBranchStatus("truth_muonMomentum",1);
  if (pfeval.flag_NCDelta){
    T_PFeval->SetBranchStatus("truth_NCDelta",1);
    T_PFeval->SetBranchStatus("truth_NprimPio",1);
  }
  if (pfeval.flag_recoprotonMomentum){
    T_PFeval->SetBranchStatus("reco_protonMomentum",1);
  }
  if (pfeval.flag_showerMomentum){
    T_PFeval->SetBranchStatus("reco_showerMomentum",1);
    T_PFeval->SetBranchStatus("reco_Nproton",1);
    T_PFeval->SetBranchStatus("truth_showerMomentum",1);
    T_PFeval->SetBranchStatus("truth_nuScatType",1);
    // oscillation formula ...
    T_PFeval->SetBranchStatus("truth_nu_momentum",1);
    T_PFeval->SetBranchStatus("neutrino_type",1);
    T_PFeval->SetBranchStatus("mcflux_ntype",1);
    T_PFeval->SetBranchStatus("mcflux_dk2gen",1);
    T_PFeval->SetBranchStatus("mcflux_gen2vtx",1);
    T_PFeval->SetBranchStatus("mcflux_ndecay",1);
  }
  if (T_PFeval->GetBranch("truth_startMomentum")){
    T_PFeval->SetBranchStatus("truth_Ntrack",1);
    T_PFeval->SetBranchStatus("truth_pdg",1); 
    T_PFeval->SetBranchStatus("truth_mother",1); 
    T_PFeval->SetBranchStatus("truth_startMomentum",1); 
  }

  WeightInfo weight;
  TTree *T_weight = (TTree*)file->Get("wcpselection/T_weight");
  weight.file_type = new std::string();
  weight.expskin_FluxUnisim= new std::vector<float>;
  weight.horncurrent_FluxUnisim= new std::vector<float>;
  weight.kminus_PrimaryHadronNormalization= new std::vector<float>;
  weight.kplus_PrimaryHadronFeynmanScaling= new std::vector<float>;
  weight.kzero_PrimaryHadronSanfordWang= new std::vector<float>;
  weight.nucleoninexsec_FluxUnisim= new std::vector<float>;
  weight.nucleonqexsec_FluxUnisim= new std::vector<float>;
  weight.nucleontotxsec_FluxUnisim= new std::vector<float>;
  weight.piminus_PrimaryHadronSWCentralSplineVariation= new std::vector<float>;
  weight.pioninexsec_FluxUnisim= new std::vector<float>;
  weight.pionqexsec_FluxUnisim= new std::vector<float>;
  weight.piontotxsec_FluxUnisim= new std::vector<float>;
  weight.piplus_PrimaryHadronSWCentralSplineVariation= new std::vector<float>;
  
  weight.All_UBGenie= new std::vector<float>;
  weight.AxFFCCQEshape_UBGenie= new std::vector<float>;
  weight.DecayAngMEC_UBGenie= new std::vector<float>;
  weight.NormCCCOH_UBGenie= new std::vector<float>;
  weight.NormNCCOH_UBGenie= new std::vector<float>;
  weight.RPA_CCQE_Reduced_UBGenie= new std::vector<float>;
  weight.RPA_CCQE_UBGenie= new std::vector<float>;
  weight.RootinoFix_UBGenie= new std::vector<float>;
  weight.ThetaDelta2NRad_UBGenie= new std::vector<float>;
  weight.Theta_Delta2Npi_UBGenie= new std::vector<float>;
  weight.TunedCentralValue_UBGenie= new std::vector<float>;
  weight.VecFFCCQEshape_UBGenie= new std::vector<float>;
  weight.XSecShape_CCMEC_UBGenie= new std::vector<float>;
  weight.splines_general_Spline= new std::vector<float>;
  weight.xsr_scc_Fa3_SCC= new std::vector<float>;
  weight.xsr_scc_Fv3_SCC= new std::vector<float>;

  weight.reinteractions_piminus_Geant4 = new std::vector<float>;
  weight.reinteractions_piplus_Geant4 = new std::vector<float>;
  weight.reinteractions_proton_Geant4 = new std::vector<float>; 
  
  TString option;
  if (rw_type == 1){
    option = "reweight";
  }else if (rw_type == 2){
    option = "reweight_cor";
  }else if (T_weight->GetBranch("expskin_FluxUnisim")){
    option = "expskin_FluxUnisim";
  }else if (T_weight->GetBranch("horncurrent_FluxUnisim")){
    option = "horncurrent_FluxUnisim";
  }else if (T_weight->GetBranch("kminus_PrimaryHadronNormalization")){
    option = "kminus_PrimaryHadronNormalization";
  }else if (T_weight->GetBranch("kplus_PrimaryHadronFeynmanScaling")){
    option = "kplus_PrimaryHadronFeynmanScaling";
  }else if (T_weight->GetBranch("kzero_PrimaryHadronSanfordWang")){
    option = "kzero_PrimaryHadronSanfordWang";
  }else if (T_weight->GetBranch("nucleoninexsec_FluxUnisim")){
    option = "nucleoninexsec_FluxUnisim";
  }else if (T_weight->GetBranch("nucleonqexsec_FluxUnisim")){
    option = "nucleonqexsec_FluxUnisim";
  }else if (T_weight->GetBranch("nucleontotxsec_FluxUnisim")){
    option = "nucleontotxsec_FluxUnisim";
  }else if (T_weight->GetBranch("piminus_PrimaryHadronSWCentralSplineVariation")){
    option = "piminus_PrimaryHadronSWCentralSplineVariation";
  }else if (T_weight->GetBranch("pioninexsec_FluxUnisim")){
    option = "pioninexsec_FluxUnisim";
  }else if (T_weight->GetBranch("pionqexsec_FluxUnisim")){
    option = "pionqexsec_FluxUnisim";
  }else if (T_weight->GetBranch("piontotxsec_FluxUnisim")){
    option = "piontotxsec_FluxUnisim";
  }else if (T_weight->GetBranch("piplus_PrimaryHadronSWCentralSplineVariation")){
    option = "piplus_PrimaryHadronSWCentralSplineVariation";
  }else if (T_weight->GetBranch("All_UBGenie")){
    option = "UBGenieFluxSmallUni";
  }else if (T_weight->GetBranch("reinteractions_piminus_Geant4")){
    option = "reinteractions_piminus_Geant4";
  }else if (T_weight->GetBranch("reinteractions_piplus_Geant4")){
    option = "reinteractions_piplus_Geant4";
  }else if (T_weight->GetBranch("reinteractions_proton_Geant4")){
    option = "reinteractions_proton_Geant4";
  }

  //std::cout << "lhagaman debug, option: " << option << "\n";
  
  set_tree_address(T_weight, weight, option);
  //std::cout << "lhagaman debug, T_eval getentries: "  << T_eval->GetEntries() << ", T_weight getentries: " << T_weight->GetEntries() << ", option: " << option << std::endl;

  std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = get_histograms(input_filename,0);

  std::set<std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > >& set_events = map_passed_events[input_filename];

  // hard-coded array of BR1gamma knob values
  std::vector<float> BR1gamma_knob_values = {-0.116749, -1.84206, -1.01624, 0.158259, -0.427051, -1.49894, 0.324356, 0.91424, -1.13741, 1.22972, 
                                              1.03936, -0.0103128, 1.82261, -0.83955, -0.474759, -0.263688, 1.91694, 0.718286, -1.40624, -0.38876, 
                                              -0.838657, -0.355183, 0.730527, 0.32331, -1.46928, -0.439985, -0.126922, 0.249811, -0.113632, 2.12674, 
                                              -0.0426383, -1.71049, -1.52058, -0.380837, -1.04232, -1.83905, 2.00862, 0.411102, 0.992786, 0.140159, 
                                              -1.01282, -1.46354, 0.194565, -1.33395, -1.3236, -1.85725, 0.492589, -0.507482, 0.906432, -0.274199, 
                                              1.65143, 0.0926213, -0.178339, 0.963408, -0.87163, -0.798143, -0.221873, 0.954312, -1.109, 0.692486, 
                                              -2.2429, 1.58851, -0.313319, 0.317757, -0.984987, -0.0646441, 0.589155, 0.883052, 0.42881, -0.540365, 
                                              0.691157, 0.586637, 0.153915, 0.245784, -0.676228, -0.424739, 1.03153, -1.41936, 0.509361, 1.77663, 
                                              0.614415, 0.882829, -0.788378, 1.58072, -0.35986, 0.451795, -0.863749, -0.720094, 0.976311, -0.179523, 
                                              -0.393457, -0.900371, -0.362396, -0.627339, -0.434928, 0.654729, 1.45166, -0.361989, -1.39039, -0.220943, 
                                              -0.514301, -0.103231, 0.684071, 0.0475307, 1.549, -0.902452, 0.037603, 0.397894, -0.5482, 1.33171, 
                                              -0.158242, -0.438244, -1.18242, 0.121432, 0.471135, 0.273454, -1.55219, 0.0335689, -0.40071, -1.22236, 
                                              0.730021, 0.417447, 0.136815, -0.138812, 0.811947, 1.53219, -0.449624, 0.172686, 1.22723, -0.179087, 
                                              0.254897, -0.039421, -0.588572, -0.458291, -0.164697, 1.39117, 0.340808, 0.285425, -0.941114, -0.27606, 
                                              -0.87664, -1.27114, 0.489467, 0.126544, 1.45356, 0.0775609, 0.445432, 0.950793, -0.107524, 1.83575, 
                                              -0.254009, -0.278065, -1.55408, -1.68702, 1.39191, 1.59957, 0.34462, -0.559965, 0.40723, 0.697773, 
                                              1.44662, -1.48115, -0.88834, 0.0447684, -0.269688, 0.963686, 1.24683, 0.320599, 1.45363, -2.24863, 
                                              -0.453747, -0.742489, 2.01879, -1.74642, -0.384873, 0.990525, 1.26116, 1.09646, -0.278814, 0.305654, 
                                              -1.11734, 0.852728, 0.00109231, 1.27902, -1.27839, 2.00059, 0.292331, 0.578787, -0.105597, 0.635094, 
                                              -0.115523, -1.3116, 0.711114, -0.0779606, -2.11692, -0.601426, -0.883515, 2.00071, -1.03561, 0.991168, 
                                              0.482889, -0.0205018, -0.342353, -0.242233, 0.37614, -1.30868, 0.0840374, -0.877406, -0.805776, -0.343845, 
                                              -0.352406, -0.205893, 0.79974, -1.4751, 0.525601, 0.00148406, -0.660796, -1.55735, 0.866562, -1.24584, 
                                              1.22088, 0.0312519, 0.524311, -0.315329, 0.557933, -0.702252, 0.331705, 0.851261, 0.79926, 0.0297268, 
                                              -0.674972, -0.259789, 0.385628, -0.314354, 2.10684, 0.0766031, 1.5018, -0.800019, 0.173723, -0.628923, 
                                              0.0699511, -1.36609, 1.24555, 1.30338, -0.929291, 0.626924, -0.14096, 0.0879994, 0.944861, -0.128607, 
                                              0.25467, 0.451079, 0.482135, -0.395157, -0.305838, 0.607631, -0.00715365, 0.600641, 1.18677, 1.28343, 
                                              0.0847713, -0.0662708, 0.0322692, 1.61935, 2.31695, 0.323979, -0.277203, 0.877043, 0.667412, -1.26465, 
                                              0.763299, -0.468031, -0.944642, -1.62988, -0.771958, -0.945209, -0.235618, 1.18782, -0.217738, -0.31688, 
                                              0.586599, 2.30366, -1.04203, 0.407962, -1.75751, 1.10416, -1.37772, -0.0154544, 0.667943, -0.789175, 
                                              -0.391488, -0.991346, -0.462262, -0.068627, 0.631519, 0.96636, 0.412728, 0.716423, -0.283235, 1.69713, 
                                              0.201634, -0.12603, -0.988541, 1.17857, 0.847212, 1.11766, -1.26282, 0.0348275, 0.583519, 0.679971, 
                                              -1.68946, -0.0346729, -0.67223, 1.09394, 0.0655047, 2.03251, -1.20923, -0.150184, -0.331123, 1.30099, 
                                              0.79657, -0.671498, 0.241091, 0.217776, 0.772156, 0.657497, -0.30889, -0.939185, 0.97638, 0.299696, 
                                              -0.711859, -0.434371, 0.255946, -0.0566167, 0.767525, -0.325056, -1.02786, -0.0802147, -0.963587, 0.586026, 
                                              0.570862, 0.984517, 0.423642, -0.0646823, -0.580852, -0.427061, -2.09155, 0.0939523, -1.32686, -1.30993, 
                                              -0.807024, 0.447053, 0.622578, 0.366325, -2.13501, -0.396309, 0.0402455, 0.539273, -0.164761, -1.64306, 
                                              -2.13555, 0.19562, 0.582938, 0.10473, -0.194025, -0.716588, 0.269256, -0.293556, 0.967054, 1.05231, 
                                              1.45958, 0.732929, 0.00287429, 0.469528, 0.150241, -0.309276, -0.536881, 0.564314, 0.0095683, -0.445847, 
                                              0.40092, -0.574748, 0.367216, -1.61667, 0.80537, 1.52171, -2.30982, -0.38894, 0.490225, -0.85252, 
                                              0.0637627, -0.193764, 0.830956, -0.546995, 1.49644, 0.618402, -0.354665, -0.757796, 0.633325, 0.0884566, 
                                              0.910275, -0.699679, 2.22951, 0.20809, 1.04308, 0.0743343, 3.10452, -0.161184, 0.144011, 0.20844, 
                                              -0.426929, 0.180594, -1.50106, 2.10624, -0.677131, 1.77591, 0.768653, 0.517504, -0.869769, -1.30806, 
                                              0.406596, 0.505189, 0.975807, -2.06136, -0.772879, -0.478636, -0.334357, -0.624104, -1.22598, -0.517347, 
                                              1.81309, -0.975408, -2.05358, -0.590852, 1.67695, 1.0599, 0.720493, 1.13906, 0.892553, 0.233791, 
                                              0.328605, 0.47809, 1.22977, -0.868212, 0.499424, 1.38521, 2.68494, -2.30381, -0.996125, -0.518615, 
                                              0.629957, -3.16407, -1.01749, -0.953748, -1.59199, -0.247813, -0.74494, -0.603797, -0.203816, 0.834678, 
                                              1.37559, -0.465014, -0.184986, 0.00804056, -1.3003, 0.158861, 1.75723, 1.01664, -0.0565817, 0.771311, 
                                              0.207813, 0.665554, 0.847766, -0.767153, 0.0322836, -0.402018, 0.857874, 1.76472, 0.692722, 0.159562, 
                                              -1.69524, 0.370932, 0.194863, 0.947485, -2.33947, 0.0918592, -0.732261, -1.37575, -1.53195, 1.45875, 
                                              -0.294938, -0.316212, 0.423585, 0.509834, 0.320823, 0.18111, 1.42213, 0.804268, -0.832184, -2.48484, 
                                              0.408772, -1.3118, -1.8105, -0.711325, 1.21753, -0.0531389, -1.3266, -1.85725, -2.60498, -2.61219, 
                                              1.05794, -1.23778, 1.15602, 1.93895, -0.811872, 1.34989, 0.0465758, -0.391398, -0.777454, 0.568834, 
                                              0.687976, -0.531665, 0.440267, 0.332448, 0.151151, 0.455383, -0.628806, 0.376317, 0.431937, -0.548751, 
                                              -1.64079, -0.511013, -0.435434, 1.57599, 1.43706, -0.55239, 1.2577, -0.461391, 0.779291, -0.875643, 
                                              -1.37192, -0.42789, 1.18077, 0.739737, 1.14166, -2.35879, 0.553115, -0.759442, 0.788734, -0.840935, 
                                              -0.0638605, 0.437462, -1.15722, 0.0845212, 0.919608, 2.11709, -0.916937, -0.237763, 1.78964, 0.563734, 
                                              -1.97224, 1.03384, 1.02426, 0.819536, 0.0693832, -0.841088, 1.16013, -0.707883, -0.514144, 1.03113, 
                                              -0.793418, 0.0608001, 1.07259, 0.557044, 0.954461, -0.106794, -0.215156, -0.109402, -1.10159, -0.668163, 
                                              -0.0824491, -0.480575, -0.363946, 0.477279, -0.709224, -0.612882, 0.38251, 1.27362, 0.463808, 1.04122, 
                                              2.47561, -2.74118, -0.293176, 1.29899, -0.424589, 1.27081, 0.053118, -0.847153, -0.206133, -1.08917};


  std::vector<int> max_lengths;
  std::vector<int> sup_lengths;
  std::map<TString, int> map_knob_length;
   
  for (size_t i=0;i!=T_eval->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);
    T_weight->GetEntry(i);

    //std::cout << i << std::endl;
    
    std::tuple<float, float, std::vector<float>, std::vector<int>, std::set<std::pair<int, float> > > event_info;
    std::get<0>(event_info) = eval.weight_cv * eval.weight_spline;
    std::get<1>(event_info) = leeweight(eval.truth_nuEnergy);

    double osc_weight = 1.0;
    bool flag_updated = false;

    for (auto it = histo_infos.begin(); it != histo_infos.end(); it++){
      TString histoname = std::get<0>(*it);
      
      //std::cout << "lhagaman debug, histoname = " << histoname << "\n";

      auto it2 = map_histoname_infos.find(histoname);
      int no = std::get<0>(it2->second);
      
      TString var_name = std::get<4>(*it);
      TString ch_name = std::get<5>(*it);
      TString add_cut = std::get<6>(*it);

      // print disabled channel names
      /*std::cout << "Printing weights disabled_ch_names:\n";
      for (const auto& entry : disabled_ch_names) {
          std::cout << entry << "\n";
      }
      std::cout << "Printing Done\n\n";*/

      auto it3 = disabled_ch_names.find(ch_name);
      if (it3 != disabled_ch_names.end()) continue;
      
      float val = get_kine_var(kine, eval, pfeval, tagger, false, var_name);

      //std::cout << "\nlhagaman checking get_cut_pass function\n";

      bool flag_pass = (get_cut_pass(ch_name, add_cut, false, eval, pfeval, tagger, kine) > 0);

      //std::cout << "lhagaman finished checking get_cut_pass function, flag_pass = " << flag_pass << "\n";
      
      if (flag_pass) {
        //std::cout << "lhagaman debug inside flag_pass\n";
        std::get<4>(event_info).insert(std::make_pair(no, val));
        if (flag_osc && is_osc_channel(ch_name) && (!flag_updated)){
          osc_weight = get_osc_weight(eval, pfeval);
          flag_updated = true;
        }
      }
    }
    // apply oscillation ...
    std::get<0>(event_info) *= osc_weight;
    //apply reweight
    double reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info());//automatically 1 if reweighting is not applied
    std::get<0>(event_info) *= reweight; 

    //if (reweight != 1) std::cout << "lhagaman debug, reweight = " << reweight << "\n";
  
    //std::cout << "lhagaman before before ifs\n";
    if (std::get<4>(event_info).size()>0){
      //std::cout << "lhagaman before ifs\n";
      if (option == "expskin_FluxUnisim"){
	std::get<2>(event_info).resize(weight.expskin_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.expskin_FluxUnisim->size());
	for (size_t j=0;j!=weight.expskin_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.expskin_FluxUnisim->at(j) - 1.0; // relative ...
	}
	//std::cout << weight.expskin_FluxUnisim->size() << std::endl;
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.expskin_FluxUnisim->size();
	else{
	  if (weight.expskin_FluxUnisim->size() > it10->second) it10->second =  weight.expskin_FluxUnisim->size();
	}
      }else if (option == "horncurrent_FluxUnisim"){
	std::get<2>(event_info).resize(weight.horncurrent_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.horncurrent_FluxUnisim->size());
	for (size_t j=0;j!= weight.horncurrent_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.horncurrent_FluxUnisim->at(j) - 1.0; // relative ...
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.horncurrent_FluxUnisim->size();
	else{
	  if (weight.horncurrent_FluxUnisim->size() > it10->second) it10->second =  weight.horncurrent_FluxUnisim->size();
	}
      }else if (option == "kminus_PrimaryHadronNormalization"){
	std::get<2>(event_info).resize(weight.kminus_PrimaryHadronNormalization->size());
	std::get<3>(event_info).push_back(weight.kminus_PrimaryHadronNormalization->size());
	for (size_t j=0;j!= weight.kminus_PrimaryHadronNormalization->size(); j++){
	  std::get<2>(event_info).at(j) = weight.kminus_PrimaryHadronNormalization->at(j) - 1.0; 
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.kminus_PrimaryHadronNormalization->size();
	else{
	  if (weight.kminus_PrimaryHadronNormalization->size() > it10->second) it10->second =  weight.kminus_PrimaryHadronNormalization->size();
	}
      }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
	std::get<2>(event_info).resize(weight.kplus_PrimaryHadronFeynmanScaling->size());
	std::get<3>(event_info).push_back(weight.kplus_PrimaryHadronFeynmanScaling->size());
	for (size_t j=0;j!=weight.kplus_PrimaryHadronFeynmanScaling->size();j++){
	  std::get<2>(event_info).at(j) = weight.kplus_PrimaryHadronFeynmanScaling->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.kplus_PrimaryHadronFeynmanScaling->size();
	else{
	  if (weight.kplus_PrimaryHadronFeynmanScaling->size() > it10->second) it10->second =  weight.kplus_PrimaryHadronFeynmanScaling->size();
	}
      }else if (option == "kzero_PrimaryHadronSanfordWang"){
	std::get<2>(event_info).resize(weight.kzero_PrimaryHadronSanfordWang->size());
	std::get<3>(event_info).push_back(weight.kzero_PrimaryHadronSanfordWang->size());
	for (size_t j=0;j!=weight.kzero_PrimaryHadronSanfordWang->size();j++){
	  std::get<2>(event_info).at(j) = weight.kzero_PrimaryHadronSanfordWang->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.kzero_PrimaryHadronSanfordWang->size();
	else{
	  if (weight.kzero_PrimaryHadronSanfordWang->size() > it10->second) it10->second =  weight.kzero_PrimaryHadronSanfordWang->size();
	}
      }else if (option == "nucleoninexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleoninexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleoninexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.nucleoninexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.nucleoninexsec_FluxUnisim->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.nucleoninexsec_FluxUnisim->size();
	else{
	  if (weight.nucleoninexsec_FluxUnisim->size() > it10->second) it10->second =  weight.nucleoninexsec_FluxUnisim->size();
	}
      }else if (option == "nucleonqexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleonqexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleonqexsec_FluxUnisim->size());
	for (size_t j=0; j!= weight.nucleonqexsec_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.nucleonqexsec_FluxUnisim->at(j) - 1.0;
	}
		auto it10 = map_knob_length.find(option);
		if (it10 == map_knob_length.end()) map_knob_length[option] = weight.nucleonqexsec_FluxUnisim->size();
	else{
	  if (weight.nucleonqexsec_FluxUnisim->size() > it10->second) it10->second =  weight.nucleonqexsec_FluxUnisim->size();
	}
      }else if (option == "nucleontotxsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.nucleontotxsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.nucleontotxsec_FluxUnisim->size());
	for (size_t j=0; j!= weight.nucleontotxsec_FluxUnisim->size(); j++){
	  std::get<2>(event_info).at(j) = weight.nucleontotxsec_FluxUnisim->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.nucleontotxsec_FluxUnisim->size();
	else{
	  if (weight.nucleontotxsec_FluxUnisim->size() > it10->second) it10->second =  weight.nucleontotxsec_FluxUnisim->size();
	}
      }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
	std::get<2>(event_info).resize(weight.piminus_PrimaryHadronSWCentralSplineVariation->size());
	std::get<3>(event_info).push_back(weight.piminus_PrimaryHadronSWCentralSplineVariation->size());
	for (size_t j=0; j!= weight.piminus_PrimaryHadronSWCentralSplineVariation->size(); j++){
	  std::get<2>(event_info).at(j) = weight.piminus_PrimaryHadronSWCentralSplineVariation->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.piminus_PrimaryHadronSWCentralSplineVariation->size();
	else{
	  if (weight.piminus_PrimaryHadronSWCentralSplineVariation->size() > it10->second) it10->second =  weight.piminus_PrimaryHadronSWCentralSplineVariation->size();
	}
      }else if (option == "pioninexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.pioninexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.pioninexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.pioninexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.pioninexsec_FluxUnisim->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.pioninexsec_FluxUnisim->size();
	else{
	  if (weight.pioninexsec_FluxUnisim->size() > it10->second) it10->second =  weight.pioninexsec_FluxUnisim->size();
	}
      }else if (option == "pionqexsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.pionqexsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.pionqexsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.pionqexsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.pionqexsec_FluxUnisim->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.pionqexsec_FluxUnisim->size();
	else{
	  if (weight.pionqexsec_FluxUnisim->size() > it10->second) it10->second =  weight.pionqexsec_FluxUnisim->size();
	}
      }else if (option == "piontotxsec_FluxUnisim"){
	std::get<2>(event_info).resize(weight.piontotxsec_FluxUnisim->size());
	std::get<3>(event_info).push_back(weight.piontotxsec_FluxUnisim->size());
	for (size_t j=0;j!=weight.piontotxsec_FluxUnisim->size();j++){
	  std::get<2>(event_info).at(j) = weight.piontotxsec_FluxUnisim->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.piontotxsec_FluxUnisim->size();
	else{
	  if (weight.piontotxsec_FluxUnisim->size() > it10->second) it10->second =  weight.piontotxsec_FluxUnisim->size();
	}
      }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
	std::get<2>(event_info).resize(weight.piplus_PrimaryHadronSWCentralSplineVariation->size());
	std::get<3>(event_info).push_back(weight.piplus_PrimaryHadronSWCentralSplineVariation->size());
	for (size_t j=0; j!= weight.piplus_PrimaryHadronSWCentralSplineVariation->size(); j++){
	  std::get<2>(event_info).at(j) = weight.piplus_PrimaryHadronSWCentralSplineVariation->at(j) - 1.0;
	}
	auto it10 = map_knob_length.find(option);
	if (it10 == map_knob_length.end()) map_knob_length[option] = weight.piplus_PrimaryHadronSWCentralSplineVariation->size();
	else{
	  if (weight.piplus_PrimaryHadronSWCentralSplineVariation->size() > it10->second) it10->second =  weight.piplus_PrimaryHadronSWCentralSplineVariation->size();
	}
      }else if (option == "reinteractions_piminus_Geant4"){
	std::get<2>(event_info).resize(weight.reinteractions_piminus_Geant4->size());
	std::get<3>(event_info).push_back(weight.reinteractions_piminus_Geant4->size());
	for (size_t j=0; j!= weight.reinteractions_piminus_Geant4->size(); j++){
	  std::get<2>(event_info).at(j) = weight.reinteractions_piminus_Geant4->at(j) - 1.0;
	}
      }else if (option == "reinteractions_piplus_Geant4"){
	std::get<2>(event_info).resize(weight.reinteractions_piplus_Geant4->size());
	std::get<3>(event_info).push_back(weight.reinteractions_piplus_Geant4->size());
	for (size_t j=0; j!= weight.reinteractions_piplus_Geant4->size(); j++){
	  std::get<2>(event_info).at(j) = weight.reinteractions_piplus_Geant4->at(j) - 1.0;
	}
      }else if (option == "reinteractions_proton_Geant4"){
	std::get<2>(event_info).resize(weight.reinteractions_proton_Geant4->size());
	std::get<3>(event_info).push_back(weight.reinteractions_proton_Geant4->size());


	if (flag_spec_weights){
	  std::vector<float> temp_vec = get_spec_weight(eval, pfeval);
	  //std::cout << temp_vec.at(0) << std::endl;
	  for (size_t j=0; j!= weight.reinteractions_proton_Geant4->size(); j++){
	    std::get<2>(event_info).at(j) = weight.reinteractions_proton_Geant4->at(j) - 1.0 + temp_vec.at(j);
	  }
	}else{
	  for (size_t j=0; j!= weight.reinteractions_proton_Geant4->size(); j++){
	    std::get<2>(event_info).at(j) = weight.reinteractions_proton_Geant4->at(j) - 1.0;
	  }
        }
	
      }else if (option == "reweight"){
        std::get<2>(event_info).resize(1000);
        std::get<3>(event_info).push_back(1000);
        if(!(flag_reweight)) reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info(true));
        for (size_t j=0;j!=1000;j++){
          if(flag_reweight){
            if (weight.weight_cv>0 && reweight!=1){
              gRandom->SetSeed(j*reweight*77777);
              double rand = gRandom->Gaus(reweight,abs(1-reweight));
              std::get<2>(event_info).at(j) = (rand-reweight)/reweight;
            }else std::get<2>(event_info).at(j) = 0;
          }else{
            gRandom->SetSeed(j*reweight*77777);
            double rand = gRandom->Gaus(1,abs(1-reweight));
            std::get<2>(event_info).at(j) = rand-1;
          }
        }
      }else if (option == "reweight_cor"){
	//std::cout << "lhagaman debug inside larger if\n";
        std::get<2>(event_info).resize(1);
        std::get<3>(event_info).push_back(1);
        //if (reweight != 1) std::cout << "lhagaman debug outside if, reweight = " << reweight << "\n";
	if(flag_reweight){
          //if (reweight != 1) std::cout << "lhagaman debug inside if, reweight = " << reweight << "\n";
          if (weight.weight_cv>0 && reweight!=1){
            std::get<2>(event_info).at(0) = (1-reweight)/reweight;
	   //std::cout << "lhagaman debug, reweight cor value: " << reweight << "\n";
	  }else{
            std::get<2>(event_info).at(0) = 0;
          }
        }else{
            reweight = get_weight("add_weight", eval, pfeval, kine, tagger, get_rw_info(true));
	          std::get<2>(event_info).at(0) = reweight-1;
        }

      }else if (option == "UBGenieFluxSmallUni"){
        int acc_no = 0;
        std::get<2>(event_info).resize(weight.All_UBGenie->size());
        std::get<3>(event_info).push_back(weight.All_UBGenie->size());
        for (size_t j=0;j!=weight.All_UBGenie->size();j++){
          
          // lhagaman modified 2024_08_01, removing BR1gamma knob
          bool remove_BR1gamma_knob = true;

          if (remove_BR1gamma_knob && pfeval.truth_NCDelta) {
            float BR1gamma_knob_value = BR1gamma_knob_values[j];
            float implied_BR1gamma_weight_linear = std::max(0.0, 1 + 0.5*BR1gamma_knob_value);
            float implied_BR1gamma_weight = implied_BR1gamma_weight_linear * implied_BR1gamma_weight_linear;
            float BR1gamma_removed_weight = weight.All_UBGenie->at(j) / implied_BR1gamma_weight;
            if (weight.weight_cv>0){
                std::get<2>(event_info).at(acc_no+j) = (BR1gamma_removed_weight - weight.weight_cv)/weight.weight_cv;
              }else{
                std::get<2>(event_info).at(acc_no+j) = 0;
            }
          } else {
            if (weight.weight_cv>0){
              std::get<2>(event_info).at(acc_no+j) = (weight.All_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
            }else{
              std::get<2>(event_info).at(acc_no+j) = 0;
            }
          }
        }

        acc_no += weight.All_UBGenie->size();

        std::get<2>(event_info).resize(acc_no + weight.AxFFCCQEshape_UBGenie->size());
        std::get<3>(event_info).push_back(weight.AxFFCCQEshape_UBGenie->size());
        for (size_t j=0; j!= weight.AxFFCCQEshape_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.AxFFCCQEshape_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.AxFFCCQEshape_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.DecayAngMEC_UBGenie->size());
        std::get<3>(event_info).push_back(weight.DecayAngMEC_UBGenie->size());
        for (size_t j=0; j!= weight.DecayAngMEC_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.DecayAngMEC_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.DecayAngMEC_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.NormCCCOH_UBGenie->size());
        std::get<3>(event_info).push_back(weight.NormCCCOH_UBGenie->size());
        for (size_t j=0; j!= weight.NormCCCOH_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.NormCCCOH_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.NormCCCOH_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.NormNCCOH_UBGenie->size());
        std::get<3>(event_info).push_back(weight.NormNCCOH_UBGenie->size());
        for (size_t j=0; j!= weight.NormNCCOH_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.NormNCCOH_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.NormNCCOH_UBGenie->size();

        if (!weight.flag_sep_28){
          std::get<2>(event_info).resize(acc_no + weight.RPA_CCQE_Reduced_UBGenie->size());
          std::get<3>(event_info).push_back(weight.RPA_CCQE_Reduced_UBGenie->size());
          for (size_t j=0; j!= weight.RPA_CCQE_Reduced_UBGenie->size(); j++){
            if (weight.weight_cv>0){
              std::get<2>(event_info).at(acc_no+j) = (weight.RPA_CCQE_Reduced_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
            }else{
              std::get<2>(event_info).at(acc_no+j) = 0;
            }
          }
          acc_no += weight.RPA_CCQE_Reduced_UBGenie->size();
        }
        
        std::get<2>(event_info).resize(acc_no + weight.RPA_CCQE_UBGenie->size());
        std::get<3>(event_info).push_back(weight.RPA_CCQE_UBGenie->size());
        for (size_t j=0; j!= weight.RPA_CCQE_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.RPA_CCQE_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.RPA_CCQE_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.ThetaDelta2NRad_UBGenie->size());
        std::get<3>(event_info).push_back(weight.ThetaDelta2NRad_UBGenie->size());
        for (size_t j=0; j!= weight.ThetaDelta2NRad_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.ThetaDelta2NRad_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.ThetaDelta2NRad_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.Theta_Delta2Npi_UBGenie->size());
        std::get<3>(event_info).push_back(weight.Theta_Delta2Npi_UBGenie->size());
        for (size_t j=0; j!= weight.Theta_Delta2Npi_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.Theta_Delta2Npi_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.Theta_Delta2Npi_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.VecFFCCQEshape_UBGenie->size());
        std::get<3>(event_info).push_back(weight.VecFFCCQEshape_UBGenie->size());
        for (size_t j=0; j!= weight.VecFFCCQEshape_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.VecFFCCQEshape_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.VecFFCCQEshape_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.XSecShape_CCMEC_UBGenie->size());
        std::get<3>(event_info).push_back(weight.XSecShape_CCMEC_UBGenie->size());
        for (size_t j=0; j!= weight.XSecShape_CCMEC_UBGenie->size(); j++){
            if (weight.weight_cv>0){
            std::get<2>(event_info).at(acc_no+j) = (weight.XSecShape_CCMEC_UBGenie->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.XSecShape_CCMEC_UBGenie->size();

        
        std::get<2>(event_info).resize(acc_no + weight.xsr_scc_Fa3_SCC->size());
        std::get<3>(event_info).push_back(weight.xsr_scc_Fa3_SCC->size());
        for (size_t j=0; j!= weight.xsr_scc_Fa3_SCC->size(); j++){
            if (weight.weight_cv>0){
              std::get<2>(event_info).at(acc_no+j) = weight.xsr_scc_Fa3_SCC->at(j)-1;//(weight.xsr_scc_Fa3_SCC->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.xsr_scc_Fa3_SCC->size();

        std::get<2>(event_info).resize(acc_no + weight.xsr_scc_Fv3_SCC->size());
        std::get<3>(event_info).push_back(weight.xsr_scc_Fv3_SCC->size());
        for (size_t j=0; j!= weight.xsr_scc_Fv3_SCC->size(); j++){
            if (weight.weight_cv>0){
              std::get<2>(event_info).at(acc_no+j) = weight.xsr_scc_Fv3_SCC->at(j)-1;//(weight.xsr_scc_Fv3_SCC->at(j) - weight.weight_cv)/weight.weight_cv;
          }else{
            std::get<2>(event_info).at(acc_no+j) = 0;
          }
        }
        acc_no += weight.xsr_scc_Fv3_SCC->size();
        
            }

            if (max_lengths.size() < std::get<3>(event_info).size()) max_lengths.resize(std::get<3>(event_info).size());
            for(size_t j = 0; j!= std::get<3>(event_info).size(); j++){
        if (max_lengths.at(j) < std::get<3>(event_info).at(j)) max_lengths.at(j) = std::get<3>(event_info).at(j);
            // if (max_length < std::get<3>(event_info).size()) max_length = std::get<3>(event_info).size();
      }
      
      set_events.insert(event_info);
    }
    
  }

  //std::cout << map_knob_length[option] << std::endl;

  if (option == "expskin_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "horncurrent_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "kminus_PrimaryHadronNormalization"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "kplus_PrimaryHadronFeynmanScaling"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "kzero_PrimaryHadronSanfordWang"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "nucleoninexsec_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "nucleonqexsec_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "nucleontotxsec_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "piminus_PrimaryHadronSWCentralSplineVariation"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "pioninexsec_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "pionqexsec_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "piontotxsec_FluxUnisim"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "piplus_PrimaryHadronSWCentralSplineVariation"){
    if (map_knob_length[option]!=0) sup_lengths.push_back(map_knob_length[option]);
    else sup_lengths.push_back(1000);
  }else if (option == "reinteractions_piminus_Geant4"){
    sup_lengths.push_back(1000);
  }else if (option == "reinteractions_piplus_Geant4"){
    sup_lengths.push_back(1000);
  }else if (option == "reinteractions_proton_Geant4"){
    sup_lengths.push_back(1000);
  }else if (option == "reweight"){
    sup_lengths.push_back(1000);
  }else if (option == "reweight_cor"){
    sup_lengths.push_back(1);
  }else if (option == "UBGenieFluxSmallUni"){
    sup_lengths.push_back(600); // all_ubgenie
    sup_lengths.push_back(1);   // AxFFCCQEshape_UBGenie-
    sup_lengths.push_back(1);   // DecayAngMEC_UBGenie
    sup_lengths.push_back(1); // NormCCCOH_UBGenie
    sup_lengths.push_back(1); // NormNCCOH_UBGenie
    if (!weight.flag_sep_28){
      sup_lengths.push_back(1); // RPA_CCQE_Reduced_UBGenie
    }
    sup_lengths.push_back(2); // RPA_CCQE_UBGenie
    sup_lengths.push_back(1); // ThetaDelta2NRad_UBGenie
    sup_lengths.push_back(1); // Theta_Delta2Npi_UBGenie
    sup_lengths.push_back(1); // VecFFCCQEshape_UBGenie
    sup_lengths.push_back(1); // XSecShape_CCMEC_UBGenie
    sup_lengths.push_back(10); // xsr_scc_Fa3_SCC
    sup_lengths.push_back(10); // xsr_scc_Fv3_SCC

  }



  
  delete file;
  return std::make_pair(max_lengths, sup_lengths);
}


