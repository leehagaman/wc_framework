#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...
#include "TCut.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "tagger.h"
#include "kine.h"
#include "eval.h"
#include "pfeval.h"

#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

namespace LEEana{
  // this is for the real data, for fake data this should be 1 ...
  double em_charge_scale = 0.95;
  //double em_charge_scale = 1.0;

  // correct reco neutrino energy ...
  double get_reco_Enu_corr(KineInfo& kine, bool flag_data);
  double get_reco_showerKE_corr(PFevalInfo& pfeval, bool flag_data);
  
  double get_reco_Eproton(KineInfo& kine);
  double get_reco_Epion(KineInfo& kine);

  double get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data, TString var_name="kine_reco_Enu");
  
  bool get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  double get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval);

  int get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  
  // generic neutrino cuts
  // TCut generic_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >15";
  bool is_generic(EvalInfo& info);
  
  // preselection cuts
  // TCut preselect_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length > 0";
  bool is_preselection(EvalInfo& info); 
  
  // nueCC cuts
  // TCut nueCC_cut = "numu_cc_flag >=0 && nue_score > 7.0";
  bool is_nueCC(TaggerInfo& tagger_info);
  bool is_loosenueCC(TaggerInfo& tagger_info);
  bool is_nueCC_1e0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_nueCC_1e0p0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_nueCC_1eNp0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);

  bool is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_1p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval); 
  bool is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);

  bool is_lowEnergy(KineInfo& kine, bool flag_data);

  bool is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  
  // numuCC cuts
  // TCut numuCC_cut = "numu_cc_flag >=0 && numu_score > 0.9";
  bool is_numuCC(TaggerInfo& tagger_info);
  bool is_loosenumuCC(TaggerInfo& tagger_info);
  bool is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data);
  bool is_numuCC_cutbased(TaggerInfo& tagger_info);
  
  // pio cuts (with and without vertex)
  // TCut pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_pi0(KineInfo& kine, bool flag_data);
  
  // must be with vertex ...
  // TCut cc_pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_cc_pi0(KineInfo& kine, bool flag_data);
  
  

 
  // NC cuts
  // TCut NC_cut = "(!cosmict_flag) && numu_score < 0.0";
  bool is_NC(TaggerInfo& tagger_info);
  bool is_NCpio_bdt(TaggerInfo& tagger_info);
  bool is_NCpio_bdt_sel(TaggerInfo& tagger_info, KineInfo& kine);
  bool is_NCdelta_bdt(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_near_bdt(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_NCdelta_presel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_NCdelta_presel_signal_blind(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  // TCut FC_cut = "match_isFC==1";
  // TCut PC_cut = "match_isFC==0";
  
  bool is_FC(EvalInfo& eval);
  
  
  // TCut truth_nueCC_inside = "abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1";
  // TCut truth_numuCC_inside = "abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1";
  bool is_truth_nueCC_inside(EvalInfo& eval);
  bool is_truth_numuCC_inside(EvalInfo& eval);
}


double LEEana::get_reco_Enu_corr(KineInfo& kine, bool flag_data){
  double reco_Enu_corr = 0;
  if (kine.kine_reco_Enu > 0){
    if (flag_data){
      for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
  	if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
  	  reco_Enu_corr +=  kine.kine_energy_particle->at(j) * em_charge_scale;
  	}else{
  	  reco_Enu_corr +=  kine.kine_energy_particle->at(j);
  	}
  	//	std::cout << "p: " << kine.kine_energy_particle->at(j) << " " << kine.kine_energy_info->at(j) << " " << kine.kine_particle_type->at(j) << " " << kine.kine_energy_included->at(j) << std::endl;
      }
      reco_Enu_corr += kine.kine_reco_add_energy;
      return reco_Enu_corr;
    }
  }
  return kine.kine_reco_Enu;
}

// added by lhagaman 2021_07_23
double LEEana::get_reco_showerKE_corr(PFevalInfo& pfeval, bool flag_data){
    if (flag_data){
	return pfeval.reco_showerKE * em_charge_scale;
    } else {
	return pfeval.reco_showerKE;
    }
}


double LEEana::get_reco_Eproton(KineInfo& kine){
  double reco_Eproton=0;
  for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
    {
      int pdgcode = kine.kine_particle_type->at(i);
      if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) // proton threshold of 35 MeV
      //if(abs(pdgcode)==2212) // no proton threshold
	reco_Eproton+=kine.kine_energy_particle->at(i);

    }
  return reco_Eproton;
}

double LEEana::get_reco_Epion(KineInfo& kine){
  double reco_Epion=0;
  for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
    {
      int pdgcode = kine.kine_particle_type->at(i);
      //if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10)  // KE threshold: 10 keV      
      if(abs(pdgcode)==211)  // no threshold
	reco_Epion+=kine.kine_energy_particle->at(i);
    }
  return reco_Epion;
}

double LEEana::get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval){
 
    double addtl_weight = 1.0;
    int use_reweight = 0;//0: no reweight, 1:reweight without accouting for constraint, 2: reweight accounting for numuCC constraint
                     //modify the variables inside the if() statment if you wish to change the reweighting
  //Calculates reweighting for the event
  if(use_reweight !=0){
    //initiates variables needed for reweighting
    double var;
    double min_var_Np;
    double max_var_Np;
    double min_var_0p;
    double max_var_0p;
    int nbins_Np;
    int nbins_0p;
    bool use_overflow;
    std::vector<double> reweight;
    std::vector<double> overflow;
  
    //sets up your desired reweighting
    if(use_reweight==1){
      var = eval.truth_energyInside;
      min_var_Np = 200;
      max_var_Np = 1100;
      min_var_0p = 100;
      max_var_0p = 800;
      nbins_Np = 15;
      nbins_0p = 15;
      use_overflow = 1;
      reweight = {0.47346, 0.513817, 0.620152, 0.711642, 0.683565, 0.599427, 0.487592, 0.368639, 0.306867, 0.316878, 0.3627, 0.382381, 0.370914, 0.339537, 0.327263,
                  2.18203, 1.99804, 1.8456, 1.85226, 1.88324, 1.8381, 1.62346, 1.321, 1.08835, 1.0664, 1.26347, 1.62921, 2.07742, 2.68451, 3.2522};
      overflow = {0.246093, 2.67171};
    } 
    if(use_reweight==2){
      var = eval.truth_energyInside;
      min_var_Np = 200;
      max_var_Np = 1400;
      min_var_0p = 100;
      max_var_0p = 800;
      nbins_Np = 24;
      nbins_0p = 14;
      use_overflow = 0;
      reweight = {0.279332, 0.388117, 0.553707, 0.679725, 0.688555, 0.617495, 0.567947, 0.560433, 0.59077, 0.617946, 0.584057, 0.519013, 0.471467, 0.458097, 0.553915, 0.735836, 0.990822, 1.26522, 1.57029, 1.73207, 1.9213, 1.99958, 2.05557, 2.18402,
                  1.2466, 1.18207, 1.16923, 1.37165, 1.51289, 1.42552, 1.11135, 0.76296, 0.574159, 0.617545, 0.864722, 1.15229, 1.64071, 2.08024};//, 2.30331, 2.67676, 3.12894, 3.86247, 4.37093, 5.51674, 7.10849, 8.05162, 7.80566, 9.49697};
      overflow = {1.0, 1.0};
    }
  
    //Finds the number of true protons for the given event
    int nTrueP = 0;
    double Eproton = 0;
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
      nTrueP+=1;
      if(pfeval.truth_startMomentum[i][3] - 0.938272 < Eproton) continue;
      Eproton = (pfeval.truth_startMomentum[i][3] - 0.938272)*1000;
    }

    //determines the reweighting bin the event falls into to assign the reweighing for the event 
    double bin_len_Np = (max_var_Np-min_var_Np)/nbins_Np;
    double bin_len_0p = (max_var_0p-min_var_0p)/nbins_0p;
    if (eval.truth_isCC==0 && pfeval.truth_NprimPio>0 && eval.match_isFC){
      if(nTrueP>0 && var>min_var_Np){
        if (var>max_var_Np && use_overflow) addtl_weight = overflow[0];
        else if(var>max_var_Np) addtl_weight = 1;
        else{
          int wbin = floor((var-min_var_Np)/bin_len_Np);
          addtl_weight = reweight[wbin];
        }
      }
      else if(var>min_var_0p){
        if (var>max_var_0p && use_overflow) addtl_weight = overflow[1];
        else if(var>max_var_0p) addtl_weight = 1;
        else{
          int wbin = floor((var-min_var_0p)/bin_len_0p)+nbins_Np;
          addtl_weight = reweight[wbin];
        }
      }
    }
  }    
	
  if (weight_name == "cv_spline"){
    //if (addtl_weight * eval.weight_cv * eval.weight_spline < 0) return 0.;	  
    return addtl_weight * eval.weight_cv * eval.weight_spline;
  }else if (weight_name == "cv_spline_cv_spline"){
    return pow(eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "unity" || weight_name == "unity_unity"){
    return 1;
  }else if (weight_name == "lee_cv_spline"){
    return (eval.weight_lee * eval.weight_cv * eval.weight_spline);
  }else if (weight_name == "lee_cv_spline_lee_cv_spline"){
    return pow(eval.weight_lee * eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "lee_cv_spline_cv_spline" || weight_name == "cv_spline_lee_cv_spline"){
    return eval.weight_lee * pow(eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "spline"){
    return eval.weight_spline;
  }else if (weight_name == "spline_spline"){
    return pow(eval.weight_spline,2);
  }else if (weight_name == "lee_spline"){
    return (eval.weight_lee * eval.weight_spline);
  }else if (weight_name == "lee_spline_lee_spline"){
    return pow(eval.weight_lee * eval.weight_spline,2);
  }else if (weight_name == "lee_spline_spline" || weight_name == "spline_lee_spline"){
    return eval.weight_lee * pow( eval.weight_spline,2);
  }else{
    std::cout <<"Unknown weights: " << weight_name << std::endl;
  }
	    
  
  return 1;
}

double LEEana::get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data , TString var_name){
  //  if (var_name == "kine_reco_Enu"){
  //  return kine.kine_reco_Enu;
  //  }else
  if (var_name == "kine_reco_Enu"){
    return get_reco_Enu_corr(kine, flag_data);
  }
  else if (var_name == "reco_showerKE"){
    return get_reco_showerKE_corr(pfeval, flag_data) * 1000.;
  }
  else if (var_name == "kine_reco_add_energy"){
    return kine.kine_reco_add_energy;
  }
  else if (var_name == "kine_reco_Eproton"){
    return get_reco_Eproton(kine);
  }
  else if (var_name == "kine_reco_Epion"){
    return get_reco_Epion(kine);
  }
  else if (var_name == "kine_pio_energy_1"){
    return kine.kine_pio_energy_1;
  }
  else if (var_name == "kine_pio_energy_2"){
    return kine.kine_pio_energy_2;
  }
  else if (var_name == "kine_pio_energy_max"){
    if(flag_data)
      return std::max(kine.kine_pio_energy_1*em_charge_scale, kine.kine_pio_energy_2*em_charge_scale);
    else
      return std::max(kine.kine_pio_energy_1, kine.kine_pio_energy_2);
  }
  else if (var_name == "kine_pio_energy_min"){
    if(flag_data)
      return std::min(kine.kine_pio_energy_1*em_charge_scale, kine.kine_pio_energy_2*em_charge_scale);
    else
      return std::min(kine.kine_pio_energy_1, kine.kine_pio_energy_2);
  }
  else if (var_name == "kine_pio_angle" || var_name == "kine_pio_costheta"){
    if (var_name == "kine_pio_angle")
      return kine.kine_pio_angle;
    else
      return TMath::Cos(kine.kine_pio_angle/180.*TMath::Pi());
  }
  else if (var_name == "thetaXZ"){
    if (pfeval.reco_muonMomentum[0]!=0 && pfeval.reco_muonMomentum[1]!=0 && pfeval.reco_muonMomentum[2]!=0)
      return (pfeval.reco_muonMomentum[0]/(sqrt(pfeval.reco_muonMomentum[0]*pfeval.reco_muonMomentum[0] + pfeval.reco_muonMomentum[2]*pfeval.reco_muonMomentum[2])));
  }
  else if (var_name == "thetaYZ"){
    if (pfeval.reco_muonMomentum[0]!=0 && pfeval.reco_muonMomentum[1]!=0 && pfeval.reco_muonMomentum[2]!=0)
      return (pfeval.reco_muonMomentum[1]/(sqrt(pfeval.reco_muonMomentum[1]*pfeval.reco_muonMomentum[1] + pfeval.reco_muonMomentum[2]*pfeval.reco_muonMomentum[2])));
  }
  else if (var_name == "match_energy"){
    return eval.match_energy;
  }else if (var_name == "pi0_energy"){
    double pi0_mass = 135;
    double alpha = fabs(kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
    return pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
  }else if (var_name == "pi0_mass"){
    if (kine.kine_pio_mass >0){
      //      TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      // TLorentzVector pio = p1 + p2;
      if (flag_data) {
	return kine.kine_pio_mass * em_charge_scale;
      }else{
	return kine.kine_pio_mass;
      }
    }else{
      return kine.kine_pio_mass;
    }
    //  }else if (var_name == "pi0_mass"){
    // return kine.kine_pio_mass;
  }else if (var_name == "nue_score"){
    return tagger.nue_score;
  }else if (var_name == "numu_score"){
    return tagger.numu_score;
  }else if (var_name == "nc_pio_score"){
    return tagger.nc_pio_score;
  }else if (var_name == "nc_delta_score"){
    return tagger.nc_delta_score;
  }else if (var_name == "shower_energy"){
    if(flag_data)
      return tagger.mip_energy*em_charge_scale;
    else
      return tagger.mip_energy;
  }else if (var_name == "shower_angle_beam"){
    return tagger.mip_angle_beam;
  }else if (var_name == "shower_angle_vertical"){
    return tagger.spt_angle_vertical;
  }else if (var_name == "shwvtx_nuvtx_dis"){
    return sqrt(pow(pfeval.reco_nuvtxX-pfeval.reco_showervtxX,2)+pow(pfeval.reco_nuvtxY-pfeval.reco_showervtxY,2)+pow(pfeval.reco_nuvtxZ-pfeval.reco_showervtxZ,2)); 
  }else if (var_name == "median_dQdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    return vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
  }else if (var_name == "median_dEdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    float median_dqdx = vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
    float alpha = 1.;
    float beta = 0.255;
    float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
    if(median_dedx<0) median_dedx = 0;
    if(median_dedx>50) median_dedx = 50;
    return median_dedx; // MeV/cm
  }else if (var_name == "reco_showervtxX"){
      return pfeval.reco_showervtxX;
  }else if (var_name == "reco_nuvtxX"){
      return pfeval.reco_nuvtxX;
  }else if (var_name == "reco_nuvtxY"){
      return pfeval.reco_nuvtxY;
  }else if (var_name == "reco_nuvtxZ"){
      return pfeval.reco_nuvtxZ;
  }else if (var_name == "reco_nuvtxU"){
    return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) - pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);
  }else if (var_name == "reco_nuvtxV"){
    return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) + pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);
  }else if (var_name == "mip_quality_n_tracks"){
      return tagger.mip_quality_n_tracks;
  }else if (var_name == "mip_quality_n_showers"){
      return tagger.mip_quality_n_showers;
  }else if (var_name == "gap_n_bad"){
      return tagger.gap_n_bad;
  }else if (var_name == "muon_KE"){
      return pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
  }else if (var_name == "reco_Emuon"){
      return pfeval.reco_muonMomentum[3]*1000; // GeV --> MeV  
  }else if (var_name == "muon_costheta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return TMath::Cos(muonMomentum.Theta());
      else
	return -2;
  }else if (var_name == "muon_theta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	//return muonMomentum.Theta(); // in radians
	return muonMomentum.Theta()*180./TMath::Pi(); // in degrees
      else
	return -1000;
  }else if (var_name == "muon_phi"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return muonMomentum.Phi()/TMath::Pi()*180.;
      else
	return -1000;
  }
  else if (var_name == "reco_Eqe_muon" || var_name == "reco_Eqe_muon_Enu_diff" || var_name == "reco_Eqe_electron" || var_name == "reco_Eqe_electron_Enu_diff"){
    // everything is in MeV
    float neutron_mass = 939.57;
    float binding_energy = 30.0;
    float muon_mass = 105.66;
    float electron_mass = 0.511;
    float proton_mass = 938.27;

    float muon_costheta;
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    if (pfeval.reco_muonMomentum[3]>0)
      muon_costheta = TMath::Cos(muonMomentum.Theta());

    float shower_costheta;
    TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);
    if (pfeval.reco_showerMomentum[3]>0)
      shower_costheta = TMath::Cos(showerMomentum.Theta());
   
    // added lhagaman 2022_04_28
    shower_costheta = TMath::Cos(tagger.mip_angle_beam/180.*TMath::Pi()); 
 
    float reco_Emuon =  pfeval.reco_muonMomentum[3]*1000.; // GeV --> MeV                                                                                                            
    float reco_Eqe_muon = 0.5 * (2*(neutron_mass-binding_energy)*reco_Emuon - (pow(neutron_mass-binding_energy,2) + pow(muon_mass,2) - pow(proton_mass,2))) / ((neutron_mass-binding_energy) - reco_Emuon + sqrt(pow(reco_Emuon,2)-pow(muon_mass,2))*muon_costheta);                              

    float reco_Eelectron= pfeval.reco_showerMomentum[3]*1000.;

    // added by lhagaman 2021_07_23, temporary, open data files don't have shower momentum variable
    // shower_costheta = TMath::Cos(tagger.mip_angle_beam/180.*TMath::Pi()); removed 2022_04_20, using open data files with showerMomentum filled
    // needs this for detsys I think
    reco_Eelectron = get_reco_showerKE_corr(pfeval, flag_data) * 1000.;


    float reco_Eqe_electron = 0.5 * (2*(neutron_mass-binding_energy)*reco_Eelectron - (pow(neutron_mass-binding_energy,2) + pow(electron_mass,2) - pow(proton_mass,2))) / ((neutron_mass-binding_energy) - reco_Eelectron + sqrt(pow(reco_Eelectron,2)-pow(electron_mass,2))*shower_costheta);

    // added lhagaman 2022_04_28
    if (shower_costheta > 0.999 || shower_costheta < -0.999) {
        reco_Eqe_electron = -1;
    }	

    if(var_name=="reco_Eqe_muon") return reco_Eqe_muon;
    else if(var_name=="reco_Eqe_electron") return reco_Eqe_electron;
    else if(var_name=="reco_Eqe_muon_Enu_diff") return reco_Eqe_muon - get_reco_Enu_corr(kine, flag_data);
    else if(var_name=="reco_Eqe_electron_Enu_diff") return reco_Eqe_electron - get_reco_Enu_corr(kine, flag_data);
  }
  else if (var_name == "proton_KE"){
      return pfeval.reco_protonMomentum[3]*1000.-938.27; // GeV--> MeV
  }else if (var_name == "proton_theta" || var_name == "proton_phi" || var_name == "proton_costheta"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      if(pfeval.reco_protonMomentum[3]>0){
	if(var_name == "proton_theta")
	  return protonMomentum.Theta()/TMath::Pi()*180.;
	else if (var_name == "proton_costheta")
	  return TMath::Cos(protonMomentum.Theta());
	else if (var_name == "proton_phi")
	  return protonMomentum.Phi()/TMath::Pi()*180.;
      }
      else
	return -1000;
  }
  else if (var_name == "shower_theta" || var_name == "shower_costheta" || var_name == "shower_costheta" || var_name == "shower_phi"){
    TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);

    if(var_name == "shower_theta")
      return tagger.mip_angle_beam;
      //return showerMomentum.Theta()/TMath::Pi()*180.;

    if(var_name == "shower_costheta"){
       
      // changed lhagaman 2022_04_20, then changed back, it seems like detvar files don't have showerMomentum saved
      return TMath::Cos(tagger.mip_angle_beam/180.*TMath::Pi());
      
      //float shower_costheta = -2.;
      //if (pfeval.reco_showerMomentum[3]>0)
      // shower_costheta = TMath::Cos(showerMomentum.Theta());

      //return shower_costheta;
    }

    if(var_name == "shower_phi"){
      if (pfeval.reco_showerMomentum[3]>0)
	return showerMomentum.Phi()/TMath::Pi()*180.;
      else
        return -1000;
    }
  }
  else if (var_name=="shower_proton_angle_sum"){
    TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
    TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);

    if(pfeval.reco_showerMomentum[3]>0 && pfeval.reco_protonMomentum[3]>0)
      return showerMomentum.Theta()/TMath::Pi()*180. + protonMomentum.Theta()/TMath::Pi()*180.;
    else 
      return -1000;
  }
  else if (var_name=="muon_proton_angle_sum"){
    TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);

    if(pfeval.reco_muonMomentum[3]>0 && pfeval.reco_protonMomentum[3]>0)
      return muonMomentum.Theta()/TMath::Pi()*180. + protonMomentum.Theta()/TMath::Pi()*180.;
    else 
      return -1000;
  }
  else if (var_name == "Ehadron"){
    /*
      // for numuCC
    if (pfeval.reco_muonMomentum[3]>0)
      return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
    else
      return -1000;
    */

    // for nueCC
    //if (pfeval.reco_showerMomentum[3]>0)
    if(flag_data)
      return get_reco_Enu_corr(kine, flag_data) - tagger.mip_energy*em_charge_scale;
    else
      return get_reco_Enu_corr(kine, flag_data) - tagger.mip_energy;
      //return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_showerMomentum[3]*1000.;
      //else
      //return -1000;

    //  }else if (var_name == "Ehadron"){
      /* Float_t Ehadron = kine.kine_reco_Enu; */
      /* for(size_t i=0; i<kine.kine_energy_particle->size(); i++) */
      /* { */
      /*     int pdgcode = kine.kine_particle_type->at(i); */
      /*     if(abs(pdgcode)==13) Ehadron = Ehadron - kine.kine_energy_particle->at(i) - 105.658; */ 
      /*     //if(abs(pdgcode)==11) Ehadron = Ehadron - kine.kine_energy_particle->at(i); */ 
      /* } */
    // return kine.kine_reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
  }else if (var_name == "Q2"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
    //  }else if (var_name == "Q2"){
    // Float_t Enu = kine.kine_reco_Enu;
    //Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    //Float_t Ehadron = Enu - Emu;
    //Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    //TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    //Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    //return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
  }else if (var_name == "x_Bjorken"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
    //  }else if (var_name == "x_Bjorken"){
    // Float_t Enu = kine.kine_reco_Enu;
    // Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    // Float_t Ehadron = Enu - Emu;
    // Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    // TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    // Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    // return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
  }else if (var_name == "N_tracks"){
      int N_tracks = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==11) continue;
          if(kine.kine_energy_particle->at(i)<10) continue;
          if(abs(pdgcode)==13 || abs(pdgcode)==211){
            N_tracks += 1;
          }
          else if(kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_tracks += 1; 
          }
      }
      return N_tracks;
  }else if (var_name == "N_other_tracks"){
	int Nothertracks = 0;
        for(size_t i=0; i<kine.kine_energy_particle->size(); i++){
                int pdgcode = kine.kine_particle_type->at(i);
                if((abs(pdgcode)==211 || abs(pdgcode)==13) && kine.kine_energy_particle->at(i)>10) Nothertracks++; // KE threshold: 10 MeV
        }
	return Nothertracks;
  }else if (var_name == "N_showers"){
      int N_showers = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)!=11) continue;
          if(kine.kine_energy_particle->at(i)>10) N_showers += 1;
      }
      return N_showers;
  }else if (var_name == "N_protons"){
      int N_protons = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)== 2212 && kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_protons += 1; 
          }
      }
      return N_protons;
  }else if (var_name == "proton_pi0_total_momentum" || var_name == "proton_pi0_invariant_mass") {
	bool debug_pf_info = 0;
	TLorentzVector max_energy_proton_momentum(-1., -1., -1., -1.);
        TLorentzVector gamma_1_momentum(-1., -1., -1., -1.);
        TLorentzVector gamma_2_momentum(-1., -1., -1., -1.);
	float max_proton_energy = 0.;
	if (debug_pf_info) {
		std::cout << "********************************* starting event  ************************\n";
		std::cout << "(100 hard coded) looping over " << pfeval.reco_Ntrack << " reco particles (" << pfeval.truth_Ntrack << " true particles)\n";
	}
	for(size_t i=0; i<100; i++) // looping over all reconstructed particles
        {
        	int pdgcode = pfeval.reco_pdg[i];
		if (debug_pf_info) std::cout << "investigating particle" << i << "//" << pfeval.reco_Ntrack << " : " << pdgcode << "\n";
          	if(abs(pdgcode)==2212 && pfeval.reco_startMomentum[i][3] > max_proton_energy){ // new max energy proton
			if (debug_pf_info) std::cout << "new max energy proton\n";
              		max_energy_proton_momentum = pfeval.reco_startMomentum[i];
			max_proton_energy = pfeval.reco_startMomentum[i][3];
          	}
		if(abs(pdgcode)==22 || abs(pdgcode)==11) { // reconstructed shower
			float shower_energy = 1000. * pfeval.reco_startMomentum[i][3];
			if (abs(shower_energy - kine.kine_pio_energy_1) / kine.kine_pio_energy_1 < 0.01) { // very close to gamma 1 energy	
				if (debug_pf_info) std::cout << "gamma 1 matched";
				if (flag_data) {
					gamma_1_momentum = TLorentzVector(em_charge_scale * pfeval.reco_startMomentum[i][0],
							em_charge_scale * pfeval.reco_startMomentum[i][1],
							em_charge_scale * pfeval.reco_startMomentum[i][2],
							em_charge_scale * pfeval.reco_startMomentum[i][3]);
				} else {
					gamma_1_momentum = pfeval.reco_startMomentum[i];
				}
			}
			if (abs(shower_energy - kine.kine_pio_energy_2) / kine.kine_pio_energy_2 < 0.01) { // very close to gamma 2 energy
				if (debug_pf_info) std::cout << "gamma 2 matched";
				if (flag_data) {
                                        gamma_2_momentum = TLorentzVector(em_charge_scale * pfeval.reco_startMomentum[i][0],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][1],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][2],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][3]);
                                } else {
                                        gamma_2_momentum = pfeval.reco_startMomentum[i];
                                }
                        }

		}
        }

	float proton_pi0_invariant_mass = -1.;
	float proton_pi0_total_momentum = -1.;
	
	float pi0_mass = 134.9768;
	float proton_mass = 938.272;

	if (max_energy_proton_momentum[3]>0 && gamma_1_momentum[3]>0 && gamma_2_momentum[3]>0) { // found proton and pi0 in PF tree
		if (debug_pf_info) std::cout << "matched proton and pi0!" << std::endl;
		TLorentzVector pi0_momentum = gamma_1_momentum + gamma_2_momentum;
		proton_pi0_invariant_mass = sqrt(pi0_mass * pi0_mass + proton_mass * proton_mass
				+ 2. * 1000. * 1000. * (max_energy_proton_momentum[3] * pi0_momentum[3]
						      - max_energy_proton_momentum[0] * pi0_momentum[0]
						      - max_energy_proton_momentum[1] * pi0_momentum[1]
						      - max_energy_proton_momentum[2] * pi0_momentum[2])); 
		proton_pi0_total_momentum = 1000. * sqrt((max_energy_proton_momentum[0] + pi0_momentum[0]) * (max_energy_proton_momentum[0] + pi0_momentum[0])
						       + (max_energy_proton_momentum[1] + pi0_momentum[1]) * (max_energy_proton_momentum[1] + pi0_momentum[1])
						       + (max_energy_proton_momentum[2] + pi0_momentum[2]) * (max_energy_proton_momentum[2] + pi0_momentum[2]));
	}

	if (debug_pf_info) std::cout << "******************* ending event *************\n";

	if (var_name == "proton_pi0_total_momentum") {
		return proton_pi0_total_momentum;
	} else if (var_name == "proton_pi0_invariant_mass") {
		return proton_pi0_invariant_mass;
	} else {
		std::cout << "No such proton-pi0 variable: " << var_name << std::endl;
	}
  }else{
    std::cout << "No such variable: " << var_name << std::endl;
    exit(EXIT_FAILURE);
  }
  return -1;
}

int LEEana::get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){
  for (auto it = map_cut_xs_bin.begin(); it != map_cut_xs_bin.end(); it++){
    TString cut_name = it->first;
    int number = it->second;

    double Emuon = pfeval.truth_muonMomentum[3]*1000; // MeV
    double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

    if (cut_file == 1){
      if (cut_name == "numuCC.inside.Enu.le.300"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400.gt.300"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 && eval.truth_nuEnergy>300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500.gt.400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 && eval.truth_nuEnergy>400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.600.gt.500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=600 && eval.truth_nuEnergy>500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.700.gt.600"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=700 && eval.truth_nuEnergy>600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.800.gt.700"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=800 && eval.truth_nuEnergy>700) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.900.gt.800"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=900 && eval.truth_nuEnergy>800) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1000.gt.900"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1000 && eval.truth_nuEnergy>900) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1100.gt.1000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1100 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1100"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1100) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1500.gt.1200"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1500 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2100.gt.1500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2100 && eval.truth_nuEnergy>1500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1400.gt.1200"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1400 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1600.gt.1400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1600 && eval.truth_nuEnergy>1400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2000.gt.1600"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2000 && eval.truth_nuEnergy>1600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2500.gt.2000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2500 && eval.truth_nuEnergy>2000) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2500) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2100"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2100) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.1500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>1500) return number;
      }else{
	         std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 2) {
      if (cut_name == "numuCC.inside.Emuon.le.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=100 && Emuon>0) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.200.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=200 && Emuon>100) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.300.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=300 && Emuon>200) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.400.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=400 && Emuon>300) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.500.gt.400"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=500 && Emuon>400) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.600.gt.500"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=600 && Emuon>500) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.700.gt.600"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=700 && Emuon>600) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.800.gt.700"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=800 && Emuon>700) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.900.gt.800"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=900 && Emuon>800) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1000.gt.900"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1000 && Emuon>900) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1200.gt.1000"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1200 && Emuon>1000) return number;
      }
      else if (cut_name == "numuCC.inside.Emuon.gt.1200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1200) return number;
      }else{
          std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 3) {
      if (cut_name == "numuCC.inside.Ehadron.le.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.200.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=200 && Ehadron>100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.300.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=300 && Ehadron>200) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.400.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=400 && Ehadron>300) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.500.gt.400"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=500 && Ehadron>400) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.600.gt.500"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=600 && Ehadron>500) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.700.gt.600"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=700 && Ehadron>600) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.800.gt.700"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=800 && Ehadron>700) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.900.gt.800"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=900 && Ehadron>800) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1000.gt.900"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1000 && Ehadron>900) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.gt.1000"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1000) return number;
      }else{
          std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }       
    }
    else if (cut_file == 4){
      if (cut_name == "numuCC.inside.Enu.le.540.gt.200"){ // recommended range: 200 - 540
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.705.gt.540"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.805.gt.705"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=805 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.920.gt.805"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>805) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1050.gt.920"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1050 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1050"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1050) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1375.gt.1200"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1375 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1570.gt.1375"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1375) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2050.gt.1570"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.4000.gt.2050"){ // recommended range: 2050 - 4000
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 5) {
      if (cut_name == "numuCC.inside.Emuon.le.226.gt.106"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.296.gt.226"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.386.gt.296"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.505.gt.386"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.577.gt.505"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.659.gt.577"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.753.gt.659"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.861.gt.753"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.984.gt.861"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1285.gt.984"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984) return number;
      }
      else if (cut_name == "numuCC.inside.Emuon.le.2506.gt.1285"){ // 1285 - 2506, only 1% > 2506 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 6) {
      if (cut_name == "numuCC.inside.Ehadron.le.100.gt.30"){ // 30 MeV - 100 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=100 && Ehadron>30) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.150.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=150 && Ehadron>100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.225.gt.150"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=225 && Ehadron>150) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.275.gt.225"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=275 && Ehadron>225) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.336.gt.275"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=336 && Ehadron>275) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.411.gt.336"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=411 && Ehadron>336) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.502.gt.411"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=502 && Ehadron>411) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.614.gt.502"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=614 && Ehadron>502) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.750.gt.614"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=750 && Ehadron>614) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1120.gt.750"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1120 && Ehadron>750) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.2500.gt.1120"){ // 1120 - 2500 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1120 && Ehadron<= 2500) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }

  }
  
  return -1;
}

bool LEEana::get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){


  float reco_Enu = get_reco_Enu_corr(kine, flag_data);
  
  double Emuon = pfeval.truth_muonMomentum[3]*1000; // MeV
  double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV
  
  bool flag_truth_inside = false; // in the active volume
  if (eval.truth_vtxX > -1 && eval.truth_vtxX <= 254.3 &&  eval.truth_vtxY >-115.0 && eval.truth_vtxY<=117.0 && eval.truth_vtxZ > 0.6 && eval.truth_vtxZ <=1036.4) flag_truth_inside = true;

  // definition of additional cuts
  std::map<std::string, bool> map_cuts_flag;
  if(is_far_sideband(kine, tagger, flag_data)) map_cuts_flag["farsideband"] = true; 
  else map_cuts_flag["farsideband"] = false; 
  
  if(is_near_sideband(kine, tagger, flag_data)) map_cuts_flag["nearsideband"] = true; 
  else map_cuts_flag["nearsideband"] = false; 
 
  if(is_nueCC(tagger)) map_cuts_flag["nueCC"] = true;
  else map_cuts_flag["nueCC"] = false;

  if(is_loosenueCC(tagger)) map_cuts_flag["loosenueCC"] = true;
  else map_cuts_flag["loosenueCC"] = false;

  if(is_loosenumuCC(tagger)) map_cuts_flag["loosenumuCC"] = true;
  else map_cuts_flag["loosenumuCC"] = false;
  
  if(is_generic(eval)) map_cuts_flag["generic"] = true;
  else map_cuts_flag["generic"] = false;

  if(eval.truth_nuEnergy <=400) map_cuts_flag["LowEintnueCC"] = true;
  else map_cuts_flag["LowEintnueCC"] = false;
  
  if (!(eval.truth_nuEnergy <=400)) map_cuts_flag["antiLowEintnueCC"] = true;
  else map_cuts_flag["antiLowEintnueCC"] = false;

  if(eval.truth_nuEnergy<=400) map_cuts_flag["LowEnu"] = true;
  else map_cuts_flag["LowEnu"] = false;
  
  if(!(eval.truth_nuEnergy<=400)) map_cuts_flag["antiLowEnu"] = true;
  else map_cuts_flag["antiLowEnu"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["badmatch"] = true;
  else map_cuts_flag["badmatch"] = false;
  
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["numuCCinFV"] = true;
  else map_cuts_flag["numuCCinFV"] = false;

  // Xs related cuts ...

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 ) map_cuts_flag["XsnumuCCinFV"] = true;
  else map_cuts_flag["XsnumuCCinFV"] = false;
  
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200) map_cuts_flag["Xs_Enu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Enu_numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon > 105.7 && Emuon<=2506) map_cuts_flag["Xs_Emu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Emu_numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron > 30 && Ehadron <=2500) map_cuts_flag["Xs_Ehad_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Ehad_numuCCinFV"] = false;

  // xs breakdown mode
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCinFV"] = true;
  else map_cuts_flag["XsecNumuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0) map_cuts_flag["XsecNC"] = true;
  else map_cuts_flag["XsecNC"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["XsecCosmic"] = true;
  else map_cuts_flag["XsecCosmic"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==1 && !(eval.truth_nuPdg==14 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200)) map_cuts_flag["XsecBkgCC"] = true;
  else map_cuts_flag["XsecBkgCC"] = false;

  // finish Xs related cuts ...
  
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["nueCCinFV"] = true;
  else map_cuts_flag["nueCCinFV"] = false;
    
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["NCinFV"] = true;
  else map_cuts_flag["NCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==0) map_cuts_flag["outFV"] = true;
  else map_cuts_flag["outFV"] = false;
      
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0) map_cuts_flag["CCpi0inFV"] = true;
  else map_cuts_flag["CCpi0inFV"] = false;
      
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio>0) map_cuts_flag["NCpi0inFV"] = true;
  else map_cuts_flag["NCpi0inFV"] = false;


  // lee hagaman adding cuts for nc delta breakdown, 2020_07_13
  // nueCCinFV, badmatch, outFV, ext, and dirt categories don't need any changing
  
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NCDelta==1) map_cuts_flag["NCDeltainFV"] = true;
  else map_cuts_flag["NCDeltainFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && pfeval.truth_NCDelta==0) map_cuts_flag["NC1Pi0inFV"] = true;
  else map_cuts_flag["NC1Pi0inFV"] = false; 
  
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_NprimPio==1 && pfeval.truth_NCDelta==0) map_cuts_flag["numuCC1Pi0inFV"] = true;
  else map_cuts_flag["numuCC1Pi0inFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_NprimPio!=1 && pfeval.truth_NCDelta==0) map_cuts_flag["numuCCotherinFV"] = true;
  else map_cuts_flag["numuCCotherinFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NprimPio!=1 && pfeval.truth_NCDelta==0) map_cuts_flag["NCotherinFV"] = true;
  else map_cuts_flag["NCotherinFV"] = false;

  // done with nc delta cuts

  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCMEC"] = true;
  else map_cuts_flag["CCMEC"] = false;
  
  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCMEC"] = true;
  else map_cuts_flag["NCMEC"] = false;
  
  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCQE"] = true;
  else map_cuts_flag["CCQE"] = false;

  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCQE"] = true;
  else map_cuts_flag["NCQE"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCRES"] = true;
  else map_cuts_flag["CCRES"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCRES"] = true;
  else map_cuts_flag["NCRES"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCDIS"] = true;
  else map_cuts_flag["CCDIS"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCDIS"] = true;
  else map_cuts_flag["NCDIS"] = false;
  
  if(pfeval.truth_nuScatType!=10 && pfeval.truth_nuScatType!=1 && pfeval.truth_nuScatType!=3 && pfeval.truth_nuScatType!=4 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["OTHER"] = true;
  else map_cuts_flag["OTHER"] = false;


  // figure out additional cuts and flag_data ...
  bool flag_add = true;
  if(add_cut == "all") flag_add = true;
  else if( (flag_data && (add_cut=="farsideband" || add_cut=="nearsideband" || add_cut=="nueCC" || add_cut=="generic" || add_cut=="loosenueCC" || add_cut=="loosenumuCC" )) || !flag_data ){ 
      std::istringstream sss(add_cut.Data());
      for(std::string line; std::getline(sss, line, '_');){
          if(map_cuts_flag.find(line)!=map_cuts_flag.end()){
              flag_add *= map_cuts_flag[line];
          }
          else{
              std::cout<<"ERROR: add_cut "<<line<<" not defined!\n";
              exit(EXIT_FAILURE);
          }
      } 
  }
  else{ 
    std::cout<<"ERROR: add_cut "<<add_cut<<" of channel "<< ch_name <<" is not assigned to sample "<<flag_data<<" [1: data; 0: mc]\n";
    std::cout<<"Please modify inc/WCPLEEANA/cuts.h\n";
    exit(EXIT_FAILURE);
  }
 
  if (!flag_add) return false;

  bool flag_generic = is_generic(eval);
  bool flag_numuCC = is_numuCC(tagger);
  bool flag_numuCC_loose = is_loosenumuCC(tagger);
  bool flag_numuCC_tight = is_numuCC_tight(tagger, pfeval);
  bool flag_numuCC_1mu0p = is_numuCC_1mu0p(tagger, kine, pfeval);
  bool flag_numuCC_lowEhad = is_numuCC_lowEhad(tagger, kine, pfeval, flag_data);
  bool flag_numuCC_cutbased = is_numuCC_cutbased(tagger);
  bool flag_nueCC = is_nueCC(tagger);
  bool flag_0p = is_0p(tagger, kine, pfeval);
  bool flag_1p = is_1p(tagger, kine, pfeval);
  bool flag_0pi = is_0pi(tagger, kine, pfeval);
  bool flag_nueCC_1e0p = is_nueCC_1e0p(tagger, kine, pfeval);
  bool flag_nueCC_1e0p0pi = is_nueCC_1e0p0pi(tagger, kine, pfeval);
  bool flag_nueCC_1eNp0pi = is_nueCC_1eNp0pi(tagger, kine, pfeval);
  bool flag_nueCC_loose = is_loosenueCC(tagger);
  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);
  bool flag_FC = is_FC(eval);
  bool flag_far = is_far_sideband(kine, tagger, flag_data);
  bool flag_near = is_near_sideband(kine, tagger, flag_data);
  bool flag_lowEnergy = is_lowEnergy(kine, flag_data);
  bool flag_ncpio_bdt = is_NCpio_bdt(tagger);
  bool flag_ncpio_bdt_sel = is_NCpio_bdt_sel(tagger, kine);
  bool flag_ncdelta_bdt = is_NCdelta_bdt(tagger, pfeval);
  bool flag_ncdelta_presel = is_NCdelta_presel(tagger, pfeval);
  bool flag_ncdelta_presel_signal_blind = is_NCdelta_presel_signal_blind(tagger, pfeval);
  bool flag_near_bdt = is_near_bdt(tagger, pfeval);

  if (ch_name == "LEE_near_nueoverlay"  || ch_name == "nueCC_near_nueoverlay" || ch_name == "nueCC_near_nueoverlay2" || ch_name == "nueCC_near_nueoverlay3"){
    if ( flag_near && flag_nueCC && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_ext" || ch_name == "BG_nueCC_near_dirt" || ch_name =="nueCC_near_bnb" || 
	    ch_name == "BG_nueCC_near_ext_2" || ch_name == "BG_nueCC_near_dirt_2" || ch_name =="nueCC_near_bnb2" || 
	    ch_name == "BG_nueCC_near_ext_3" || ch_name == "BG_nueCC_near_dirt_3" || ch_name =="nueCC_near_bnb3"){
    if (flag_near && flag_nueCC) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_overlay" || ch_name == "BG_nueCC_near_overlay2" || ch_name == "BG_nueCC_near_overlay3"){
    if ( flag_near && flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_far_nueoverlay"  || ch_name == "nueCC_far_nueoverlay" || ch_name == "nueCC_far_nueoverlay2" || ch_name == "nueCC_far_nueoverlay3"){
    if ( flag_far && flag_nueCC && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_ext" || ch_name == "BG_nueCC_far_dirt" || ch_name =="nueCC_far_bnb" || 
	    ch_name == "BG_nueCC_far_ext_2" || ch_name == "BG_nueCC_far_dirt_2" || ch_name =="nueCC_far_bnb2" || 
	    ch_name == "BG_nueCC_far_ext_3" || ch_name == "BG_nueCC_far_dirt_3" || ch_name =="nueCC_far_bnb3"){
    if (flag_far && flag_nueCC) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_overlay" || ch_name == "BG_nueCC_far_overlay2" || ch_name == "BG_nueCC_far_overlay3"){
    if ( flag_far && flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_nearfar_nueoverlay"  || ch_name == "nueCC_nearfar_nueoverlay" || ch_name == "nueCC_nearfar_nueoverlay2" || ch_name == "nueCC_nearfar_nueoverlay3"){
    if ( (flag_near || flag_far) && flag_nueCC && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_nearfar_ext" || ch_name == "BG_nueCC_nearfar_dirt" || ch_name =="nueCC_nearfar_bnb" || 
	    ch_name == "BG_nueCC_nearfar_ext_2" || ch_name == "BG_nueCC_nearfar_dirt_2" || ch_name =="nueCC_nearfar_bnb2" || 
	    ch_name == "BG_nueCC_nearfar_ext_3" || ch_name == "BG_nueCC_nearfar_dirt_3" || ch_name =="nueCC_nearfar_bnb3"){
    if ((flag_near || flag_far) && flag_nueCC) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_nearfar_overlay" || ch_name == "BG_nueCC_nearfar_overlay2" || ch_name == "BG_nueCC_nearfar_overlay3"){
    if ( (flag_near || flag_far) && flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;




  }else if (ch_name == "LEE_midenergy_nueoverlay"  || ch_name == "nueCC_midenergy_nueoverlay" || ch_name == "nueCC_midenergy_nueoverlay2" || ch_name == "nueCC_midenergy_nueoverlay3"){
    if ( (flag_far && reco_Enu<1400) && flag_nueCC && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_midenergy_ext" || ch_name == "BG_nueCC_midenergy_dirt" || ch_name =="nueCC_midenergy_bnb" || 
	    ch_name == "BG_nueCC_midenergy_ext_2" || ch_name == "BG_nueCC_midenergy_dirt_2" || ch_name =="nueCC_midenergy_bnb2" || 
	    ch_name == "BG_nueCC_midenergy_ext_3" || ch_name == "BG_nueCC_midenergy_dirt_3" || ch_name =="nueCC_midenergy_bnb3"){
    if ((flag_far && reco_Enu<1400) && flag_nueCC) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_midenergy_overlay" || ch_name == "BG_nueCC_midenergy_overlay2" || ch_name == "BG_nueCC_midenergy_overlay3"){
    if ( (flag_far && reco_Enu<1400) && flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_midenergy_1e0p_nueoverlay"  || ch_name == "nueCC_midenergy_1e0p_nueoverlay" || ch_name == "nueCC_midenergy_1e0p_nueoverlay2" || ch_name == "nueCC_midenergy_1e0p_nueoverlay3"){
    if ( (flag_far && reco_Enu<1400) && flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_midenergy_1e0p_ext" || ch_name == "BG_nueCC_midenergy_1e0p_dirt" || ch_name =="nueCC_midenergy_1e0p_bnb" || 
	    ch_name == "BG_nueCC_midenergy_1e0p_ext_2" || ch_name == "BG_nueCC_midenergy_1e0p_dirt_2" || ch_name =="nueCC_midenergy_1e0p_bnb2" || 
	    ch_name == "BG_nueCC_midenergy_1e0p_ext_3" || ch_name == "BG_nueCC_midenergy_1e0p_dirt_3" || ch_name =="nueCC_midenergy_1e0p_bnb3"){
    if ((flag_far && reco_Enu<1400) && flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_midenergy_1e0p_overlay" || ch_name == "BG_nueCC_midenergy_1e0p_overlay2" || ch_name == "BG_nueCC_midenergy_1e0p_overlay3"){
    if ( (flag_far && reco_Enu<1400) && flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_midenergy_1eNp_nueoverlay"  || ch_name == "nueCC_midenergy_1eNp_nueoverlay" || ch_name == "nueCC_midenergy_1eNp_nueoverlay2" || ch_name == "nueCC_midenergy_1eNp_nueoverlay3"){
    if ( (flag_far && reco_Enu<1400) && flag_nueCC && !flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_midenergy_1eNp_ext" || ch_name == "BG_nueCC_midenergy_1eNp_dirt" || ch_name =="nueCC_midenergy_1eNp_bnb" || 
	    ch_name == "BG_nueCC_midenergy_1eNp_ext_2" || ch_name == "BG_nueCC_midenergy_1eNp_dirt_2" || ch_name =="nueCC_midenergy_1eNp_bnb2" || 
	    ch_name == "BG_nueCC_midenergy_1eNp_ext_3" || ch_name == "BG_nueCC_midenergy_1eNp_dirt_3" || ch_name =="nueCC_midenergy_1eNp_bnb3"){
    if ((flag_far && reco_Enu<1400) && flag_nueCC && !flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_midenergy_1eNp_overlay" || ch_name == "BG_nueCC_midenergy_1eNp_overlay2" || ch_name == "BG_nueCC_midenergy_1eNp_overlay3"){
    if ( (flag_far && reco_Enu<1400) && flag_nueCC && !flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;



  }else if (ch_name == "LEE_near_1e0p_nueoverlay"  || ch_name == "nueCC_near_1e0p_nueoverlay" || ch_name == "nueCC_near_1e0p_nueoverlay2" || ch_name == "nueCC_near_1e0p_nueoverlay3"){
    if ( flag_near && flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_1e0p_ext" || ch_name == "BG_nueCC_near_1e0p_dirt" || ch_name =="nueCC_near_1e0p_bnb" || 
	    ch_name == "BG_nueCC_near_1e0p_ext_2" || ch_name == "BG_nueCC_near_1e0p_dirt_2" || ch_name =="nueCC_near_1e0p_bnb2" || 
	    ch_name == "BG_nueCC_near_1e0p_ext_3" || ch_name == "BG_nueCC_near_1e0p_dirt_3" || ch_name =="nueCC_near_1e0p_bnb3"){
    if (flag_near && flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_1e0p_overlay" || ch_name == "BG_nueCC_near_1e0p_overlay2" || ch_name == "BG_nueCC_near_1e0p_overlay3"){
    if ( flag_near && flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_far_1e0p_nueoverlay"  || ch_name == "nueCC_far_1e0p_nueoverlay" || ch_name == "nueCC_far_1e0p_nueoverlay2" || ch_name == "nueCC_far_1e0p_nueoverlay3"){
    if ( flag_far && flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_1e0p_ext" || ch_name == "BG_nueCC_far_1e0p_dirt" || ch_name =="nueCC_far_1e0p_bnb" || 
	    ch_name == "BG_nueCC_far_1e0p_ext_2" || ch_name == "BG_nueCC_far_1e0p_dirt_2" || ch_name =="nueCC_far_1e0p_bnb2" || 
	    ch_name == "BG_nueCC_far_1e0p_ext_3" || ch_name == "BG_nueCC_far_1e0p_dirt_3" || ch_name =="nueCC_far_1e0p_bnb3"){
    if (flag_far && flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_1e0p_overlay" || ch_name == "BG_nueCC_far_1e0p_overlay2" || ch_name == "BG_nueCC_far_1e0p_overlay3"){
    if ( flag_far && flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_nearfar_1e0p_nueoverlay"  || ch_name == "nueCC_nearfar_1e0p_nueoverlay" || ch_name == "nueCC_nearfar_1e0p_nueoverlay2" || ch_name == "nueCC_nearfar_1e0p_nueoverlay3"){
    if ( (flag_near || flag_far) && flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_nearfar_1e0p_ext" || ch_name == "BG_nueCC_nearfar_1e0p_dirt" || ch_name =="nueCC_nearfar_1e0p_bnb" || 
	    ch_name == "BG_nueCC_nearfar_1e0p_ext_2" || ch_name == "BG_nueCC_nearfar_1e0p_dirt_2" || ch_name =="nueCC_nearfar_1e0p_bnb2" || 
	    ch_name == "BG_nueCC_nearfar_1e0p_ext_3" || ch_name == "BG_nueCC_nearfar_1e0p_dirt_3" || ch_name =="nueCC_nearfar_1e0p_bnb3"){
    if ((flag_near || flag_far) && flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_nearfar_1e0p_overlay" || ch_name == "BG_nueCC_nearfar_1e0p_overlay2" || ch_name == "BG_nueCC_nearfar_1e0p_overlay3"){
    if ( (flag_near || flag_far) && flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;



  }else if (ch_name == "LEE_near_1eNp_nueoverlay"  || ch_name == "nueCC_near_1eNp_nueoverlay" || ch_name == "nueCC_near_1eNp_nueoverlay2" || ch_name == "nueCC_near_1eNp_nueoverlay3"){
    if ( flag_near && flag_nueCC && !flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_1eNp_ext" || ch_name == "BG_nueCC_near_1eNp_dirt" || ch_name =="nueCC_near_1eNp_bnb" || 
	    ch_name == "BG_nueCC_near_1eNp_ext_2" || ch_name == "BG_nueCC_near_1eNp_dirt_2" || ch_name =="nueCC_near_1eNp_bnb2" || 
	    ch_name == "BG_nueCC_near_1eNp_ext_3" || ch_name == "BG_nueCC_near_1eNp_dirt_3" || ch_name =="nueCC_near_1eNp_bnb3"){
    if (flag_near && flag_nueCC && !flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_1eNp_overlay" || ch_name == "BG_nueCC_near_1eNp_overlay2" || ch_name == "BG_nueCC_near_1eNp_overlay3"){
    if ( flag_near && flag_nueCC && !flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_far_1eNp_nueoverlay"  || ch_name == "nueCC_far_1eNp_nueoverlay" || ch_name == "nueCC_far_1eNp_nueoverlay2" || ch_name == "nueCC_far_1eNp_nueoverlay3"){
    if ( flag_far && flag_nueCC && !flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_1eNp_ext" || ch_name == "BG_nueCC_far_1eNp_dirt" || ch_name =="nueCC_far_1eNp_bnb" || 
	    ch_name == "BG_nueCC_far_1eNp_ext_2" || ch_name == "BG_nueCC_far_1eNp_dirt_2" || ch_name =="nueCC_far_1eNp_bnb2" || 
	    ch_name == "BG_nueCC_far_1eNp_ext_3" || ch_name == "BG_nueCC_far_1eNp_dirt_3" || ch_name =="nueCC_far_1eNp_bnb3"){
    if (flag_far && flag_nueCC && !flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_1eNp_overlay" || ch_name == "BG_nueCC_far_1eNp_overlay2" || ch_name == "BG_nueCC_far_1eNp_overlay3"){
    if ( flag_far && flag_nueCC && !flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;

  }else if (ch_name == "LEE_nearfar_1eNp_nueoverlay"  || ch_name == "nueCC_nearfar_1eNp_nueoverlay" || ch_name == "nueCC_nearfar_1eNp_nueoverlay2" || ch_name == "nueCC_nearfar_1eNp_nueoverlay3"){
    if ( (flag_near || flag_far) && flag_nueCC && !flag_nueCC_1e0p && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_nearfar_1eNp_ext" || ch_name == "BG_nueCC_nearfar_1eNp_dirt" || ch_name =="nueCC_nearfar_1eNp_bnb" || 
	    ch_name == "BG_nueCC_nearfar_1eNp_ext_2" || ch_name == "BG_nueCC_nearfar_1eNp_dirt_2" || ch_name =="nueCC_nearfar_1eNp_bnb2" || 
	    ch_name == "BG_nueCC_nearfar_1eNp_ext_3" || ch_name == "BG_nueCC_nearfar_1eNp_dirt_3" || ch_name =="nueCC_nearfar_1eNp_bnb3"){
    if ((flag_near || flag_far) && flag_nueCC && !flag_nueCC_1e0p) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_nearfar_1eNp_overlay" || ch_name == "BG_nueCC_nearfar_1eNp_overlay2" || ch_name == "BG_nueCC_nearfar_1eNp_overlay3"){
    if ( (flag_near || flag_far) && flag_nueCC && !flag_nueCC_1e0p && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;



  }else if (ch_name == "LEE_near_FC_nueoverlay"  || ch_name == "nueCC_near_FC_nueoverlay" || ch_name == "nueCC_near_FC_nueoverlay2" || ch_name == "nueCC_near_FC_nueoverlay3"){
    if ( flag_near && flag_nueCC && flag_FC && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_FC_ext" || ch_name == "BG_nueCC_near_FC_dirt" || ch_name =="nueCC_near_FC_bnb" || 
	    ch_name == "BG_nueCC_near_FC_ext_2" || ch_name == "BG_nueCC_near_FC_dirt_2" || ch_name =="nueCC_near_FC_bnb2" || 
	    ch_name == "BG_nueCC_near_FC_ext_3" || ch_name == "BG_nueCC_near_FC_dirt_3" || ch_name =="nueCC_near_FC_bnb3"){
    //nueCC FC
    if (flag_near && flag_nueCC && flag_FC) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_FC_overlay" || ch_name == "BG_nueCC_near_FC_overlay2" || ch_name == "BG_nueCC_near_FC_overlay3"){
    if ( flag_near && flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "LEE_near_PC_nueoverlay" || ch_name == "nueCC_near_PC_nueoverlay" || ch_name == "nueCC_near_PC_nueoverlay2" || ch_name == "nueCC_near_PC_nueoverlay3" ){
    // nueCC PC
    if ( flag_near && flag_nueCC && (!flag_FC) && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_PC_ext" || ch_name == "BG_nueCC_near_PC_dirt" || ch_name == "nueCC_near_PC_bnb" ||
	    ch_name == "BG_nueCC_near_PC_ext_2" || ch_name == "BG_nueCC_near_PC_dirt_2" || ch_name == "nueCC_near_PC_bnb2" ||
	    ch_name == "BG_nueCC_near_PC_ext_3" || ch_name == "BG_nueCC_near_PC_dirt_3" || ch_name == "nueCC_near_PC_bnb3"
	    ){
    // nueCC PC
    if ( flag_near && flag_nueCC && (!flag_FC)) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_near_PC_overlay" ||ch_name == "BG_nueCC_near_PC_overlay2" || ch_name == "BG_nueCC_near_PC_overlay3"){
    if ( flag_near && flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;


  }else if (ch_name == "LEE_far_FC_nueoverlay"  || ch_name == "nueCC_far_FC_nueoverlay" || ch_name == "nueCC_far_FC_nueoverlay2" || ch_name == "nueCC_far_FC_nueoverlay3"){
    if ( flag_far && flag_nueCC && flag_FC && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_FC_ext" || ch_name == "BG_nueCC_far_FC_dirt" || ch_name =="nueCC_far_FC_bnb" || 
	    ch_name == "BG_nueCC_far_FC_ext_2" || ch_name == "BG_nueCC_far_FC_dirt_2" || ch_name =="nueCC_far_FC_bnb2" || 
	    ch_name == "BG_nueCC_far_FC_ext_3" || ch_name == "BG_nueCC_far_FC_dirt_3" || ch_name =="nueCC_far_FC_bnb3"){
    //nueCC FC
    if (flag_far && flag_nueCC && flag_FC) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_FC_overlay" || ch_name == "BG_nueCC_far_FC_overlay2" || ch_name == "BG_nueCC_far_FC_overlay3"){
    if ( flag_far && flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;
  }else if (ch_name == "LEE_far_PC_nueoverlay" || ch_name == "nueCC_far_PC_nueoverlay" || ch_name == "nueCC_far_PC_nueoverlay2" || ch_name == "nueCC_far_PC_nueoverlay3" ){
    // nueCC PC
    if ( flag_far && flag_nueCC && (!flag_FC) && flag_truth_inside) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_PC_ext" || ch_name == "BG_nueCC_far_PC_dirt" || ch_name == "nueCC_far_PC_bnb" ||
	    ch_name == "BG_nueCC_far_PC_ext_2" || ch_name == "BG_nueCC_far_PC_dirt_2" || ch_name == "nueCC_far_PC_bnb2" ||
	    ch_name == "BG_nueCC_far_PC_ext_3" || ch_name == "BG_nueCC_far_PC_dirt_3" || ch_name == "nueCC_far_PC_bnb3"
	    ){
    // nueCC PC
    if ( flag_far && flag_nueCC && (!flag_FC)) return true; 
    else return false;
  }else if (ch_name == "BG_nueCC_far_PC_overlay" ||ch_name == "BG_nueCC_far_PC_overlay2" || ch_name == "BG_nueCC_far_PC_overlay3"){
    if ( flag_far && flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; 
    else return false;


  }else if (ch_name == "LEE_FC_nueoverlay"  || ch_name == "nueCC_FC_nueoverlay" || ch_name == "nueCC_FC_nueoverlay2" || ch_name == "nueCC_FC_nueoverlay3"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_ext" || ch_name == "BG_nueCC_FC_dirt" || ch_name =="nueCC_FC_bnb" || 
	    ch_name == "BG_nueCC_FC_ext_2" || ch_name == "BG_nueCC_FC_dirt_2" || ch_name =="nueCC_FC_bnb2" || 
	    ch_name == "BG_nueCC_FC_ext_3" || ch_name == "BG_nueCC_FC_dirt_3" || ch_name =="nueCC_FC_bnb3"){
    //nueCC FC
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_overlay" || ch_name == "BG_nueCC_FC_overlay2" || ch_name == "BG_nueCC_FC_overlay3"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "LEE_PC_nueoverlay" || ch_name == "nueCC_PC_nueoverlay" || ch_name == "nueCC_PC_nueoverlay2" || ch_name == "nueCC_PC_nueoverlay3" ){
    // nueCC PC
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_ext" || ch_name == "BG_nueCC_PC_dirt" || ch_name == "nueCC_PC_bnb" ||
	    ch_name == "BG_nueCC_PC_ext_2" || ch_name == "BG_nueCC_PC_dirt_2" || ch_name == "nueCC_PC_bnb2" ||
	    ch_name == "BG_nueCC_PC_ext_3" || ch_name == "BG_nueCC_PC_dirt_3" || ch_name == "nueCC_PC_bnb3"
	    ){
    // nueCC PC
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_overlay" ||ch_name == "BG_nueCC_PC_overlay2" || ch_name == "BG_nueCC_PC_overlay3"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;

  }else if (ch_name == "numuCC_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_nopi0_nonueCC_FC_bnb" ||
	    ch_name == "numuCC_nopi0_nonueCC_FC_overlay2" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext_2" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt_2" || ch_name == "numuCC_nopi0_nonueCC_FC_bnb2" ||
	    ch_name == "numuCC_nopi0_nonueCC_FC_overlay3" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext_3" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt_3" || ch_name == "numuCC_nopi0_nonueCC_FC_bnb3"
	    ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_nopi0_nonueCC_PC_bnb" ||
	    ch_name == "numuCC_nopi0_nonueCC_PC_overlay2" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext_2" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt_2" || ch_name == "numuCC_nopi0_nonueCC_PC_bnb2" ||
	    ch_name == "numuCC_nopi0_nonueCC_PC_overlay3" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext_3" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt_3" || ch_name == "numuCC_nopi0_nonueCC_PC_bnb3" 
){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_FC_overlay" || ch_name =="BG_CCpi0_nonueCC_FC_ext" || ch_name == "BG_CCpi0_nonueCC_FC_dirt" || ch_name == "CCpi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi02_nonueCC_FC_overlay" || ch_name =="BG_CCpi02_nonueCC_FC_ext" || ch_name == "BG_CCpi02_nonueCC_FC_dirt" || ch_name == "CCpi02_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi03_nonueCC_FC_overlay" || ch_name =="BG_CCpi03_nonueCC_FC_ext" || ch_name == "BG_CCpi03_nonueCC_FC_dirt" || ch_name == "CCpi03_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_PC_overlay" || ch_name == "BG_CCpi0_nonueCC_PC_ext" || ch_name == "BG_CCpi0_nonueCC_PC_dirt" || ch_name == "CCpi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi02_nonueCC_PC_overlay" || ch_name == "BG_CCpi02_nonueCC_PC_ext" || ch_name == "BG_CCpi02_nonueCC_PC_dirt" || ch_name == "CCpi02_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi03_nonueCC_PC_overlay" || ch_name == "BG_CCpi03_nonueCC_PC_ext" || ch_name == "BG_CCpi03_nonueCC_PC_dirt" || ch_name == "CCpi03_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_overlay" || ch_name == "BG_NCpi0_nonueCC_ext" || ch_name == "BG_NCpi0_nonueCC_dirt" || ch_name == "NCpi0_nonueCC_bnb"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi02_nonueCC_overlay" || ch_name == "BG_NCpi02_nonueCC_ext" || ch_name == "BG_NCpi02_nonueCC_dirt" || ch_name == "NCpi02_nonueCC_bnb"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi03_nonueCC_overlay" || ch_name == "BG_NCpi03_nonueCC_ext" || ch_name == "BG_NCpi03_nonueCC_dirt" || ch_name == "NCpi03_nonueCC_bnb"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_FC_overlay" || ch_name == "BG_NCpi0_nonueCC_FC_ext" || ch_name == "BG_NCpi0_nonueCC_FC_dirt" || ch_name == "NCpi0_nonueCC_FC_bnb"){
    if (flag_NC && flag_pi0 && flag_FC && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi02_nonueCC_FC_overlay" || ch_name == "BG_NCpi02_nonueCC_FC_ext" || ch_name == "BG_NCpi02_nonueCC_FC_dirt" || ch_name == "NCpi02_nonueCC_FC_bnb"){
    if (flag_NC && flag_pi0 && flag_FC && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi03_nonueCC_FC_overlay" || ch_name == "BG_NCpi03_nonueCC_FC_ext" || ch_name == "BG_NCpi03_nonueCC_FC_dirt" || ch_name == "NCpi03_nonueCC_FC_bnb"){
    if (flag_NC && flag_pi0 && flag_FC && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_PC_overlay" || ch_name == "BG_NCpi0_nonueCC_PC_ext" || ch_name == "BG_NCpi0_nonueCC_PC_dirt" || ch_name == "NCpi0_nonueCC_PC_bnb"){
    if (flag_NC && flag_pi0 && (!flag_FC) && (!flag_nueCC) ) return true;
    //if (flag_NC && flag_pi0 && (!flag_0p) && (!flag_FC) && (!flag_nueCC) ) return true; // hack, Np
    else return false;
  }else if (ch_name == "NCpi02_nonueCC_PC_overlay" || ch_name == "BG_NCpi02_nonueCC_PC_ext" || ch_name == "BG_NCpi02_nonueCC_PC_dirt" || ch_name == "NCpi02_nonueCC_PC_bnb"){
    if (flag_NC && flag_pi0 && (!flag_FC) && (!flag_nueCC) ) return true;
    //if (flag_NC && flag_pi0 && (!flag_0p) && (!flag_FC) && (!flag_nueCC) ) return true; // hack, Np
    else return false;
  }else if (ch_name == "NCpi03_nonueCC_PC_overlay" || ch_name == "BG_NCpi03_nonueCC_PC_ext" || ch_name == "BG_NCpi03_nonueCC_PC_dirt" || ch_name == "NCpi03_nonueCC_PC_bnb"){
    if (flag_NC && flag_pi0 && (!flag_FC) && (!flag_nueCC) ) return true;
    //if (flag_NC && flag_pi0 && (!flag_0p) && (!flag_FC) && (!flag_nueCC) ) return true; // hack, Np
    else return false;
  }else if (ch_name == "nueCC_bnb" || ch_name == "nueCC_nueoverlay"){   // side band ...
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay" || ch_name == "nueCC_bnb") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb" || ch_name == "all_but_nueCC_overlay" || ch_name == "all_but_nueCC_ext" || ch_name == "all_but_nueCC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay" || ch_name != "all_but_nueCC_overlay") return true;
    else return false;
  }else if (ch_name == "nueCC_bnb1" || ch_name == "nueCC_nueoverlay1"){
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay1" || ch_name == "nueCC_bnb1") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb1" || ch_name == "all_but_nueCC_overlay1" || ch_name == "all_but_nueCC_ext1" || ch_name == "all_but_nueCC_dirt1"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay1" || ch_name != "all_but_nueCC_overlay1") return true;
    else return false;
  }

  // dividing nueCC into 0p and Np, for BNB 
  else if(ch_name == "nueCC_lowEhad_FC_nueoverlay"){
    if (flag_nueCC_1e0p && flag_FC && flag_truth_inside) return true;
    //if (flag_nueCC_1e0p && flag_near && flag_FC && flag_truth_inside) return true; // near sideband hack

    else return false;
  }
  else if (ch_name == "BG_nueCC_lowEhad_FC_ext" || ch_name =="BG_nueCC_lowEhad_FC_dirt" || ch_name == "nueCC_lowEhad_FC_bnb" || ch_name == "nueCC_lowEhad_FC_numi"){
    if (flag_nueCC_1e0p && flag_FC) return true;
    //if (flag_nueCC_1e0p && flag_near && flag_FC) return true; // near sideband hack

    else return false;
  }
  else if (ch_name == "BG_nueCC_lowEhad_FC_overlay"){
    if (flag_nueCC_1e0p && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    //if (flag_nueCC_1e0p && flag_FC && flag_near && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; // near sideband hack
    else return false;
  }
  else if(ch_name == "nueCC_lowEhad_PC_nueoverlay"){
    if (flag_nueCC_1e0p && (!flag_FC) && flag_truth_inside) return true;
    //if (flag_nueCC_1e0p && flag_near && (!flag_FC) && flag_truth_inside) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_lowEhad_PC_ext" || ch_name =="BG_nueCC_lowEhad_PC_dirt" || ch_name == "nueCC_lowEhad_PC_bnb"|| ch_name == "nueCC_lowEhad_PC_numi"){
    if (flag_nueCC_1e0p && (!flag_FC) ) return true;
    //if (flag_nueCC_1e0p && flag_near && (!flag_FC) ) return true; // near sideband hack

    else return false;
  }
  else if (ch_name == "BG_nueCC_lowEhad_PC_overlay"){
    if (flag_nueCC_1e0p && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true;
    //if (flag_nueCC_1e0p && flag_near && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true; // near sideband hack

    else return false;
  }
  else if(ch_name == "nueCC_highEhad_FC_nueoverlay"){
    if (flag_nueCC && (!flag_nueCC_1e0p) && flag_FC && flag_truth_inside ) return true;
    //if (flag_nueCC && flag_near && (!flag_nueCC_1e0p) && flag_FC && flag_truth_inside ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_highEhad_FC_ext" || ch_name =="BG_nueCC_highEhad_FC_dirt" || ch_name == "nueCC_highEhad_FC_bnb" || ch_name == "nueCC_highEhad_FC_numi"){
    if (flag_nueCC && (!flag_nueCC_1e0p) && flag_FC ) return true;
    //if (flag_nueCC && flag_near && (!flag_nueCC_1e0p) && flag_FC ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_highEhad_FC_overlay"){
    if (flag_nueCC && (!flag_nueCC_1e0p) && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    //if (flag_nueCC && flag_near && (!flag_nueCC_1e0p) && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; // near sideband hack
    else return false;
  }
  else if(ch_name == "nueCC_highEhad_PC_nueoverlay"){
    if (flag_nueCC && (!flag_nueCC_1e0p) && (!flag_FC) && flag_truth_inside ) return true;
    //if (flag_nueCC && flag_near && (!flag_nueCC_1e0p) && (!flag_FC) && flag_truth_inside ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_highEhad_PC_ext" || ch_name =="BG_nueCC_highEhad_PC_dirt" || ch_name == "nueCC_highEhad_PC_bnb" || ch_name == "nueCC_highEhad_PC_numi"){
    if (flag_nueCC && (!flag_nueCC_1e0p) && (!flag_FC) ) return true;
    //if (flag_nueCC && flag_near && (!flag_nueCC_1e0p) && (!flag_FC) ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_highEhad_PC_overlay"){
    if (flag_nueCC && (!flag_nueCC_1e0p) && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true;
    //if (flag_nueCC && flag_near && (!flag_nueCC_1e0p) && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true; // near sideband hack
    else return false;
  }
  ////////////////////////
  // 1e0/Np0pi selection
  else if(ch_name == "nueCC_1e0p0pi_FC_nueoverlay"){
    //if (flag_nueCC_1e0p0pi && flag_FC && flag_truth_inside) return true;
    if (flag_nueCC_1e0p0pi && (flag_near || flag_far) && flag_FC && flag_truth_inside) return true; // near sideband hack

    else return false;
  }
  else if (ch_name == "BG_nueCC_1e0p0pi_FC_ext" || ch_name =="BG_nueCC_1e0p0pi_FC_dirt" || ch_name == "nueCC_1e0p0pi_FC_bnb" || ch_name == "nueCC_1e0p0pi_FC_numi"){
    //if (flag_nueCC_1e0p0pi && flag_FC) return true;
    if (flag_nueCC_1e0p0pi && (flag_near || flag_far) && flag_FC) return true; // near sideband hack

    else return false;
  }
  else if (ch_name == "BG_nueCC_1e0p0pi_FC_overlay"){
    //if (flag_nueCC_1e0p0pi && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    if (flag_nueCC_1e0p0pi && flag_FC && (flag_near || flag_far) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; // near sideband hack
    else return false;
  }
  else if(ch_name == "nueCC_1e0p0pi_PC_nueoverlay"){
    //if (flag_nueCC_1e0p0pi && (!flag_FC) && flag_truth_inside) return true;
    if (flag_nueCC_1e0p0pi && (flag_near || flag_far) && (!flag_FC) && flag_truth_inside) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_1e0p0pi_PC_ext" || ch_name =="BG_nueCC_1e0p0pi_PC_dirt" || ch_name == "nueCC_1e0p0pi_PC_bnb"|| ch_name == "nueCC_1e0p0pi_PC_numi"){
    //if (flag_nueCC_1e0p0pi && (!flag_FC) ) return true;
    if (flag_nueCC_1e0p0pi && (flag_near || flag_far) && (!flag_FC) ) return true; // near sideband hack

    else return false;
  }
  else if (ch_name == "BG_nueCC_1e0p0pi_PC_overlay"){
    //if (flag_nueCC_1e0p0pi && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true;
    if (flag_nueCC_1e0p0pi && (flag_near || flag_far) && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true; // near sideband hack

    else return false;
  }
  else if(ch_name == "nueCC_1eNp0pi_FC_nueoverlay"){
    //if (flag_nueCC && flag_nueCC_1eNp0pi && flag_FC && flag_truth_inside ) return true;
    if (flag_nueCC && (flag_near || flag_far) && flag_nueCC_1eNp0pi && flag_FC && flag_truth_inside ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_1eNp0pi_FC_ext" || ch_name =="BG_nueCC_1eNp0pi_FC_dirt" || ch_name == "nueCC_1eNp0pi_FC_bnb" || ch_name == "nueCC_1eNp0pi_FC_numi"){
    //if (flag_nueCC && flag_nueCC_1eNp0pi && flag_FC ) return true;
    if (flag_nueCC && (flag_near || flag_far) && flag_nueCC_1eNp0pi && flag_FC ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_1eNp0pi_FC_overlay"){
    //if (flag_nueCC && flag_nueCC_1eNp0pi && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    if (flag_nueCC && (flag_near || flag_far) && flag_nueCC_1eNp0pi && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true; // near sideband hack
    else return false;
  }
  else if(ch_name == "nueCC_1eNp0pi_PC_nueoverlay"){
    //if (flag_nueCC && flag_nueCC_1eNp0pi && (!flag_FC) && flag_truth_inside ) return true;
    if (flag_nueCC && (flag_near || flag_far) && flag_nueCC_1eNp0pi && (!flag_FC) && flag_truth_inside ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_1eNp0pi_PC_ext" || ch_name =="BG_nueCC_1eNp0pi_PC_dirt" || ch_name == "nueCC_1eNp0pi_PC_bnb" || ch_name == "nueCC_1eNp0pi_PC_numi"){
    //if (flag_nueCC && flag_nueCC_1eNp0pi && (!flag_FC) ) return true;
    if (flag_nueCC && (flag_near || flag_far) && flag_nueCC_1eNp0pi && (!flag_FC) ) return true; // near sideband hack
    else return false;
  }
  else if (ch_name == "BG_nueCC_1eNp0pi_PC_overlay"){
    //if (flag_nueCC && flag_nueCC_1eNp0pi && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true;
    if (flag_nueCC && (flag_near || flag_far) && flag_nueCC_1eNp0pi && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) ) return true; // near sideband hack
    else return false;
  }


  // far sideband nue score plot
  else if (ch_name == "testA_bnb" || ch_name == "testA_nueoverlay"){
    //if (flag_truth_inside && flag_nueCC_loose && flag_far && ch_name == "testA_nueoverlay" || (flag_nueCC_loose && flag_far && ch_name == "testA_bnb")) return true; // loose nueCC
    if ( flag_FC && reco_Enu<=600 && flag_truth_inside && flag_nueCC_loose && ch_name == "testA_nueoverlay" || (flag_FC && reco_Enu<=600 && flag_nueCC_loose && ch_name == "testA_bnb")) return true; // loose nueCC FC open data low E
    //if ( reco_Enu<=500 && flag_truth_inside && flag_nueCC_loose && ch_name == "testA_nueoverlay" || (reco_Enu<=500 && flag_nueCC_loose && ch_name == "testA_bnb")) return true; // loose nueCC hack
    //if (flag_FC && flag_truth_inside && flag_nueCC_loose && flag_far && ch_name == "testA_nueoverlay" || (flag_FC && flag_nueCC_loose && flag_far && ch_name == "testA_bnb")) return true; // loose nueCC FC
    //if ( (flag_numuCC_loose && flag_far && ch_name == "testA_bnb")) return true; // loose numuCC
    //if ( (flag_FC && flag_numuCC_loose && flag_far && ch_name == "testA_bnb")) return true; // loose numuCC FC
    else return false;
  }
  else if (ch_name == "testA_overlay" || ch_name == "testA_ext" || ch_name == "testA_dirt"){
    //if ( (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_nueCC_loose && flag_far && ch_name == "testA_overlay") || (flag_nueCC_loose && flag_far && ch_name != "testA_overlay") ) return true; // loose nueCC
    if ( (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_nueCC_loose && flag_FC && reco_Enu<=600 && ch_name == "testA_overlay") || (flag_nueCC_loose && flag_FC && reco_Enu<=600&& ch_name != "testA_overlay") ) return true; // loose nueCC FC open data low E
    //if ( (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_nueCC_loose && reco_Enu<=500 && ch_name == "testA_overlay") || (flag_nueCC_loose && reco_Enu<=500 && ch_name != "testA_overlay") ) return true; // loose nueCC hack
    //if ( flag_FC &&(!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_nueCC_loose && flag_far && ch_name == "testA_overlay") || (flag_FC && flag_nueCC_loose && flag_far && ch_name != "testA_overlay") ) return true; // loose nueCC FC
    //if ( (flag_numuCC_loose && (!flag_nueCC) && flag_far && ch_name == "testA_overlay") || (flag_numuCC_loose && flag_far && ch_name != "testA_overlay") ) return true; // loose numuCC
    //if ( (reco_Enu<=500 && flag_numuCC_loose && (!flag_nueCC) && flag_far && ch_name == "testA_overlay") || (reco_Enu<=500 && flag_numuCC_loose && flag_far && ch_name != "testA_overlay") ) return true; // loose numuCC hack
    //if ( (flag_FC && flag_numuCC_loose && (!flag_nueCC) && flag_far && ch_name == "testA_overlay") || (flag_FC && flag_numuCC_loose && flag_far && ch_name != "testA_overlay") ) return true; // loose numuCC FC
    else return false;
  }

  else if (ch_name == "testB_bnb" || ch_name == "testB_nueoverlay"){
    //if (!flag_FC && flag_truth_inside && flag_nueCC_loose && flag_far && ch_name == "testB_nueoverlay" || (!flag_FC && flag_nueCC_loose && flag_far && ch_name == "testB_bnb")) return true; // loose nueCC PC
    if (flag_FC && reco_Enu<=600 && flag_truth_inside && flag_nueCC_loose && ch_name == "testB_nueoverlay" || (flag_FC && reco_Enu<=600 && flag_nueCC_loose && ch_name == "testB_bnb")) return true; // loose nueCC FC open data low E
    //if ( (!flag_FC && flag_numuCC_loose && flag_far && ch_name == "testB_bnb")) return true; // loose numuCC PC
    else return false;
  }
  else if (ch_name == "testB_overlay" || ch_name == "testB_ext" || ch_name == "testB_dirt"){
    //if ( flag_numuCC && (!flag_nueCC) && (!flag_cc_pi0) && flag_far ) return true; // numu score case
    //if ( !flag_FC && (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_nueCC_loose && flag_far && ch_name == "testB_overlay") || (!flag_FC && flag_nueCC_loose && flag_far && ch_name != "testB_overlay") ) return true; // loose nueCC PC
    if ( flag_FC && reco_Enu<=600 && (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_nueCC_loose && ch_name == "testB_overlay") || (flag_FC && reco_Enu<=600 && flag_nueCC_loose && ch_name != "testB_overlay") ) return true; // loose nueCC FC open data low E
    //if ( (!flag_FC && flag_numuCC_loose && (!flag_nueCC) && flag_far && ch_name == "testB_overlay") || (!flag_FC && flag_numuCC_loose && flag_far && ch_name != "testB_overlay") ) return true; // loose numuCC PC
    else return false;
  }

  else if (ch_name == "testC_bnb" || ch_name == "testC_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testC_nueoverlay" || ch_name == "testC_bnb") return true;
    else return false;
  }else if (ch_name == "testC_overlay" || ch_name == "testC_ext" || ch_name == "testC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testC_overlay" || ch_name != "testC_overlay") return true;
    else return false;
  }else if (ch_name == "testD_bnb" || ch_name == "testD_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testD_nueoverlay" || ch_name == "testD_bnb") return true;
    else return false;
  }else if (ch_name == "testD_overlay" || ch_name == "testD_ext" || ch_name == "testD_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testD_overlay" || ch_name != "testD_overlay") return true;
    else return false;
 // Janet's requests: <600 MeV numuCC PC, FC for three variables = 6 obs channels 
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;

  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;

  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_FC_bnb"){
    //if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    if (reco_Enu<500 && flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_PC_bnb"){
    //if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    if (reco_Enu<500 && flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;


  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;


  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_bnb"){
    //if (flag_numuCC_tight && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    if (reco_Enu<600 && flag_numuCC_tight && !flag_numuCC_1mu0p  && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC  && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_bnb"){
    //if (flag_numuCC_tight && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    if ((reco_Enu>600 && reco_Enu<=800) && flag_numuCC_tight && !flag_numuCC_1mu0p && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC  && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra3_nopi0_nonueCC_overlay" || ch_name == "BG_numuCC_extra3_nopi0_nonueCC_ext" || ch_name =="BG_numuCC_extra3_nopi0_nonueCC_dirt" || ch_name == "numuCC_extra3_nopi0_nonueCC_bnb"){
    if (reco_Enu>800 && flag_numuCC_tight && !flag_numuCC_1mu0p && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_numuCC_tight && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC  && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra4_nopi0_nonueCC_overlay" || ch_name == "BG_numuCC_extra4_nopi0_nonueCC_ext" || ch_name =="BG_numuCC_extra4_nopi0_nonueCC_dirt" || ch_name == "numuCC_extra4_nopi0_nonueCC_bnb"){
    if ((reco_Enu>600 && reco_Enu<1400) && flag_numuCC_tight && !flag_numuCC_1mu0p && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_numuCC_tight && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_far && flag_numuCC  && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;



  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu<500) return true; // temp
    else return false;
  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu<500) return true; // temp
    else return false;

  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu<500) return true; // temp
    else return false;
  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    //if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu<500) return true; // temp
    else return false;

  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

 // Mike Shaevitz >800 MeV nueCC PC+FC 1 obs channel
  }else if (ch_name == "nueCC_extra_nueoverlay"){
    if (flag_nueCC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext" || ch_name == "BG_nueCC_extra_dirt" || ch_name =="nueCC_extra_bnb"){
    if (flag_nueCC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
 // cut-based numuCC FC/PC 2 obs channels   
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_cutbased && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_cutbased && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

 // generic selection nu PC+FC 1 obs channel   
  }else if (ch_name == "generic_nu_overlay" || ch_name == "BG_generic_nu_ext" || ch_name =="BG_generic_nu_dirt" || ch_name == "generic_nu_bnb"){
    if (flag_generic) return true;
    else return false;
 // numuCC selection PC+FC 1 obs channel   
  }else if (ch_name == "numuCC_overlay" || ch_name == "BG_numuCC_ext" || ch_name =="BG_numuCC_dirt" || ch_name == "numuCC_bnb"){
    if (flag_numuCC) return true;
    else return false;
 // cutbased numuCC selection PC+FC 1 obs channel   
  }else if (ch_name == "numuCC_cutbased_overlay" || ch_name == "BG_numuCC_cutbased_ext" || ch_name =="BG_numuCC_cutbased_dirt" || ch_name == "numuCC_cutbased_bnb"){
    if (flag_numuCC_cutbased) return true;
    else return false;
 // nueCC 3 variables: n_trakcs, n_showers, gap_n_bad, FC/PC x3 = 6 channels; 4 additional channels 
  }else if (ch_name == "nueCC2_FC_nueoverlay"){
    //if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    if (flag_near && flag_nueCC && flag_FC && flag_truth_inside) return true;// hack
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_ext" || ch_name == "BG_nueCC2_FC_dirt" || ch_name =="nueCC2_FC_bnb"){
    //if (flag_nueCC && flag_FC) return true;
    if (flag_near && flag_nueCC && flag_FC) return true;//hack
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_overlay"){
    //if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    if (flag_near && flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;//hack
    else return false;
  }else if (ch_name == "nueCC2_PC_nueoverlay" ){
    //if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    if (flag_near && flag_nueCC && (!flag_FC) && flag_truth_inside) return true;//hack
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_ext" || ch_name == "BG_nueCC2_PC_dirt" || ch_name == "nueCC2_PC_bnb"){
    //if (flag_nueCC && (!flag_FC)) return true;
    if (flag_near && flag_nueCC && (!flag_FC)) return true;//hack
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_overlay"){
    //if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    if (flag_near && flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;//hack
    else return false;
  }else if (ch_name == "nueCC3_FC_nueoverlay"){
    //if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    if (flag_near && flag_nueCC && flag_FC && flag_truth_inside) return true; //hack
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_ext" || ch_name == "BG_nueCC3_FC_dirt" || ch_name =="nueCC3_FC_bnb"){
    //if (flag_nueCC && flag_FC) return true;
    if (flag_near && flag_nueCC && flag_FC) return true;//hack
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_overlay"){
    //if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    if (flag_near && flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;//hack
    else return false;
  }else if (ch_name == "nueCC3_PC_nueoverlay" ){
    //if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    if (flag_near && flag_nueCC && (!flag_FC) && flag_truth_inside) return true;//hack
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_ext" || ch_name == "BG_nueCC3_PC_dirt" || ch_name == "nueCC3_PC_bnb"){
    //if (flag_nueCC && (!flag_FC)) return true;
    if (flag_near && flag_nueCC && (!flag_FC)) return true;//hack
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_overlay"){
    //if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    if (flag_near && flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;//hack
    else return false;
    // add some cuts for Xs related cases ...
  }else if (ch_name == "numuCC_FC_bnb" || ch_name == "BG_numuCC_FC_ext" || ch_name == "BG_numuCC_FC_dirt"
	    || ch_name == "numuCC1_FC_bnb" || ch_name == "BG_numuCC1_FC_ext" || ch_name == "BG_numuCC1_FC_dirt"
	    || ch_name == "numuCC2_FC_bnb" || ch_name == "BG_numuCC2_FC_ext" || ch_name == "BG_numuCC2_FC_dirt"
	    ){
    if (flag_numuCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb" || ch_name == "BG_numuCC_PC_ext" || ch_name == "BG_numuCC_PC_dirt"
	    || ch_name == "numuCC1_PC_bnb" || ch_name == "BG_numuCC1_PC_ext" || ch_name == "BG_numuCC1_PC_dirt"
	    || ch_name == "numuCC2_PC_bnb" || ch_name == "BG_numuCC2_PC_ext" || ch_name == "BG_numuCC2_PC_dirt"
	    ){
    //if (flag_FC) return true; // quick hack test ...
    if (flag_numuCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC_background_PC_overlay"
  	    || ch_name == "numuCC1_signal_FC_overlay" || ch_name == "numuCC1_signal_PC_overlay" || ch_name == "numuCC1_background_FC_overlay" || ch_name == "numuCC1_background_PC_overlay"
  	    || ch_name == "numuCC2_signal_FC_overlay" || ch_name == "numuCC2_signal_PC_overlay" || ch_name == "numuCC2_background_FC_overlay" || ch_name == "numuCC2_background_PC_overlay"
  	    ){
    if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC1_signal_FC_overlay" || ch_name == "numuCC2_signal_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC1_signal_PC_overlay" || ch_name == "numuCC2_signal_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC1_background_FC_overlay" || ch_name == "numuCC2_background_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_PC_overlay" || ch_name == "numuCC1_background_PC_overlay" || ch_name == "numuCC2_background_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }
    return false;
  } else if (ch_name == "numuCC_signal_Enu_FC_overlay" || ch_name == "numuCC_signal_Enu_PC_overlay" || ch_name == "numuCC_background_Enu_FC_overlay" || ch_name == "numuCC_background_Enu_PC_overlay"
	     || ch_name == "numuCC1_signal_Enu_FC_overlay" || ch_name == "numuCC1_signal_Enu_PC_overlay" || ch_name == "numuCC1_background_Enu_FC_overlay" || ch_name == "numuCC1_background_Enu_PC_overlay"
	     || ch_name == "numuCC2_signal_Enu_FC_overlay" || ch_name == "numuCC2_signal_Enu_PC_overlay" || ch_name == "numuCC2_background_Enu_FC_overlay" || ch_name == "numuCC2_background_Enu_PC_overlay"
	    ){
    if (ch_name == "numuCC_signal_Enu_FC_overlay" || ch_name == "numuCC1_signal_Enu_FC_overlay" || ch_name == "numuCC2_signal_Enu_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Enu_PC_overlay" || ch_name == "numuCC1_signal_Enu_PC_overlay" || ch_name == "numuCC2_signal_Enu_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Enu_FC_overlay" || ch_name == "numuCC1_background_Enu_FC_overlay" || ch_name == "numuCC2_background_Enu_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Enu_PC_overlay" || ch_name == "numuCC1_background_Enu_PC_overlay" || ch_name == "numuCC2_background_Enu_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }
    return false;
    
  } else if (ch_name == "numuCC_signal_Emu_FC_overlay" || ch_name == "numuCC_signal_Emu_PC_overlay" || ch_name == "numuCC_background_Emu_FC_overlay" || ch_name == "numuCC_background_Emu_PC_overlay"
	     || ch_name == "numuCC1_signal_Emu_FC_overlay" || ch_name == "numuCC1_signal_Emu_PC_overlay" || ch_name == "numuCC1_background_Emu_FC_overlay" || ch_name == "numuCC1_background_Emu_PC_overlay"
	     || ch_name == "numuCC2_signal_Emu_FC_overlay" || ch_name == "numuCC2_signal_Emu_PC_overlay" || ch_name == "numuCC2_background_Emu_FC_overlay" || ch_name == "numuCC2_background_Emu_PC_overlay"
	     ){
    if (ch_name == "numuCC_signal_Emu_FC_overlay" || ch_name == "numuCC1_signal_Emu_FC_overlay" || ch_name == "numuCC2_signal_Emu_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Emu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Emu_PC_overlay" || ch_name == "numuCC1_signal_Emu_PC_overlay" || ch_name == "numuCC2_signal_Emu_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Emu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Emu_FC_overlay" || ch_name == "numuCC1_background_Emu_FC_overlay" || ch_name == "numuCC2_background_Emu_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Emu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Emu_PC_overlay" || ch_name == "numuCC1_background_Emu_PC_overlay" || ch_name == "numuCC2_background_Emu_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Emu_numuCCinFV"])) return true;
    }
    return false;

  } else if (ch_name == "numuCC_signal_Ehad_FC_overlay" || ch_name == "numuCC_signal_Ehad_PC_overlay" || ch_name == "numuCC_background_Ehad_FC_overlay" || ch_name == "numuCC_background_Ehad_PC_overlay"
	     || ch_name == "numuCC1_signal_Ehad_FC_overlay" || ch_name == "numuCC1_signal_Ehad_PC_overlay" || ch_name == "numuCC1_background_Ehad_FC_overlay" || ch_name == "numuCC1_background_Ehad_PC_overlay"
	     || ch_name == "numuCC2_signal_Ehad_FC_overlay" || ch_name == "numuCC2_signal_Ehad_PC_overlay" || ch_name == "numuCC2_background_Ehad_FC_overlay" || ch_name == "numuCC2_background_Ehad_PC_overlay"
	    ){
    if (ch_name == "numuCC_signal_Ehad_FC_overlay" || ch_name == "numuCC1_signal_Ehad_FC_overlay" || ch_name == "numuCC2_signal_Ehad_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Ehad_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Ehad_PC_overlay" || ch_name == "numuCC1_signal_Ehad_PC_overlay" || ch_name == "numuCC2_signal_Ehad_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Ehad_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Ehad_FC_overlay" || ch_name == "numuCC1_background_Ehad_FC_overlay" || ch_name == "numuCC2_background_Ehad_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Ehad_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Ehad_PC_overlay" || ch_name == "numuCC1_background_Ehad_PC_overlay" || ch_name == "numuCC2_background_Ehad_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Ehad_numuCCinFV"])) return true;
    }
    return false;

    
  }else if (ch_name == "numuCC_FC_bnb_L800MeV" || ch_name == "BG_numuCC_FC_ext_L800MeV" || ch_name == "BG_numuCC_FC_dirt_L800MeV"
	    || ch_name == "numuCC1_FC_bnb_L800MeV" || ch_name == "BG_numuCC1_FC_ext_L800MeV" || ch_name == "BG_numuCC1_FC_dirt_L800MeV"
	    || ch_name == "numuCC2_FC_bnb_L800MeV" || ch_name == "BG_numuCC2_FC_ext_L800MeV" || ch_name == "BG_numuCC2_FC_dirt_L800MeV"
	    ){
    if (flag_numuCC && flag_FC && reco_Enu<800) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb_L800MeV" || ch_name == "BG_numuCC_PC_ext_L800MeV" || ch_name == "BG_numuCC_PC_dirt_L800MeV"
	    || ch_name == "numuCC1_PC_bnb_L800MeV" || ch_name == "BG_numuCC1_PC_ext_L800MeV" || ch_name == "BG_numuCC1_PC_dirt_L800MeV"
	    || ch_name == "numuCC2_PC_bnb_L800MeV" || ch_name == "BG_numuCC2_PC_ext_L800MeV" || ch_name == "BG_numuCC2_PC_dirt_L800MeV"
	    ){
    if (flag_numuCC && (!flag_FC) && reco_Enu<800) return true;
    else return false;
  }else if (ch_name == "numuCC_FC_overlay_L800MeV" || ch_name == "numuCC_PC_overlay_L800MeV" 
	    || ch_name == "numuCC1_FC_overlay_L800MeV" || ch_name == "numuCC1_PC_overlay_L800MeV" 
	    || ch_name == "numuCC2_FC_overlay_L800MeV" || ch_name == "numuCC2_PC_overlay_L800MeV"   ){
    if (ch_name == "numuCC_FC_overlay_L800MeV" || ch_name == "numuCC1_FC_overlay_L800MeV" || ch_name == "numuCC2_FC_overlay_L800MeV"){
      if (flag_numuCC && flag_FC && reco_Enu<800) return true;
    }else if (ch_name == "numuCC_PC_overlay_L800MeV" || ch_name == "numuCC1_PC_overlay_L800MeV" || ch_name == "numuCC2_PC_overlay_L800MeV" ){
      if (flag_numuCC && (!flag_FC) && reco_Enu<800) return true;
    }
    return false;
  }else if (ch_name == "numuCC_FC_overlay" || ch_name == "numuCC_PC_overlay" 
	    || ch_name == "numuCC1_FC_overlay" || ch_name == "numuCC1_PC_overlay" 
	    || ch_name == "numuCC2_FC_overlay" || ch_name == "numuCC2_PC_overlay"   ){
    if (ch_name == "numuCC_FC_overlay" || ch_name == "numuCC1_FC_overlay" || ch_name == "numuCC2_FC_overlay"){
      if (flag_numuCC && flag_FC ) return true;
    }else if (ch_name == "numuCC_PC_overlay" || ch_name == "numuCC1_PC_overlay" || ch_name == "numuCC2_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) ) return true;
    }
    return false;

    }else if (ch_name == "nc_pio_energy_FC" || ch_name == "nc_pio_score_FC"
	    || ch_name == "nc_pio_energy_FC_ncpio_overlay" || ch_name == "nc_pio_score_FC_ncpio_overlay"
	    || ch_name == "nc_pio_energy_FC_ncdelta_overlay" || ch_name == "nc_pio_score_FC_ncdelta_overlay"
	    || ch_name == "nc_pio_energy_FC_overlay" || ch_name == "nc_pio_score_FC_overlay"
	    || ch_name == "nc_pio_energy_FC_ext" || ch_name == "nc_pio_score_FC_ext"
	    || ch_name == "nc_pio_energy_FC_dirt" || ch_name == "nc_pio_score_FC_dirt"
	    ){
    if (ch_name == "nc_pio_energy_FC" 
	|| ch_name == "nc_pio_energy_FC_ext" 
	|| ch_name == "nc_pio_energy_FC_dirt" ){
	    if (flag_ncpio_bdt && flag_FC) return true;
    }else if (ch_name == "nc_pio_energy_FC_ncpio_overlay" ){
	    if (flag_ncpio_bdt && flag_FC && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0 
				    && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true; 
    }else if (ch_name == "nc_pio_energy_FC_ncdelta_overlay" ){
	    if (flag_ncpio_bdt && flag_FC && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_pio_energy_FC_overlay" ){
	    if (flag_ncpio_bdt && flag_FC && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
					    && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)))
			    && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
	    //if (flag_ncpio_bdt && flag_FC) return true;
    }else if (ch_name == "nc_pio_score_FC"
		    || ch_name == "nc_pio_score_FC_ext"
		    || ch_name == "nc_pio_score_FC_dirt"){
	    if (flag_FC) return true;
    }else if (ch_name == "nc_pio_score_FC_ncpio_overlay"){ 
	    if (flag_FC && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0 
				    && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true; 
    }else if (ch_name == "nc_pio_score_FC_ncdelta_overlay"){
	    if (flag_FC && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_pio_score_FC_overlay"){
	    if (flag_FC && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
					    && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
			    && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
	    //  if (flag_FC) return true;
    }

    return false;


    // added by lhagaman, 2021_07_19

    }else if (ch_name == "fc_0p_score_excess_bin" || ch_name == "fc_0p_score_excess_bin_overlay" || ch_name == "fc_0p_score_excess_bin_ext" || ch_name == "fc_0p_score_excess_bin_dirt") {
	    if (flag_FC && flag_0p && (-11.39 < tagger.nc_delta_score) && (tagger.nc_delta_score < -9.39)) return true;
    }else if (ch_name == "fc_0p_cut_excess_bin" || ch_name == "fc_0p_cut_excess_bin_overlay" || ch_name == "fc_0p_cut_excess_bin_ext" || ch_name == "fc_0p_cut_excess_bin_dirt") {

	    int Nothertracks = 0;
	    for(size_t i=0; i<kine.kine_energy_particle->size(); i++){
		    int pdgcode = kine.kine_particle_type->at(i);
		    if((abs(pdgcode)==211 || abs(pdgcode)==13) && kine.kine_energy_particle->at(i)>10) Nothertracks++; // KE threshold: 10 MeV
	    }
	    //std::cout << "num other tracks: " << Nothertracks << "0p: " << flag_0p << "numu_score: " << tagger.numu_score << "kine_reco_Enu: " << get_reco_Enu_corr(kine, flag_data) << std::endl;

	    if (flag_FC && flag_0p && (0.9 < tagger.numu_score) && (Nothertracks == 1) && (300. < get_reco_Enu_corr(kine, flag_data)) && (get_reco_Enu_corr(kine, flag_data) < 900.)) return true;
	    // first part is getting the names for the 0p_Np cov_input.txt file, second part is adding a bunch of nc delta names in case you want to plot a bunch of things at once	
    }else if (ch_name == "numuCC_overlay_fc_0p" || ch_name == "numuCC_bnb_fc_0p" || ch_name == "BG_numuCC_ext_fc_0p" || ch_name == "BG_numuCC_dirt_fc_0p") {
	    if (flag_FC && flag_0p && tagger.numu_score > 0.9) return true;
    
    }else if (ch_name == "nc_delta_0p_01_true_Np_nc_delta_overlay" || ch_name == "nc_delta_0p_01_true_0p_nc_delta_overlay"){
            int nTrueP = 0;
	    if (flag_FC && flag_ncdelta_bdt && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
    		for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      			if(pfeval.truth_pdg[i] != 2212) continue;
      			if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
      			nTrueP+=1;
		}
		return ((nTrueP==0 && ch_name=="nc_delta_0p_01_true_0p_nc_delta_overlay") || (nTrueP>0 && ch_name=="nc_delta_0p_01_true_Np_nc_delta_overlay"));
	    }
	    return false;
   
    }else if (ch_name == "nc_delta_Np_01_true_0p_nc_delta_overlay" || ch_name == "nc_delta_Np_01_true_Np_nc_delta_overlay"){
            int nTrueP = 0;
            if (flag_FC && flag_ncdelta_bdt && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                for(size_t i=0; i<pfeval.truth_Ntrack; i++){
                        if(pfeval.truth_pdg[i] != 2212) continue;
                        if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
                        nTrueP+=1;
                }
                return ((nTrueP==0 && ch_name=="nc_delta_Np_01_true_0p_nc_delta_overlay") || (nTrueP>0 && ch_name=="nc_delta_Np_01_true_Np_nc_delta_overlay"));
            }
            return false;


    }else if (ch_name == "nc_delta_0p_01" || ch_name == "nc_delta_0p_02" || ch_name == "nc_delta_0p_03" || ch_name == "nc_delta_0p_04"
		    || ch_name == "nc_delta_0p_05" || ch_name == "nc_delta_0p_06" || ch_name == "nc_delta_0p_07" || ch_name == "nc_delta_0p_08"
		    || ch_name == "nc_delta_0p_09" || ch_name == "nc_delta_0p_10" || ch_name == "nc_delta_0p_11" || ch_name == "nc_delta_0p_12"
		    || ch_name == "nc_delta_0p_13" || ch_name == "nc_delta_0p_14" || ch_name == "nc_delta_0p_15" || ch_name == "nc_delta_0p_16"
		    || ch_name == "nc_delta_0p_17" || ch_name == "nc_delta_0p_18" || ch_name == "nc_delta_0p_19" || ch_name == "nc_delta_0p_20"){
	    if (flag_FC && flag_ncdelta_bdt && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01" || ch_name == "nc_delta_Np_02" || ch_name == "nc_delta_Np_03" || ch_name == "nc_delta_Np_04"
            || ch_name == "nc_delta_Np_05" || ch_name == "nc_delta_Np_06" || ch_name == "nc_delta_Np_07" || ch_name == "nc_delta_Np_08"
            || ch_name == "nc_delta_Np_09" || ch_name == "nc_delta_Np_10" || ch_name == "nc_delta_Np_11" || ch_name == "nc_delta_Np_12"
            || ch_name == "nc_delta_Np_13" || ch_name == "nc_delta_Np_14" || ch_name == "nc_delta_Np_15" || ch_name == "nc_delta_Np_16"
            || ch_name == "nc_delta_Np_17" || ch_name == "nc_delta_Np_18" || ch_name == "nc_delta_Np_19" || ch_name == "nc_delta_Np_20"){
                  if (flag_FC && flag_ncdelta_bdt && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_ext" || ch_name == "nc_delta_0p_02_ext" || ch_name == "nc_delta_0p_03_ext" || ch_name == "nc_delta_0p_04_ext"
            || ch_name == "nc_delta_0p_05_ext" || ch_name == "nc_delta_0p_06_ext" || ch_name == "nc_delta_0p_07_ext" || ch_name == "nc_delta_0p_08_ext"
            || ch_name == "nc_delta_0p_09_ext" || ch_name == "nc_delta_0p_10_ext" || ch_name == "nc_delta_0p_11_ext" || ch_name == "nc_delta_0p_12_ext"
            || ch_name == "nc_delta_0p_13_ext" || ch_name == "nc_delta_0p_14_ext" || ch_name == "nc_delta_0p_15_ext" || ch_name == "nc_delta_0p_16_ext"
            || ch_name == "nc_delta_0p_17_ext" || ch_name == "nc_delta_0p_18_ext" || ch_name == "nc_delta_0p_19_ext" || ch_name == "nc_delta_0p_20_ext"){
                  if (flag_FC && flag_ncdelta_bdt && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_ext" || ch_name == "nc_delta_Np_02_ext" || ch_name == "nc_delta_Np_03_ext" || ch_name == "nc_delta_Np_04_ext"
            || ch_name == "nc_delta_Np_05_ext" || ch_name == "nc_delta_Np_06_ext" || ch_name == "nc_delta_Np_07_ext" || ch_name == "nc_delta_Np_08_ext"
            || ch_name == "nc_delta_Np_09_ext" || ch_name == "nc_delta_Np_10_ext" || ch_name == "nc_delta_Np_11_ext" || ch_name == "nc_delta_Np_12_ext"
            || ch_name == "nc_delta_Np_13_ext" || ch_name == "nc_delta_Np_14_ext" || ch_name == "nc_delta_Np_15_ext" || ch_name == "nc_delta_Np_16_ext"
            || ch_name == "nc_delta_Np_17_ext" || ch_name == "nc_delta_Np_18_ext" || ch_name == "nc_delta_Np_19_ext" || ch_name == "nc_delta_Np_20_ext"){
                  if (flag_FC && flag_ncdelta_bdt && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_dirt" || ch_name == "nc_delta_0p_02_dirt" || ch_name == "nc_delta_0p_03_dirt" || ch_name == "nc_delta_0p_04_dirt"
            || ch_name == "nc_delta_0p_05_dirt" || ch_name == "nc_delta_0p_06_dirt" || ch_name == "nc_delta_0p_07_dirt" || ch_name == "nc_delta_0p_08_dirt"
            || ch_name == "nc_delta_0p_09_dirt" || ch_name == "nc_delta_0p_10_dirt" || ch_name == "nc_delta_0p_11_dirt" || ch_name == "nc_delta_0p_12_dirt"
            || ch_name == "nc_delta_0p_13_dirt" || ch_name == "nc_delta_0p_14_dirt" || ch_name == "nc_delta_0p_15_dirt" || ch_name == "nc_delta_0p_16_dirt"
            || ch_name == "nc_delta_0p_17_dirt" || ch_name == "nc_delta_0p_18_dirt" || ch_name == "nc_delta_0p_19_dirt" || ch_name == "nc_delta_0p_20_dirt"){
                  if (flag_FC && flag_ncdelta_bdt && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_dirt" || ch_name == "nc_delta_Np_02_dirt" || ch_name == "nc_delta_Np_03_dirt" || ch_name == "nc_delta_Np_04_dirt"
            || ch_name == "nc_delta_Np_05_dirt" || ch_name == "nc_delta_Np_06_dirt" || ch_name == "nc_delta_Np_07_dirt" || ch_name == "nc_delta_Np_08_dirt"
            || ch_name == "nc_delta_Np_09_dirt" || ch_name == "nc_delta_Np_10_dirt" || ch_name == "nc_delta_Np_11_dirt" || ch_name == "nc_delta_Np_12_dirt"
            || ch_name == "nc_delta_Np_13_dirt" || ch_name == "nc_delta_Np_14_dirt" || ch_name == "nc_delta_Np_15_dirt" || ch_name == "nc_delta_Np_16_dirt"
            || ch_name == "nc_delta_Np_17_dirt" || ch_name == "nc_delta_Np_18_dirt" || ch_name == "nc_delta_Np_19_dirt" || ch_name == "nc_delta_Np_20_dirt"){
                  if (flag_FC && flag_ncdelta_bdt && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_nc_delta_overlay" || ch_name == "nc_delta_0p_01_nc_delta_overlay_add" || ch_name == "nc_delta_0p_02_nc_delta_overlay" || ch_name == "nc_delta_0p_03_nc_delta_overlay" || ch_name == "nc_delta_0p_04_nc_delta_overlay"
            || ch_name == "nc_delta_0p_05_nc_delta_overlay" || ch_name == "nc_delta_0p_06_nc_delta_overlay" || ch_name == "nc_delta_0p_07_nc_delta_overlay" || ch_name == "nc_delta_0p_08_nc_delta_overlay"
            || ch_name == "nc_delta_0p_09_nc_delta_overlay" || ch_name == "nc_delta_0p_10_nc_delta_overlay" || ch_name == "nc_delta_0p_11_nc_delta_overlay" || ch_name == "nc_delta_0p_12_nc_delta_overlay"
            || ch_name == "nc_delta_0p_13_nc_delta_overlay" || ch_name == "nc_delta_0p_14_nc_delta_overlay" || ch_name == "nc_delta_0p_15_nc_delta_overlay" || ch_name == "nc_delta_0p_16_nc_delta_overlay"
            || ch_name == "nc_delta_0p_17_nc_delta_overlay" || ch_name == "nc_delta_0p_18_nc_delta_overlay" || ch_name == "nc_delta_0p_19_nc_delta_overlay" || ch_name == "nc_delta_0p_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_nc_delta_overlay" || ch_name == "nc_delta_Np_01_nc_delta_overlay_add" || ch_name == "nc_delta_Np_02_nc_delta_overlay" || ch_name == "nc_delta_Np_03_nc_delta_overlay" || ch_name == "nc_delta_Np_04_nc_delta_overlay"
            || ch_name == "nc_delta_Np_05_nc_delta_overlay" || ch_name == "nc_delta_Np_06_nc_delta_overlay" || ch_name == "nc_delta_Np_07_nc_delta_overlay" || ch_name == "nc_delta_Np_08_nc_delta_overlay"
            || ch_name == "nc_delta_Np_09_nc_delta_overlay" || ch_name == "nc_delta_Np_10_nc_delta_overlay" || ch_name == "nc_delta_Np_11_nc_delta_overlay" || ch_name == "nc_delta_Np_12_nc_delta_overlay"
            || ch_name == "nc_delta_Np_13_nc_delta_overlay" || ch_name == "nc_delta_Np_14_nc_delta_overlay" || ch_name == "nc_delta_Np_15_nc_delta_overlay" || ch_name == "nc_delta_Np_16_nc_delta_overlay"
            || ch_name == "nc_delta_Np_17_nc_delta_overlay" || ch_name == "nc_delta_Np_18_nc_delta_overlay" || ch_name == "nc_delta_Np_19_nc_delta_overlay" || ch_name == "nc_delta_Np_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_nc_pi0_overlay" || ch_name == "nc_delta_0p_02_nc_pi0_overlay" || ch_name == "nc_delta_0p_03_nc_pi0_overlay" || ch_name == "nc_delta_0p_04_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_05_nc_pi0_overlay" || ch_name == "nc_delta_0p_06_nc_pi0_overlay" || ch_name == "nc_delta_0p_07_nc_pi0_overlay" || ch_name == "nc_delta_0p_08_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_09_nc_pi0_overlay" || ch_name == "nc_delta_0p_10_nc_pi0_overlay" || ch_name == "nc_delta_0p_11_nc_pi0_overlay" || ch_name == "nc_delta_0p_12_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_13_nc_pi0_overlay" || ch_name == "nc_delta_0p_14_nc_pi0_overlay" || ch_name == "nc_delta_0p_15_nc_pi0_overlay" || ch_name == "nc_delta_0p_16_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_17_nc_pi0_overlay" || ch_name == "nc_delta_0p_18_nc_pi0_overlay" || ch_name == "nc_delta_0p_19_nc_pi0_overlay" || ch_name == "nc_delta_0p_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_nc_pi0_overlay" || ch_name == "nc_delta_Np_02_nc_pi0_overlay" || ch_name == "nc_delta_Np_03_nc_pi0_overlay" || ch_name == "nc_delta_Np_04_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_05_nc_pi0_overlay" || ch_name == "nc_delta_Np_06_nc_pi0_overlay" || ch_name == "nc_delta_Np_07_nc_pi0_overlay" || ch_name == "nc_delta_Np_08_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_09_nc_pi0_overlay" || ch_name == "nc_delta_Np_10_nc_pi0_overlay" || ch_name == "nc_delta_Np_11_nc_pi0_overlay" || ch_name == "nc_delta_Np_12_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_13_nc_pi0_overlay" || ch_name == "nc_delta_Np_14_nc_pi0_overlay" || ch_name == "nc_delta_Np_15_nc_pi0_overlay" || ch_name == "nc_delta_Np_16_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_17_nc_pi0_overlay" || ch_name == "nc_delta_Np_18_nc_pi0_overlay" || ch_name == "nc_delta_Np_19_nc_pi0_overlay" || ch_name == "nc_delta_Np_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                        && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_overlay" || ch_name == "nc_delta_0p_02_overlay" || ch_name == "nc_delta_0p_03_overlay" || ch_name == "nc_delta_0p_04_overlay"
            || ch_name == "nc_delta_0p_05_overlay" || ch_name == "nc_delta_0p_06_overlay" || ch_name == "nc_delta_0p_07_overlay" || ch_name == "nc_delta_0p_08_overlay"
            || ch_name == "nc_delta_0p_09_overlay" || ch_name == "nc_delta_0p_10_overlay" || ch_name == "nc_delta_0p_11_overlay" || ch_name == "nc_delta_0p_12_overlay"
            || ch_name == "nc_delta_0p_13_overlay" || ch_name == "nc_delta_0p_14_overlay" || ch_name == "nc_delta_0p_15_overlay" || ch_name == "nc_delta_0p_16_overlay"
            || ch_name == "nc_delta_0p_17_overlay" || ch_name == "nc_delta_0p_18_overlay" || ch_name == "nc_delta_0p_19_overlay" || ch_name == "nc_delta_0p_20_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                               && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_overlay" || ch_name == "nc_delta_Np_02_overlay" || ch_name == "nc_delta_Np_03_overlay" || ch_name == "nc_delta_Np_04_overlay"
            || ch_name == "nc_delta_Np_05_overlay" || ch_name == "nc_delta_Np_06_overlay" || ch_name == "nc_delta_Np_07_overlay" || ch_name == "nc_delta_Np_08_overlay"
            || ch_name == "nc_delta_Np_09_overlay" || ch_name == "nc_delta_Np_10_overlay" || ch_name == "nc_delta_Np_11_overlay" || ch_name == "nc_delta_Np_12_overlay"
            || ch_name == "nc_delta_Np_13_overlay" || ch_name == "nc_delta_Np_14_overlay" || ch_name == "nc_delta_Np_15_overlay" || ch_name == "nc_delta_Np_16_overlay"
            || ch_name == "nc_delta_Np_17_overlay" || ch_name == "nc_delta_Np_18_overlay" || ch_name == "nc_delta_Np_19_overlay" || ch_name == "nc_delta_Np_20_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false; // start added 2022_04_20
    }else if (ch_name == "near_bdt_0p_01" || ch_name == "near_bdt_0p_02" || ch_name == "near_bdt_0p_03" || ch_name == "near_bdt_0p_04"
                    || ch_name == "near_bdt_0p_05" || ch_name == "near_bdt_0p_06" || ch_name == "near_bdt_0p_07" || ch_name == "near_bdt_0p_08"
                    || ch_name == "near_bdt_0p_09" || ch_name == "near_bdt_0p_10" || ch_name == "near_bdt_0p_11" || ch_name == "near_bdt_0p_12"
                    || ch_name == "near_bdt_0p_13" || ch_name == "near_bdt_0p_14" || ch_name == "near_bdt_0p_15" || ch_name == "near_bdt_0p_16"
                    || ch_name == "near_bdt_0p_17" || ch_name == "near_bdt_0p_18" || ch_name == "near_bdt_0p_19" || ch_name == "near_bdt_0p_20"){
            if (flag_FC && flag_near_bdt && flag_0p) return true;
                  return false;
    }else if (ch_name == "near_bdt_Np_01" || ch_name == "near_bdt_Np_02" || ch_name == "near_bdt_Np_03" || ch_name == "near_bdt_Np_04"
            || ch_name == "near_bdt_Np_05" || ch_name == "near_bdt_Np_06" || ch_name == "near_bdt_Np_07" || ch_name == "near_bdt_Np_08"
            || ch_name == "near_bdt_Np_09" || ch_name == "near_bdt_Np_10" || ch_name == "near_bdt_Np_11" || ch_name == "near_bdt_Np_12"
            || ch_name == "near_bdt_Np_13" || ch_name == "near_bdt_Np_14" || ch_name == "near_bdt_Np_15" || ch_name == "near_bdt_Np_16"
            || ch_name == "near_bdt_Np_17" || ch_name == "near_bdt_Np_18" || ch_name == "near_bdt_Np_19" || ch_name == "near_bdt_Np_20"){
                  if (flag_FC && flag_near_bdt && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "near_bdt_0p_01_ext" || ch_name == "near_bdt_0p_02_ext" || ch_name == "near_bdt_0p_03_ext" || ch_name == "near_bdt_0p_04_ext"
            || ch_name == "near_bdt_0p_05_ext" || ch_name == "near_bdt_0p_06_ext" || ch_name == "near_bdt_0p_07_ext" || ch_name == "near_bdt_0p_08_ext"
            || ch_name == "near_bdt_0p_09_ext" || ch_name == "near_bdt_0p_10_ext" || ch_name == "near_bdt_0p_11_ext" || ch_name == "near_bdt_0p_12_ext"
            || ch_name == "near_bdt_0p_13_ext" || ch_name == "near_bdt_0p_14_ext" || ch_name == "near_bdt_0p_15_ext" || ch_name == "near_bdt_0p_16_ext"
            || ch_name == "near_bdt_0p_17_ext" || ch_name == "near_bdt_0p_18_ext" || ch_name == "near_bdt_0p_19_ext" || ch_name == "near_bdt_0p_20_ext"){
                  if (flag_FC && flag_near_bdt && flag_0p) return true;
                  return false;
    }else if (ch_name == "near_bdt_Np_01_ext" || ch_name == "near_bdt_Np_02_ext" || ch_name == "near_bdt_Np_03_ext" || ch_name == "near_bdt_Np_04_ext"
            || ch_name == "near_bdt_Np_05_ext" || ch_name == "near_bdt_Np_06_ext" || ch_name == "near_bdt_Np_07_ext" || ch_name == "near_bdt_Np_08_ext"
            || ch_name == "near_bdt_Np_09_ext" || ch_name == "near_bdt_Np_10_ext" || ch_name == "near_bdt_Np_11_ext" || ch_name == "near_bdt_Np_12_ext"
            || ch_name == "near_bdt_Np_13_ext" || ch_name == "near_bdt_Np_14_ext" || ch_name == "near_bdt_Np_15_ext" || ch_name == "near_bdt_Np_16_ext"
            || ch_name == "near_bdt_Np_17_ext" || ch_name == "near_bdt_Np_18_ext" || ch_name == "near_bdt_Np_19_ext" || ch_name == "near_bdt_Np_20_ext"){
                  if (flag_FC && flag_near_bdt && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "near_bdt_0p_01_dirt" || ch_name == "near_bdt_0p_02_dirt" || ch_name == "near_bdt_0p_03_dirt" || ch_name == "near_bdt_0p_04_dirt"
            || ch_name == "near_bdt_0p_05_dirt" || ch_name == "near_bdt_0p_06_dirt" || ch_name == "near_bdt_0p_07_dirt" || ch_name == "near_bdt_0p_08_dirt"
            || ch_name == "near_bdt_0p_09_dirt" || ch_name == "near_bdt_0p_10_dirt" || ch_name == "near_bdt_0p_11_dirt" || ch_name == "near_bdt_0p_12_dirt"
            || ch_name == "near_bdt_0p_13_dirt" || ch_name == "near_bdt_0p_14_dirt" || ch_name == "near_bdt_0p_15_dirt" || ch_name == "near_bdt_0p_16_dirt"
            || ch_name == "near_bdt_0p_17_dirt" || ch_name == "near_bdt_0p_18_dirt" || ch_name == "near_bdt_0p_19_dirt" || ch_name == "near_bdt_0p_20_dirt"){
                  if (flag_FC && flag_near_bdt && flag_0p) return true;
                  return false;
    }else if (ch_name == "near_bdt_Np_01_dirt" || ch_name == "near_bdt_Np_02_dirt" || ch_name == "near_bdt_Np_03_dirt" || ch_name == "near_bdt_Np_04_dirt"
            || ch_name == "near_bdt_Np_05_dirt" || ch_name == "near_bdt_Np_06_dirt" || ch_name == "near_bdt_Np_07_dirt" || ch_name == "near_bdt_Np_08_dirt"
            || ch_name == "near_bdt_Np_09_dirt" || ch_name == "near_bdt_Np_10_dirt" || ch_name == "near_bdt_Np_11_dirt" || ch_name == "near_bdt_Np_12_dirt"
            || ch_name == "near_bdt_Np_13_dirt" || ch_name == "near_bdt_Np_14_dirt" || ch_name == "near_bdt_Np_15_dirt" || ch_name == "near_bdt_Np_16_dirt"
            || ch_name == "near_bdt_Np_17_dirt" || ch_name == "near_bdt_Np_18_dirt" || ch_name == "near_bdt_Np_19_dirt" || ch_name == "near_bdt_Np_20_dirt"){
                  if (flag_FC && flag_near_bdt && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "near_bdt_0p_01_nc_delta_overlay" || ch_name == "near_bdt_0p_01_nc_delta_overlay_add" || ch_name == "near_bdt_0p_02_nc_delta_overlay" || ch_name == "near_bdt_0p_03_nc_delta_overlay" || ch_name == "near_bdt_0p_04_nc_delta_overlay"
            || ch_name == "near_bdt_0p_05_nc_delta_overlay" || ch_name == "near_bdt_0p_06_nc_delta_overlay" || ch_name == "near_bdt_0p_07_nc_delta_overlay" || ch_name == "near_bdt_0p_08_nc_delta_overlay"
            || ch_name == "near_bdt_0p_09_nc_delta_overlay" || ch_name == "near_bdt_0p_10_nc_delta_overlay" || ch_name == "near_bdt_0p_11_nc_delta_overlay" || ch_name == "near_bdt_0p_12_nc_delta_overlay"
            || ch_name == "near_bdt_0p_13_nc_delta_overlay" || ch_name == "near_bdt_0p_14_nc_delta_overlay" || ch_name == "near_bdt_0p_15_nc_delta_overlay" || ch_name == "near_bdt_0p_16_nc_delta_overlay"
            || ch_name == "near_bdt_0p_17_nc_delta_overlay" || ch_name == "near_bdt_0p_18_nc_delta_overlay" || ch_name == "near_bdt_0p_19_nc_delta_overlay" || ch_name == "near_bdt_0p_20_nc_delta_overlay"){
                  if (flag_FC && flag_near_bdt && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "near_bdt_Np_01_nc_delta_overlay" || ch_name == "near_bdt_Np_01_nc_delta_overlay_add" || ch_name == "near_bdt_Np_02_nc_delta_overlay" || ch_name == "near_bdt_Np_03_nc_delta_overlay" || ch_name == "near_bdt_Np_04_nc_delta_overlay"
            || ch_name == "near_bdt_Np_05_nc_delta_overlay" || ch_name == "near_bdt_Np_06_nc_delta_overlay" || ch_name == "near_bdt_Np_07_nc_delta_overlay" || ch_name == "near_bdt_Np_08_nc_delta_overlay"
            || ch_name == "near_bdt_Np_09_nc_delta_overlay" || ch_name == "near_bdt_Np_10_nc_delta_overlay" || ch_name == "near_bdt_Np_11_nc_delta_overlay" || ch_name == "near_bdt_Np_12_nc_delta_overlay"
            || ch_name == "near_bdt_Np_13_nc_delta_overlay" || ch_name == "near_bdt_Np_14_nc_delta_overlay" || ch_name == "near_bdt_Np_15_nc_delta_overlay" || ch_name == "near_bdt_Np_16_nc_delta_overlay"
            || ch_name == "near_bdt_Np_17_nc_delta_overlay" || ch_name == "near_bdt_Np_18_nc_delta_overlay" || ch_name == "near_bdt_Np_19_nc_delta_overlay" || ch_name == "near_bdt_Np_20_nc_delta_overlay"){
                  if (flag_FC && flag_near_bdt && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "near_bdt_0p_01_nc_pi0_overlay" || ch_name == "near_bdt_0p_02_nc_pi0_overlay" || ch_name == "near_bdt_0p_03_nc_pi0_overlay" || ch_name == "near_bdt_0p_04_nc_pi0_overlay"
            || ch_name == "near_bdt_0p_05_nc_pi0_overlay" || ch_name == "near_bdt_0p_06_nc_pi0_overlay" || ch_name == "near_bdt_0p_07_nc_pi0_overlay" || ch_name == "near_bdt_0p_08_nc_pi0_overlay"
            || ch_name == "near_bdt_0p_09_nc_pi0_overlay" || ch_name == "near_bdt_0p_10_nc_pi0_overlay" || ch_name == "near_bdt_0p_11_nc_pi0_overlay" || ch_name == "near_bdt_0p_12_nc_pi0_overlay"
            || ch_name == "near_bdt_0p_13_nc_pi0_overlay" || ch_name == "near_bdt_0p_14_nc_pi0_overlay" || ch_name == "near_bdt_0p_15_nc_pi0_overlay" || ch_name == "near_bdt_0p_16_nc_pi0_overlay"
            || ch_name == "near_bdt_0p_17_nc_pi0_overlay" || ch_name == "near_bdt_0p_18_nc_pi0_overlay" || ch_name == "near_bdt_0p_19_nc_pi0_overlay" || ch_name == "near_bdt_0p_20_nc_pi0_overlay"){
                  if (flag_FC && flag_near_bdt && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "near_bdt_Np_01_nc_pi0_overlay" || ch_name == "near_bdt_Np_02_nc_pi0_overlay" || ch_name == "near_bdt_Np_03_nc_pi0_overlay" || ch_name == "near_bdt_Np_04_nc_pi0_overlay"
            || ch_name == "near_bdt_Np_05_nc_pi0_overlay" || ch_name == "near_bdt_Np_06_nc_pi0_overlay" || ch_name == "near_bdt_Np_07_nc_pi0_overlay" || ch_name == "near_bdt_Np_08_nc_pi0_overlay"
            || ch_name == "near_bdt_Np_09_nc_pi0_overlay" || ch_name == "near_bdt_Np_10_nc_pi0_overlay" || ch_name == "near_bdt_Np_11_nc_pi0_overlay" || ch_name == "near_bdt_Np_12_nc_pi0_overlay"
            || ch_name == "near_bdt_Np_13_nc_pi0_overlay" || ch_name == "near_bdt_Np_14_nc_pi0_overlay" || ch_name == "near_bdt_Np_15_nc_pi0_overlay" || ch_name == "near_bdt_Np_16_nc_pi0_overlay"
            || ch_name == "near_bdt_Np_17_nc_pi0_overlay" || ch_name == "near_bdt_Np_18_nc_pi0_overlay" || ch_name == "near_bdt_Np_19_nc_pi0_overlay" || ch_name == "near_bdt_Np_20_nc_pi0_overlay"){
                  if (flag_FC && flag_near_bdt && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                        && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "near_bdt_0p_01_overlay" || ch_name == "near_bdt_0p_02_overlay" || ch_name == "near_bdt_0p_03_overlay" || ch_name == "near_bdt_0p_04_overlay"
            || ch_name == "near_bdt_0p_05_overlay" || ch_name == "near_bdt_0p_06_overlay" || ch_name == "near_bdt_0p_07_overlay" || ch_name == "near_bdt_0p_08_overlay"
            || ch_name == "near_bdt_0p_09_overlay" || ch_name == "near_bdt_0p_10_overlay" || ch_name == "near_bdt_0p_11_overlay" || ch_name == "near_bdt_0p_12_overlay"
            || ch_name == "near_bdt_0p_13_overlay" || ch_name == "near_bdt_0p_14_overlay" || ch_name == "near_bdt_0p_15_overlay" || ch_name == "near_bdt_0p_16_overlay"
            || ch_name == "near_bdt_0p_17_overlay" || ch_name == "near_bdt_0p_18_overlay" || ch_name == "near_bdt_0p_19_overlay" || ch_name == "near_bdt_0p_20_overlay"){
                  if (flag_FC && flag_near_bdt && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                               && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "near_bdt_Np_01_overlay" || ch_name == "near_bdt_Np_02_overlay" || ch_name == "near_bdt_Np_03_overlay" || ch_name == "near_bdt_Np_04_overlay"
            || ch_name == "near_bdt_Np_05_overlay" || ch_name == "near_bdt_Np_06_overlay" || ch_name == "near_bdt_Np_07_overlay" || ch_name == "near_bdt_Np_08_overlay"
            || ch_name == "near_bdt_Np_09_overlay" || ch_name == "near_bdt_Np_10_overlay" || ch_name == "near_bdt_Np_11_overlay" || ch_name == "near_bdt_Np_12_overlay"
            || ch_name == "near_bdt_Np_13_overlay" || ch_name == "near_bdt_Np_14_overlay" || ch_name == "near_bdt_Np_15_overlay" || ch_name == "near_bdt_Np_16_overlay"
            || ch_name == "near_bdt_Np_17_overlay" || ch_name == "near_bdt_Np_18_overlay" || ch_name == "near_bdt_Np_19_overlay" || ch_name == "near_bdt_Np_20_overlay"){
                  if (flag_FC && flag_near_bdt && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false; // end added 2022_04_20
    }else if (ch_name == "nc_delta_0p_01_overlay_entire" || ch_name == "nc_delta_0p_02_overlay_entire" || ch_name == "nc_delta_0p_03_overlay_entire" || ch_name == "nc_delta_0p_04_overlay_entire"){
                  if (flag_FC && flag_ncdelta_bdt && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_overlay_entire" || ch_name == "nc_delta_Np_02_overlay_entire" || ch_name == "nc_delta_Np_03_overlay_entire" || ch_name == "nc_delta_Np_04_overlay_entire"){
                  if (flag_FC && flag_ncdelta_bdt && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01" || ch_name == "nc_delta_Xp_02" || ch_name == "nc_delta_Xp_03" || ch_name == "nc_delta_Xp_04"
            || ch_name == "nc_delta_Xp_05" || ch_name == "nc_delta_Xp_06" || ch_name == "nc_delta_Xp_07" || ch_name == "nc_delta_Xp_08"
            || ch_name == "nc_delta_Xp_09" || ch_name == "nc_delta_Xp_10" || ch_name == "nc_delta_Xp_11" || ch_name == "nc_delta_Xp_12"
            || ch_name == "nc_delta_Xp_13" || ch_name == "nc_delta_Xp_14" || ch_name == "nc_delta_Xp_15" || ch_name == "nc_delta_Xp_16"
            || ch_name == "nc_delta_Xp_17" || ch_name == "nc_delta_Xp_18" || ch_name == "nc_delta_Xp_19" || ch_name == "nc_delta_Xp_20"){
                  if (flag_FC && flag_ncdelta_bdt) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_ext" || ch_name == "nc_delta_Xp_02_ext" || ch_name == "nc_delta_Xp_03_ext" || ch_name == "nc_delta_Xp_04_ext"
            || ch_name == "nc_delta_Xp_05_ext" || ch_name == "nc_delta_Xp_06_ext" || ch_name == "nc_delta_Xp_07_ext" || ch_name == "nc_delta_Xp_08_ext"
            || ch_name == "nc_delta_Xp_09_ext" || ch_name == "nc_delta_Xp_10_ext" || ch_name == "nc_delta_Xp_11_ext" || ch_name == "nc_delta_Xp_12_ext"
            || ch_name == "nc_delta_Xp_13_ext" || ch_name == "nc_delta_Xp_14_ext" || ch_name == "nc_delta_Xp_15_ext" || ch_name == "nc_delta_Xp_16_ext"
            || ch_name == "nc_delta_Xp_17_ext" || ch_name == "nc_delta_Xp_18_ext" || ch_name == "nc_delta_Xp_19_ext" || ch_name == "nc_delta_Xp_20_ext"){
                  if (flag_FC && flag_ncdelta_bdt) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_dirt" || ch_name == "nc_delta_Xp_02_dirt" || ch_name == "nc_delta_Xp_03_dirt" || ch_name == "nc_delta_Xp_04_dirt"
            || ch_name == "nc_delta_Xp_05_dirt" || ch_name == "nc_delta_Xp_06_dirt" || ch_name == "nc_delta_Xp_07_dirt" || ch_name == "nc_delta_Xp_08_dirt"
            || ch_name == "nc_delta_Xp_09_dirt" || ch_name == "nc_delta_Xp_10_dirt" || ch_name == "nc_delta_Xp_11_dirt" || ch_name == "nc_delta_Xp_12_dirt"
            || ch_name == "nc_delta_Xp_13_dirt" || ch_name == "nc_delta_Xp_14_dirt" || ch_name == "nc_delta_Xp_15_dirt" || ch_name == "nc_delta_Xp_16_dirt"
            || ch_name == "nc_delta_Xp_17_dirt" || ch_name == "nc_delta_Xp_18_dirt" || ch_name == "nc_delta_Xp_19_dirt" || ch_name == "nc_delta_Xp_20_dirt"){
                  if (flag_FC && flag_ncdelta_bdt) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_nc_delta_overlay" || ch_name == "nc_delta_Xp_02_nc_delta_overlay" || ch_name == "nc_delta_Xp_03_nc_delta_overlay" || ch_name == "nc_delta_Xp_04_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_05_nc_delta_overlay" || ch_name == "nc_delta_Xp_06_nc_delta_overlay" || ch_name == "nc_delta_Xp_07_nc_delta_overlay" || ch_name == "nc_delta_Xp_08_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_09_nc_delta_overlay" || ch_name == "nc_delta_Xp_10_nc_delta_overlay" || ch_name == "nc_delta_Xp_11_nc_delta_overlay" || ch_name == "nc_delta_Xp_12_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_13_nc_delta_overlay" || ch_name == "nc_delta_Xp_14_nc_delta_overlay" || ch_name == "nc_delta_Xp_15_nc_delta_overlay" || ch_name == "nc_delta_Xp_16_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_17_nc_delta_overlay" || ch_name == "nc_delta_Xp_18_nc_delta_overlay" || ch_name == "nc_delta_Xp_19_nc_delta_overlay" || ch_name == "nc_delta_Xp_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_nc_pi0_overlay" || ch_name == "nc_delta_Xp_02_nc_pi0_overlay" || ch_name == "nc_delta_Xp_03_nc_pi0_overlay" || ch_name == "nc_delta_Xp_04_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_05_nc_pi0_overlay" || ch_name == "nc_delta_Xp_06_nc_pi0_overlay" || ch_name == "nc_delta_Xp_07_nc_pi0_overlay" || ch_name == "nc_delta_Xp_08_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_09_nc_pi0_overlay" || ch_name == "nc_delta_Xp_10_nc_pi0_overlay" || ch_name == "nc_delta_Xp_11_nc_pi0_overlay" || ch_name == "nc_delta_Xp_12_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_13_nc_pi0_overlay" || ch_name == "nc_delta_Xp_14_nc_pi0_overlay" || ch_name == "nc_delta_Xp_15_nc_pi0_overlay" || ch_name == "nc_delta_Xp_16_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_17_nc_pi0_overlay" || ch_name == "nc_delta_Xp_18_nc_pi0_overlay" || ch_name == "nc_delta_Xp_19_nc_pi0_overlay" || ch_name == "nc_delta_Xp_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_overlay" || ch_name == "nc_delta_Xp_02_overlay" || ch_name == "nc_delta_Xp_03_overlay" || ch_name == "nc_delta_Xp_04_overlay"
            || ch_name == "nc_delta_Xp_05_overlay" || ch_name == "nc_delta_Xp_06_overlay" || ch_name == "nc_delta_Xp_07_overlay" || ch_name == "nc_delta_Xp_08_overlay"
            || ch_name == "nc_delta_Xp_09_overlay" || ch_name == "nc_delta_Xp_10_overlay" || ch_name == "nc_delta_Xp_11_overlay" || ch_name == "nc_delta_Xp_12_overlay"
            || ch_name == "nc_delta_Xp_13_overlay" || ch_name == "nc_delta_Xp_14_overlay" || ch_name == "nc_delta_Xp_15_overlay" || ch_name == "nc_delta_Xp_16_overlay"
            || ch_name == "nc_delta_Xp_17_overlay" || ch_name == "nc_delta_Xp_18_overlay" || ch_name == "nc_delta_Xp_19_overlay" || ch_name == "nc_delta_Xp_20_overlay"){
                  if (flag_FC && flag_ncdelta_bdt && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p" || ch_name == "nc_pi0_2_0p"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np" || ch_name == "nc_pi0_2_Np"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp" || ch_name == "nc_pi0_2_Xp"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_ext" || ch_name == "nc_pi0_2_0p_ext"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_ext" || ch_name == "nc_pi0_2_Np_ext"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_ext" || ch_name == "nc_pi0_2_Xp_ext"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_dirt" || ch_name == "nc_pi0_2_0p_dirt"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_dirt" || ch_name == "nc_pi0_2_Np_dirt"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_dirt" || ch_name == "nc_pi0_2_Xp_dirt"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_nc_delta_overlay" || ch_name == "nc_pi0_0p_nc_delta_overlay_add" || ch_name == "nc_pi0_2_0p_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_nc_delta_overlay" || ch_name == "nc_pi0_Np_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Np_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_nc_delta_overlay" || ch_name == "nc_pi0_Xp_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Xp_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;

    }else if (ch_name == "nc_pi0_Np_true_Np_nc_delta_overlay" || ch_name == "nc_pi0_Np_true_0p_nc_delta_overlay"){
            int nTrueP = 0;
            if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                for(size_t i=0; i<pfeval.truth_Ntrack; i++){
                        if(pfeval.truth_pdg[i] != 2212) continue;
                        if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
                        nTrueP+=1;
                }
                return ((nTrueP==0 && ch_name=="nc_pi0_Np_true_0p_nc_delta_overlay") || (nTrueP>0 && ch_name=="nc_pi0_Np_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "nc_pi0_0p_true_Np_nc_delta_overlay" || ch_name == "nc_pi0_0p_true_0p_nc_delta_overlay"){
            int nTrueP = 0;
            if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                for(size_t i=0; i<pfeval.truth_Ntrack; i++){
                        if(pfeval.truth_pdg[i] != 2212) continue;
                        if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
                        nTrueP+=1;
                }
                return ((nTrueP==0 && ch_name=="nc_pi0_0p_true_0p_nc_delta_overlay") || (nTrueP>0 && ch_name=="nc_pi0_0p_true_Np_nc_delta_overlay"));
            }
            return false;



    }else if (ch_name == "nc_pi0_0p_nc_pi0_overlay" || ch_name == "nc_pi0_2_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_nc_pi0_overlay" || ch_name == "nc_pi0_2_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_nc_pi0_overlay" || ch_name == "nc_pi0_2_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_overlay" || ch_name == "nc_pi0_2_0p_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_overlay" || ch_name == "nc_pi0_2_Np_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_overlay" || ch_name == "nc_pi0_2_Xp_overlay"){
                  if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
   

    }else if (ch_name == "nc_pi0_0p_old" || ch_name == "nc_pi0_2_0p_old"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_old" || ch_name == "nc_pi0_2_Np_old"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_old" || ch_name == "nc_pi0_2_Xp_old"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_old_ext" || ch_name == "nc_pi0_2_0p_old_ext"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_old_ext" || ch_name == "nc_pi0_2_Np_old_ext"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_old_ext" || ch_name == "nc_pi0_2_Xp_old_ext"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_old_dirt" || ch_name == "nc_pi0_2_0p_old_dirt"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_old_dirt" || ch_name == "nc_pi0_2_Np_old_dirt"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_old_dirt" || ch_name == "nc_pi0_2_Xp_old_dirt"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_old_nc_delta_overlay" || ch_name == "nc_pi0_0p_old_nc_delta_overlay_add" || ch_name == "nc_pi0_2_0p_old_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_old_nc_delta_overlay" || ch_name == "nc_pi0_Np_old_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Np_old_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_old_nc_delta_overlay" || ch_name == "nc_pi0_Xp_old_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Xp_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_old_nc_pi0_overlay" || ch_name == "nc_pi0_2_0p_old_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_old_nc_pi0_overlay" || ch_name == "nc_pi0_2_Np_old_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_old_nc_pi0_overlay" || ch_name == "nc_pi0_2_Xp_old_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_old_overlay" || ch_name == "nc_pi0_2_0p_old_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_old_overlay" || ch_name == "nc_pi0_2_Np_old_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_old_overlay" || ch_name == "nc_pi0_2_Xp_old_overlay"){
                  if (flag_FC && flag_ncpio_bdt && (!flag_ncdelta_bdt) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;



    }else if (ch_name == "cc_pi0_0p" || ch_name == "cc_pi0_2_0p"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np" || ch_name == "cc_pi0_2_Np"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp" || ch_name == "cc_pi0_2_Xp"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_ext" || ch_name == "cc_pi0_2_0p_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_ext" || ch_name == "cc_pi0_2_Np_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_ext" || ch_name == "cc_pi0_2_Xp_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_dirt" || ch_name == "cc_pi0_2_0p_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_dirt" || ch_name == "cc_pi0_2_Np_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_dirt" || ch_name == "cc_pi0_2_Xp_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_nc_delta_overlay" || ch_name == "cc_pi0_0p_nc_delta_overlay_add" || ch_name == "cc_pi0_2_0p_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_nc_delta_overlay" || ch_name == "cc_pi0_Np_nc_delta_overlay_add" || ch_name == "cc_pi0_2_Np_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_nc_delta_overlay" || ch_name == "cc_pi0_Xp_nc_delta_overlay_add" || ch_name == "cc_pi0_2_Xp_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_nc_pi0_overlay" || ch_name == "cc_pi0_2_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_nc_pi0_overlay" || ch_name == "cc_pi0_2_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_nc_pi0_overlay" || ch_name == "cc_pi0_2_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_overlay" || ch_name == "cc_pi0_2_0p_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_overlay" || ch_name == "cc_pi0_2_Np_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_overlay" || ch_name == "cc_pi0_2_Xp_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_0p"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0)  && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_nc_delta_overlay" || ch_name == "numuCC_0p_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_nc_delta_overlay" || ch_name == "numuCC_Np_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_nc_delta_overlay" || ch_name == "numuCC_Xp_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0)  && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    
    // "old" here means that CC Pi0 events haven't been removed. 
    }else if (ch_name == "numuCC_0p_old"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_old"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_old"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_old_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_old_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_old_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_old_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_old_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_old_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_old_nc_delta_overlay" || ch_name == "numuCC_0p_old_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_old_nc_delta_overlay" || ch_name == "numuCC_Np_old_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_old_nc_delta_overlay" || ch_name == "numuCC_Xp_old_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;

    }else if (ch_name == "numuCC_Np_old_true_Np_nc_delta_overlay" || ch_name == "numuCC_Np_old_true_0p_nc_delta_overlay"){
            int nTrueP = 0;
            if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                for(size_t i=0; i<pfeval.truth_Ntrack; i++){
                        if(pfeval.truth_pdg[i] != 2212) continue;
                        if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
                        nTrueP+=1;
                }
                return ((nTrueP==0 && ch_name=="numuCC_Np_old_true_0p_nc_delta_overlay") || (nTrueP>0 && ch_name=="numuCC_Np_old_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "numuCC_0p_old_true_Np_nc_delta_overlay" || ch_name == "numuCC_0p_old_true_0p_nc_delta_overlay"){
            int nTrueP = 0;
            if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                for(size_t i=0; i<pfeval.truth_Ntrack; i++){
                        if(pfeval.truth_pdg[i] != 2212) continue;
                        if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
                        nTrueP+=1;
                }
                return ((nTrueP==0 && ch_name=="numuCC_0p_old_true_0p_nc_delta_overlay") || (nTrueP>0 && ch_name=="numuCC_0p_old_true_Np_nc_delta_overlay"));
            }
            return false;




    }else if (ch_name == "numuCC_0p_old_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_old_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_old_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_old_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_old_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_old_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt) && (!flag_ncdelta_bdt) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
 
    }else if (ch_name == "numuCC_0p_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && (!flag_cc_pi0)) return true;
                  return false;

    }else if (ch_name == "nc_pi0_1p" || ch_name == "nc_pi0_1p_ext" || ch_name == "nc_pi0_1p_dirt" || ch_name == "nc_pi0_1p_overlay_entire" || ch_name == "nc_pi0_2_1p" || ch_name == "nc_pi0_2_1p_ext" || ch_name == "nc_pi0_2_1p_dirt" || ch_name == "nc_pi0_2_1p_overlay_entire"){
	    if (flag_FC && flag_ncpio_bdt_sel && (!flag_ncdelta_bdt) && flag_1p) return true;
                  return false;

    }else if (ch_name == "cc_pi0_1p" || ch_name == "cc_pi0_1p_ext" || ch_name == "cc_pi0_1p_dirt" || ch_name == "cc_pi0_1p_overlay_entire" || ch_name == "cc_pi0_2_1p" || ch_name == "cc_pi0_2_1p_ext" || ch_name == "cc_pi0_2_1p_dirt" || ch_name == "cc_pi0_2_1p_overlay_entire"){
            if (flag_FC && flag_cc_pi0 && (!flag_ncpio_bdt_sel) && (!flag_ncdelta_bdt) && flag_1p) return true;
                  return false;



// end added by lhagaman

  }else if (ch_name == "nc_delta_energy_FC_0p" || ch_name == "nc_delta_score_FC_0p"
	    || ch_name == "nc_delta_energy_FC_0p_ncpio_overlay" || ch_name == "nc_delta_score_FC_0p_ncpio_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_ncdelta_overlay" || ch_name == "nc_delta_score_FC_0p_ncdelta_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_overlay" || ch_name == "nc_delta_score_FC_0p_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_ext" || ch_name == "nc_delta_score_FC_0p_ext"
	    || ch_name == "nc_delta_energy_FC_0p_dirt" || ch_name == "nc_delta_score_FC_0p_dirt" || ch_name == "nc_delta_score_FC_0p_overlay_entire"
            || ch_name == "nc_delta_score_signal_blind_FC_0p" || ch_name == "nc_delta_score_signal_blind_FC_0p_ext"|| ch_name == "nc_delta_score_signal_blind_FC_0p_dirt"
            || ch_name == "nc_delta_score_signal_blind_FC_0p_ncpio_overlay" || ch_name == "nc_delta_score_signal_blind_FC_0p_ncdelta_overlay" || ch_name == "nc_delta_score_signal_blind_FC_0p_overlay"
            || ch_name == "nc_delta_score_signal_blind_FC_0p_overlay_entire"){

    if (ch_name == "nc_delta_energy_FC_0p" ||  ch_name == "nc_delta_energy_FC_0p_ext" || ch_name == "nc_delta_energy_FC_0p_dirt"){
      if (flag_FC && flag_ncdelta_bdt && flag_0p) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_ncpio_overlay"){
      if (flag_FC && flag_ncdelta_bdt && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0 
						     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_ncdelta_overlay"){
      if (flag_FC && flag_ncdelta_bdt && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_overlay"){
      if (flag_FC && flag_ncdelta_bdt && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC && flag_ncdelta_bdt && flag_0p) return true;
    }else if (ch_name == "nc_delta_score_FC_0p" ||  ch_name == "nc_delta_score_FC_0p_ext" || ch_name == "nc_delta_score_FC_0p_dirt"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel && pfeval.reco_showerKE > 0) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_ncpio_overlay"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0 
				  && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_ncdelta_overlay"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_overlay"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel && pfeval.reco_showerKE > 0 && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
      					       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_overlay_entire"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel && pfeval.reco_showerKE > 0) return true;  //      if (flag_FC  && flag_0p) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_0p" ||  ch_name == "nc_delta_score_signal_blind_FC_0p_ext" || ch_name == "nc_delta_score_signal_blind_FC_0p_dirt"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_0p_ncpio_overlay"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                  && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_0p_ncdelta_overlay"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_0p_overlay"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0 && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                               && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_0p_overlay_entire"){
      if (flag_FC  && flag_0p && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0) return true;  //      if (flag_FC  && flag_0p) return true;}else{
    
    return false;}
  }else if (ch_name == "nc_delta_energy_FC_Np" || ch_name == "nc_delta_score_FC_Np"
	    || ch_name == "nc_delta_energy_FC_Np_ncpio_overlay" || ch_name == "nc_delta_score_FC_Np_ncpio_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_ncdelta_overlay" || ch_name == "nc_delta_score_FC_Np_ncdelta_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_overlay" || ch_name == "nc_delta_score_FC_Np_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_ext" || ch_name == "nc_delta_score_FC_Np_ext"
	    || ch_name == "nc_delta_energy_FC_Np_dirt" || ch_name == "nc_delta_score_FC_Np_dirt" || ch_name == "nc_delta_score_FC_Np_overlay_entire"
            || ch_name == "nc_delta_score_signal_blind_FC_Np" || ch_name == "nc_delta_score_signal_blind_FC_Np_ext"|| ch_name == "nc_delta_score_signal_blind_FC_Np_dirt"
            || ch_name == "nc_delta_score_signal_blind_FC_Np_ncpio_overlay" || ch_name == "nc_delta_score_signal_blind_FC_Np_ncdelta_overlay" || ch_name == "nc_delta_score_signal_blind_FC_Np_overlay"
            || ch_name == "nc_delta_score_signal_blind_FC_Np_overlay_entire"){

    if (ch_name == "nc_delta_energy_FC_Np" ||  ch_name == "nc_delta_energy_FC_Np_ext" || ch_name == "nc_delta_energy_FC_Np_dirt"){
      if (flag_FC && flag_ncdelta_bdt && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_ncpio_overlay"){
      if (flag_FC && flag_ncdelta_bdt && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
      					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_ncdelta_overlay"){
      if (flag_FC && flag_ncdelta_bdt && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_overlay"){
      if (flag_FC && flag_ncdelta_bdt && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC && flag_ncdelta_bdt && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_score_FC_Np" ||  ch_name == "nc_delta_score_FC_Np_ext" || ch_name == "nc_delta_score_FC_Np_dirt"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel && pfeval.reco_showerKE > 0) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_ncpio_overlay"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
      					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_ncdelta_overlay"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true; 
    }else if (ch_name == "nc_delta_score_FC_Np_overlay"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel && pfeval.reco_showerKE > 0 && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      // if (flag_FC  && (!flag_0p)) return true;
    } else if (ch_name == "nc_delta_score_FC_Np_overlay_entire"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel && pfeval.reco_showerKE > 0) return true;  //      if (flag_FC  && flag_0p) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_Np" ||  ch_name == "nc_delta_score_signal_blind_FC_Np_ext" || ch_name == "nc_delta_score_signal_blind_FC_Np_dirt"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_Np_ncpio_overlay"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                  && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_Np_ncdelta_overlay"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0 && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_Np_overlay"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0 && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio>0
                                               && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_signal_blind_FC_Np_overlay_entire"){
      if (flag_FC  && (!flag_0p) && flag_ncdelta_presel_signal_blind && pfeval.reco_showerKE > 0) return true;  //      if (flag_FC  && flag_0p) return true;}else{
    }else{
      return false;
    }
  }else{
    std::cout << "Not sure what cut: " << ch_name << std::endl;
  }
  
  return false;
}


bool LEEana::is_lowEnergy(KineInfo& kine, bool flag_data){
  bool flag = false;

  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  
  if (reco_Enu<=500 ) flag = true;
  return flag;

}

bool LEEana::is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;

  bool flag_numuCC = is_numuCC(tagger);
  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);

  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  
  if ((reco_Enu>=800 && tagger.nue_score >=0) ||
      (tagger.nue_score<=0 && (flag_numuCC || (flag_pi0 && flag_NC) ))) flag = true;
  return flag;
}

bool LEEana::is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  
  if (reco_Enu < 800 && tagger.nue_score>0 && (reco_Enu>=600 || tagger.nue_score<=7)) flag = true;
  
  return flag ;
}

bool LEEana::is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  if (reco_Enu < 600 && tagger.nue_score>7) flag = true;
  return flag;
}




bool LEEana::is_truth_nueCC_inside(EvalInfo& eval){
  bool flag = false;

  if (fabs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;
  
  return flag;
}

bool LEEana::is_truth_numuCC_inside(EvalInfo& eval){
   bool flag = false;

  if (fabs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;
  
  return flag;
}



bool LEEana::is_FC(EvalInfo& eval){
  if (eval.match_isFC){
    return true;
  }else{
    return false;
  }
}

bool LEEana::is_cc_pi0(KineInfo& kine, bool flag_data){

  bool flag = false;
  
  if (flag_data){
    if (kine.kine_pio_mass>0){
      //     TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      // TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      //TLorentzVector pio = p1 + p2;
      //pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;
      
      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300)
	flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
      flag = true;
  }
  
  return flag;
}


bool LEEana::is_pi0(KineInfo& kine, bool flag_data){
  bool flag = false;

  if (flag_data){
    if (kine.kine_pio_mass>0){
      //      TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      // TLorentzVector pio = p1 + p2;
      // pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;
      
      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300)
	flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
      flag = true;
  }
  
  return flag;
}


bool LEEana::is_NCpio_bdt_sel(TaggerInfo& tagger_info, KineInfo& kine){ // added lhagaman 2022_06_10 in order to add Giacomo's extra energy cuts
  bool flag = false;
  if (tagger_info.nc_pio_score > 1.816 && tagger_info.numu_cc_flag >=0 && kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.) flag = true;
  return flag;
}

bool LEEana::is_NCpio_bdt(TaggerInfo& tagger_info){
  bool flag = false;
  // if (tagger_info.nc_pio_score > 1.68 && tagger_info.numu_cc_flag >=0) flag = true;
  if (tagger_info.nc_pio_score > 1.816 && tagger_info.numu_cc_flag >=0) flag = true;
  return flag;
}

// modified by lhagaman, 2021_07_28, used to not have the reco_showerKE cut here
bool LEEana::is_NCdelta_bdt(TaggerInfo& tagger_info, PFevalInfo& pfeval){
	bool flag = false;
  if (tagger_info.nc_delta_score > 2.61 && tagger_info.numu_cc_flag >=0 && pfeval.reco_showerKE > 0) flag = true;
  return flag;
}

// added lhagaman, 2022_04_20
bool LEEana::is_near_bdt(TaggerInfo& tagger_info, PFevalInfo& pfeval){
        bool flag = false;
  if (-1.39 < tagger_info.nc_delta_score && tagger_info.nc_delta_score < 2.61 && tagger_info.numu_cc_flag >=0 && pfeval.reco_showerKE > 0) flag = true;
  return flag;
}

bool LEEana::is_NCdelta_presel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
        bool flag = false;
  if (tagger_info.numu_cc_flag >=0 && pfeval.reco_showerKE > 0) flag = true;
  return flag;
}

// added lhagaman, 2022_04_20
bool LEEana::is_NCdelta_presel_signal_blind(TaggerInfo& tagger_info, PFevalInfo& pfeval){
        bool flag = false;
  if (tagger_info.numu_cc_flag >=0 && pfeval.reco_showerKE > 0 && tagger_info.nc_delta_score < 2.61) flag = true;
  return flag;
}

bool LEEana::is_NC(TaggerInfo& tagger_info){
  bool flag = false;
  if ((!tagger_info.cosmict_flag) && tagger_info.numu_score < 0)
    flag = true;
  
  return flag;
}


bool LEEana::is_numuCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
    flag = true;
  
  return flag;
}

bool LEEana::is_loosenumuCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > -10.0)
    flag = true;
  
  return flag;
}

bool LEEana::is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0)
    flag = true;
  
  return flag;
}

bool LEEana::is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
    // 1 lepton <=1 proton 0 charged pion                                                                
    // 1 lepton guaranteed by numu cc flag                       
    // using pi0 flag to remove pi0 component in channel definition                      
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
	int pdgcode = kine.kine_particle_type->at(i);
	if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
	if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
    if(Nproton==0) flag = true;
  }

  return flag;
}

bool LEEana::is_1p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
    // 1 lepton <=1 proton 0 charged pion                                                                
    // 1 lepton guaranteed by numu cc flag                       
    // using pi0 flag to remove pi0 component in channel definition                      
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
        int pdgcode = kine.kine_particle_type->at(i);
        if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
        if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
    if(Nproton==1) flag = true;
  }

  return flag;
}



bool LEEana::is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
    // 1 lepton <=1 proton 0 charged pion                                                                
    // 1 lepton guaranteed by numu cc flag                       
    // using pi0 flag to remove pi0 component in channel definition                      
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
	int pdgcode = kine.kine_particle_type->at(i);
	if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
	if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
    if(Npion==0) flag = true;
  }

  return flag;
}

bool LEEana::is_nueCC_1e0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.nue_score > 7.0){
    // 1 lepton <=1 proton 0 charged pion                                                                                                                                          
    // 1 lepton guaranteed by numu cc flag                                                                                                                                         
    // using pi0 flag to remove pi0 component in channel definition                                                                                                                
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
	int pdgcode = kine.kine_particle_type->at(i);
	if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?                                                                  
	if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV                                                                              
      }
    if(Nproton==0) flag = true;
  }

  return flag;
}

bool LEEana::is_nueCC_1e0p0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.nue_score > 7.0){
    // 1 lepton <=1 proton 0 charged pion                                                                                                                                          
    // 1 lepton guaranteed by numu cc flag                                                                                                                                         
    // using pi0 flag to remove pi0 component in channel definition                                                                                                                
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
	int pdgcode = kine.kine_particle_type->at(i);
	if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?                                                                  
	if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV                                                                              
      }
    if(Nproton==0 && Npion==0) flag = true;
  }

  return flag;
}

bool LEEana::is_nueCC_1eNp0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.nue_score > 7.0){
    // 1 lepton <=1 proton 0 charged pion                                                                                                                                          
    // 1 lepton guaranteed by numu cc flag                                                                                                                                         
    // using pi0 flag to remove pi0 component in channel definition                                                                                                                
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
	int pdgcode = kine.kine_particle_type->at(i);
	if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?                                                                  
	if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV                                                                              
      }
    if(Nproton>0 && Npion==0) flag = true;
  }

  return flag;
}

bool LEEana::is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;
  
  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){ 
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm? 
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV 
      }
      if(Nproton==0) flag = true;
  } 
  
  return flag;
}


bool LEEana::is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data){
    bool flag = false;
   
    if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){

      double reco_Enu = get_reco_Enu_corr(kine, flag_data);
      
      Float_t Ehadron = reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
      if(Ehadron<200) // MeV
        {
	  flag = true;
        }
    }
    return flag;
}

bool LEEana::is_numuCC_cutbased(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag==1 && tagger_info.cosmict_flag==0) 
    flag = true;
  
  return flag;
}


bool LEEana::is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;
  // default 7.0
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    //  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score <= 7.0 && tagger_info.nue_score > 0)
    flag = true;
  
  return flag;
}

bool LEEana::is_loosenueCC(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 0.0)
    flag = true;
  
  return flag;
}

bool LEEana::is_generic(EvalInfo& eval){
  // not very useful for the main analysis
  bool flag = is_preselection(eval);

  flag = flag && (eval.stm_clusterlength > 15);
  return flag;
}

bool LEEana::is_preselection(EvalInfo& eval){
  bool flag = false;

  // match code ...
  int tmp_match_found = eval.match_found;
  if (eval.is_match_found_int){
    tmp_match_found = eval.match_found_asInt;
  }

  if (tmp_match_found == 1 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) flag = true;
  
  
  return flag;
}


#endif
