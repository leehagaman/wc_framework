namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes

  TString spectra_file = "./data_inputs/merge.root";
  TString flux_Xs_directory = "./data_inputs/hist_rootfiles/XsFlux/";
  TString detector_directory = "./data_inputs/hist_rootfiles/DetVar/";
  TString mc_directory = "./data_inputs/mc_stat/";

  //////////
  
  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 19;
 
  int syst_cov_mc_stat_begin = 7;// files in mc_directory above
  //int syst_cov_mc_stat_end   = 99;
  int syst_cov_mc_stat_end   = 7; // this version should be for no LEE

  
  int channels_observation = 0;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse
                               // NOTE: This value is not used in the lastest version

  ///////////////////////////////
  
  int array_LEE_ch[1] = {0}; // for 1d strength. Pay attention to the setting.
                             // Confict will happen if both the 1d and 2d(defined below) work.
                             // Three options: 1d works, 2d works, or neither works

  //int array_LEE_ch[12] = {2, 5, 8, 11, 14, 17, 3, 6, 9, 12, 15, 18};
  
  /////////////////////////////// separated Np and 0p strengthes: 2d strengthes

  //int array_LEE_Np_ch[6] = {2, 5, 8, 11, 14, 17};// element value "0" will not be set to the LEE_ch
  //int array_LEE_0p_ch[6] = {3, 6, 9, 12, 15, 18};// element value "0" will not be set to the LEE_ch

  //int array_LEE_Np_ch[8] = {2, 5, 8, 11, 14, 17, 20, 23};// element value "0" will not be set to the LEE_ch
  //int array_LEE_0p_ch[8] = {3, 6, 9, 12, 15, 18, 21, 24};// element value "0" will not be set to the LEE_ch

  // Just doing signal scaling for the signal channels!
  int array_LEE_Np_ch[4] = {2, 5, 8, 11};// element value "0" will not be set to the LEE_ch
  int array_LEE_0p_ch[4] = {3, 6, 9, 12};// element value "0" will not be set to the LEE_ch

  /*
  Uncollapsed Bin #        Collapsed Bin #    cov_input.txt block index (channel #)    description
  1-2,                     1-2                1,                                       WC 1gNp background
  3-4,                     1-2                2,                                       WC 1gNp True Np signal
  5-6,                     1-2                3,                                       WC 1gNp True 0p signal
  7-8,                     3-4                4,                                       WC 1g0p background
  9-10,                    3-4                5,                                       WC 1g0p True Np signal
  11-12,                   3-4                6,                                       WC 1g0p True 0p signal
  13-14,                   5-6                7,                                       Pandora 1g1p background
  15-16,                   5-6                8,                                       Pandora 1g1p True Np signal
  17-18,                   5-6                9,                                       Pandora 1g1p True 0p signal
  19-20,                   6-7                10,                                      Pandora 1g0p background
  21-22,                   6-7                11,                                      Pandora 1g0p True Np signal
  23-24,                   6-7                12,                                      Pandora 1g0p True 0p signal
  25-40,                   8-23               13,                                      WC NC Pi0 Np background
  41-56,                   8-23               14,                                      WC NC Pi0 Np True Np signal
  57-72,                   8-23               15,                                      WC NC Pi0 Np True 0p signal
  73-88,                   24-39              16,                                      WC NC Pi0 0p background
  89-104,                  24-39              17,                                      WC NC Pi0 0p True Np signal
  105-120,                 24-39              18,                                      WC NC Pi0 0p True 0p signal
  121-136,                 40-55              19,                                      WC numuCC Np background
  137-152,                 40-55              20,                                      WC numuCC Np True Np signal
  153-168,                 40-55              21,                                      WC numuCC Np True 0p signal
  169-184,                 40-55              22,                                      WC numuCC 0p background
  185-200,                 40-55              23,                                      WC numuCC 0p True Np signal
  201-216,                 40-55              24,                                      WC numuCC 0p True 0p signal
  */

  // Warning: need to change Lee_test->array_no_stat_bins = new int[4];	in read_TLee_v20.cxx if we change the length here
  int array_no_stat_bins[4] = {1, 3, 5, 7}; // in reco space, these are the explicitly empty overflow bins where we want no uncertainties
  int num_no_stat_bins = 4;
  
  //int array_LEE_Np_ch[1] = {0};
  //int array_LEE_0p_ch[1] = {0};
  
  ///////////////////////////////
  
  /// some places may need to be changed when use different file-formats
  /// void TLee::Set_Spectra_MatrixCov()
  /// (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  /// (*) map_Lee_ch                     -----> tag LEE channels
  /// (*) map_detectorfile_str           -----> detector files
   
  ////////// display graphics flag

  bool flag_display_graphics = 0;
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs    = 1;
  bool flag_syst_detector   = 1;
  bool flag_syst_additional = 0;
  bool flag_syst_mc_stat    = 0;
  bool flag_syst_mc_data_stat_cor = 0; // adds off-diagonal elements to the data and pred stat cor matrices

  bool flag_syst_reweight        = 0;
  bool flag_syst_reweight_cor    = 0;

  //double Lee_strength_for_outputfile_covariance_matrix = 1;
  double Lee_strength_for_outputfile_covariance_matrix = 0;

  double Lee_Np_strength_for_outputfile_covariance_matrix = 0;
  double Lee_0p_strength_for_outputfile_covariance_matrix = 0;
  
  bool flag_plotting_systematics   = 1; // TEMPORARY
  
  ////////// goodness of fit
  
  double Lee_strength_for_GoF         = 0;
  
  double Lee_Np_strength_for_GoF      = 0;
  double Lee_0p_strength_for_GoF      = 0;

  bool flag_GoF_output2file_default_0 = 0;

  bool flag_lookelsewhere             = 0;
  
  bool flag_both_numuCC            = 0;// 1
  bool flag_CCpi0_FC_by_numuCC     = 0;// 2
  bool flag_CCpi0_PC_by_numuCC     = 0;// 3
  bool flag_NCpi0_by_numuCC        = 0;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = 0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = 0;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = 0;// 7
  bool flag_nueCC_FC_by_all        = 0;// 8

  ////////// Lee strength fitting -- data

  bool flag_Lee_strength_data = 0;
  bool flag_Lee_scan_data     = 0;

  bool flag_GOF = 0;
}
