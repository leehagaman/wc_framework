#!/bin/bash

# nue
root -l -b -q 'calculate_num.C(0.2, 0.54,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")' 
root -l -b -q 'calculate_num.C(0.54, 0.705,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")' 
root -l -b -q 'calculate_num.C(0.705, 0.805,  1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'
root -l -b -q 'calculate_num.C(0.805, 0.920,  1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'
root -l -b -q 'calculate_num.C(0.920, 1.05,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'
root -l -b -q 'calculate_num.C(1.05, 1.2,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'
root -l -b -q 'calculate_num.C(1.2, 1.375,    1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'
root -l -b -q 'calculate_num.C(1.375, 1.57,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'
root -l -b -q 'calculate_num.C(1.57, 2.05,    1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'
root -l -b -q 'calculate_num.C(2.05, 4.0,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_nue_total")'

# numu (for numuCC no pi0 channel)
root -l -b -q 'calculate_num.C(0.2, 0.54,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")' 
root -l -b -q 'calculate_num.C(0.54, 0.705,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")' 
root -l -b -q 'calculate_num.C(0.705, 0.805,  1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(0.805, 0.920,  1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(0.920, 1.05,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.05, 1.2,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.2, 1.375,    1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.375, 1.57,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.57, 2.05,    1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(2.05, 4.0,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'


# numu (for numuCC pi0 channel)
root -l -b -q 'calculate_num.C(0.275, 4.,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")' 
root -l -b -q 'calculate_num.C(0.54, 0.705,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")' 
root -l -b -q 'calculate_num.C(0.705, 0.805,  1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(0.805, 0.920,  1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(0.920, 1.05,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.05, 1.2,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.2, 1.375,    1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.375, 1.57,   1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(1.57, 2.05,    1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'
root -l -b -q 'calculate_num.C(2.05, 4.0,     1000, 6.516748e+20, 0, 0, "gh_averaged_bnb_5GeV_flux.root", "gh_averaged_numu_total")'

