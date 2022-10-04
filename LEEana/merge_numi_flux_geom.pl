#!/usr/bin/perl

system("./bin/applyNuMIGeomtryWeights -i/data1/xqian/MicroBooNE/processed_checkout_rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root -o/data1/xqian/MicroBooNE/processed_checkout_rootfiles/prodgenie_numi_intrinsic_nue_overlay_run1/nucleoninexsec_FluxUnisim.root -mFHC -w/home/xqian/wire-cell/wcp-uboone-bdt/scripts/NuMI_Geometry_Weights_Histograms.root &");


system("./bin/applyNuMIGeomtryWeights -i/data1/xqian/MicroBooNE/processed_checkout_rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run3.root -o/data1/xqian/MicroBooNE/processed_checkout_rootfiles/prodgenie_numi_intrinsic_nue_overlay_run3/nucleoninexsec_FluxUnisim.root -mRHC -w/home/xqian/wire-cell/wcp-uboone-bdt/scripts/NuMI_Geometry_Weights_Histograms.root &");

system("./bin/applyNuMIGeomtryWeights -i/data1/xqian/MicroBooNE/processed_checkout_rootfiles/checkout_prodgenie_numi_nu_overlay_run1.root -o/data1/xqian/MicroBooNE/processed_checkout_rootfiles/prodgenie_numi_nu_overlay_run1/nucleoninexsec_FluxUnisim.root -mFHC -w/home/xqian/wire-cell/wcp-uboone-bdt/scripts/NuMI_Geometry_Weights_Histograms.root &");

system("./bin/applyNuMIGeomtryWeights -i/data1/xqian/MicroBooNE/processed_checkout_rootfiles/checkout_prodgenie_numi_nu_overlay_run3.root -o/data1/xqian/MicroBooNE/processed_checkout_rootfiles/prodgenie_numi_nu_overlay_run3/nucleoninexsec_FluxUnisim.root -mRHC -w/home/xqian/wire-cell/wcp-uboone-bdt/scripts/NuMI_Geometry_Weights_Histograms.root &");
