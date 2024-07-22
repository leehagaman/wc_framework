echo making folders...

mkdir two_d_rootfiles/new_folder_RENAME
mkdir two_d_rootfiles/new_folder_RENAME/wc_only
mkdir two_d_rootfiles/new_folder_RENAME/glee_only
mkdir two_d_rootfiles/new_folder_RENAME/joint

echo doing wc_only...

rm fc_files/*.root
cp fc_files/wc_only/* fc_files
root -l -q cal_CL.cc
mv file_map_CL.root two_d_rootfiles/new_folder_RENAME/wc_only
root -l -q cal_CL_Asimov.cc 
mv file_map_CL_Asimov.root two_d_rootfiles/new_folder_RENAME/wc_only
root -l -q cal_CL_Wilks.cc 
mv file_map_CL_Wilks.root two_d_rootfiles/new_folder_RENAME/wc_only
root -l -q cal_CL_Asimov_Wilks.cc 
mv file_map_CL_Asimov_Wilks.root two_d_rootfiles/new_folder_RENAME/wc_only
root -l -q cal_CL_Mixed_Asimov_Data.cc
mv file_map_CL_mixed_Asimov_data.root two_d_rootfiles/new_folder_RENAME/wc_only
root -l -q cal_CL_Mixed_Asimov_Data_Wilks.cc
mv file_map_CL_mixed_Asimov_data_Wilks.root two_d_rootfiles/new_folder_RENAME/wc_only

echo doing glee_only...

rm fc_files/*.root
cp fc_files/glee_only/* fc_files
root -l -q cal_CL.cc 
mv file_map_CL.root two_d_rootfiles/new_folder_RENAME/glee_only
root -l -q cal_CL_Asimov.cc 
mv file_map_CL_Asimov.root two_d_rootfiles/new_folder_RENAME/glee_only
root -l -q cal_CL_Wilks.cc 
mv file_map_CL_Wilks.root two_d_rootfiles/new_folder_RENAME/glee_only
root -l -q cal_CL_Asimov_Wilks.cc 
mv file_map_CL_Asimov_Wilks.root two_d_rootfiles/new_folder_RENAME/glee_only
root -l -q cal_CL_Mixed_Asimov_Data.cc
mv file_map_CL_mixed_Asimov_data.root two_d_rootfiles/new_folder_RENAME/glee_only
root -l -q cal_CL_Mixed_Asimov_Data_Wilks.cc
mv file_map_CL_mixed_Asimov_data_Wilks.root two_d_rootfiles/new_folder_RENAME/glee_only

echo doing joint...

rm fc_files/*.root
cp fc_files/joint/* fc_files
root -l -q cal_CL.cc 
mv file_map_CL.root two_d_rootfiles/new_folder_RENAME/joint
root -l -q cal_CL_Asimov.cc 
mv file_map_CL_Asimov.root two_d_rootfiles/new_folder_RENAME/joint
root -l -q cal_CL_Wilks.cc 
mv file_map_CL_Wilks.root two_d_rootfiles/new_folder_RENAME/joint
root -l -q cal_CL_Asimov_Wilks.cc 
mv file_map_CL_Asimov_Wilks.root two_d_rootfiles/new_folder_RENAME/joint
root -l -q cal_CL_Mixed_Asimov_Data.cc
mv file_map_CL_mixed_Asimov_data.root two_d_rootfiles/new_folder_RENAME/joint
root -l -q cal_CL_Mixed_Asimov_Data_Wilks.cc
mv file_map_CL_mixed_Asimov_data_Wilks.root two_d_rootfiles/new_folder_RENAME/joint

echo done

