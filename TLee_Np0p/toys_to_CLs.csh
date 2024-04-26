echo making folders...

mkdir two_d_rootfiles/new_folder_RENAME
mkdir two_d_rootfiles/new_folder_RENAME/wc_only
mkdir two_d_rootfiles/new_folder_RENAME/glee_only
mkdir two_d_rootfiles/new_folder_RENAME/joint

echo doing wc_only...

cp fc_files/wc_only/* fc_files
root -b -l cal_CL.cc < /dev/null
mv file_map_CL.root two_d_rootfiles/new_folder_RENAME/wc_only
root -b -l cal_CL_Asimov.cc < /dev/null
mv file_map_CL_Asimov.root two_d_rootfiles/new_folder_RENAME/wc_only
root -b -l cal_CL_Wilks.cc < /dev/null
mv file_map_CL_Wilks.root two_d_rootfiles/new_folder_RENAME/wc_only
root -b -l cal_CL_Asimov_Wilks.cc < /dev/null
mv file_map_CL_Asimov_Wilks.root two_d_rootfiles/new_folder_RENAME/wc_only

echo doing glee_only...

cp fc_files/glee_only/* fc_files
root -b -l cal_CL.cc < /dev/null
mv file_map_CL.root two_d_rootfiles/new_folder_RENAME/glee_only
root -b -l cal_CL_Asimov.cc < /dev/null
mv file_map_CL_Asimov.root two_d_rootfiles/new_folder_RENAME/glee_only
root -b -l cal_CL_Wilks.cc < /dev/null
mv file_map_CL_Wilks.root two_d_rootfiles/new_folder_RENAME/glee_only
root -b -l cal_CL_Asimov_Wilks.cc < /dev/null
mv file_map_CL_Asimov_Wilks.root two_d_rootfiles/new_folder_RENAME/glee_only

echo doing joint...

cp fc_files/joint/* fc_files
root -b -l cal_CL.cc < /dev/null
mv file_map_CL.root two_d_rootfiles/new_folder_RENAME/joint
root -b -l cal_CL_Asimov.cc < /dev/null
mv file_map_CL_Asimov.root two_d_rootfiles/new_folder_RENAME/joint
root -b -l cal_CL_Wilks.cc < /dev/null
mv file_map_CL_Wilks.root two_d_rootfiles/new_folder_RENAME/joint
root -b -l cal_CL_Asimov_Wilks.cc < /dev/null
mv file_map_CL_Asimov_Wilks.root two_d_rootfiles/new_folder_RENAME/joint

echo done

