echo making folders...

mkdir text_CL_outputs/new_folder_RENAME

echo doing wc_only...
rm input_fc_files/*.root
cp input_fc_files/wc_only/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/wc_only_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/wc_only_asimov.txt

echo doing glee_only...
rm input_fc_files/*.root
cp input_fc_files/glee_only/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_only_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_only_asimov.txt

echo doing joint...
rm input_fc_files/*.root
cp input_fc_files/joint/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/joint_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/joint_asimov.txt

echo doing glee_1g1p...
rm input_fc_files/*.root
cp input_fc_files/glee_1g1p_only/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_1g1p_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_1g1p_asimov.txt

echo doing glee_1g0p...
rm input_fc_files/*.root
cp input_fc_files/glee_1g0p_only/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_1g0p_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_1g0p_asimov.txt

echo doing wc_1gNp...
rm input_fc_files/*.root
cp input_fc_files/wc_1gNp_only/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/wc_1gNp_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/wc_1gNp_asimov.txt

echo doing wc_1g0p...
rm input_fc_files/*.root
cp input_fc_files/wc_1g0p_only/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/wc_1g0p_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/wc_1g0p_asimov.txt

echo done



