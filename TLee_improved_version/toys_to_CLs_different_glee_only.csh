echo making folders...

mkdir text_CL_outputs/new_folder_RENAME

echo doing glee_only...
rm input_fc_files/*.root
cp input_fc_files/glee_only/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_only_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_only_asimov.txt

echo doing glee_only_tweaked...
rm input_fc_files/*.root
cp input_fc_files/glee_only_tweaked/* input_fc_files
cp input_fc_files/file_data_data.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_only_tweaked_data.txt
cp input_fc_files/file_data_asimov.root input_fc_files/file_data.root
root -l -q plot_stat_FC.cc+ >> text_CL_outputs/new_folder_RENAME/glee_only_tweaked_asimov.txt

echo done



