nohup ./read_TLee_v20 -f -2 > log_files/log_-2.txt &
sleep 1
nohup ./read_TLee_v20 -f -1 > log_files/log_-1.txt &
sleep 1
nohup ./read_TLee_v20 -f 0 > log_files/log_00.txt &


