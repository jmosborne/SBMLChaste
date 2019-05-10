#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.

cd ../../../

initial_sim=0;
offset=0;
num_runs=1;
num_sweeps=1;

for (( i=0 ; i<${num_sweeps} ; i++))
do
	start_sim=`expr $i \* $num_runs + $initial_sim + $offset`;
	offset_num_runs=`expr $num_runs-$offset`;
	echo $start_sim
    	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
    	# ">" directs std::cout to the file.
    	# "2>&1" directs std::cerr to the same place.
    	# "&" on the end lets the script carry on and not wait until this has finished.
		nice -20 projects/TanCollaboration/build/optimised/TestSweepCryptWithSbmlSrnModelRunner -run_index $start_sim -num_runs $offset_num_runs -is_tan true -is_gamma_1 false > projects/TanCollaboration/test/output/TanCryptRun_${i}_Output.txt 2>&1 &
		nice -20 projects/TanCollaboration/build/optimised/TestSweepCryptWithSbmlSrnModelRunner -run_index $start_sim -num_runs $offset_num_runs -is_tan false -is_gamma_1 true > projects/TanCollaboration/test/output/VLGamma1CryptRun_${i}_Output.txt 2>&1 &
done

cd projects/TanCollaboration/test

echo "Jobs submitted"
