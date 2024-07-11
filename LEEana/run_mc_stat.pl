#!/usr/bin/perl

for (my $i=0;$i != 100; $i++){
    my $lee_strength = 0.15 * $i;
    # check if lee_strength is 0, then add 0.00001
    if ($lee_strength == 0){
        $lee_strength = 0.00001;
    }
    if ($i % 30 == 29){
	system("./bin/merge_hist -r0 -e2 -l$lee_strength > ./mc_stat/$i\.log");
    }else{
	system("./bin/merge_hist -r0 -e2 -l$lee_strength > ./mc_stat/$i\.log &");
    }
    
}
