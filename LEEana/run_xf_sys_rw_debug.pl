#!/usr/bin/perl
my $num1 = scalar(@ARGV);

for (my $i=1;$i<20;$i++){
    if ($i% 9 == 8){
	if ($num1==0){
	    system("./bin/xf_cov_matrix -r$i > xf_debug_$i.txt");
	}else{
	    print "Oscillation! \n"; 
	    system("./bin/xf_cov_matrix -r$i -o1");
	}
    }else{
	if ($num1 ==0){
	    system("./bin/xf_cov_matrix -r$i > xf_debug_$i.txt &");
	}else{
	    print "Oscillation! \n"; 
	    system("./bin/xf_cov_matrix -r$i -o1&");
	}
    }
}
