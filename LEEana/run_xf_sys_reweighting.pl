#!/usr/bin/perl
my $num1 = scalar(@ARGV);

for (my $i=1;$i<20;$i++){
    if ($i% 7 == 6){ # end the current simultaneous statements, don't use the ampersand
	if ($num1==0){
	    system("./bin/xf_cov_matrix -r$i ");
	}else{
	    print "Oscillation! \n"; 
	    system("./bin/xf_cov_matrix -r$i -o1");
	}
    }else{
	if ($num1 ==0){
	    system("./bin/xf_cov_matrix -r$i &");
	}else{
	    print "Oscillation! \n"; 
	    system("./bin/xf_cov_matrix -r$i -o1&");
	}
    }
}
