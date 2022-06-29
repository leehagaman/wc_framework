#!/usr/bin/perl
my $num1 = scalar(@ARGV);

for (my $i=1;$i<=10;$i++){
    if ($i == 5){ # this is dE/dx, we skip this one (otherwise causes a harmless but distracting seg fault)
        $i++;
    }
    if ($num1 ==0){
	system("./bin/det_cov_matrix -r$i &");
    }else{
	print "Oscillation! \n"; 
	system("./bin/det_cov_matrix -r$i -o1 &");
    }
}
