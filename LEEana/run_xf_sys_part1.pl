#!/usr/bin/perl

for (my $i=1;$i<11;$i++){
    if ($i == 10){
	system("./bin/xf_cov_matrix -r$i ");
    }else{
	system("./bin/xf_cov_matrix -r$i &");
    }
}
