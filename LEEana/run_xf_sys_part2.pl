#!/usr/bin/perl

for (my $i=11;$i<18;$i++){
    if ($i == 17){
	system("./bin/xf_cov_matrix -r$i ");
    }else{
	system("./bin/xf_cov_matrix -r$i &");
    }
}
