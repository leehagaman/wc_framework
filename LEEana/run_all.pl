#!/usr/bin/perl


# convert_histo.pl
open(infile,"./configurations/cv_input.txt");

my $num1 = scalar(@ARGV);
if ($num1 !=0){
    $num1 = $ARGV[0];
    print "$num1\n";
}
my $num = 0;
while(<infile>){
    @temp = split(/\s+/,$_);

    if ($temp[0] == -1) {
        last;
    }

    if ($temp[0] !=-1 && $temp[0] != "\#file"){
        print "$temp[3]\n";
        if ($num %12 == 11){
            if ($num1 == 0){
                system("./bin/convert_checkout_hist $temp[3] $temp[4]");
            }elsif ($num1==1){
                system("./bin/convert_checkout_hist_xs $temp[3] $temp[4]");
            }elsif ($num1==2){
                system("./bin/convert_checkout_hist $temp[3] $temp[4] -o1");
            }
        }else{
            if ($num1 == 0){
                system("./bin/convert_checkout_hist $temp[3] $temp[4]&");
            }elsif ($num1 == 1){
                system("./bin/convert_checkout_hist_xs $temp[3] $temp[4]&");
            }elsif ($num1 == 2){
                system("./bin/convert_checkout_hist $temp[3] $temp[4] -o1&");
            }
        }
    }
    $num ++;
};

# merge_hist steps

system("./bin/merge_hist -r0 -l0 -e2 > ./mc_stat/0.log");
system("./bin/merge_hist -r0 -l0");

# run_xf_sys.pl
my $num1 = scalar(@ARGV);

for (my $i=1;$i<18;$i++){
    if ($i == 6 or $i == 12 or $i == 17){
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

# TLee step
system("./bin/read_TLee_v20");

# plotting step
system("./bin/plot_hist -r0 -l0 -e3 -t3 -sfile_collapsed_covariance_matrix.root");





