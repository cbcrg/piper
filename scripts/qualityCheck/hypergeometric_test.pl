#!/soft/bin/perl -w
use Pod::Usage;
use Getopt::Long;


my $help                   = 0;
my ($all,$count1,$count2,$jointcount);

# parsing parameters
my $result = GetOptions ("all=s"        => \$all,
			 "count1=s"     => \$count1,
			 "count2=s"     => \$count2,
			 "jointcount=s" => \$jointcount,
                         "help"         => \$help);
$help = 1 unless ($result);

if ((! $help) && (! defined $all)) {
  warn "Error! -all parameter is missing.";
  $help = 1;
}
if ((! $help) && (! defined $count1)) {
  warn "Error! -count1 parameter is missing.";
  $help = 1;
}
if ((! $help) && (! defined $count2)) {
  warn "Error! -count2 parameter is missing.";
  $help = 1;
}
if ((! $help) && (! defined $jointcount)) {
  warn "Error! -jointcount parameter is missing.";
  $help = 1;
}
pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1) if ($help);



my ($p2,$z,$ex)=statistics();
print "expected: $ex  z-score: $z p-value: $p2\n";




##FUNCTIONS
sub statistics{

    $logPdep=logfac(3)+logfac($count1-$jointcount)+logfac($count2-$jointcount)+logfac($jointcount)+logfac($all+$jointcount-$count1-$count2)-logfac($all+3);
    $logPind=logfac($count1)+logfac($all-$count1)+logfac($count2)+logfac($all-$count2)-2*logfac($all+1);
    $zscore=exp($logPdep-$logPind);

    $EXP=($count1+1)*($count2+1)/($all+2);

    $iEXP=int($EXP);

    my $pvalue=0;
    $rjointcount = int($jointcount+0.5);
    $rcount1 = int($count1+0.5);
    $rcount2 = int($count2+0.5);
    $min=$rcount2;
    if ($rcount2>$rcount1){
	$min=$rcount1;
    }
    $start=$rcount1+$rcount2-$all;
    if ($start<0){
	$start=0;
    }
    if ($rjointcount>$iEXP){
	for ($w=$rjointcount;$w<=$min;$w++){
	    $pvalue+=exp(lgam($rcount1+1)-lgam($rcount1-$w+1)-lgam($w+1)+lgam($all-$rcount1+1)-lgam($all-$rcount1-$rcount2+$w+1)-lgam($rcount2-$w+1)-lgam($all+1)+lgam($all-$rcount2+1)+lgam($rcount2+1));

	}
    }
    else{

	for ($w=$rjointcount;$w>=$start;$w--){
	    $pvalue+=exp(lgam($rcount1+1)-lgam($rcount1-$w+1)-lgam($w+1)+lgam($all-$rcount1+1)-lgam($all-$rcount1-$rcount2+$w+1)-lgam($rcount2-$w+1)-lgam($all+1)+lgam($all-$rcount2+1)+lgam($rcount2+1));

	}
    }
    return ($pvalue,$zscore,$iEXP);
}

sub lgam {
    my ($x) = @_;
    my $res;
    if($x == 0){
        $x = 1;
    }
    $res = ($x + 0.5) * log($x + 5.5) - $x - 5.5;
    $res += log(2.5066282746310005 *
               (1.000000000190015 +
                76.18009172947146 / ($x + 1) -
                86.50532032941677 / ($x + 2) +
                24.01409824083091 / ($x + 3) -
                1.231739572460166 / ($x + 4) +
                0.001208650973866179 / ($x + 5) -
                0.000005395239384953 / ($x + 6)
                ) / $x);
    return $res;
}

sub logfac{
    $result=lgam($_[0]+1);
    return $result;
}



=head1 NAME
hypergeometric_test.pl - it performs an hypergeometric test


=head1 SYNOPSIS

perl hypergeometric_test.pl  -all integer -count1 integer -count2 integer -jointcount integer [-help]

=head1 DESCRIPTION

=item * It accepts 4 numbers that represent the size of all the set (-all), the size of une subgroup (-count1) and the size of another one (-count2)

=item * It assess the statistical significance of the overlap found among the 2 sugroups (-jointcount)

=item * It returns 3 values, an e_value, a z score and a p_value (which should correspond to the probability in a CUMULATIVE hypergeometric distribution)

=item * In other words you have a population, like and urn (all) and you extract a certain number of balls (count1) without reputting the balls inside after each estraction

=item * Given a certain number of withe balls in the urn (count2) you wanna know what is the probability of finding a certain number of withe balls in your extraction (-jointcount)

=item * You can compare the results with this online service:   http://easycalculation.com/statistics/hypergeometric-distribution.php

=item * For more informations ask Ionas Erb!

=head1 OPTIONS

=item B<-all> I<integer>

integer

=item B<-count1> I<integer>

integer

=item B<-count2> I<integer>

integer

=item B<-jointcount> I<integer>

integer


Help

=cut
