$test_file="/u/group/hallc_ana/xem/pass3/replay/radcor/working/TEST/5.766_18.0_2.out";

open(IN, $test_file) or die "couldn't open ${test_file}\n";
while(<IN>)
{
  $line=$_;
  chomp($line);
  $line=~ s/^ +//;
  ($e1,$rest)=split (/ +/,$line,2);
  if ($e1==5.766)
  {
    ($e2,$thc, $qsq, $x, $dis, $qe, $sigma, $sigrad, $rest)=split(/ +/,$rest,9);
    $e2=~ s/0+$//;
    $DATA{$e2}{'qsq'}=$qsq;
    $DATA{$e2}{'x'}=$x;
    $DATA{$e2}{'sigma'}=$sigma;
    
  }
}
close(IN);

$collect=1;
$total=0;
$ave=0;
$ave_e=0;
$first=1;
open (OUT, ">resolution_out.dat") or die "couldn't open output file\n";
foreach $ep (sort {$a<=>$b} (keys %DATA))
{
    
  $qsq=$DATA{$ep}{'qsq'};
  $x=$DATA{$ep}{'x'};
  $sigma=$DATA{$ep}{'sigma'};
  if($first)
  {
    $bins=int(0.001*$ep/0.0005);
    if ($bins % 2)
    {
      $bins++;
    }
    $size=0.0005*$bins;
    $estart=$ep;
    $ep_stop=$ep+$size;
    $first=0;
    $ave_e=$ep;
    $ave=$sigma;
    $total++;
    print "starting with ${ep} going for ${size}\n";
  }
  if ($collect)
  {
    if($ep<=$ep_stop)
    {
      $ave_e=$ave_e+$ep;
      $ave=$ave+$sigma;
      $total++;
    }
    else
    {
      $ave_e=$ave_e+$ep;
      $ave=$ave+$sigma;
      $total++;
      $ave_e=$ave_e/$total;
      $ave=$ave/$total;
#      $ave_e=int($ave_e*10000);
#      $mod=$ave_e % 5;
#      print "mod of ${ave_e} is ${mod}\n";
#      $ave_e=$ave_e/10000;
#      $eref=int($ave_e/5);
#      $closest_e=$eref*0.0005+0.00025;
      $closest_bin=int((($ave_e-$e_start)+0.5*0.0005)/0.0005);
      $closest_e=$e_start+$closest_bin*0.0005;
      if ($closest_e==$ave_e) 
      {
	$var=1;
      }
      elsif($closest_e<$ave_e) 
      {
	$var=2;
      }
      else
      {
	$var=3;
      }
#      if ($mod>3)
#      {
#	$closest_e=$closest_e+0.0005;
#      }
#      $ave_e=$ave_e/10000;
      $sigma_ref=$DATA{$closest_e}{'sigma'};
      $x=$DATA{$closest_e}{'x'};
      print OUT "${ave_e} ${ave} ${closest_e} ${sigma_ref} ${x} ${var}\n";
      $collect=1;
      $total=0;
      $ave=0;
      $ave_e=0;
      $first=1;
      
    }
    
  }

}
close (OUT);
