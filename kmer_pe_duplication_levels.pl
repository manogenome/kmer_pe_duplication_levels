# Script        : kmer_pe_duplication_levels.pl
# Input Format  : Gzipped read1 and read2 files
# Usage         : perl ~/kmer_pe_duplication_levels.pl <R1.gz> <R2.gz> <kmer_length> <output_prefix>
# Author & Date : Manojkumar Sumathiselvaraju & July 30 2013
# Contact       : manojsbiotech [at] gmail [dot] com

$T = localtime();
print "$T Analysis Begun\n";

$trim_len = $ARGV[2];
open(F1,"gunzip -c $ARGV[0] |");
open(F2,"gunzip -c $ARGV[1] |");
while(<F1>){
        $p1_H1 = $_; $p1_S = <F1>;
        $p1_H2 = <F1>; $p1_Q = <F1>;
        $p2_H1 = <F2>; $p2_S = <F2>;
        $p2_H2 = <F2>; $p2_Q = <F2>;

        chomp ($p1_S, $p2_S);
        $p1_S = substr($p1_S,0,$trim_len);
        $p2_S = substr($p2_S,0,$trim_len);
        $uniq{$p1_S}{$p2_S} += 1;
        $pair2{$p2_S} += 1;
}
close F1; close F2;

$T = localtime();
print "$T Reading Files Complete\n";

$max = 0;
foreach $k1 (keys %uniq){
        foreach $k2 (keys %{$uniq{$k1}}){
                $d{$uniq{$k1}{$k2}} += 1;
                $r1 += $uniq{$k1}{$k2};
        }
        $max = $r1 if ($max < $r1);
        $read1{$r1} += 1; $r1 = ();
}

%uniq = (); 

foreach $k (keys %pair2){
        $read2{$pair2{$k}} += 1;
}

%pair2 = ();

$max = $max*10;
for($i=10; $i<=$max; $i*=10){
        $j = $i*10;
        $k = 51 if($i == 10);
        $k = $i+1 if($i > 10);
        if($i > 1){ 
                $s{$k}=$j;
        }
}

open(O,">$ARGV[2]\_$ARGV[3]\_Pair_PE_Duplicates.txt");
print O "Duplication_Levels,Unique_Pairs,All_Pairs,Unique_Fraction,All_Fraction\n";
$df = cal_duplicates(\%d,\%s);
close O;

$T = localtime();
print "$T $ARGV[2] $ARGV[3] PE Duplication = $df\n";

open(O,">$ARGV[2]\_$ARGV[3]\_R1_PE_Duplicates.txt");
print O "Duplication_Levels,Unique_Reads,All_Reads,Unique_Fraction,All_Fraction\n";
$df = cal_duplicates(\%read1,\%s);
close O;

$T = localtime();
print "$T $ARGV[2] $ARGV[3] R1 Duplication = $df\n";

open(O,">$ARGV[2]\_$ARGV[3]\_R2_PE_Duplicates.txt");
print O "Duplication_Levels,Unique_Reads,All_Reads,Unique_Fraction,All_Fraction\n";
$df = cal_duplicates(\%read2,\%s);
close O;

$T = localtime();
print "$T $ARGV[2] $ARGV[3] R2 Duplication = $df\n";

sub cal_duplicates () {
        my ($h1_ref, $h2_ref) = @_; %d2 = (); $uc = (); $dc = (); $df = ();
        %h1 = %$h1_ref; %h2 = %$h2_ref;
        for($i=1; $i<51; $i+=1){
                $d2{$i} = 0;
        }
        foreach $k1 (sort keys %h1){
                if($k1<51){
                        $d2{$k1} = $h1{$k1}*$k1;
                        $d3{$k1} = $h1{$k1};
                        $uc += $h1{$k1};
                        $dc += $h1{$k1}*$k1;
                }
                else{
                        foreach $k2 (keys %h2){
                                if (($k2<=$k1) and ($k1<=$h2{$k2})){
                                        $d2{"$k2:$h2{$k2}"} += $h1{$k1}*$k1;
                                        $d3{"$k2:$h2{$k2}"} += $h1{$k1};
                                        $uc += $h1{$k1};
                                        $dc += $h1{$k1}*$k1;
                                }
                        }
                }
        }
        $flag = 1; $x = 0; $y = 0;
        foreach $k (sort {$a <=> $b} keys %d2){
                $x = $d2{$k}/$dc; $y = $d3{$k}/$uc;
                print O "$k,$d3{$k},$d2{$k},$y,$x\n";
                $df += $x if($flag > 1);
                $flag += 1; $x = 0; $y = 0;
        }
        $flag = 0; return($df);
}
