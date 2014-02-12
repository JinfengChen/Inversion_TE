#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
GetOptions (\%opt,"headers:s", "region:s" ,"project:s","help");


my $help=<<USAGE;
perl ../../bin/sub_draw_files.pl --region ../../bin/Region.txt
--region:
Regions to draw, other regions will be skiped.
Name	Start1	End1	Start2	End2
OG	161178	167039	265109	273845		
OS	114737	120734	232590	240318
OS.evolved	114737	121164	232620	241638
HEG4.update	138420	142951	258982	264015
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{project} ||="test";

### common value
my @feature=("OG_10_2247279_2720175","OS_Chr10_2616367_3045941","OS_Chr10_2616367_3045941.evolved","HEG4_chr10_2404092_2868093.update");

###Regions
my $region_draw=read_regions($opt{region});
`mkdir Draw_region`;

for(my $i=0;$i<@feature;$i++){
   my $genegff=parseGFF("$feature[$i].gene.gff", $region_draw->{$feature[$i]});
   my $tegff=parseTEGFF("$feature[$i].te.gff", $region_draw->{$feature[$i]});
   my $subfa=getfasta("$feature[$i].fasta", $region_draw->{$feature[$i]});
} 

`cp $Bin/act/*.pl Draw_region`;
print "Run in Draw_region:\nperl runblast2seq.pl\nperl run2act.pl\n";

#####
sub parseGFF
{
my ($gff, $region)=@_;
my $ofile = './Draw_region/'.basename($gff);
print "$ofile\n";
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
#print "$region->[0]\t$region->[1]\t$region->[2]\t$region->[3]\n";
open OUT, ">$ofile" or die "$!"; 
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[3] >= $region->[0] and $unit[4] <= $region->[1]){
       $unit[3] = $unit[3] - $region->[0];
       $unit[4] = $unit[4] - $region->[0];
    }elsif($unit[3] >= $region->[2] and $unit[4] <= $region->[3]){
       $unit[3] = $unit[3] - $region->[2] + $region->[1]-$region->[0]+1001;
       $unit[4] = $unit[4] - $region->[2] + $region->[1]-$region->[0]+1001;
    }else{
       next;
    }
    my $line = join("\t", @unit);
    print OUT "$line\n";
    if ($unit[2]=~/mRNA/){
        #print "GENE: $unit[3]\t$unit[4]\t$unit[8]\n";

        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$line\n";
        $hash{$id}=$record;
        
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$line\n";
    }

}
close IN;
close OUT;
return \%hash;
}

#####
sub parseTEGFF
{
my ($gff, $region)=@_;
my $ofile = './Draw_region/'.basename($gff);
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open OUT, ">$ofile" or die "$!";
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[3] >= $region->[0] and $unit[4] <= $region->[1]){
       $unit[3] = $unit[3] - $region->[0];
       $unit[4] = $unit[4] - $region->[0];
    }elsif($unit[3] >= $region->[2] and $unit[4] <= $region->[3]){
       $unit[3] = $unit[3] - $region->[2] + $region->[1]-$region->[0]+1001;
       $unit[4] = $unit[4] - $region->[2] + $region->[1]-$region->[0]+1001;
    }else{
       next;
    }
    my $line = join("\t", @unit);
    print OUT "$line\n";
    if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
        $id=$1;
        $hash{$id}="$line";
    }

}
close IN;
close OUT;
return \%hash;
}


####
# Name    Start1  End1    Start2  End2
# OG      161178  167039  265109  273845          
# OS      114737  120734  232590  240318
# OS.evolved      114737  121164  232620  241638
# HEG4.update     138420  142951  258982  264015

sub read_regions{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $length = $unit[2] - $unit[1] + 1 + $unit[4] - $unit[3] + 1;
    $hash{$unit[0]}=[$unit[1],$unit[2],$unit[3],$unit[4], $length];
}
close IN;
return \%hash;
}



sub getfasta
{
$/=">";
my %hash;
my ($file, $region)=@_;
my $ofile='Draw_region/'.basename($file);
open OUT, ">$ofile" or die "$!";
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    my $part1 = substr($seq, $region->[0], $region->[1]-$region->[0]+1);
    my $part2 = substr($seq, $region->[2], $region->[3]-$region->[2]+1);
    my $gap   = "N"x1001;    
    my $subseq= $part1.$gap.$part2;
    print OUT ">$head\n$subseq\n";
}
close IN;
close OUT;
$/="\n";
return \%hash;
}
