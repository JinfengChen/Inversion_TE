#!/usr/bin/perl
use Getopt::Long;
use SVG;

GetOptions (\%opt,"headers:s","project:s","help");


my $help=<<USAGE;
Draw N way using ACT files. Modify the head in \@feature.
perl $0 --headers --project 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{project} ||="test";


our $height=600;
our $width=800;
my $svg=SVG->new(width=>$width,height=>$height);

### common value
#my $rate  =100000/600; ### kb/width
my @feature=("OG_10_2247279_2720175","OS_Chr10_2616367_3045941","OS_Chr10_2616367_3045941.evolved","HEG4_chr10_2404092_2868093.update");
#my @feature=split(",",$opt{headers});
my $maxlen=maxlen(@feature);
my $rate  =($maxlen+0.2*$maxlen)/600; ### kb/width
print "MAX: $maxlen\nRATE: $rate\n";

###draw chromosome feature
my $firstx=200; ### x of first chromosome
my $firsty=95; ### y of first chromosome
my $yinterval=100;  ### interval between very chromosome, adjust with the number of chromosome to draw

for(my $i=0;$i<@feature;$i++){
   my $genegff=parseGFF("$feature[$i].gene.gff");
   my $width=getfastalen("$feature[$i].fasta")/$rate;
   my $xx1=$firstx;
   my $xx2=$xx1+$width;
   my $yy1=$firsty+$i*$yinterval;
   $svg=drawXgene($svg,$genegff,$rate,$xx1,$xx2,$yy1,$feature[$i]);
   if (-f "$feature[$i].te.gff"){
      my $tegff=parseTEGFF("$feature[$i].te.gff");
      $svg=drawXTE($svg,$tegff,$rate,$xx1,$yy1);
   }
} 


###Draw chromosome links
my $bary;
for(my $i=1;$i<@feature;$i++){
   my $act=$feature[$i-1]."VS".$feature[$i]."4ACT";
   my $h1=($i-1)*$yinterval+$firsty+15;
   my $h2=$i*$yinterval+$firsty-15;
   $bary=$h2+60;
   drawlinkACT($act,$rate,$h1,$h2);
}
##

###Draw scale bar
my $barx1=550;
my $barx2=$barx1+5000/$rate;
my $scalebar=$svg->line(
     x1=>$barx1,y1=>$bary,
     x2=>$barx2,y2=>$bary,
     style=>{stroke=>'black'}
   );
my $scale=$svg->text(
     x=>$barx2+5,y=>$bary,
     style=>{
         'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
     }
)->cdata("5 Kb");

$svg=legend($svg,$barx1-150,$bary-5,"GENE","black");
$svg=legend($svg,$barx1-100,$bary-5,"LTR","red");
$svg=legend($svg,$barx1-50,$bary-5,"DNA","blue");

my $outfile="$opt{project}.svg";
writesvg($outfile,$svg);

#############

sub legend{
my ($svg,$x,$y,$name,$color)=@_;
my $xtext=$x+18; 
 $svg->rectangle(
            x=>$x,y=>$y,
            width=>10,height=>10,
            style=>{
                fill=>$color
            }
 );
 $svg->text(
            x=>$xtext,y=>$y+6,
            style=>{
               'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
            }
 )->cdata($name);
 
return $svg;
}



sub maxlen
{
my (@header)=@_;
my $max=0;
for(my $i=0;$i<@header;$i++){
   my $fa=$header[$i].".fasta";
   my $len=getfastalen($fa);
   $max= $len > $max ? $len : $max;
}
return $max;
}


#############
sub parseACT
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split(" ",$_);
   push (@com,[$unit[2],$unit[3],$unit[5],$unit[6]]);
}
close IN;
return \@com;
}

#############
sub parsecoord
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   next unless ($_=~/^\s*\d+/);
   my @unit=split(" ",$_);
   push (@com,[$unit[1],$unit[2],$unit[4],$unit[5]]);
}
close IN;
return \@com;
}


#############
sub parseblastm8
{
my ($file)=@_;
my @com;
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   push (@com,[$unit[6],$unit[7],$unit[8],$unit[9]]);
}
close IN;
return \@com;
}

sub drawlinkACT
{
my ($act,$rate,$h1,$h2)=@_;
my $identity=50;
my $lencut=50;
my @links;
open IN, "$act" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_ eq "");
    my @array=split(" ",$_);
    if ($array[1] >= $identity and $array[3]-$array[2] > $lencut ){
        push (@links,"$array[2]\t$array[3]\t$array[5]\t$array[6]");
    }
}
close IN;

foreach (@links){
    my @unit=split("\t",$_);
    my $qleft=$unit[0]/$rate+200;
    my $qright=$unit[1]/$rate+200;
    my $tleft=$unit[2]/$rate+200;
    my $tright=$unit[3]/$rate+200;
    my $color;
    if ($tright < $tleft){
         $color='red';
    }else{
         $color='#778899';
    }

    my $qheight=$h1;
    my $theight=$h2;
    my $xv=[$qleft,$qright,$tright,$tleft];
    my $yv=[$qheight,$qheight,$theight,$theight];
    my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
    my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                        #fill=>'#FFDAB9' 
                     }
              );
}
}

#######drawlink blastm8
sub drawlinkm8
{
my ($m8,$rate)=@_;
my $identity=50;
my $lencut=50;
my @links;
open IN, "$m8" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_ eq "");
    my @array=split(" ",$_);
    if ($array[2] >= $identity and $array[7]-$array[6] > $lencut ){
        push (@links,"$array[6]\t$array[7]\t$array[8]\t$array[9]");
    }
}
close IN;

foreach (@links){
    my @unit=split("\t",$_);
    my $qleft=$unit[0]/$rate+50;
    my $qright=$unit[1]/$rate+50;
    my $tleft=$unit[2]/$rate+50;
    my $tright=$unit[3]/$rate+50;
    my $color;
    if ($tright < $tleft){
         $color='red';
    }else{
         $color='#778899';
    }

    my $qheight=110;
    my $theight=290;
    my $xv=[$qleft,$qright,$tright,$tleft];
    my $yv=[$qheight,$qheight,$theight,$theight];
    my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
    my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                        #fill=>'#FFDAB9'
                     }
              );
}
}



sub drawdotplot
{
my ($svg,$compare,$rate,$x,$y)=@_;
foreach my $match (@$compare){
   my $x1=$match->[0]/$rate+$x;
   my $y1=$y-$match->[2]/$rate;
   my $x2=$match->[1]/$rate+$x;
   my $y2=$y-$match->[3]/$rate;
   #print "$x1\t$x2\t$y1\t$y2\n";
   my $line=$svg->line(
     x1=>$x1,y1=>$y1,
     x2=>$x2,y2=>$y2,
     style=>{stroke=>'black'}
   );   
}
return $svg;
}


sub drawXgene
{
my ($svg,$refgenegff,$rate,$x1,$x2,$y,$head)=@_;
print "Start:$x1\tEnd:$x2\n";
my $strandline=$svg->line(
     x1=>$x1,y1=>$y,
     x2=>$x2,y2=>$y,
     style=>{stroke=>'black'}
);
my $title=$svg->text(
     x=>$x1-190,y=>$y
)->cdata("$head");

foreach my $g (keys %$refgenegff){
    my @line=split("\n",$refgenegff->{$g});
    my @pos;
    my $strand;
    foreach my $e (@line){
        my @unit=split("\t",$e);
        if ($unit[2] eq "mRNA"){
           $strand=$unit[6];
        }else{
           push (@pos,[$unit[3],$unit[4]]);
        }
    }
    @pos=sort {$a->[0] <=> $b->[1]} @pos;
    my $gstart=$pos[0][0]/$rate+$x1;
    my $gend  =$pos[$#pos][1]/$rate+$x1; 
    if ($strand eq "+"){
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-50,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$g");
=cut
       my $geneline=$svg->line(
          x1=>$gstart,y1=>$y-6,
          x2=>$gend,y2=>$y-6,
          style=>{stroke=>'black'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
           my $start=$e->[0]/$rate+$x1;
           my $elen =($e->[1]-$e->[0]+1)/$rate;
           my $exony=$y-10;
           my $exon=$svg->rectangle(
              x=>$start, y=>$exony,
              width=>$elen,height=>8,
              style=>{
                fill=>'black'
              }
           );
       }   
    }else{
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-10,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$g");
=cut
       my $geneline=$svg->line(
          x1=>$gstart,y1=>$y+6,
          x2=>$gend,y2=>$y+6,
          style=>{stroke=>'black'}
       );
       foreach my $e (sort {$a->[0] <=> $b->[1]} @pos){
           my $start=$e->[0]/$rate+$x1;
           my $elen =($e->[1]-$e->[0]+1)/$rate;
           my $exony=$y+2;
           my $exon=$svg->rectangle(
              x=>$start, y=>$exony,
              width=>$elen,height=>8,
              style=>{
                fill=>'black'
              }
           );
       }
    }
}
return $svg;
}
#####
sub drawXTE
{
my ($svg,$reftegff,$rate,$x1,$y)=@_;
foreach my $te (keys %$reftegff){
    my @line=split("\t",$reftegff->{$te});
    my $gstart=$line[3]/$rate+$x1;
    my $gend  =$line[4]/$rate+$x1;
    my $strand=$line[6];
    my $type=$1 if ($line[8]=~/Class=(.*?);/);
        my $color="gray";
    if ($type=~/DNA/){
       $color="blue";
    }elsif($type=~/LTR/){
       $color="red";
    }
    $type=~s/DNA\///;
    $type=~s/LTR\///;
    #print "$te\t$y\t$gstart\t$gend\t$strand\t$type\n";
    if ($strand eq "+"){
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y-14,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$type");
=cut 
             my $qleft =$gstart;
              my $qright=$gend;
              my $tright=$gend;
              my $tleft =$gstart;
              my $qheight=$y-2;
              my $theight=$y-10;
              my $xv=[$qleft,$qright,$tright,$tleft];
              my $yv=[$qheight,$qheight,$theight,$theight];
              my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
              my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                     }
              );

    }else{
=pod
       my $geneid=$svg->text(
          x=>$gstart,y=>$y+14,
          style=>{
             'font-size'=>'50%','text-anchor'=>'start','font-weight'=>'100'
          }
       )->cdata("$type");
=cut
              my $qleft =$gstart;
              my $qright=$gend;
              my $tright=$gend;
              my $tleft =$gstart;
              my $qheight=$y+2;
              my $theight=$y+10;
              my $xv=[$qleft,$qright,$tright,$tleft];
              my $yv=[$qheight,$qheight,$theight,$theight];
              my $points =$svg->get_path(
                     x=>$xv,y=>$yv,
                     -type=>'polyline',
                     -closed=>'true'
              );
              my $tag=$svg->polyline(
                     %$points,
                     style=>{
                        fill=>$color
                     }
              );      
    }
}
return $svg;
}


###
sub getfastalen
{
$/=">";
my %hash;
my $len;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    $len=length $seq;
}
$/="\n";
return $len;
}


####
sub drawbox
{
my ($svg,$x1,$x2,$y1,$y2)=@_;
my $hline=$svg->line(
     x1=>$x1,y1=>$y2,
     x2=>$x2,y2=>$y2,
     style=>{stroke=>'black'}
);
my $vline=$svg->line(
     x1=>$x1,y1=>$y1,
     x2=>$x1,y2=>$y2,
     style=>{stroke=>'black'}
);
return $svg;
}


#########
sub drawxaxis
{
my ($svg,$x1,$x2,$y,$min,$max,$step,$rate)=@_;
#print "$x1\t$x2\t$y\n";
my $xaxis=$svg->line(
     x1=>$x1,y1=>$y,
     x2=>$x2,y2=>$y,
     style=>{stroke=>'black'}
);

my $bottomline=$svg->line(
     x1=>$x2,y1=>$y,
     x2=>$x2,y2=>$y+5,
     style=>{stroke=>'black'}
);
my $tail =int ($max/1000);
=pod
my $bottomtext=$svg->text(
     x=>$x2,y=>$y+20,
     style=>{
         'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
     }
)->cdata("$tail kb");
=cut

for(my $i=$min;$i<$max;$i=$i+$step){
     my $tempx=$x1+($i-$min+1)/$rate;
     #print "$tempx\t$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$tempx,y1=>$y,
         x2=>$tempx,y2=>$y+5,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000);
     my $text=$svg->text(
         x=>$tempx+3,y=>$y+20,
         style=>{
             'font-size'=>'70%','text-anchor'=>'end','font-weight'=>'100'
         }
     )->cdata($tempi);
}
return $svg;
}

#########
sub drawyaxis
{
my ($svg,$y1,$y2,$x,$min,$max,$step,$rate)=@_;
#print "$x1\t$x2\t$y\n";
my $yaxis=$svg->line(
     x1=>$x,y1=>$y1,
     x2=>$x,y2=>$y2,
     style=>{stroke=>'black'}
);

my $bottomline=$svg->line(
     x1=>$x,y1=>$y2,
     x2=>$x+5,y2=>$y2,
     style=>{stroke=>'black'}
);
my $tail =int ($max/1000);
=pod
my $bottomtext=$svg->text(
     x=>$x,y=>$y2+20,
     style=>{
         'font-size'=>'70%','text-anchor'=>'start','font-weight'=>'100'
     }
)->cdata("$tail kb");
=cut
for(my $i=$min;$i<=$max;$i=$i+$step){
     my $tempy=$y1-($i-$min+1)/$rate;
     #print "$tempy\t$min\t$step\t$i\t$rate\n";
     #print "$tempx\t$y\n";
     my $line=$svg->line(
         x1=>$x,y1=>$tempy,
         x2=>$x+5,y2=>$tempy,
         style=>{stroke=>'black'}
     );
     my $tempi=int ($i/1000);
     my $text=$svg->text(
         x=>$x+20,y=>$tempy+3,
         style=>{
             'font-size'=>'70%','text-anchor'=>'end','font-weight'=>'100'
         }
     )->cdata($tempi);
}
return $svg;
}
#####
sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;
return \%hash;
}

#####
sub parseTEGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
        $id=$1;
        $hash{$id}="$_";
    }

}
close IN;
return \%hash;
}


################################### sub for write svg to file
sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/rhome/cjinfeng/software/tools/draw/svg2xxx_release/svg2xxx $file -t pdf -m 700";
}


