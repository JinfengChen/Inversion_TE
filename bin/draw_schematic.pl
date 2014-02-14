#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use SVG;

GetOptions (\%opt,"project:s","help");


my $help=<<USAGE;
perl draw_schematic.pl --project test

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}



my $svg=SVG->new(width=>800,height=>600);

$svg = Pre_inversion();
$svg = Inversion_intrachr();
my $outfile="$opt{project}.svg";
writesvg($outfile,$svg);

################
#################
sub Pre_inversion{
my $xstart=100; ## x start of line
my $ystart=180; ## y start of line
my $interval=120; ## interval between different strain
my $width=30; ## width of chr
#### line
$svg = line($xstart, $ystart, $xstart+575, $ystart, $svg);
#### Ditto
$svg = transposon_minus($xstart+20,$ystart-$width/2,45,$width,' ','Ditto','lightblue',$svg);
#### mPing1
$svg = transposon_plus($xstart+90,$ystart-$width/2,45,$width,'TAA','mPing','orange',$svg);
#### inversion
$svg = transposon_plus($xstart+165,$ystart-$width/2,150,$width,' ','Inversion','gray',$svg);
#### mPing2
$svg = transposon_minus($xstart+400,$ystart-$width/2,45,$width,'TTA','mPing','orange',$svg);
#### Fragment
$svg = transposon_plus($xstart+460,$ystart-$width/2,15,$width,' ',' ', 'green',$svg);
#### mPing3
$svg = transposon_minus($xstart+500,$ystart-$width/2,45,$width,'TAA','mPing','orange',$svg);

return $svg;
}



#################
sub Inversion_intrachr{
my $xstart=100; ## x start of line
my $ystart=300; ## y start of line
my $interval=120; ## interval between different strain
my $width=30; ## width of chr
#### line
$svg = line($xstart, $ystart, $xstart+275, $ystart, $svg);
#### inversion
$svg = curve_chr($xstart+300, $ystart-$width/2, $xstart+300, $ystart+$interval+$width/2, $width, 'Inversion', $svg);
#### Ditto
$svg = transposon_minus($xstart+20,$ystart-$width/2,45,$width,' ','Ditto','lightblue',$svg);
#### mPing1
$svg = transposon_plus($xstart+90,$ystart-$width/2,45,$width,'TAA','mPing','orange',$svg);
#### start part of inversion
$svg = transposon($xstart+165,$ystart-$width/2,140,$width,' ',' ','gray',$svg);

#### line 
$svg = line($xstart, $ystart+$interval, $xstart + 250, $ystart+$interval, $svg);
#### mPing3
$svg = transposon_plus($xstart+90,$ystart+$interval-$width/2,45,$width,'TTA','mPing','orange',$svg);
#### Fragment
$svg = transposon_minus($xstart+170,$ystart+$interval-$width/2,15,$width,' ',' ', 'green',$svg);
#### mPing2
$svg = transposon_plus($xstart+210,$ystart+$interval-$width/2,45,$width,'TAA','mPing','orange',$svg);
#### end part of inversion
$svg = transposon_minus($xstart+290,$ystart+$interval-$width/2,12,$width,' ',' ', 'gray',$svg);

return $svg;
}



##################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}

sub line{
my ($x1, $y1, $x2, $y2, $svg)=@_;
my $l = $svg->line(
       x1 => $x1, y1 => $y1,
       x2 => $x2, y2 => $y2,
       style => {
           fill => 'black',
           stroke => 'black'
       }
   );

return $svg;
}




###
#draw curved chromosome
#x1, y1 is the left up postion
#x2, y2 is the left down position
#widht is the width of chromosome bar
#svg is the object to draw on
###
sub curve_chr{
my ($x1,$y1,$x2,$y2,$width,$title,$svg)=@_;
###exterior arc for inversion
my $rx1=($y2-$y1)/2;
my $ry1=$rx1;
###interior arc for inversion
my $x3=$x1;
my $x4=$x2;
my $y3=$y1+$width;
my $y4=$y2-$width;
my $rx2=($y4-$y3)/2;
my $ry2=$rx2;
###
my $string = "M$x1,$y1 A$rx1,$ry1 0 0,1 $x2,$y2 L$x4,$y4 A$rx2,$ry2 0 0,0 $x3,$y3 L$x1,$y1";
#my $string = "M$x1,$y1 A$rx1,$ry1 0 0,1 $x2,$y2 L$x4,$y4 ";
#print $string,"\n";

my $tag = $svg->path(
        d => $string,
        style => {
            'fill' => 'gray',
            #'fill'   => 'green',
        }
    );

my $tx = $x1 + $rx2 + $width/3;
my $ty = $y1 + $ry1;

print "$tx\t$ty\n";
my $g = $svg->group(
       transform => "translate($tx, $ty) rotate(90)"
    );
my $name = $g->text(
       'text-anchor' => 'middle'
   )->cdata($title);


return $svg;
}

##############
#Draw box of transposon or gene
##############
sub transposon{
my ($x1,$y1,$width,$height,$tsd,$title,$color,$svg)=@_;
my $x_upl=$x1;
my $x_upr=$x1+$width;
my $x_downl=$x1;
my $x_downr=$x1+$width;
my $y_upl=$y1;
my $y_upr=$y1;
my $y_downl=$y1+$height;
my $y_downr=$y1+$height;
my $xv = [$x_upl, $x_upr, $x_downr, $x_downl];
my $yv = [$y_upl, $y_upr, $y_downr, $y_downl];

my $points = $svg->get_path(
       x => $xv,
       y => $yv,
       -type => 'path',
       -closed => 'true'
   );
my $tag = $svg->path(
       %$points,
       style=>{
          fill => $color
       }
   );

my $tx = $x1 + $width/2;
my $ty = $y1 + $height*3/5;
my $name = $svg->text(
       x => $tx, y => $ty,
       'text-anchor' => 'middle'
   )->cdata($title);

my $tsd_x1 = $x1 - 25;
my $tsd_x2 = $x1 + $width;
my $tsd_y  = $y1 + $height/2;

my $tsd_left = $svg->text(
       x => $tsd_x1, y => $tsd_y,
       'text-anchor' => 'start'
   )->cdata($tsd);
my $tsd_right = $svg->text(
       x => $tsd_x2, y => $tsd_y,
       'text-anchor' => 'start'
   )->cdata($tsd);

return $svg;
}

###############
#Draw box of transposon or gene on plus strand
##############
sub transposon_plus{
my ($x1,$y1,$width,$height,$tsd,$title,$color,$svg)=@_;
my $x_upl=$x1;
my $x_upr=$x1+$width;
my $x_downl=$x1;
my $x_downr=$x1+$width;
my $x_top=$x1+$width+5;
my $y_upl=$y1;
my $y_upr=$y1;
my $y_downl=$y1+$height;
my $y_downr=$y1+$height;
my $y_top=$y1+$height/2;
my $xv = [$x_upl, $x_upr, $x_top, $x_downr, $x_downl];
my $yv = [$y_upl, $y_upr, $y_top, $y_downr, $y_downl];
print "$sv->[0]\n$yv->[0]\n";

my $points = $svg->get_path(
       x => $xv,
       y => $yv,
       -type => 'path',
       -closed => 'true'
   );
my $tag = $svg->path(
       %$points,
       style=>{
          fill => $color
       }
   );

my $tx = $x1 + $width/2;
my $ty = $y1 + $height*3/5;
my $name = $svg->text(
       x => $tx, y => $ty,
       'text-anchor' => 'middle'
   )->cdata($title);

my $tsd_x1 = $x1 - 25;
my $tsd_x2 = $x1 + $width + 5;
my $tsd_y  = $y1 + $height/2;

my $tsd_left = $svg->text(
       x => $tsd_x1, y => $tsd_y,
       'text-anchor' => 'start'
   )->cdata($tsd);
my $tsd_right = $svg->text(
       x => $tsd_x2, y => $tsd_y,
       'text-anchor' => 'start'
   )->cdata($tsd);

return $svg;
}

#############
#Draw box of transposon or gene on minus strand
##############
sub transposon_minus{
my ($x1,$y1,$width,$height,$tsd,$title,$color,$svg)=@_;
my $x_upl=$x1;
my $x_upr=$x1+$width;
my $x_downl=$x1;
my $x_downr=$x1+$width;
my $x_top=$x1-5;
my $y_upl=$y1;
my $y_upr=$y1;
my $y_downl=$y1+$height;
my $y_downr=$y1+$height;
my $y_top=$y1+$height/2;
my $xv = [$x_top, $x_upl, $x_upr, $x_downr, $x_downl];
my $yv = [$y_top, $y_upl, $y_upr, $y_downr, $y_downl];
#print "$sv->[0]\n$yv->[0]\n";

my $points = $svg->get_path(
       x => $xv,
       y => $yv,
       -type => 'path',
       -closed => 'true'
   );
my $tag = $svg->path(
       %$points,
       style=>{
          fill => $color
       }
   );

my $tx = $x1 + $width/2;
my $ty = $y1 + $height*3/5;
my $name = $svg->text(
       x => $tx, y => $ty,
       'text-anchor' => 'middle'
   )->cdata($title);

my $tsd_x1 = $x1 - 30;
my $tsd_x2 = $x1 + $width;
my $tsd_y  = $y1 + $height/2;

my $tsd_left = $svg->text(
       x => $tsd_x1, y => $tsd_y,
       'text-anchor' => 'start'
   )->cdata($tsd);
my $tsd_right = $svg->text(
       x => $tsd_x2, y => $tsd_y,
       'text-anchor' => 'start'
   )->cdata($tsd);

return $svg;
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


 
