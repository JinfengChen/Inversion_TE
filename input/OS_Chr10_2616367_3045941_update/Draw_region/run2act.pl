use FindBin qw($Bin $Script);

opendir DIR, "./" or die "can not open my dir";
foreach my $file (readdir DIR){
    if ($file=~/(.*)\.blast/){
	   system "perl $Bin/format_blastn.pl $1";
	   system "perl $Bin/toACT.pl $1";
       print "$file\n";
       push (@seq,$1);
    }
}
close DIR;










