perl ../../bin/sub_draw_files.pl --region ../../bin/Region.txt
cd Draw_region
perl runblast2seq.pl
perl run2act.pl
rm *.fasta.n* *.temp *.blast
perl ../../../bin/drawRegionNway.pl --project inversion
