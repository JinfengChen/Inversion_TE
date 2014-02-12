python Insert_TE.py --fasta ../input/OS_Chr10_2616367_3045941_update/OS_Chr10_2616367_3045941.fasta --gff ../input/OS_Chr10_2616367_3045941_update/OS_Chr10_2616367_3045941.gene.gff --TE ../input/OS_Chr10_2616367_3045941_update/OS_Chr10_2616367_3045941.te.gff --insertion TEinsertion.table

cd ../input/OS_Chr10_2616367_3045941_update
perl ../../bin/sub_draw_files.pl --region ../../bin/Region.txt
cd Draw_region
perl runblast2seq.pl
perl run2act.pl
perl ../../../bin/drawRegionNway.pl --project inversion

