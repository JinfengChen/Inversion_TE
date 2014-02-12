#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def usage():
    test="name"
    message='''
python Insert_TE.py --fasta ../input/OS_Chr10_2616367_3045941_update/OS_Chr10_2616367_3045941.fasta --gff ../input/OS_Chr10_2616367_3045941_update/OS_Chr10_2616367_3045941.gene.gff --TE ../input/OS_Chr10_2616367_3045941_update/OS_Chr10_2616367_3045941.te.gff --insertion TEinsertion.table
Insert TE sequence into sequence, update gff and create new act files.
--fasta: fasta sequence to patch gaps
--gff: gene gff
--TE: te gff
--insertion: need to start from the largest position, because once the file is modified it will not influence the position of smaller one
TE	Chr	insertion_site	strand	TSD	class	sequence
MPING	Chr10	237878	-	TAA	DNA/Harbinger	GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTTTCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGTCCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAACTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGTTTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATTGTGACTGGCC

    '''
    print message


'''

'''
def update_act(path, prefix):
    embl  = 'perl ./act/GFF2embl.pl -gff %s -embl %s -fasta %s\n' % (path + '/' + prefix + '.evolved.gene.gff', path + '/' + prefix + '.evolved.gene.embl' , path + '/' + prefix + '.evolved.fasta')
    embl += 'perl ./act/gffrepeat2embl.pl -repeat %s -embl %s -title %s\n' % (path + '/' + prefix + '.evolved.te.gff', path + '/' + prefix + '.evolved.gene.embl', path + '/' + prefix + '.evolved')
    embl += 'mv %s %s' % (path + '/' + prefix + '.evolved.merge', path + '/' + prefix + '.evolved.embl')    
    #print embl
    os.system(embl)    

    act   = 'cp %s ./\n' % (path + '/' + '*.fasta')  
    act  += 'perl ./act/runblast2seq.pl\n'
    act  += 'perl ./act/run2act.pl\n'
    act  += 'rm *.fasta *.blast *.temp *.nhr *.nin *.nsq\n'
    act  += 'mv *4ACT %s' % (path)
    #print act
    os.system(act)
 
'''
chr10   maker   gene    9314    13560   .       +       .       ID=HEG4_Os10g12220;Name=queuine tRNA-ribosyltransferase,putative,expressed;
chr10   maker   mRNA    9314    13560   .       +       .       ID=HEG4_Os10g12220.1;Parent=HEG4_Os10g12220;
'''
def update_gff(gff, update, path, update_file):
    ofile = open (path + '/' + update_file, 'w')
    with open (gff, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if int(unit[3]) > update[0]:
                    unit[3] = str(int(unit[3]) + update[1])
                    unit[4] = str(int(unit[4]) + update[1])
                    newline = '\t'.join(unit)
                    print >> ofile, newline
                else:
                    print >> ofile, line   
    ofile.close()


def write_seq(seqid, seq, ofile, path):
    ofile = open (path + '/' + ofile , 'w')
    newrecord = SeqRecord(Seq(seq), id=seqid, description="")
    SeqIO.write(newrecord, ofile, "fasta")
    ofile.close()

'''
get fasta sequence in file
'''
def fasta_seq(fastafile):
    for record in SeqIO.parse(fastafile,"fasta"):
        data = [str(record.id), str(record.seq)]
    return data



'''
chr10   RepeatMasker    Transposon      262078  262507  4057    -       .	ID=TE242864;Target=MPING 1 430;Class=DNA/Harbinger;PercDiv=0.0;PercDel=0.0;PercIns=0.0;
'''
def insertion_gff(gff, update, path, update_file, insertion):
    flag = 0
    ofile = open (path + '/' + update_file, 'w')
    with open (gff, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if int(unit[3]) > update[0]:
                    if flag == 0:
                        print >> ofile, insertion
                        flag = 1        
                    unit[3] = str(int(unit[3]) + update[1])
                    unit[4] = str(int(unit[4]) + update[1])
                    newline = '\t'.join(unit)
                    print >> ofile, newline
                else:
                    print >> ofile, line   
    ofile.close()




'''
te:
MPING   Chr10 237878  -       TAA     DNA/Harbinger   GGCCAGTCACAA

'''
def insert_TE(ref, te, path, gene, TE):
    insertions = 0
    fasta = ref[1]
    fa = ''
    with open (te, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                insertions += 1 
                unit = re.split(r'\t',line)
                update = [int(unit[2]), len(unit[6])+len(unit[4])]
                update_inf = update_fasta(fasta, ref[0]+'.evolved', update, unit[6], unit[4], 'TE_' + str(insertions) + '.fasta', './')
                update_gff(gene, update_inf, './', 'TE_' + str(insertions) + '.gene.gff')
                update_gff(TE, update_inf, './', 'TE_' + str(insertions) + '.TE.gff')
                insert_line = '%s\tRepeatMasker\tTransposon\t%s\t%s\t4057\t%s\t.\tID=%sTarget=%s 1 %s;Class=%s;PercDiv=0.0;PercDel=0.0    ;PercIns=0.0;' % (unit[1], unit[2], str(int(unit[2])+len(unit[6])), unit[3], 'Add'+str(insertions), unit[0], str(len(unit[6])), unit[5])
                insertion_gff(TE, update_inf, './', 'TE_' + str(insertions) + '.TE.IN.gff', insert_line) 
                temp = fasta_seq('TE_' + str(insertions) + '.fasta')
                fasta= temp[1]
                gene = 'TE_' + str(insertions) + '.gene.gff'
                TE   = 'TE_' + str(insertions) + '.TE.IN.gff'
                fa   = 'TE_' + str(insertions) + '.fasta'
    mv  = 'mv %s %s\n' % (fa, path + '/' + ref[0] + '.evolved.fasta') 
    mv += 'mv %s %s\n' % (gene, path + '/' + ref[0] + '.evolved.gene.gff')
    mv += 'mv %s %s\n' % (TE, path + '/' + ref[0] + '.evolved.te.gff')
    mv += 'rm TE_*'
    #print mv
    os.system(mv)
    
    update_act(path, ref[0])

'''
fasta: fasta id and sequence need to update, [id, seq]
seqid: sequence id of update_fasta
update: [0] is the boundary for TE insertion, [1] is the length that add into sequence TE+TSD
repeat: repeat sequence
TSD: TSD sequence
update_fasta: new update fasta
path: where to write the update_fasta
'''
def update_fasta(fasta, seqid, update, repeat, TSD, update_fasta, path):
        before_TE = fasta[0:update[0]]
        after_TE  = fasta[update[0]:]
        new_fasta = before_TE + repeat + TSD + after_TE
        boundary   = update[0] ### after this boundary, the gene gff are need to update
        update_len = update[1] ### after boundary, the gene gff need to add update_len

        '''update files'''
        write_seq(seqid, new_fasta, update_fasta, path)
        return [boundary, update_len] 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta')
    parser.add_argument('-g', '--gff')
    parser.add_argument('-T', '--TE')
    parser.add_argument('-i', '--insertion')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.insertion) > 0
    except:
        usage()
        sys.exit(2)

    '''get sequence'''
    ref = fasta_seq(args.fasta)

    '''files and paths'''
    directory  = os.path.abspath(os.path.dirname(args.fasta))

    '''insert TE, update files'''    
    update = insert_TE(ref, args.insertion, directory, args.gff, args.TE)
    #update_gff(args.gff, update, directory, ref_prefix + '.gene.update.gff')
    #update_act(directory + '/' + ref_prefix + '.update.fasta', directory, ref_prefix)

if __name__ == '__main__':
    main()

