### script to edit and align mtDNA sequence data ###

# import modules
import os
import Bio
from Bio import SeqIO
import glob

# assemble mtDNA data from raw reads using mitofinder (run on cluster)
os.chdir('/media/burke/bigMac/ethan/syma_cleaned/')
for reads in glob.iglob('*R1.clean.fq'):
    R1 = glob.iglob('*R1.clean.fq'):
    R2 = glob.iglob('*R2.clean.fq'):
    id = filepath.split('_')[0]
    dir_command = 'mkdir ' + id
    os.system(dir_command)
    cd_command = 'cd ' + id + '/'
    os.system(cd_command)
    mitofinder = '/media/burke/bigMac/ethan/MitoFinder/mitofinder '
    mt_command = mitofinder + '-j ' + id + '-1 ' + R1 + '-2 ' + R2 + ' -r /media/burke/bigMac/ethan/t_sanctus.gb -p 24'
    os.system(mt_command)
    back_command = 'cd /media/burke/bigMac/ethan/mtDNA/'
    os.system(back_command)

# concatenate sequences
os.chdir('/Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/')
for filepath in glob.iglob('*.fasta'):
    con_seq = ""
    id = filepath.split('_')[0]
    output_filename = id + "_concat.fasta"
    output_handle = open(output_filename, "w")
    for rec in SeqIO.parse(filepath, "fasta"):
        con_seq += rec.seq
        con_seq_h = ">"+id+"\n"+con_seq+"\n"
        fasta = open(output_filename, "w+")
        print("Now writing " + id)
        fasta.write(str(con_seq_h))
        fasta.close()

# extract ND2
os.chdir('/Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/')
for filepath in glob.iglob('*.fasta'):
    id = filepath.split('_')[0]
    wanted_ids = id + "@ND2"
    output_filename = id + "_ND2.fasta"
    count = 0
    total = 0
    output_handle = open(output_filename, "w")
    for record in SeqIO.parse(filepath, "fasta"):
        total = total + 1
        if record.id in wanted_ids:
            count = count + 1
            SeqIO.write(record, output_handle, "fasta")
    output_handle.close()
    print(str(count) + " records selected out of " + str(total))

# make new directories
os.system('mkdir /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/concatenated/')
os.system('mv *concat* /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/concatenated/')
os.system('mkdir /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/ND2')
os.system('mv *ND2* /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/ND2')

# reconcatenate full mtDNA assemblies for maaft
os.system('cat /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/concatenated/*.fasta > syma_mtDNA_full.fasta')

# run maaft for full alignment
maaft_params = 'mafft --auto '
input = '/Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/concatenated/syma_mtDNA_full.fasta '
output = '> /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/concatenated/syma_mtDNA_full_msa.fasta'
os.system(maaft_params + input + output)

# reconcatenate ND2 assemblies for maaft
os.system('cat /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/ND2/*.fasta > syma_ND2_full.fasta')

# run maaft for ND2 alone
maaft_params = 'mafft --auto '
input = '/Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/concatenated/syma_ND2_full.fasta '
output = '> /Users/ethanlinck/Dropbox/syma_speciation/mtDNA/final/ND2/syma_ND2_full_msa.fasta'
os.system(maaft_params + input + output)
