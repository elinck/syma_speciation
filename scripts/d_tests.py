import sys
import random
import os
import csv

os.chdir("/Users/ethanlinck/Dropbox/syma_speciation/data/")

# read alignment
alignment = []
with open("syma_nuc_d3.fasta", 'r') as aln_file: #read in banana alignment
	for line in aln_file:
		alignment.append(str.strip(line))

torotoro_sequence = ''.join(alignment[alignment.index('>torotoro')+1:alignment.index('>megarhyncha')]) #Get torotoro sequence
megarhyncha_sequence = ''.join(alignment[alignment.index('>megarhyncha')+1:alignment.index('>ochracea')]) #Get megarhyncha sequence
ochracea_sequence = ''.join(alignment[alignment.index('>ochracea')+1:alignment.index('>senegalensis')]) #get ochracea sequence
senegalensis_sequence = ''.join(alignment[alignment.index('>senegalensis')+1:8]) #get ochracea sequence


#Function to calculcate D
def calculate_D(sp1_genome, sp2_genome, sp3_genome, sp4_genome):

	ABBA_count, BABA_count = 0, 0

	for i in range(len(sp1_genome)):
		if sp2_genome[i] == sp3_genome[i] and sp2_genome[i] != sp4_genome[i] and sp4_genome[i] == sp1_genome[i]: #ABBA sites
			ABBA_count += 1
		elif sp2_genome[i] == sp4_genome[i] and sp3_genome[i] != sp4_genome[i]	and sp3_genome[i] == sp1_genome[i]: #BABA sites
			BABA_count += 1

	D = (float(ABBA_count) - float(BABA_count))/(float(ABBA_count) + float(BABA_count)) #D statistic

	return D


#Function to estimate variance of D and D3 using block jackknife
def D_block_bootstrap(sp1_genome, sp2_genome, sp3_genome, sp4_genome, n_replicates, filename):

	D_estimates = []

	for i in range(n_replicates): #for each bootstrap replicate

		sp1_bootstrapped, sp2_bootstrapped, sp3_bootstrapped, sp4_bootstrapped = [], [], [], [] #to store bootstrapped genomes

		for j in range(100): #for each resampling window

			sampling_index = random.randint(0, len(sp1_genome)) #index to slice

			#sample windows
			sp1_bootstrapped.append(sp1_genome[sampling_index:sampling_index+10000])
			sp2_bootstrapped.append(sp2_genome[sampling_index:sampling_index+10000])
			sp3_bootstrapped.append(sp3_genome[sampling_index:sampling_index+10000])
			sp4_bootstrapped.append(sp4_genome[sampling_index:sampling_index+10000])

		#concatenate resampling windows
		sp1_bootstrapped = ''.join(sp1_bootstrapped)
		sp2_bootstrapped = ''.join(sp2_bootstrapped)
		sp3_bootstrapped = ''.join(sp3_bootstrapped)
		sp4_bootstrapped = ''.join(sp4_bootstrapped)

		#estimate statistics
		D_estimates.append(calculate_D(sp1_bootstrapped, sp2_bootstrapped, sp3_bootstrapped, sp4_bootstrapped))

		# write to csv
		with open(filename, 'wb+') as f:
    			writer = csv.writer(f)
    			for item in D_estimates:
        			writer.writerow([item])

	#estimate means
	mean_D = sum(D_estimates)/float(len(D_estimates))

	return mean_D

# write results from all three possible species tree topologies to file
print(D_block_bootstrap(ochracea_sequence, megarhyncha_sequence, torotoro_sequence, senegalensis_sequence, 200, "/Users/ethanlinck/Dropbox/syma_speciation/data/d_d3_och_meg.csv"))
print(D_block_bootstrap(ochracea_sequence, torotoro_sequence, megarhyncha_sequence, senegalensis_sequence, 200, "/Users/ethanlinck/Dropbox/syma_speciation/data/d_d3_och_tor.csv"))
print(D_block_bootstrap(megarhyncha_sequence, torotoro_sequence, ochracea_sequence, senegalensis_sequence, 200, "/Users/ethanlinck/Dropbox/syma_speciation/data/d_d3__tor.csv"))
