#!/usr/bin/python

# http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/

from Bio.SeqIO import parse
from Bio.Seq import Seq
from numpy import nan
import sys

seqFile = sys.argv[1] # "Mytilus_edulis#solo200.cds"
ingroup = sys.argv[2] # "Mytilus_edulis"
outgroup = sys.argv[3] # "provincialis"


def polyM(x):
	# x = alignement[ingroup]
	res = {} # returned output containing syno / nonSyno lengths, piS and piN
	LSyn = 0.0 # number of synonymous and non synonymous positions
	LNonSyn = 0.0
	seq1 = x[0]
	nCodons = len(seq1) / 3
	codons = {} # codons[codon#][list of codons over N individuals]
	cnt = -1
	for i in range(nCodons): # loop over codons
		cnt += 1
		codons[cnt] = []
		pos = range(i*3, (i*3)+3, 1)
		for j in x:
			tmp = "".join(j[pos[0]:(pos[2]+1)])
			codons[cnt].append(tmp)
	N = len(codons[0]) # number of individuals
	# retainedCodons = list of retained codons; nTranslatedElements = list of number of alternative translated elements;  nElements = list of number of alternative elements/codons
	retainedCodons = []
	nElements = []
	nTranslatedElements = []
	cnt = -1
	for i in codons:
		cnt += 1
		test = 1
		nElements.append(len(set(codons[i]))) # number of alternative codons
		nTranslatedElements.append(len(set([ Seq(j).translate()[0] for j in codons[i] ]))) # number of alternative aa
		for j in range(N):
#version1			if "N" in codons[i][j]: # reject codon positions if one individual has a N
			if codons[i][j] not in codonTable: # version2: reject codon positions if codon not found in codonTable
				test = 0
		if test == 1: # a codon is retained if no individual has a "N"
			for k in range(N):
				LSyn += codonTable[codons[i][j]]['nS']
				LNonSyn += codonTable[codons[i][j]]['nN']
			retainedCodons.append(cnt)
	LSyn /= N
	LNonSyn /= N
	piS = 0.0
	piN = 0.0
	if len(retainedCodons) >= 10:
		for i in retainedCodons:
			if nElements[i] == 2 and nTranslatedElements[i] == 1:
				tmp = list(set(codons[i]))
				codon1 = tmp[0]
				codon2 = tmp[1]
				piS += (codons[i].count(codon1) * codons[i].count(codon2))/(N * (N-1) / 2.0)
			if nElements[i] == 2 and nTranslatedElements[i] == 2:
				tmp = list(set(codons[i]))
				codon1 = tmp[0]
				codon2 = tmp[1]
				piN += (codons[i].count(codon1) * codons[i].count(codon2))/(N * (N-1) / 2.0)
		res["LSyno"] = LSyn
		res["LNonSyno"] = LNonSyn
		res["piS"] = piS/LSyn
		res["piN"] = piN/LNonSyn
	else:
		res["LSyno"] = nan
		res["LNonSyno"] = nan
		res["piS"] = nan
		res["piN"] = nan
	res["max_N_indiv"] = N
	return(res)


# bases
bases = ["A", "T", "G", "C"]


# table codons
codonTable = {}
for pos1 in bases:
	for pos2 in bases:
		for pos3 in bases:
			codon = pos1 + pos2 + pos3
			aa = Seq(codon).translate()
			if aa != "*":
				codonTable[codon] = {}
				codonTable[codon]["aa"] = aa[0]
				codonTable[codon]["nS"] = 0.0
				codonTable[codon]["nN"] = 0.0
				for i in bases:
					if i != pos1:
						if Seq(pos1 + pos2 + pos3).translate() == Seq(i + pos2 + pos3).translate():
							codonTable[codon]["nS"] += 1
						else:
							codonTable[codon]["nN"] += 1
					if i != pos2:
						if Seq(pos1 + pos2 + pos3).translate() == Seq(pos1 + i + pos3).translate():
							codonTable[codon]["nS"] += 1
						else:
							codonTable[codon]["nN"] += 1
					if i != pos3:
						if Seq(pos1 + pos2 + pos3).translate() == Seq(pos1 + pos2 + i).translate():
							codonTable[codon]["nS"] += 1
						else:
							codonTable[codon]["nN"] += 1
				codonTable[codon]["nS"] /= 3.0
				codonTable[codon]["nN"] /= 3.0


# record sequences with 'alignement'
alignement = {}
input = parse(seqFile, "fasta")
for i in input:
	tmp = i.id.split("|")
	sequence = tmp[0]
	species = tmp[1]
	if sequence not in alignement:
		alignement[sequence] = {}
	if species not in alignement[sequence]:
		alignement[sequence][species] = []
	L = len(i.seq)
	L = L - L%3
	alignement[sequence][species].append(i.seq[0:L])
input.close()


# treatment of sequences
LSyn, LNonSyn, piS, piN, N, loci = [], [], [], [], [], []


print("analyse")
for i in alignement.keys():
	res = polyM(alignement[i][ingroup])
	LSyn.append(res['LSyno'])
	LNonSyn.append(res['LNonSyno'])
	piS.append(res['piS'])
	piN.append(res['piN'])
	N.append(res['max_N_indiv'])
	loci.append(i)


output = "Loci\tnInd\tLsyn\tLasyn\tpiS\tpiN\n"
for i in range(len(N)):
	output += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(loci[i], N[i], LSyn[i], LNonSyn[i], piS[i], piN[i])


outfile = open("test_output_pNpSdNdS", "w")
outfile.write(output)
outfile.close()


