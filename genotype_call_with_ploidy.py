#!/usr/bin/env python
# input vcfs after GenotypeGVCFs, and after splitting multiallelic sites.

import gzip
import numpy as np
import sys
import pandas as pd
import itertools

genes = ["KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL5A;KIR2DL5B","KIR2DS1","KIR2DS2","KIR2DS3;KIR2DS5","KIR2DS4","KIR3DL1","KIR3DL2","KIR3DL3","KIR3DS1","KIR2DL5A;KIR2DL5B","KIR2DS3;KIR2DS5","KIR2DL4"]
refs = ["KIR2DL1_001","KIR2DL2_0010101","KIR2DL3_0010101","KIR2DL5A","KIR2DS1_001","KIR2DS2_0010101","KIR2DS3","KIR2DS4_0010101","KIR3DL1_0010101","KIR3DL2_0010101","KIR3DL3_00101","KIR3DS1_010","KIR2DL5A_0010101","KIR2DS3_00101","KIR2DL4_00101"]

GENE_REF = dict(zip(genes, refs))

gene_var = {}
vars = []
gene_TypeDetPos = {}
var_TypeDetDose = {}
gene_header = {}
var_dosage = {}

with open("./data/"+gene+".difpos.all.txt","rt") as file:
	for line in file:
		if line.find("#") == -1:
			line = line.rstrip().split("\t")
			gene_TypeDetPos.setdefault(gene,[]).append([line[0],line[1],line[2]])
			varname = [gene, line[0],line[1],line[2]]
			var_TypeDetDose[tuple(varname)] = line[3:]
		else:
			gene_header[gene] = line.rstrip().split("\t")


## input real genotype and dosage info
for method in callmethod:
	with gzip.open(per_sample_vcf_file,"rt") as file:
		for line in file:
			if line.find("#") == -1:
				line = line.rstrip().split("\t")
				if line[5] != "LowQual":
					dose = line[9].split(":")[0].split("/")
					dosage = 0
					for d in dose:
						if d != ".":
							dosage += int(d)
						else:
							dosage = "NA"
							break
					if dosage != "NA":
						var_dosage[tuple([line[0],line[1],line[3],line[4]])] = dosage
						gene_var.setdefault(line[0],[]).append([line[1],line[3],line[4]])
	
	ploidy = len(dose)
	
	## output common variants for each gene (ref, alt and number) for each sample
	dosage_header = "#GENE\tPOS\tREF\tALT\t" + sample
	gene_commonvar ={}
	
	this_ref = GENE_REF[gene]
	OUT_dosage = open(dosage_file,"w")
	OUT_ref = open(reference_file,"w")
	vars_in_vcf = gene_var[gene]
	vars_in_typedet = gene_TypeDetPos[gene]
	commonvar = []
	commonvar_with_flag = {}
	
	for var in vars_in_vcf:
		if var in vars_in_typedet:
			commonvar.append(var)
			commonvar_with_flag[tuple(var)] = "SAME"
		else:
			for typevar in vars_in_typedet:
				if var[0] == typevar[0] and var[1] == typevar[2] and var[2] == typevar[1]:
					commonvar.append(var)
					commonvar_with_flag[tuple(var)] = "OPP"
				elif (var[0] == typevar[0] and var[1] == typevar[1] and var[2] == "<NON_REF>") or (var[0] == typevar[0] and var[1] == typevar[1] and var[2] == "."):
					newvar = [var[0],var[1],typevar[2]]
					commonvar.append(newvar)
					commonvar_with_flag[tuple(newvar)] = "SAME"
					var_dosage[tuple([gene]+newvar)] = var_dosage[tuple([gene]+var)]
				elif (var[0] == typevar[0] and var[1] == typevar[2] and var[2] == "<NON_REF>") or (var[0] == typevar[0] and var[1] == typevar[2] and var[2] == "."):
					newvar = [var[0],var[1],typevar[1]]
					commonvar.append(newvar)
					commonvar_with_flag[tuple(newvar)] = "OPP"
					var_dosage[tuple([gene]+newvar)] = var_dosage[tuple([gene]+var)]
	
	print(dosage_header, file = OUT_dosage)
	this_gene_header = "\t".join(gene_header[gene][0:3])+"\t"+this_ref+"\t"+"\t".join(gene_header[gene][3:])
	print(this_gene_header, file = OUT_ref)
	for var in commonvar:
		varname = [gene]
		varname += var
		out1 = "\t".join(varname)
		out2 = var_dosage[tuple(varname)]
		out = out1 + "\t" + str(out2)
		print(out, file = OUT_dosage)
		if commonvar_with_flag[tuple(var)] == "SAME":
			out1 = "\t".join(var)
			out1 += "\t0"
			out2 = "\t".join(var_TypeDetDose[tuple(varname)])
			out = out1 + "\t" + out2
			print(out, file = OUT_ref)
		elif commonvar_with_flag[tuple(var)] == "OPP":
			out = var[0] + "\t" + var[2] + "\t" + var[1] + "\t1"
			varname = (gene,var[0],var[2],var[1])
			for dose in var_TypeDetDose[varname]:
				if dose == "1":
					out += "\t0"
				elif dose == "0":
					out += "\t1"
				else:
					print("ERROR")
			print(out, file = OUT_ref)
	
	OUT_dosage.close()
	OUT_ref.close()
	
	# read ref and dosage, and remove NA, and return possible alleles
	geno = pd.read_table(dosage_file, header=0)
	geno = np.array(geno.iloc[:,4])
	type = pd.read_table(reference_file, header=0)
	type = type.T[3:]
	type_dose = {}
	for i in range(len(type.index)):
		type_dose[type.index[i]] = np.array(type.iloc[i],dtype="int")
	
	# check for ambiguous alleles
	uniq_type_dose = {}
	allelepool = []
	allelepool_type = []
	for this_type in type_dose:
		count = 0
		if len(allelepool) >= 1:
			for i in range(len(allelepool)):
				if np.array_equal(type_dose[this_type], allelepool[i]):
					pooledtype = allelepool_type[i]
					del uniq_type_dose[pooledtype]
					allelepool_type[i] = pooledtype + "-" + this_type
					uniq_type_dose[allelepool_type[i]] = type_dose[this_type]
					count += 1
		if count == 0:
			uniq_type_dose[this_type] = type_dose[this_type]
			allelepool.append(type_dose[this_type])
			allelepool_type.append(this_type)
	
	
	# make combinations according to ploidy
	uniq_type_dose_list = []
	uniq_type_name_list = []
	uniq_type_shortname_list = []
	count = 0
	for key in uniq_type_dose:
		uniq_type_name_list.append(key)
		uniq_type_dose_list.append(uniq_type_dose[key])
		uniq_type_shortname_list.append(gene+"::"+str(count))
		count += 1
	ALL_COMB_name = list(itertools.combinations_with_replacement(uniq_type_name_list,ploidy))
	ALL_COMB_dose = list(itertools.combinations_with_replacement(uniq_type_dose_list,ploidy))
	ALL_COMB_shortname = list(itertools.combinations_with_replacement(uniq_type_shortname_list,ploidy))
	ploid_name = []
	ploid_dose = []
	for i in range(len(ALL_COMB_name)):
		name = "/".join(ALL_COMB_name[i])
		shortname = "/".join(ALL_COMB_shortname[i])
		dose = np.sum(np.array(ALL_COMB_dose[i]),axis=0)
		ploid_name.append(name)
		ploid_dose.append(dose)
	
	## check for ambiguous ploidies
	uniq_ploid_name = []
	uniq_ploid_dose = []
	tmp_dose_name = {}
	
	for i in range(len(ploid_dose)):
		dose = ploid_dose[i]
		name = ploid_name[i]
		tmp_dose_name.setdefault(tuple(dose),[]).append(name)
	
	for dose in tmp_dose_name:
		uniq_ploid_dose.append(np.array(dose))
		uniq_ploid_name.append("-or-".join(tmp_dose_name[dose]))
	
	
	## determine ploidy-combination for each sample
	sample_dose = np.array(geno,dtype="float") # float is for nans
	count = 0
	closeness = []
	
	for j in range(len(uniq_ploid_dose)):
		if np.allclose(sample_dose,uniq_ploid_dose[j]):
			SAMPLE_GENETYPE = uniq_ploid_name[j]
			count += 1
		else:
			diff = abs(sample_dose - uniq_ploid_dose[j])
			tmp = sum(diff)
			closeness.append(tmp)
	
	if count == 0:
		min_index = [h for h, x in enumerate(closeness) if x == min(closeness)]
		CLOSETYPE = []
		for id in min_index:
			CLOSETYPE.append(uniq_ploid_name[id])
			closeval = str(closeness[id])
		closetype_out = "Close_to_" + "-OR-".join(CLOSETYPE) + "[score=" + closeval + "]"
		SAMPLE_GENETYPE = closetype_out
	
	OUT = open(output_allele_file,"w")
	out = "\t".join([sample,gene,SAMPLE_GENETYPE])
	print(out, file=OUT)
	OUT.close()

