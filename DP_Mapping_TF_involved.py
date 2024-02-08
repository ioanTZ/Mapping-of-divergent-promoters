import os
import re
import pprint
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# import pyranges as pr


# processing of the dictionaries for the TFs
def print_dict(dict1,dict2, count1,count2):
	print("tf\tcontrol\tdisease")
	deleted = []

	for item in dict1.keys():
		dict1[item] = (100*dict1[item])/count1
		if(dict1[item] < 1):
			deleted.append(item)

	for item in deleted:
		del dict1[item]

	deleted = []	
	for item in dict2.keys():
		dict2[item] = (100*dict2[item])/count2
		if(item in dict1.keys()):
			continue
		elif(dict2[item] < 1):
			deleted.append(item)

	for item in deleted:
		del dict2[item]

	for item1, amount1 in dict1.items():
		if(item1 in dict2.keys()):
			print("{}\t{:.2f}%\t{:.2f}%".format(item1,amount1, dict2[item1]))
		
			
		
		else:
			print("{}\t{:.2f}%\t0".format(item1,amount1))			
	

merge_distance = 1000 

ref_genome = "/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Divergent_ML/Annotation_sorted.bed"
ref_gen_changed = "/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Divergent_ML/Annotation_changed.bed"
out_sorted = "/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Divergent_ML/Annot_changed_sorted.bed"


#change the values in Annotation_sorted based on strand, sort them again, write in new file
# f = open(ref_gen_changed, "w")
# with open(ref_genome,"r") as ref:
#   for line in ref:
		
#       l = line.split("\t")

#       FivePSite = l[1]
#       ThreePSite = l[2]
#       strand = l[5]

#       if(strand == "+"):
#           ThreePSite = int(FivePSite) + 1
#       elif(strand == "-"):
#           FivePSite = int(ThreePSite) - 1

#       newLine = [l[0], str(FivePSite), str(ThreePSite),l[3],l[4],l[5],l[6]]
#       newLine = ("\t").join(newLine)
#       f.write(newLine)

# ref.close()
# f.close()

# cmd = f"sort -k 1,1 -k2,2n  {ref_gen_changed} > {out_sorted}"
# os.system(cmd)




control_path = r'/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Divergent_ML/Tissue_TSS/1.DeepTSS_outFiles/Control/final_files/'
res1 = []

for path in os.listdir(control_path):
	if os.path.isfile(os.path.join(control_path, path)):
		res1.append(path)
		


disease_path = r'/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Divergent_ML/Tissue_TSS/1.DeepTSS_outFiles/Disease/final_files/'
res2 = []

for path in os.listdir(disease_path):
	if os.path.isfile(os.path.join(disease_path, path)):
		res2.append(path)
	

#CLOSEST 
control_out_path = r'/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Control_out/'
disease_out_path = r'/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Disease_out/'

# removing old files from the following directories
cmd = f'rm {control_out_path + "*"}'
os.system(cmd)
cmd = f'rm {disease_out_path + "*"}'
os.system(cmd)

# control
for file in res1:
	to_annotate = f"{control_path}{file}"
	cmd = f"/usr/bin/bedtools closest -a {to_annotate} -b {out_sorted} -s -d -t first > {os.path.dirname(control_out_path)}/Control.{file}"
	os.system(cmd)

#disease
for file in res2:
	to_annotate = f"{disease_path}{file}"                                                                                              												
	cmd = f"/usr/bin/bedtools closest -a {to_annotate} -b {out_sorted} -s -d -t first > {os.path.dirname(disease_out_path)}/Disease.{file}"
	os.system(cmd)



#finding the divergent promoters in each of these files 

res1 = []
out_path = r'/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/'
for path in os.listdir(control_out_path):
	if os.path.isfile(os.path.join(control_out_path, path)):
		res1.append(path)

res2 = []

for path in os.listdir(disease_out_path):
  if os.path.isfile(os.path.join(disease_out_path, path)):
      res2.append(path)


#MERGE

for file in res1:
	to_annotate = f"{control_out_path}{file}"
	cmd = f"/usr/bin/bedtools merge -i {to_annotate} -d {merge_distance} -c 1,12,13,14 -o count,collapse,collapse,collapse > {os.path.dirname(out_path)}/Merge.{merge_distance}.{file}"
	os.system(cmd)

for file in res2:
	to_annotate = f"{disease_out_path}{file}"
	cmd = f"/usr/bin/bedtools merge -i {to_annotate} -d {merge_distance} -c 1,12,13,14 -o count,collapse,collapse,collapse > {os.path.dirname(out_path)}/Merge.{merge_distance}.{file}"
	os.system(cmd)


#clean up the files

res1=[]
for path in os.listdir(out_path):
	if os.path.isfile(os.path.join(out_path, path)):
		res1.append(path)



for file in res1:
	div_prom = 0
	buffer = ''
	# print(file)
	with open(out_path + file,"r+") as ref:
		for line in ref:
			l = line.split("\t")
			if(l[3] != 1):
				gene_names = l[4].split(",")
				for name in gene_names:
					j=0
					count = 0
					while j< len(gene_names):
						if name == gene_names[j]:
							count += 1
							if count > 1:
								del gene_names[j]
						j+=1
				if(len(gene_names) != 1 and gene_names[0] != gene_names[1]):
					l[4] = ",".join(gene_names)
					
					strands = l[5].split(",")
					interactions = l[6].split(",")

					for strand in strands:
						j = 0
						count = 0
						while j<len(strands):
							if strand == strands[j]:
								count += 1
								if count > 1:
									del strands[j]
									del interactions[j]
							j+=1
					
					if(len(strands)!=1):
							
						l[5] = ",".join(strands)
						l[6] = ",".join(interactions)
						if(not("\n" in l[6])):
							l[6] += "\n"
						
						newLine = l[0:7]
						newLine = "\t".join(newLine)
						buffer += newLine 
						div_prom += 1
					else:
						continue	
				else:
					continue
			else:
				continue
				
	ref.close()
	ref = open(out_path + file,"w")
	ref.write(buffer)
	ref.close()
		

#count the interactions
for file in res1:
	interactions = {}
	div_prom = 0
	with open(out_path + file,"r") as ref:
		for line in ref:
			l = line.split("\t")
			div_prom += 1
			types = l[6].split(",")
			types[1] = types[1].replace("\n","")

			form1 = types[0] + "," + types[1]
			form2 = types[1] + "," + types[0]

			keys = interactions.keys()
			if(form1 in keys):
				interactions[form1] += 1
			elif(form2 in keys):
				interactions[form2] += 1
			else:
				interactions.update({form1 : 1})
									
	
	print(file)
	sorted_interactions = sorted(interactions.items(), key = lambda x:x[1], reverse = True)
	interactions = dict(sorted_interactions)
	print_dict(interactions, div_prom)	
	print("\n")		
	
	ref.close()			

#enlongation of the divergent promoters

for file in res1:

	with open(out_path + file, "r") as ref:
		buffer = ''
		for line in ref:

			l = line.split("\t")
			dist = int(l[2]) - int(l[1])
			strand = l[5].split(",")
			
			if(dist > 2000):
				amount = dist - 2000
				extra1 = int((3*amount)/4)
				extra2 = amount-extra1
				if(strand[0] == "+"):
					l[1] = str(int(l[1]) + extra1)
					l[2] = str(int(l[2]) - extra2)
				else:
					l[1] = str(int(l[1]) + extra2)
					l[2] = str(int(l[2]) - extra1)	

			else:
				amount = 2000 - dist
				extra1 = int((3*amount)/4)
				extra2 = amount-extra1
				
				if(strand[0] == "+"):
					l[1] = str(int(l[1]) - extra1)
					l[2] = str(int(l[2]) + extra2)
				else:
					l[1] = str(int(l[1]) - extra2)
					l[2] = str(int(l[2]) + extra1)	

				
			newLine = "\t".join(l)
			buffer += newLine 
		
		ref.close()
		ref = open(out_path + file,"w")
		ref.write(buffer)
		ref.close()

#build the fasta files
# res1 already has all the files we need to make the fasta files
fasta = "/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Divergent_ML/hg38.fa"

div_prom_path = out_path

for i in res1:
	file = f"{div_prom_path}{i}"
	buffer = i.split(".")
	new_name = buffer[2] + "." + buffer[3]
	cmd = f"/usr/bin/bedtools getfasta -s -fi {fasta} -bed {file} -fo {out_path}/fasta/{new_name}.fasta"
	os.system(cmd)


#get the involved transcription factors
fasta_path = r'/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/fasta/'
res = []
for path in os.listdir(fasta_path):
	if os.path.isfile(os.path.join(fasta_path, path)):
		print(path)
		res.append(path)

all_motifs = "/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/Divergent_ML/Human_727Motifs_meme_Jaspar2022.txt"

#FIMO
# getting the transcription factors
for i in res:
	file = f"{fasta_path}{i}"
	print(file)
	new_name = i.replace(".fasta","")
	new_name = new_name.replace(".","_")
	print(new_name)
	cmd = f"fimo --parse-genomic-coord --verbosity 1 --thresh 0.00001 --oc {out_path}/{new_name} {all_motifs} {file}"
	os.system(cmd)


#plot for each tissue with the TFs in control and disease state
tissues = ["Bladder","Blood","Bone","Brain","Breast","Cervix","Esophagus","Gall-Bladder","Kidney","Liver","Lung","Ovary","Pancreas","Prostate","Stomach","Uterus"]

for tissue in tissues: 
	out_path = f"/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/FIMO/Control_{tissue}/fimo.gff"

	tfs1 = {}
	div_prom1 = 0

	with open(out_path, "r") as ref:
			buffer = ''
			for line in ref:
				if("##gff" in line):
					continue
				else:
					div_prom1 += 1

					l = line.split(".\t")
					
					tf = l[1].split(";")
					tf = tf[0].replace("Name=","")
					tf = tf.split("_")
					tf = tf[0]

					keys = tfs1.keys()
					if(tf in keys):
						tfs1[tf] += 1
					else:
						tfs1.update({tf : 1})


	sorted_tfs1 = sorted(tfs1.items(), key = lambda x:x[1], reverse = True)
	tfs1 = dict(sorted_tfs1)
		

	out_path = f"/home/ioanna/Downloads/Divergent_ML-20221210T092443Z-001/FIMO/Disease_{tissue}/fimo.gff"

	tfs2 = {}
	div_prom2 = 0

	with open(out_path, "r") as ref:
			buffer = ''
			for line in ref:
				if("##gff" in line):
					continue
				else:
					div_prom2 += 1

					l = line.split(".\t")
					
					tf = l[1].split(";")
					tf = tf[0].replace("Name=","")
					tf = tf.split("_")
					tf = tf[0]

					keys = tfs2.keys()
					if(tf in keys):
						tfs2[tf] += 1
					else:
						tfs2.update({tf : 1})


	sorted_tfs2 = sorted(tfs2.items(), key = lambda x:x[1], reverse = True)
	tfs2 = dict(sorted_tfs2)
	print_dict(tfs1,tfs2,div_prom1, div_prom2)



	mylist = [tfs1, tfs2]
	df = pd.DataFrame.from_records(mylist)
	
	df = df.fillna(0)
	
	print(df)
	df = df.T
	print(df)
	df.columns = ['Control','Disease']
	print(df)


	plt.rcParams["figure.figsize"] = [6.50, 6.50]
	plt.rcParams["figure.autolayout"] = True

	# Array for horizontal bar's position
	ind = list(np.arange(0,df.shape[0]))
	ind = np.array(ind)
	print(ind)

	# Bar's width
	width = 0.4

	fig, ax = plt.subplots()

	# Horizontal bar plot
	ax.barh(ind, np.array(df.loc[:,"Control"]), width, color='orange', label='Control')
	ax.barh(ind + width, np.array(df.loc[:,"Disease"]), width, color='blue', label='Disease')

	# Set Y-axis ticks and ticklabels
	ax.set(yticks=ind + width, yticklabels=np.array(df.index.values),
	ylim=[2*width - 1, len(ind)])


	for i in ax.patches:
	    plt.text(i.get_width()+0.1, i.get_y()+0.05,
	             str(round((i.get_width()), 3)),
	             fontsize = 6, fontweight ='bold',
	             color ='grey')

	# Legend at the upper right corner
	ax.legend(loc='upper right')
	plt.title(tissue)
	# Display the plot
	plt.show()

