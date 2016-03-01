__author__ = 'junheeyun'

from sys import argv
from FILE_READER import FILE_READ
from REF_GENE_PROCESS import REF_GENE

#f1 = FILE_READ("newPanel_design.edit.txt","\t")
f1 = FILE_READ(argv[1],"\t")
REF = REF_GENE("data/refGene.exonAA.edit.fa","data/refGene.exonNuc.edit.fa")
f2 = FILE_READ("data/refgene_list","\t")
f3 = FILE_READ("data/rep_refGene.txt","\t")

#rw = open("test_input.result.txt","w")
rw1 = open(argv[1].replace(".txt",".result.txt"),"w")


input_file = f1.returning()
refgene_ids = f2.returning_selected(4,5)
rep_refgene_ids = f3.returning_selected(0,1)


def proce(aam):
	temp = aam[0]
	switch = 0
	for x in aam:
		if x.isdigit()==True:
			if switch < 2:
				temp = temp + x
		else:
			switch += 1
	return temp

def selecting(comp, comp_arr, ext_el): ######word, array, extracted position

	result = []
	selected_comps = [comp in x for x in comp_arr]
	if not selected_comps:
		print comp + " is not contained in array"

	selected_comp_true = [i for i, x in enumerate(selected_comps) if x==True]

	for true_el in selected_comp_true:
		true_el = int(true_el)
		result.append(comp_arr[true_el][ext_el])

	result = filter(None, result)
	result = list(set(result))

	return result

for a in input_file:
	#print str(a)
	gene = a[0]
	if gene=="March10":
		gene = "MARCH10"
	elif gene=="MLL3":
		gene = "KMT2C"
	elif gene=="LST3":
		gene = "KMT2C"

	aa = proce(a[1])
	#refs = []

	refs = selecting(gene, refgene_ids, 1)
	reps = selecting(gene, rep_refgene_ids, 1)

	total_refs = refs+reps
	total_refs = list(set(total_refs))
	output_result = REF.sources(gene, aa, total_refs)


	empty_count = 0

	for o in output_result:
		if "WRONG_REF" not in o:
			empty_count+=1
			rw1.write("\t".join(o)+"\n")

	if empty_count == 0:
		rw1.write(gene+"\t"+a[1]+"\t"+"PLEASE CHECK THE CODE\n")

rw1.close()







