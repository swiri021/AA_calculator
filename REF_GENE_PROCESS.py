__author__ = 'junheeyun'
from FILE_READER import FILE_READ

class REF_GENE:

	def exons_processing(self, ref, refgroup, cds_refgroup):

		#print ref
		#print cds_refgroup

		result_pos_arr = []
		result_seq = ''
		cds_result_seq = ''
		chromosome = ''

		#### r_name : >NM_032291_hg19_1_25 3 0 1 chr1:67000042-67000051+
		#### r_seq : MME

		total_exon = 0
		rev_flag = 1

		total_exon_check = []
		temp_ag = []

		for r in refgroup:
			r_name = r[0]
			r_seq = r[1]

			r_arr1 = r_name.split(" ")
			name_info = r_arr1[0]
			origin_pos = str(r_arr1[4])

			pos = origin_pos.split(":")
			chromosome = pos[0]

			ppos = pos[1].split("-")
			start = int(ppos[0])
			end = str(ppos[1])
			if end.find("+") > -1:
				end = end[:len(end)-1]
				rev_flag = 0
			end = int(end)

			r_arr2 = name_info.split("_")
			exon_numb = int(r_arr2[3])
			total_exon = int(r_arr2[4])
			total_exon_check.append(total_exon)

			cds_seq = self.matching_cds(name_info+":"+origin_pos, cds_refgroup, rev_flag)
			temp_ag.append([exon_numb,chromosome,start,end,r_seq, total_exon, cds_seq])

		set_texon = list(set(total_exon_check))

		if len(set_texon)>1:
			correct_total = max(set_texon)
			total_exon = correct_total
			temp_ag = self.chr_exon_clear(temp_ag,correct_total)

		#print temp_ag


		counter = 1

		while counter <= total_exon:
			for ex_count in temp_ag:
				chromosome = ex_count[1] ####change chromosome

				if counter == int(ex_count[0]):
					t_seq = ex_count[4]
					result_seq = result_seq + t_seq
					cds_result_seq = cds_result_seq + ex_count[6]
					#result_seq = result_seq + t_seq[::-1]

					for vv in range(int(ex_count[3])-int(ex_count[2])+1):
						if rev_flag==0:
							result_pos_arr.append(int(ex_count[2])+vv)
						else:
							result_pos_arr.append(int(ex_count[3])-vv)
			counter+=1

		result_seq = result_seq.strip()
		cds_result_seq = cds_result_seq.strip()

		return [result_seq,result_pos_arr,rev_flag, chromosome, cds_result_seq]


	def matching_cds(self,tag, cds_dat, rev):
		for c in cds_dat:
			c_name = c[0]
			c_seq = c[1]

			c_arr1 = c_name.split(" ")
			pos = str(c_arr1[4])
			cname_info= c_arr1[0]

			tag_comp = cname_info+":"+pos

			if tag == tag_comp:
				if rev==0:
					return c_seq
				else:
				#print tag +"------"+tag_comp
					return c_seq[::-1]



	def calculating(self, aa, group_seq, gene_name, ref_name):

		aminos = group_seq[0]
		pos_group = group_seq[1]
		flag = int(group_seq[2])
		chrom = group_seq[3]

		cds = group_seq[4]

		input_amino = aa[0]
		input_amino_numb = int(aa[1:])
		array_numb = input_amino_numb - 1


		pos_result = []

		try:
			if aminos[array_numb] == input_amino:
				if flag==0:

					pos_result = [gene_name, aa, ref_name, chrom,str(pos_group[array_numb * 3]), str(pos_group[array_numb * 3]+1), str(pos_group[array_numb * 3]+2)]
					cds_result = self.cds_calculating(array_numb * 3, cds)

					pos_result = pos_result + cds_result

				else:

					pos_result = [gene_name, aa, ref_name,chrom,str(pos_group[array_numb * 3]-2), str(pos_group[array_numb * 3]-1), str(pos_group[array_numb * 3])]
					cds_result = self.cds_calculating(array_numb * 3, cds)

					pos_result = pos_result + cds_result


			else:
				pos_result = [gene_name, aa, ref_name,chrom,"WRONG_REF","WRONG_REF","WRONG_REF","WRONG_REF","WRONG_REF","WRONG_REF"]

			return pos_result

		except IndexError:
			return [gene_name, aa, ref_name,chrom,"WRONG_REF","WRONG_REF","WRONG_REF","WRONG_REF","WRONG_REF","WRONG_REF"]


	def cds_calculating(self, numbering, dat):
		return [dat[numbering]+str(numbering+1), dat[numbering+1]+str(numbering+2), dat[numbering+2]+str(numbering+3)]



	def chr_exon_clear(self, arr, filter_stand):
		temp_cleared = []
		cleared_arr = []
		for x in arr:
			if len(x[1].split("_")) == 1: ###chromosome chr5_dbb_
				temp_cleared.append(x)

		for y in temp_cleared:
			if filter_stand==y[5]:
				cleared_arr.append(y)

		return cleared_arr

	def ext_exongroup(self, ref_id, base_dat):
		res =[]
		bool_dat = [ref_id in d[0] for d in base_dat]
		index_dat = [i for i, x in enumerate(bool_dat) if x==True]

		for ind in index_dat:
			ind = int(ind)
			res.append([base_dat[ind][0], base_dat[ind+1][0]])
		return res

	def sources(self, gene_symbol, aa, refseq_arr):
		self.gene_symbol = gene_symbol
		self.refseq_arr = refseq_arr
		refs_result = []

		for ref in self.refseq_arr:
			ref =  ref+"_hg19"

			ref_group = self.ext_exongroup(ref,self.file_data)
			cds_ref_group = self.ext_exongroup(ref,self.cds_file_data)

			ref_group_seq = self.exons_processing(ref,ref_group,cds_ref_group)
			refs_result.append(self.calculating(aa, ref_group_seq, gene_symbol, ref))

		return refs_result


	def __init__(self, file_name, cds_file_name):
		f1 = FILE_READ(file_name,"\t")
		self.file_data = f1.returning()

		f2 = FILE_READ(cds_file_name,"\t")
		self.cds_file_data = f2.returning()


	def __del__(self):
		print "DONE"