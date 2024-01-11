#!/usr/bin/python
import sys
import statistics as stats

samples = [ line.strip().strip('"') for line in open("UDN_samples.txt", 'r')]

vntr_vcf = open("./UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.RARE.vcf", 'w')
with open("/oak/stanford/groups/euan/projects/tannerj/UDN/updated_results/Fast_001/SV_calls/vamos/Fast_001.GRCh38.vntr.vamos.vcf", 'r') as tmp_vcf:
	for line in tmp_vcf:
		if line.startswith("##"):
			print(line.strip(), file=vntr_vcf)
		elif line.startswith("#"):
			print('##INFO=<ID=MEAN_RLEN,Number=1,Type=String,Description="Mean Repeat Length">', file=vntr_vcf)
			print('##INFO=<ID=SD_RLEN,Number=1,Type=String,Description="Repeat Length standard deviation">', file=vntr_vcf)
			print('##INFO=<ID=OUTLIER_COUNT,Number=1,Type=String,Description="Count per VNTR of extreme VNTR outliers">', file=vntr_vcf)
			print('##FORMAT=<ID=RN,Number=1,Type=String,Description="Repeat Number">', file=vntr_vcf)
			print('##FORMAT=<ID=MND,Number=1,Type=String,Description="Repeat mean neighbor distance (distance to k nearest neighbors)">', file=vntr_vcf)
			print('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples), file=vntr_vcf)
		else:
			break

current_id=""
prev_id=""
vntr_dict={}
sample=""
vntr_id=""
repeat_length=""
chrom=""
pos=""
RU=""
mean="NA"
sd="NA"
outlier_count=0
with open("./UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.csv", 'r') as vntr_csv:
	for i,line in enumerate(vntr_csv):
		if i == 0: continue
		toks = line.strip().split(",")
		current_id = toks[0]
		if i > 1 and current_id != prev_id:
			out = [chrom,pos,prev_id,'N','<VNTR>','.','PASS']
			out.append("RU=%s;SVTYPE=VNTR;MEAN_RLEN=%s;SD_RLEN=%s;OUTLIER_COUNT=%d" % (RU,mean,sd, outlier_count))
			out.append("GT:RN:MND")
			for s in samples:
				if s not in vntr_dict:
					out.append(":".join(['./.','NA','NA']))
				elif "H2" not in vntr_dict[s]:
					out.append(":".join(vntr_dict[s]["H1"]))
				elif vntr_dict[s]["H1"][0] == "1" and vntr_dict[s]["H2"][0] == "1":
					out.append(":".join(["1/1",",".join([vntr_dict[s]["H1"][1],vntr_dict[s]["H2"][1]]),",".join([vntr_dict[s]["H1"][2],vntr_dict[s]["H2"][2]])]))
				elif vntr_dict[s]["H1"][0] == "1":
					out.append(":".join(["0/1",",".join([vntr_dict[s]["H2"][1],vntr_dict[s]["H1"][1]]),",".join([vntr_dict[s]["H2"][2],vntr_dict[s]["H1"][2]])]))
				elif vntr_dict[s]["H2"][0] == "1":
					out.append(":".join(["0/1",",".join([vntr_dict[s]["H1"][1],vntr_dict[s]["H2"][1]]),",".join([vntr_dict[s]["H1"][2],vntr_dict[s]["H2"][2]])]))
				else: 
					out.append(":".join(["0/0",",".join([vntr_dict[s]["H1"][1],vntr_dict[s]["H2"][1]]),",".join([vntr_dict[s]["H1"][2],vntr_dict[s]["H2"][2]])]))
			if outlier_count > 0:
				print('\t'.join(out), file=vntr_vcf)
			vntr_dict={}
			outlier_count=0
		chrom,pos,RU=current_id.split(':')
		sample=toks[1]
		haplo=toks[2]
		repeat_length=int(toks[3])
		mnd=toks[4]
		if mnd=='': mnd=0.0
		else: mnd=float(mnd)
		outlier=int(toks[5] == "TRUE")
		mean=float(toks[7]) if toks[7] != '' else 'NA'
		sd=float(toks[8]) if toks[8] != '' else 'NA'
		if outlier: outlier_count+=1
		if not sample in vntr_dict: 
			vntr_dict[sample] = {}
		vntr_dict[sample][haplo] = list(map(str,(outlier, repeat_length, mnd)))
		prev_id=current_id

out = [chrom,pos,prev_id,'N','<VNTR>','.','PASS']
out.append("RU=%s;SVTYPE=VNTR;MEAN_RLEN=%s;SD_RLEN=%s;OUTLIER_COUNT=%d" % (RU,mean,sd, outlier_count))
out.append("GT:RN:MND")
for s in samples:
	if s not in vntr_dict:
		out.append(":".join(['./.','NA','NA']))
	elif "H2" not in vntr_dict[s]:
		out.append(":".join(vntr_dict[s]["H1"]))
	elif vntr_dict[s]["H1"][0] == "1" and vntr_dict[s]["H2"][0] == "1":
		out.append(":".join(["1/1",",".join([vntr_dict[s]["H1"][1],vntr_dict[s]["H2"][1]]),",".join([vntr_dict[s]["H1"][2],vntr_dict[s]["H2"][2]])]))
	elif vntr_dict[s]["H1"][0] == "1":
		out.append(":".join(["0/1",",".join([vntr_dict[s]["H2"][1],vntr_dict[s]["H1"][1]]),",".join([vntr_dict[s]["H2"][2],vntr_dict[s]["H1"][2]])]))
	elif vntr_dict[s]["H2"][0] == "1":
		out.append(":".join(["0/1",",".join([vntr_dict[s]["H1"][1],vntr_dict[s]["H2"][1]]),",".join([vntr_dict[s]["H1"][2],vntr_dict[s]["H2"][2]])]))
	else: 
		out.append(":".join(["0/0",",".join([vntr_dict[s]["H1"][1],vntr_dict[s]["H2"][1]]),",".join([vntr_dict[s]["H1"][2],vntr_dict[s]["H2"][2]])]))
if outlier_count > 0:
    print('\t'.join(out), file=vntr_vcf)
