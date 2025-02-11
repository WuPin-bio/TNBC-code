# script to calculate CCF and define subclonal/clonal mutations
#1. ......................................generate input data......................................
import os
import re
import subprocess
import gzip
final_cell=[]
file=open('./cell-list.txt','r')
lines=file.readlines()
for line in lines:
        data=line.rstrip()
        if data not in final_cell:
                final_cell.append(data)
print('total '+str(len(final_cell))+' cell with WES and WGS!')

chr_len={}
file=open('./chrom_len.tsv','r')
lines=file.readlines()
for line in lines:
        data=line.rstrip().split('\t')
        chr_len.update({data[0]:data[1]})#chr chr_len

#..read fiels......
#vcf giles
patientname='P01'
path1='./1-WES-mutation'
#bam files
path2='./0-Bam-files'
#varbin result
path3='./varbin-result'

ll_header=['mutation_id','ref_counts','var_counts','normal_cn','minor_cn','major_cn']
#read total mutation list
cell_list={}#cell:[bam file, varbin result]
all_mut_list=[]
file=open(path1+'/WES-mutation.maf','r')
lines=file.readlines()
for line in lines[1:]:
        data=line.rstrip().split('\t')
        mut='_'.join(data[0:4])+'_'+data[6]+'_'+data[7] #gene_chr_start_end_ref_alt
        cell=data[-1]
        if (cell.split('_')[0]==patientname) and (cell in final_cell):
                if cell not in cell_list.keys():
                        cell_list.update({cell:[]})
                if mut not in all_mut_list:
                        all_mut_list.append(mut)
print('total '+str(len(cell_list.keys()))+' cells in patient '+patientname)
print('total '+str(len(all_mut_list))+' mutations in patient '+patientname)
#require paired bam files & varbin data
bam1_list=os.listdir(path2)
for filename in bam1_list:
        cell=filename.split('.')[0]
        if cell in cell_list.keys():
                if filename[-4:]=='.bam':
                        cell_list[cell].append(filename)
bam1_list=os.listdir(path3)
for filename in bam1_list:
        cell=filename.split('.')[0]
        if cell in cell_list.keys():
                name=cell+'.hg38.100k.k100.nobad.varbin.data.txt'
                if name not in cell_list[cell]:
                        cell_list[cell].append(name)
for key in cell_list.keys():
        if len(cell_list[key])!=2:
                print('error')
                print(key)
                print(cell_list[key])
#read counts
for filename in cell_list.keys():
        two_base=[]
        bam1=cell_list[filename][0]
        mut_list={}
        for mut in all_mut_list:
                if mut not in mut_list.keys():
                        mut_list.update({mut:[]})
                else:
                        print('error')
        for mut in mut_list.keys():
                print(mut)
                if ('-' in mut.split('_')[-1]) or ('-' in mut.split('_')[-2]):
                        #a indel mutation
                        if '-'==mut.split('_')[-1]:#deletion
                                print('deletion')
                                region=str(mut.split('_')[1])+':'+str(int(mut.split('_')[2])-1)+'-'+str(int(mut.split('_')[2])-1)
                                print('region: '+str(region))
                                del_len=0-(int(mut.split('_')[3])-int(mut.split('_')[2])+1)
                                print('del_len: '+str(del_len))
                                input_file=path2+'/'+bam1
                                (status, output) = subprocess.getstatusoutput('samtools mpileup -r %s -f Homo_sapiens_assembly38.fasta --ff DUP -B -Q 13 %s' % (region,input_file))
                                print(output)
                                if int(status)==0 and len(output.split('\n'))==2:
                                        info=output.split('\n')[-1]
                                        data3=info[:-1].split('\t')
                                        all_count=data3[3]
                                        del_num=len(data3[4].split(str(del_len)))-1
                                        ref_count=int(all_count)-del_num
                                        mut_list[mut].append(str(ref_count))
                                        mut_list[mut].append(str(del_num))
                                        print('ref counts: '+str(ref_count))
                                        print('alt counts: '+str(del_num))
                                else:
                                        print('indel error3')
                                        mut_list[mut].append('0')
                                        mut_list[mut].append('0')
                        else:
                                #insertion
                                print('insertion')
                                region=str(mut.split('_')[1])+':'+str(int(mut.split('_')[2]))+'-'+str(int(mut.split('_')[2]))
                                print('region: '+str(region))
                                insert_len=len(mut.split('_')[5])
                                print('insert_len: '+str(insert_len))
                                input_file=path2+'/'+bam1
                                (status, output) = subprocess.getstatusoutput('samtools mpileup -r %s -f Homo_sapiens_assembly38.fasta --ff DUP -B -Q 13 %s' % (region,input_file))
                                print(output)
                                if int(status)==0 and len(output.split('\n'))==2:
                                        info=output.split('\n')[-1]
                                        data3=info[:-1].split('\t')
                                        all_count=data3[3]
                                        del_num=len(data3[4].split('+'+str(insert_len)))-1
                                        ref_count=int(all_count)-del_num
                                        mut_list[mut].append(str(ref_count))
                                        mut_list[mut].append(str(del_num))
                                        print('ref counts: '+str(ref_count))
                                        print('alt counts: '+str(del_num))
                                else:
                                        print('indel error3')
                                        mut_list[mut].append('0')
                                        mut_list[mut].append('0')

                else:
                        #point mutation
                        if len(mut.split('_')[-1])==1 and len(mut.split('_')[-2])==1:
                                #single base
                                print('SNVs')
                                region=str(mut.split('_')[1])+':'+str(int(mut.split('_')[2]))+'-'+str(int(mut.split('_')[2]))
                                input_file=path2+'/'+bam1
                                (status, output) = subprocess.getstatusoutput('samtools mpileup -r %s -f Homo_sapiens_assembly38.fasta --ff DUP -B -Q 13 %s' % (region,input_file))
                                #print(status)
                                print(output)
                                if int(status)==0 and len(output.split('\n'))==2:
                                        info=output.split('\n')[-1]
                                        data3=info[:-1].split('\t')
                                        all_count=data3[3]
                                        map_result=data3[4]
                                        mismatch=re.findall(r'[ATCGNatcgn]',map_result)
                                        indel_insert=re.findall(r'\+[0-9]+[ATCGNatcgn]+',map_result)
                                        indel_del=re.findall(r'\-[0-9]+[ATCGNatcgn]+',map_result)
                                        start_alt=re.findall(r'\^+[ATCGNatcgn]',map_result)
                                        indel_base=[]
                                        for indel1 in indel_insert:
                                                length=int(re.findall(r"\d+\.?\d*", indel1)[0])
                                                if length<=len(indel1)-1-len(str(length)):
                                                        for i in range(1+len(str(length)),length+1+len(str(length))):
                                                                indel_base.append(re.findall(r'[ATCGNatcgn]',indel1)[i-1-len(str(length))])

                                                else:
                                                        for i in range(1+len(str(length)),len(indel1)):
                                                                indel_base.append(re.findall(r'[ATCGNatcgn]',indel1)[i-1-len(str(length))])
                                        for indel1 in indel_del:
                                                length=int(re.findall(r"\d+\.?\d*", indel1)[0])
                                                if length<=len(indel1)-1-len(str(length)):
                                                        for i in range(1+len(str(length)),length+1+len(str(length))):
                                                                indel_base.append(re.findall(r'[ATCGNatcgn]',indel1)[i-1-len(str(length))])
                                                else:
                                                        for i in range(1+len(str(length)),len(indel1)):
                                                                indel_base.append(re.findall(r'[ATCGNatcgn]',indel1)[i-1-len(str(length))])
                                        start=[]
                                        for indel1 in start_alt:
                                                for indel_item in re.findall(r'[ATCGNatcgn]',indel1):
                                                        start.append(indel_item)
                                        specific_base=[]
                                        for base in mismatch:
                                                for i in range(0,len(base)):
                                                        specific_base.append(base[i])
                                        for base in start:
                                                if base in specific_base:
                                                        specific_base.remove(base)
                                        for base in indel_base:
                                                specific_base.remove(base)
                                        specific_count=0
                                        for base in specific_base:
                                                if base.upper()==mut.split('_')[-1] or base.lower()==mut.split('_')[-1] :
                                                        specific_count=specific_count+1
                                        ref_count=int(all_count)-specific_count
                                        mut_list[mut].append(str(ref_count))
                                        mut_list[mut].append(str(specific_count))
                                        print('ref counts: '+str(ref_count))
                                        print('alt counts: '+str(specific_count))
                                else:
                                        print('snv error3')
                                        mut_list[mut].append('0')
                                        mut_list[mut].append('0')
                        else:
                                #two base or more base mutation.
                                two_base.append(1)
        print('two base mutation: '+str(len(two_base)))

        #extract CN info
        varbin_file=cell_list[filename][1]
        dic_varbin={}
        file=open(path3+'/'+varbin_file,'r')
        lines=file.readlines()
        for i in range(1,len(lines)-1):
                region_start='.'
                region_end='.'
                line1=lines[i].rstrip().split('\t')
                chr1='.'
                if line1[0] not in ['23','24']:
                        chr1='chr'+line1[0]
                else:
                        if line1[0]=='23':
                                chr1='chrX'
                        if line1[0]=='24':
                                chr1='chrY'
                line2=lines[i+1].rstrip().split('\t')
                chr2='.'
                if line2[0] not in ['23','24']:
                        chr2='chr'+line2[0]
                else:
                        if line2[0]=='23':
                                chr2='chrX'
                        if line2[0]=='24':
                                chr2='chrY'
                if chr1==chr2:
                        region_start=str(line1[1])
                        region_end=str(int(line2[1])-1)
                else:
                        region_start=str(line1[1])
                        region_end=str(chr_len[chr1])
                CN=round(float(line1[-1])*2)
                if region_start!='.' and region_end!='.':
                        region_id=chr1+'_'+region_start+'_'+region_end
                        if region_id not in dic_varbin:
                                dic_varbin.update({region_id:CN})
                        else:
                                print('error_cn')
                else:
                        print('error_cn')
        line=lines[-1].rstrip().split('\t')
        chr3='.'
        region_start='.'
        region_end='.'
        CN=round(float(line[-1])*2)
        if line[0] not in ['23','24']:
                chr3='chr'+line[0]
        else:
                if line[0]=='23':
                        chr3='chrX'
                if line[0]=='24':
                        chr3='chrY'
        if chr3!='.':
                region_start=str(line[1])
                region_end=str(chr_len[chr3])
                if region_start!='.' and region_end!='.':
                        region_id=chr3+'_'+region_start+'_'+region_end
                        if region_id not in dic_varbin:
                                dic_varbin.update({region_id:CN})
                        else:
                                print('error_cn')
                else:
                        print('error_cn')
        else:
                print('error_cn')

        for mut in mut_list.keys():
                chr_id=mut.split('_')[1]
                start=mut.split('_')[2]
                end=mut.split('_')[3]
                CN_value='.'
                for key in dic_varbin.keys():
                        if chr_id==key.split('_')[0] and int(key.split('_')[1])<=int(start)<=int(end)<=int(key.split('_')[2]):
                                CN_value=dic_varbin[key]
                if CN_value!='.':
                        mut_list[mut].append('2')
                        mut_list[mut].append('0')
                        mut_list[mut].append(str(CN_value))
                else:
                        print('error: without CN info')


        for mut in mut_list.keys():
                if len(mut_list[mut])!=5:
                        print(mut+': error-without 5 info')


        out=open('./6-pyclone/1-input/'+filename+'.txt','w')
        out.write('\t'.join(ll_header))
        for mut in mut_list.keys():
                if len(mut_list[mut])==5:
                        out.write('\n'+mut+'\t'+'\t'.join(mut_list[mut]))
        out.close()
        print(filename+': ok')

#2. ......................................run pyclone......................................
PyClone run_analysis_pipeline --in_files \
Example_A_5.txt \
Example_A_9.txt \
Example_BC_14.txt \
Example_BC_9.txt \
--working_dir Example \
--prior total_copy_number \
--plot_file_format pdf \
--tumour_contents \
1 \
0.8 \
0.7 \
0.99 


#3. ......................................calculate merged CCF......................................
import re
import os
#depth info
dic_depth={}
path='./6-pyclone/1-input'
filelist=os.listdir(path)
for filename in filelist:
        cell=filename.split('.')[0]
        if cell not in dic_depth.keys():
                dic_depth.update({cell:[]})#[mut,depth]
        file=open(path+'/'+filename,'r')
        lines=file.readlines()
        for line in lines[1:]:
                data=line.rstrip().split('\t')
                depth=int(data[1])+int(data[2])
                dic_depth[cell].append([data[0],depth])

dic_info={}
dic_mutation={}#mutation list per patient
path='./6-pyclone/3-pyclone-output'
dirlist=os.listdir(path)
for patient in dirlist:
        if patient not in dic_mutation.keys():
                dic_mutation.update({patient:[]})
        file=open(path+'/'+patient+'/tables/loci.tsv','r')
        lines=file.readlines()
        for line in lines[1:]:
                data=line.rstrip().split('\t')
                mut=data[0]
                cell=data[1].split('.')[0]
                CCF=data[3]
                if cell not in dic_info.keys():
                        dic_info.update({cell:[]})
                depth='.'
                for item in dic_depth[cell]:
                        if item[0]==mut:
                                depth=item[1]
                if depth!='.':
                        ll=[mut,CCF,depth]
                        dic_info[cell].append(ll)
                else:
                        print('without depth info!')
                if mut not in dic_mutation[patient]:
                        dic_mutation[patient].append(mut)

for patient in dirlist:
        out=open('./6-pyclone/4-merged-CCF/'+patient+'-merged.txt','w')
        tumor_list=[]
        for cell in dic_info.keys():
                if cell.split('_')[0]==patient:
                        tumor='_'.join(cell.split('_')[0:2])
                        if tumor not in tumor_list:
                                tumor_list.append(tumor)
        header=['mutation_mergedCCF']
        for tumor in tumor_list:
                header.append(tumor)
        out.write('\t'.join(header))
        for mutation in dic_mutation[patient]:
                header=[mutation]
                for tumor in tumor_list:
                        num1=0 #sum of CCF*depth
                        num2=0 #sum of depth
                        for cell in dic_info.keys():
                                if '_'.join(cell.split('_')[0:2])==tumor:
                                        for item in dic_info[cell]:
                                                if item[0]==mutation:
                                                        CCF=float(item[1])
                                                        depth=float(item[2])
                                                        num1=num1+(CCF*depth)
                                                        num2=num2+depth
                        merged_CCF=num1/num2
                        if merged_CCF<1:
                                header.append(str(merged_CCF))
                        else:
                                header.append('1')
                out.write('\n'+'\t'.join(header))
        out.close()








# script to calculate OR value
import re
import os
import scipy.stats as stats
import subprocess

driver_gene=[]
file=open('./Breast_cancer_driver_genes.txt','r')
lines=file.readlines()
for line in lines[1:]:
        data=line.rstrip().split('\t')
        if data[0] not in driver_gene:
                driver_gene.append(data[0])
print('total driver genes in breast cancer: '+str(len(driver_gene)))

path='./1-WES-mutation'
func_list=['Frame_Shift_Del','In_Frame_Del','Missense_Mutation','Nonsense_Mutation','Splice_Site']
dic_sample={}
file=open(path+'/WES-mutation.maf','r')
lines=file.readlines()
for line in lines[1:]:
        data=line.rstrip().split('\t')
        gene=data[0]
        func=data[4]
        sample=data[-1]
        if ('_BC' in sample) or ('_A' in sample):
                if sample not in dic_sample.keys():
                        dic_sample.update({sample:[]})
                if (func in func_list) and (gene in driver_gene):
                        if gene not in dic_sample[sample]:
                                dic_sample[sample].append(gene)

dic_patient={}
for key in dic_sample.keys():
        patient_id=key.split('_')[0]
        tumor_id=key.split('_')[0]+'_'+key.split('_')[1]
        if patient_id not in dic_patient.keys():
                dic_patient.update({patient_id:[[],[]]})
        if '_BC' in tumor_id:
                if tumor_id not in dic_patient[patient_id][0]:
                        dic_patient[patient_id][0].append(tumor_id)
        if '_A' in tumor_id:
                if tumor_id not in dic_patient[patient_id][1]:
                        dic_patient[patient_id][1].append(tumor_id)

for patient in dic_patient.keys():
        if len(dic_patient[patient][0])==1 and len(dic_patient[patient][1])==1:
                print(patient)
                #all mutated driver gene
                mutated_driver_gene=[]
                for key in dic_sample.keys():
                        if key.split('_')[0] ==patient:
                                for gene in dic_sample[key]:
                                        if gene not in mutated_driver_gene:
                                                mutated_driver_gene.append(gene)

                #calculate OR for each mutated driver genes (per patient)
                out=open('./10-OR-value/1-mutated-driver-gene/'+patient+'-OR-for-mutatedDriverGene.txt','w')
                header=['mutated_driver_gene','N.cell with mut in PT','N.cell without mut in PT','N.cell with mut in A','N.cell without mut in A','Pvalue','oddRatio']
                out.write('\t'.join(header)+'\n')
                for gene in mutated_driver_gene:
                        header=[gene]
                        PT=dic_patient[patient][0][0]
                        MT=dic_patient[patient][1][0]
                        value1=0 #number of cell cluster with mutated driver gene in PT
                        value2=0 #number of cell cluster with mutated driver gene in MT
                        value3=0 #number of cell cluster without mutated driver gene in PT
                        value4=0 #number of cell cluster without mutated driver gene in MT
                        for sample in dic_sample.keys():
                                if '_'.join(sample.split('_')[0:2])==PT:
                                        if gene in dic_sample[sample]:
                                                value1=value1+1
                                        else:
                                                value3=value3+1
                                if '_'.join(sample.split('_')[0:2])==MT:
                                        if gene in dic_sample[sample]:
                                                value2=value2+1
                                        else:
                                                value4=value4+1
                        data=[value1,value3,value2,value4]
                        print(data)
                        r_script_path='./fisher-test.R'
                        r_command = ["Rscript", r_script_path]
                        r_command.extend(map(str, data))
                        result=subprocess.run(r_command, stdout=subprocess.PIPE)
                        print(result.stdout)
                        #print(result.stdout.decode())
                        #print(result.stdout.decode().split('\n')[0].split(' '))
                        ll=[]
                        for item in result.stdout.decode().split('\n')[0].split(' ')[1:]:
                                if item!='':
                                        ll.append(item)
                        p_value=ll[0]
                        odd_ratio=ll[1]
                        print('pvalue: '+str(p_value))
                        print('odd_ratio: '+str(odd_ratio))
                        header.append(str(value1))
                        header.append(str(value3))
                        header.append(str(value2))
                        header.append(str(value4))
                        header.append(p_value)
                        header.append(odd_ratio)
                        if len(header)==7:
                                out.write('\t'.join(header)+'\n')
                        else:
                                print('error')
                out.close()
