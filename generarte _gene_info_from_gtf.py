import random
import os
import linecache
import sys
import errno
import re
#path of hg19 chromsomes(file per each) you can download them from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/ but be cation here we assume each in this script each chromosome file contins only one line for the whole sequence while the files of this sie contains fatsa format
#or you can run the script "chrs_from_genome.py" to create the chromomsme files
chrs_path = "C:\Users\QCRI-L002\Desktop\ucsc_hg19_chrs\\"

#path of the hg19 gtf file and you can download it from https://usegalaxy.org/library_common/download_dataset_from_folder?library_id=4ab3a886a95d362e&show_deleted=False&cntrller=library&use_panels=False&id=05ccf9e20303b392
gtf_path = "C:\Users\QCRI-L002\Desktop\hg19_genes.gtf"
#path of the gene list. We can exract it from the htseq-count file
gene_path = "C:\Users\QCRI-L002\Desktop\genes_from_htseq.txt"


#a list of the number of genes usch that each element is an array of exons for this gene
group_by_gene={}

#an array of genes from htseq-count file 
genes_htseq = []

#an aray of genes grom gtf
genes_gtf = []
#a list of the length of all exons for a gene
exons_len_for_a_gene ={}
#a list of the length of all exons for a gene
exons_for_a_gene ={}
#a list of all chromosomes
chromosomes = {}
#a list of the gc content of all exons for a gene
exons_gc_for_a_gene ={}


#read files from the specific path
def list_files(path):
    # returns a list of names (with extension, without full path) of all files 
    # in folder path
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(name)
    return files

#return the gc content for a specific DNA string
def GC(s):
    No_c=0
    No_g=0
    No_a=0
    No_t=0
    for i in range(len(s)):
        if s[i]=='A' or s[i]=='a':
            No_a = No_a+1
        elif s[i]=='C' or s[i]=='c':
            No_c = No_c+1
        elif s[i]=='G' or s[i]=='g':
            No_g = No_g+1
        elif s[i]=='T' or s[i]=='t':
            No_t = No_t+1
    No_gc = float(No_g+No_c)
    No_all = float(No_g+No_c+No_t+No_a)
    return (No_gc/No_all)

#read the gtf file and map each exon to its gene
def read_gene_list_file():
    f= open(gene_path,'r')
    for line in f:
        s=line.strip()
        genes_htseq.append(s)
    f.close()
        
def read_ucsc_gtf_file():
    f= open(gtf_path,'r')
    for line in f:
        s=line.strip()
        if(s.find("exon")>-1):
            s=s.split()
            vc ={}
            vc["chr"]=s[0]
            vc["start"]=int(s[3])
            vc["end"]=int(s[4])
            vc["direction"]=s[6]
            
            try:
                index = s.index("gene_id")
                temp = s[index+1].split("\"")
                try:
                    le = len(group_by_gene[temp[1]])
                    group_by_gene[temp[1]].append(vc)
                except:
                    genes_gtf.append(temp[1])
                    group_by_gene[temp[1]]=[]
                    group_by_gene[temp[1]].append(vc)
            except:
                try:
                    index = s.index("gene_name")
                    temp = s[index+1].split("\"")
                    try:
                        le = len(group_by_gene[temp[1]])
                        group_by_gene[temp[1]].append(vc)
                    except:
                        genes_gtf.append(temp[1])
                        group_by_gene[temp[1]]=[]
                        group_by_gene[temp[1]].append(vc)
                except:
                    print "Error: The line \n"+line+"does not conatin any gene_id or gene_name\n"+s[8]
    f.close()

# In this function we will reomve the repeated exons that are coming from reading the gtf file beacuse it contains the same exon multiple time one for each transcript
def remove_repeated_exons_for_each_gene():
    for x in genes_gtf:
        temp = []
        for y in group_by_gene[x]:
            if y not in temp:
                temp.append(y)
        group_by_gene[x]=temp

def calculate_the_length_of_all_exons_in_gene():
    for x in genes_gtf:
        exons_len_for_a_gene[x]=0
        for y in group_by_gene[x]:
            exons_len_for_a_gene[x]=exons_len_for_a_gene[x]+(y["end"]-y["start"]+1)
def concantenate_the_all_exons_in_a_gene():
    for x in genes_gtf:
        exons_for_a_gene[x]=""
        for y in group_by_gene[x]:
            exons_for_a_gene[x]=exons_for_a_gene[x]+chromosomes[y["chr"]][y["start"]-1:y["end"]]
        if(exons_len_for_a_gene[x]==len(exons_for_a_gene[x])):
           pass
        else:
           print "error in equality of gene length for gen  "+x+"\n"
def calculate_the_length_of_all_exons_in_gene():
    for x in genes_gtf:
        exons_len_for_a_gene[x]=0
        for y in group_by_gene[x]:
            exons_len_for_a_gene[x]=exons_len_for_a_gene[x]+(y["end"]-y["start"]+1)

def calculate_the_gc_content_of_all_exons_in_gene():
    for x in genes_gtf:
        exons_gc_for_a_gene[x]=GC(exons_for_a_gene[x])
                
    

#file the chromosomes list with the sequence of the chromosomes, such that each chromosome in the path of chromosomes and exist in the file in one line
def raed_chromosomes():
    list_chrs = list_files(chrs_path)
    for x in list_chrs:
        name=x[:x.find(".txt")]
        #print name
        f_1=open(chrs_path+x,'r')
        seq=f_1.read()
        chromosomes[name]=seq
        #print chromsomes[name][1:5]
        f_1.close()
def write_gene_info(fName):
    f=open(fName,'w')
    f.write("gene_name\texons_length\texons_GC_content\tNo_of_exons\n")
    for x in genes_htseq:
        f.write(x+"\t"+str(exons_len_for_a_gene[x])+"\t"+str(exons_gc_for_a_gene[x])+"\t"+str(len(group_by_gene[x]))+"\n")
    f.close()
        
if __name__ == "__main__":
    read_gene_list_file();
    read_ucsc_gtf_file();
    remove_repeated_exons_for_each_gene();
    calculate_the_length_of_all_exons_in_gene();
    raed_chromosomes();
    concantenate_the_all_exons_in_a_gene();
    calculate_the_gc_content_of_all_exons_in_gene();
    write_gene_info("gene_Info.csv")
       
    """
    print len(group_by_gene)
    print len(genes_htseq)
    print len(genes_gtf)
    print len(group_by_gene["BRCA1"])
    remove_repeated_exons_for_each_gene()
    print len(group_by_gene["BRCA1"])
    """
    
            





