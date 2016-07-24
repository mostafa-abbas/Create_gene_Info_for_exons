#This script write the chromsomes in files such thst each file contains only one line; given the genome file in fatsa format
import os
genome_path = "ucsc.hg19.fasta"
chrs_directory_name = "ucsc_hg19_chrs"
if __name__ == "__main__":
    if not os.path.exists(chrs_directory_name):
        os.makedirs(chrs_directory_name)
    file = open(genome_path,'r')
    k=0
    f=open("test.txt",'w')
    for line in file:
        s=line.replace("\n",'')
        if(s.find(">")>-1):
            f.close()
            f=open(chrs_directory_name+"\\"+s[1:]+".txt",'w')
            print s
        else:
            f.write(s)
            pass
    f.close()
    file.close()
    


    
