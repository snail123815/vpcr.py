import sys
import re
import glob
from Bio import SeqIO
import subprocess


#Theo Allnutt, 2016
#virtual PCR to find amplicons in genome assemblies
#usage:
#python vpcr.py primers.txt folder_of_assemblies_/ hits.outfile lengths.outfile 20000 2
#primer file is in format: name<tab>forward<tab>reverse<return>
#with new line for each primer pair
# hits is the name of file with summary of all hits found
# lengths is a tab table of all lengths found
#20000=max length of amplicon
#2=number of mismatches
#requires tntblast to be installed in ~/bin/tntblast-2.01/bin/tntblast or edit line 58 to change location

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))
				  

folder = sys.argv[2] 

filelist=glob.glob(folder+"/*")

filelist.sort(key=tokenize)

primers = sys.argv[1]

g=open(sys.argv[3],'w') #sequence hits fasta

h=open(sys.argv[4],'w') #sequence hit lengths summary

maxlen = sys.argv[5]

mismatches = sys.argv[6]

#get mlva names
primerfile=open(primers,'r')
primernames=[]
for i in primerfile:
	k=i.split("\t")
	primernames.append(k[0])

primernames.sort()
h.write("\t"+"\t".join(str(p)for p in primernames)+"\n")


for f in filelist:
	data={}
	
	p1= subprocess.Popen("~/bin/tntblast-2.01/bin/tntblast  -i %s -d %s -o temp.fasta -e 45 -m 1 -L T -l %s --primer-clamp 3 --max-mismatch %s " %(primers,f,maxlen,mismatches), shell=True).wait() 

#output format= fastafile num hits /n mlvaname hitdesc hitlength	hitseq
	fastaname=f.split("/")[-1]
	c=0
	t=open('temp.fasta','r')
	#g.write(fastaname+"\n")
	h.write(fastaname+"\t")
	for x in SeqIO.parse(t,'fasta'):
		line=x.description.split(" ")
		id=line[-1].rstrip("\n")
		c=c+1
		g.write(">"+fastaname+"-"+str(c)+"\t"+str(id)+" "+str(line[:-1])+" "+str(len(x.seq))+"\n"+str(x.seq)+"\n")
	
		if id not in data.keys():
			data[id]=[]
			data[id].append(len(x.seq))
		else:
			data[id].append(len(x.seq))
	print data
	print primernames
	for i in primernames:
		if i in data.keys():
			h.write(" / ".join(str(p)for p in data[i])+"\t")
		else:
			h.write("0\t")
	h.write("\n")
	
p2=subprocess.Popen("rm temp.fasta",shell=True).wait

	
	
