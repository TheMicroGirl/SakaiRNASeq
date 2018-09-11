import re

ChromosomeLocation = "Chromosome.fasta"
pOsak1Location = "pOSAk1.fasta"
pO157Location = "pO157.fasta"

#Formatting the E.coli O157:H7 Chromosome
sakaiChromosome = open(ChromosomeLocation,'r')
scFile= open('SakaiChromosomeHeader.fasta', 'w')
for scline in sakaiChromosome:
	if scline[0]==">":
		s = ">"+re.search('ECs(\d{4})', scline).group(0)+"\n" #search function, look for  ECs which is followed by exaclty 4 digits, () defines a group
		scFile.write(s)
	else: scFile.write(scline)
print "Chromosome is done"
sakaiChromosome.close()
scFile.close()

#Formatting the pOSAK1 plasmid
pOsak1File = open(pOsak1Location,'r')
poFile = open('pOsak1Header.fasta', 'w')
for poLine in pOsak1File:
	if poLine[0]==">":
		try:
			p = ">"+re.search('gene=([a-zA-Z]{4})', poLine).group(1)+"\n"			#search function, look for  gene = and then take only letters in the ()
			poFile.write(p)
		except:
			poFile.write(poLine)
	else: poFile.write(poLine)
print "pOSAK1 is done"
pOsak1File.close()
poFile.close()


#Formatting the pO157 plasmid
pO157Location  = open(pO157Location,'r')
plasmidOFile= open('pO157Header.fasta', 'w')
for plasmidOLine in pO157Location:
	if plasmidOLine[0]==">":
		try:
			o = ">"+re.search('gene=([a-zA-Z]{4})', plasmidOLine).group(1)+"\n" #search function, look for gene = and then take only letters in the ()
			plasmidOFile.write(o)
		except:
			plasmidOFile.write(plasmidOLine)
	else: plasmidOFile.write(plasmidOLine)
print "pO157 is done"
pO157Location.close()
plasmidOFile.close()

#Combine the 3 files into an index file which can be used to create a blast database
openChromosome= open('SakaiChromosomeHeader.fasta', 'r')
openPO157= open('pO157Header.fasta', 'r')
openPOsak1= open('pOsak1Header.fasta', 'r')

blastDatabaseTable = open('E_coli_Sakai_Genome.fasta','w') 

for write1 in openChromosome:
		blastDatabaseTable.write(write1)

for write2 in openPO157:
		blastDatabaseTable.write(write2)

for write3 in openPOsak1:
		blastDatabaseTable.write(write3)

openChromosome.close()
openPO157.close()
openPOsak1.close()
blastDatabaseTable.close()

print "Shit is done. So done."



