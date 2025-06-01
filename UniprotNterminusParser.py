import re
import csv
import sys
import os
import string
import traceback
import argparse
import Bio
import itertools
from Bio.Seq import Seq

print("USAGE:  In file, -f, --fin, Uniprot file to scan.  The script scans all Uniprot entries provided, which can include the entire database or sections, for the N-terminal starting sequence specified (don't just type M).  It also counts DNA, RNA and Nucleotide interactions by default output as a grid, although many proteins are not annotated.  If desired, any GO, Interpro or Pfam identifiers for each protein identified can be generated. Specify name header for output files with the -o or --fout option, 2 or 5 files will be generated based on preferences, --reg, -r the N terminal amino acids searched for in capital 1 letter code, must include all starting with M, in blocks of 10 tab spaced based on Uniprot fasta format, if --disc or -d, yes or Y you will have pfam, GO and Interprot discriptors for each protein in seperate files.  You have to specify --disc, -d as yes or no, it's a simple text parser with switches.")
arg_parser = argparse.ArgumentParser( description = "Convert vcf file to fasta." )
arg_parser.add_argument( "--fin", "-f", type=str, required=True )
arg_parser.add_argument( "--fout", "-o", type=str, required=True)
arg_parser.add_argument( "--reg", "-r", type=str, required=True)
arg_parser.add_argument( "--disc", "-d", type=str, required=True)
arguments = arg_parser.parse_args()

filename = arguments.fin
filename2 = (arguments.fout+".fasta")
filename3 = (arguments.fout+"_Pfamtype.tsv")
filename4 = (arguments.fout+"_Interprotype.tsv")
filename5 = (arguments.fout+"_counts.tsv")
filename6 = (arguments.fout+"_GOtype.tsv")
Reg = arguments.reg
Disc = arguments.disc
fileout = []
fileout2 = []
fileout3 = []
fileout4 = []
fileout5 = []


fileout = open(filename2, 'a')

fileout4 = open(filename5, 'a')

templist = []
templist2 = []
templist3 = []
templist4 = []
templist5 = []
templist6 = []
templist7 = []
templist8 = []
templist9 = []
templist10 = []
templist11 = []

x1=">"
x2="GO;"
x3="CC "
x5="AC "
x6="GN "
x7="DR "
x8="Pfam"
x9="InterPro"
x10="SQ  "
x11=("     "+Reg)
x12="//"

tl2 = []
tl3 = []
tl4 = []
tl5 = []
tl6 = []
tl7 = []
tl9 = []
tl10 = []
tl11 = []

tl2s = []
tl3s = []
tl4s = []
tl5s = []
tl6s = []
tl7s = []
tl9s = []
tl10s = []
tl11s = []

c1s = []
c2s = []
c3s = []
c4s = []
c5s = []
c6s = []
c7s = []
c8s = []
c9s = []

DNA="DNAbindingentries.tsv"
RNA="RNAbindingentries.tsv"
NUC="Nucleotidebindingentries.tsv"
DNA2="DNAInterpro.tsv"
RNA2="RNAInterpro.tsv"
NUC2="NucleotideInterpro.tsv"
DNA3="DNAInterpro2.tsv"
RNA3="RNAInterpro2.tsv"
NUC3="NucleotideInterpro2.tsv"
with open(filename, 'r') as filen:
	for linen, line3 in enumerate(filen):

		if line3.find(x5) != -1:
			line3_1 = line3.split("   ")
			p11 = line3_1[1]
			p111 = p11.replace(';','')
			p122 = p111.replace('\n',' ')
			p123 = p122.split(" ")
			p1 = p123[0]
		if line3.find(x6) != -1: 
			if line3.find("NAME=") != -1:
				lin3_1 = line3.split("   ")
				p222 = lin3_1[1]
				p22 = p222.split("=")
				p223 = p22[1]
				p233 = p223.replace('\n',' ')
				p2 = p233.replace(';','')
			else:
				lin3_1 = line3.split("   ")
				p22 = line3_1[1]
				p2 = p22.replace('\n',' ')
		if Disc == "yes" or Disc == "Y":
			c1 = 0
			c2 = 0
			c3 = 0
			c4 = 0
			c5 = 0
			c6 = 0
			c7 = 0
			c8 = 0
			c9 = 0
			if line3.find(x7) != -1:
				if line3.find(x2) != -1:
					
					linen4 = linen
					for linen4, linena in enumerate(filen):
						if linena.find(x2) != -1:
							filep7 = open(DNA2, 'r');
							
							for linep7 in filep7:
								linep7_split = linep7.split("\t")
								pf7 = linep7_split[6]
								pf7r = pf7.replace('\n','')
								pf7h = linep7_split[1]
								pf7_split = pf7r.split(",")
								for i in range(0,len(pf7_split)):	
									i1 = pf7_split[i]
									if linena.find(i1) != -1:
										tl9 = (p1,str(i1),str(pf7h))
										tl9_l = '\t'.join(str (v) for v in tl9)
										tl9s.append(tl9_l)
				
									else:
										continue
					
							filep7.close()
							filep8 = open(RNA2, 'r');
							for linep8 in filep8:
								linep8_split = linep8.split("\t")
								pf8 = linep8_split[6]
								pf8r = pf8.replace('\n','')
								pf8h = linep8_split[1]
								pf8_split = pf8r.split(",")
								for j in range(0,len(pf8_split)):
							
									j2 = pf8_split[j]

									if linena.find(j2) != -1:
										tl10 = (p1,str(j2),str(pf8h))
										tl10_l = '\t'.join(str (v) for v in tl10)
										tl10s.append(tl10_l)

									else:
										continue
						
							filep8.close()
							filep9 = open(NUC2, 'r');
							for linep9 in filep9:
								linep9_split = linep9.split("\t")
								pf9 = linep9_split[6]
								pf9r = pf9.replace('\n','')
								pf9h = linep9_split[1]
								pf9_split = pf9r.split(",")
								for iii in range(0,len(pf9_split)):
									iii3 = pf9_split[iii]
									if linena.find(iii3) != -1:
										tl11 = (p1,str(iii3),str(pf9h))
										tl11_l = '\t'.join(str (v) for v in tl10)
										tl11s.append(tl11_l)

									else:
										continue
							filep9.close()
							
												
						else:
							
							break
					else:
						continue
				if line3.find(x9) != -1:
					
					linen6 = (linen-1)			
					for linen6, linenc in enumerate(filen):
						if linenc.find(x9) != -1:
							filep4 = open(DNA3, 'r');
							
							for linep4 in filep4:
								linep4_split = linep4.split("\t")
								pf4 = linep4_split[0]
								pf44 = linep4_split[1]
								if linenc.find(pf4) != -1:
									tl5 = (p1,str(pf4),str(pf44))
									tl5_l = '\t'.join(str (v) for v in tl5)
									tl5s.append(tl5_l)
									
								else:
									continue
							filep4.close()
							filep5 = open(RNA3, 'r');
							for linep5 in filep5:
								linep5_split = linep5.split("\t")
								pf5 = linep5_split[0]
								pf55 = linep5_split[1]
								if linenc.find(pf5) != -1:
									tl6 = (p1,str(pf5),str(pf55))
									tl6_l = '\t'.join(str (v) for v in tl6)
									tl6s.append(tl6_l)

								else:
									continue
							filep5.close()
							filep6 = open(NUC3, 'r');
							for linep6 in filep6:
								linep6_split = linep6.split("\t")
								pf6 = linep6_split[0]
								pf66 = linep6_split[1]
								if linenc.find(pf6) != -1:
									tl7 = (p1,str(pf6),str(pf66))
									tl7_l = '\t'.join(str (v) for v in tl7)
									tl7s.append(tl7_l)


								else:
									continue
							filep6.close()
							
						else:
							
							break
					
					else:
						continue
						
						
				if line3.find(x8) != -1:
					
					linen5 = linen		
					for linen5, linenb in enumerate(filen):
						if linenb.find(x8) != -1:
							filep1 = open(DNA, 'r');
							
							for linep1 in filep1:
								linep1_split = linep1.split("\t")
								pf1 = linep1_split[0]
								pf11 = linep1_split[1]
								if linenb.find(pf1) != -1:
									tl2 = (p1,str(pf1),str(pf11))
									tl2_l = '\t'.join(str (v) for v in tl2)
									tl2s.append(tl2_l)
									
												
								else:
									continue
							filep1.close()
							filep2 = open(RNA, 'r');
							for linep2 in filep2:
								linep2_split = linep2.split("\t")
								pf2 = linep2_split[0]
								pf22 = linep2_split[1]
								if linenb.find(pf2) != -1:
									tl3 = (p1,str(pf2),str(pf22))
									tl3_l = '\t'.join(str (v) for v in tl3)
									tl3s.append(tl3_l)
									
												
								else:
									continue
							filep2.close()
							filep3 = open(NUC, 'r');
							for linep3 in filep3:
								linep3_split = linep3.split("\t")
								pf3 = linep3_split[0]
								pf33 = linep3_split[1]
								if linenb.find(pf3) != -1:
									tl4 = (p1,str(pf3),str(pf33))
									tl4_l = '\t'.join(str (v) for v in tl4)
									tl4s.append(tl4_l)
									
												
								else:
										continue
							filep3.close()
							
							
						else:
							
							break
					else:
						continue			
				else:
					continue
					
			if line3.find(x10) != -1:
				frstln = filen.readline(linen+1)
				firstln2 = filen.readline(linen+2)
				linenad1 = linen - 1
				if frstln.find(x11) != -1:
					#templist8.append(str(c123)+str(c456)+str(c789))	
					templist.append("> "+p1+"   "+"GName= "+p2+" "+line3)
					templist.append(frstln)
					templist.append(firstln2)
					if (len(tl2s) != 0):
						tl2_ll = '\n'.join(str (v) for v in tl2s)
						templist2.append(tl2_ll)
						c1 = 1
					if (len(tl3s) != 0):
						tl3_ll = '\n'.join(str (v) for v in tl3s)
						templist3.append(tl3_ll)
						c2 = 1
					if (len(tl4s) != 0):
						tl4_ll = '\n'.join(str (v) for v in tl4s)
						templist4.append(tl4_ll)
						c3 = 1
					if (len(tl5s) != 0):
						tl5_ll = '\n'.join(str (v) for v in tl5s)
						templist5.append(tl5_ll)
						c4 = 1
					if (len(tl6s)!= 0):
						tl6_ll = '\n'.join(str (v) for v in tl6s)
						templist6.append(tl6_ll)
						c5 = 1
					if (len(tl7s) != 0):
						tl7_ll = '\n'.join(str (v) for v in tl7s)
						templist7.append(tl7_ll)
						c6 = 1
					if (len(tl9s) != 0):
						c7 = 1
						tl9_ll = '\n'.join(str (v) for v in tl9s)
						templist9.append(tl9_ll)
					if (len(tl10s) != 0):
						c8 = 1
						tl10_ll = '\n'.join(str (v) for v in tl10s)
						templist10.append(tl10_ll)
					if (len(tl11s) != 0):
						c9 = 1
						tl11_ll = '\n'.join(str (v) for v in tl11s)
						templist11.append(tl11_ll)
					if firstln2.find(x12) == -1:
						for linenad1, lineg in enumerate(filen):
							if lineg.find(x12) == -1:
								templist.append(lineg)
							elif lineg.find(x12) != -1:
								templist8.append(p1+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(c3)+"\t"+str(c4)+"\t"+str(c5)+"\t"+str(c6)+"\t"+str(c7)+"\t"+str(c8)+"\t"+str(c9)+"\n")
								break
					if firstln2.find(x12) != -1:
						templist8.append(p1+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(c3)+"\t"+str(c4)+"\t"+str(c5)+"\t"+str(c6)+"\t"+str(c7)+"\t"+str(c8)+"\t"+str(c9)+"\n")
						continue
					else:
						continue
				else:
					continue
			if line3.find(x12) != -1:
				tl2 = []
				tl3 = []
				tl4 = []
				tl5 = []
				tl6 = []
				tl7 = []
				tl9 = []
				tl10 = []
				tl11 = []
				tl2s = []
				tl3s = []
				tl4s = []
				tl5s = []
				tl6s = []
				tl7s = []
				tl9s = []
				tl10s = []
				tl11s = []
						
			else:
				continue
							

		if Disc == "no" or Disc == "N":	
			c1 = 0
			c2 = 0
			c3 = 0
			c4 = 0
			c5 = 0
			c6 = 0
			c7 = 0
			c8 = 0
			c9 = 0
			if line3.find(x2) != -1:
				
				filep4 = open(DNA2, 'r');
				for linep4 in filep4:
					linep4_split = linep4.split("\t")
					pf4 = linep4_split[6]
					pf4r = pf4.replace('\n','')
					pf4rr = pf4r.replace(' ','')
					pf44 = pf4rr.split(",")
					for k in range(0,len(pf44)):
						k1 = pf44[k]
						if line3.find(k1) != -1:
							c7=1
							c7s.append(str(c7))
							
							
						else:
							continue
					
				filep4.close()
				filep5 = open(RNA2, 'r');
				for linep5 in filep5:
					linep5_split = linep5.split("\t")
					pf5 = linep5_split[6]
					pf5r = pf5.replace('\n','')
					pf5rr = pf5r.replace(' ','')
					pf55 = pf5rr.split(",")
					for kk in range(0,len(pf55)):
						kk2 = pf55[kk]
						if line3.find(kk2) != -1:
							c8=1
							c8s.append(str(c8))
							
							
						else:
							continue
					
				filep5.close()
				filep6 = open(NUC2, 'r');
				for linep6 in filep6:
					linep6_split = linep6.split("\t")
					pf6 = linep6_split[6]
					pf6r = pf6.replace('\n','')
					pf6rr = pf6r.replace(' ','')
					pf66 = pf6r.split(",")
					for kkk in range(0,len(pf66)):
						kkk3 = pf66[kkk]
						if line3.find(kkk3) != -1:
							c9=1
							c9s.append(str(c9))
							
							
						else:
							continue
					
				filep6.close()	

				#c7s.append(str(c7))
				#c8s.append(str(c8))
				#c9s.append(str(c9))
			if line3.find(x9) != -1:
				
				filep4 = open(DNA2, 'r');
				for linep4 in filep4:
					linep4_split = linep4.split("\t")
					pf4 = linep4_split[0]
					if line3.find(pf4) != -1:
						c4=1
						c4s.append(str(c4))
						
					else:
						continue
				filep4.close()
				filep5 = open(RNA2, 'r');
				for linep5 in filep5:
					linep5_split = linep5.split("\t")
					pf5 = linep5_split[0]
					if line3.find(pf5) != -1:
						c5=1
						c5s.append(str(c5))
						
					else:
						continue
				filep5.close()
				filep6 = open(NUC2, 'r');
				
				for linep6 in filep6:
					linep6_split = linep6.split("\t")
					pf6 = linep6_split[0]
					if line3.find(pf6) != -1:
						c6=1
						c6s.append(str(c6))
						
					else:
						continue
				filep6.close()	
				#c4s.append(str(c4))
				#c5s.append(str(c5))
				#c6s.append(str(c6))
			if line3.find(x8) != -1:
				
				filep1 = open(DNA, 'r');
				for linep1 in filep1:
					linep1_split = linep1.split("\t")
					pf1 = linep1_split[0]
					if line3.find(pf1) != -1:
						c1=1
						c1s.append(str(c1))
						
					else:
						continue
				filep1.close()
				filep2 = open(RNA, 'r');
				for linep2 in filep2:
					linep2_split = linep2.split("\t")
					pf2 = linep2_split[0]
					if line3.find(pf2) != -1:
						c2=1
						c2s.append(str(c2))
						
					else:
						continue
				filep2.close()
				filep3 = open(NUC, 'r');
				for linep3 in filep3:
					linep3_split = linep3.split("\t")
					pf3 = linep3_split[0]
					if line3.find(pf3) != -1:
						c3=1
						c3s.append(str(c3))
						
					else:
						continue
				filep3.close()
				#c1s.append(str(c1))
				#c2s.append(str(c2))
				#c3s.append(str(c3))

			if line3.find(x10) != -1:
				c1 = 0
				c2 = 0
				c3 = 0
				c4 = 0
				c5 = 0
				c6 = 0
				c7 = 0
				c8 = 0
				c9 = 0
				c1ss = str(c1s)
				c2ss = str(c2s)
				c3ss = str(c3s)
				c4ss = str(c4s)
				c5ss = str(c5s)
				c6ss = str(c6s)
				c7ss = str(c7s)
				c8ss = str(c8s)
				c9ss = str(c9s)
				qt2 = "1"
				frstln = filen.readline(linen+1)
				firstln2 = filen.readline(linen+2)
				linenad1 = linen-1
				if frstln.find(x11) != -1:
					#templist8.append(p1+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(c3)+"\t"+str(c4)+"\t"+str(c5)+"\t"+str(c6)+"\t"+str(c7)+"\t"+str(c8)+"\t"+str(c9)+"\n")
					#templist8.append(p1+"\t"+str(c123l)+str(c456)+str(c789))
					templist.append("> "+p1+"   "+"GName= "+p2+" "+line3)
					templist.append(frstln)
					templist.append(firstln2)
					if c1ss.find(qt2) != -1:
						c1 = 1
					if c2ss.find(qt2) != -1:
						c2 = 1
					if c3ss.find(qt2) != -1:
						c3 = 1
					if c4ss.find(qt2) != -1:
						c4 = 1
					if c5ss.find(qt2) != -1:
						c5 = 1
					if c6ss.find(qt2) != -1:
						c6 = 1
					if c7ss.find(qt2) != -1:
						c7 = 1
					if c8ss.find(qt2) != -1:
						c8 = 1
					if c9ss.find(qt2) != -1:
						c9 = 1
					if firstln2.find(x12) == -1:
						for linenad1, lineg in enumerate(filen):
							if lineg.find(x12) == -1:
								templist.append(lineg)
							elif lineg.find(x12) != -1:
								templist8.append(p1+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(c3)+"\t"+str(c4)+"\t"+str(c5)+"\t"+str(c6)+"\t"+str(c7)+"\t"+str(c8)+"\t"+str(c9)+"\n")
								break
					if firstln2.find(x12) != -1:
						templist8.append(p1+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(c3)+"\t"+str(c4)+"\t"+str(c5)+"\t"+str(c6)+"\t"+str(c7)+"\t"+str(c8)+"\t"+str(c9)+"\n")
						continue
					else:

						continue 			
					
					
				else:
					continue
			if line3.find(x12) != -1:
				c1s = []
				c2s = []
				c3s = []
				c4s = []
				c5s = []
				c6s = []
				c7s = []
				c8s = []
				c9s = []
				c1ss = []
				c2ss = []
				c3ss = []
				c4ss = []
				c5ss = []
				c6ss = []
				c7ss = []
				c8ss = []
				c9ss = []
				
			
			else:
				continue
		else:
			continue
					
fileout.write(''.join(templist))
fileout4.write('Protein ID\tPfam DNA count\tPfam RNA count\tPfam Nucleotide count\tIntPro DNA count\tIntPro RNA count\tIntPro Nucleotide count\tGo DNA\tGO RNA\tGo Nucleotide\n')
fileout4.write(''.join(templist8))
fileout4.close()
#fileout2.write(''.join(templist2))
fileout.close()
if Disc == "yes" or Disc == "Y":
	fileout2 = open(filename3, 'a')
	fileout2.write('ProteinID\tPfam ID DNA\t Pfam Discriptor\n')
	fileout2.write('\n'.join(templist2))
	fileout2.write('\nProteinID\tPfam ID RNA\t Pfam Discriptor\n')
	fileout2.write('\n'.join(templist3))
	fileout2.write('\nProteinID\tPfam ID nucleotide\t Pfam Discriptor\n')
	fileout2.write('\n'.join(templist4))
	fileout2.close()
	fileout3 = open(filename4, 'a')	
	fileout3.write('ProteinID\tIntPro DNA\tInterpro Discriptor\n')
	fileout3.write('\n'.join(templist5))
	fileout3.write('\nProteinID\tIntPro RNA\tInterpro Discriptor\n')
	fileout3.write('\n'.join(templist6))
	fileout3.write('\nProteinID\tIntPro Nucleotide\tInterpro Discriptor\n')
	fileout3.write('\n'.join(templist7))		
	fileout3.close()
	fileout5 = open(filename6, 'a')	
	fileout5.write('ProteinID\tGO DNA\tGO Discriptor\n')
	fileout5.write('\n'.join(templist9))
	fileout5.write('\nProteinID\tGO RNA\tGO Discriptor\n')
	fileout5.write('\n'.join(templist10))
	fileout5.write('\nProteinID\tGO Nucleotide\tGO Discriptor\n')
	fileout5.write('\n'.join(templist11))	
	fileout5.close()
#
