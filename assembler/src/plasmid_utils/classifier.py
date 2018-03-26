

#Laplasian Smoothing
def lap_smooth(table):
    ### adds 1 to each element of the table. Don't forget to add # of rows to denominator
    for i in table:#
      #  if len(i)!=3: print(i)
        i[1]=int(i[1])+1# = [int(i)+1 for i in i[1:]]
        i[2]=int(i[2])+1
    return table;



def naive_bayes(input_list):
    # Takes list of HMMS,returns list: [Plasmid/Chromosome, PPlasm,Pchrom]
    #Constants
    plasm_genes=1041995
    chrom_genes=2990937



# Open table
    with open("/Nancy/mrayko/jgi_10k/pfam_count_total_fin_sorted10.txt", 'r') as infile:
        table=infile.readlines()

# create table - each row contain one feature. 
    table = [i.split() for i in table]
    for i in table:
        i=[x for x in i if x] 

# keep first 3 rows
    table = [i[:3] for i in table]


# Add probabilities

#P (plasm) of each hmm - # of occurences on plamids div. by (total number of genes in plasmids + num. of hmms)



    table=lap_smooth(table)

    for i in table:
        i.append (int(i[1])/(chrom_genes+len(table)))	
        i.append (int(i[2])/(plasm_genes+len(table)))

    plasm=1
    chrom=1
    for i in input_list:
        for j in table:
            if i==j[0]: 
                chrom=chrom*j[3]
                plasm=plasm*j[4]
    res = plasm/chrom
    if res > 1: return "NBC_Plasmid"  #["Plasmid", plasm, chrom]
    else: return "NBC_Chromosome" # ["Chromosome", plasm, chrom]

   
#print ("P plasmid = ", plasm) 
#print ("P chromosome = ", chrom)





# Apply classifier



   



def main():
#lap_smooth(table)
#	print (table[:3])
	input_list=("MCM_OB", "Phage_integrase", "Phage_integrase")
	print (naive_bayes(input_list))


#	return

if __name__ == "__main__":
    # execute only if run as a script
    main()
#Сколько всего генов в пламзидах и хромосомах?
#print (434/0.000144770685574)
#print (110/0.000104607027865)
#print (int(1)/2990937)
#print (plasm_genes+len(table))
#print (int(1)/(plasm_genes+len(table)))
#27     i.append (int(i[1])/2990937.0)
#28     i.append (int(i[2])/1041995.0)
