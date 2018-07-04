from sklearn.naive_bayes import MultinomialNB


from math import log
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


    table=lap_smooth(table)
    for i in table:
    	for j in i[1:3]:
    		if j==0: print (j)

    for i in table:
        i.append (float(i[1])/(chrom_genes+len(table)))	
        i.append (float(i[2])/(plasm_genes+len(table)))

    plasm=0
    chrom=0
    for i in input_list:
        for j in table:
            if i==j[0]:
                chrom=chrom+log(j[3])
                plasm=plasm+log(j[4])
                
    res = plasm - chrom
    if res > 0: return "Plasmid" + "  " +  str(plasm) + " " + str(chrom)
    else: return str("Chromosome" + " " +  str(plasm) + " " +  str(chrom))


def create_vector_pfams(hmms): # list of hmm lists

    with open("/Nancy/mrayko/chromosomal_removal_test/verification_test/pfam_names.list", "r") as infile:
        pfam_list=infile.readlines()
    pfam_list = [i.strip() for i in pfam_list] 


    vector=[len(pfam_list)*[0]]*len(hmms)
    for i in range(0,len(hmms)): # take each hit
       for j in hmms[i]:
        #print (j)
        hit_index = pfam_list.index(j)
        vector[i][hit_index]+=1
    return vector


def scikit_multNB (input_list):
    import pickle
    # load it again
#    with open('/Nancy/mrayko/PlasmidVerify/scikit_nbc/plasmids_and_chunks/my_dumped_classifier_chunks.pkl', 'rb') as fid:
    with open('/Nancy/mrayko/PlasmidVerify/scikit_nbc/plasmids_and_chunks/my_dumped_classifier_chunks.pkl', 'rb') as fid:
        clf = pickle.load(fid)

    a=create_vector_pfams([input_list])
#    print (clf.predict(a))
 #
#     print (clf.predict_proba(a))

#    return str(clf.predict(a)), str(clf.predict_proba(a))
    return clf.predict(a)[0] + "\t" +  "\t".join(map(str,clf.predict_log_proba(a)[0])) 
#print(clf.predict(a))
#print(clf.predict_proba(a))




# Apply classifier

def main():
#lap_smooth(table)
#	print (table[:3])
	#input_list=("MCM_OB", "Phage_integrase", "Phage_integrase")

    #t="DUF11 CARDB TRAP_beta DUF2393"
#    t="Y1_Tnp Y1_Tnp"
 #   t=("HTH_17", "DUF5447", "Phage_Coat_B", "DUF2523", "Zot", "Phage_integrase", "PhdYeFM_antitox","ParE_toxin")
    t = "Terminase_2 Resolvase HTH_17 Phage_AlpA VirE DUF3924 HTH_39"
    input_list=(t.split(" "))
   
    print (naive_bayes(input_list))
    print (scikit_multNB(input_list))
    return






if __name__ == "__main__":
    # execute only if run as a script
    main()
