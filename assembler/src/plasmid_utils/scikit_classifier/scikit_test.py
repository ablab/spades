import sys
import numpy as np
from sklearn.naive_bayes import MultinomialNB
from sklearn.model_selection import train_test_split
import collections

#plasmid_file =  "plasmids_1000.tblout" # "all_refseq_plasmids_PfamA_001_tblout"   #sys.argv[1] "plasmids_1000.tblout" 
#chromosome_file =  "chromosomes_10000.tblout" #  "total_filtered.fna_tblout"  #sys.argv[2]  "chromosomes_1000.tblout"

with open("pfam_names.list", "r") as infile:
        pfam_list=infile.readlines()
pfam_list = [i.strip() for i in pfam_list] 




def get_annot_from_tblout(tblout_pfam):
# Takes tblout file, returns list of contig names + hmms (e.g. ['NODE_4_length_6662_cov_295.444039_cutoff_5', 'ABC2_membrane', 'HTH_3', 'Rep_3', 'Relaxase', 'MobC'])    
    with open(tblout_pfam, "r") as infile:
        tblout_pfam=infile.readlines()

    tblout_pfam = [i.split() for i in tblout_pfam] 
    contigs = collections.OrderedDict()
    for i in range(3, len(tblout_pfam)-10):
        name = tblout_pfam[i][2].rsplit("_", maxsplit=1)[0]
        if (name not in contigs) and (float(tblout_pfam[i][4]) < 0.01)  : # if there's the SRR in protein ID and e-value < 0.01 : 
            contigs[name]=[tblout_pfam[i][0]]
           # print (contigs)
        else:
            if tblout_pfam[i][2] != tblout_pfam[i-1][2] and float(tblout_pfam[i][4]) < 0.01:
                contigs[name].append(tblout_pfam[i][0])
    
    return contigs
            



def create_vector_pfams(hmms): # list of hmm lists
   vector=[len(pfam_list)*[0]]*len(hmms)
   for i in range(0,len(hmms)): # take each hit
       for j in hmms[i]:
       	#print (j)
        hit_index = pfam_list.index(j)
        vector[i][hit_index]+=1
   return vector


    
print ("Extracting hmms from tblout file...")


#plasmid_train = get_annot_from_tblout(plasmid_file)
#chromosome_train = get_annot_from_tblout(chromosome_file)

#print (chromosome_train)


#input_list=(t.split(" "))

p1="Resolvase Sulfate_transp Usp KAP_NTPase HTH_Tnp_1 Resolvase AlbA_2 MobA_MobL TraC DUF3847 RepA_N DUF536 G_glu_transpept YitT_membrane Resolvase Resolvase"
p2="DUF334 HTH_3 Resolvase"
p3="Mob_Pre MFS_1 Rep_trans Mob_Pre"
p4="DUF536 Pkinase ABC_membrane PadR DUF1700 DUF4969 Resolvase Beta-lactamase2 Peptidase_M56 Penicillinase_R Resolvase rve AAA_22 ACP_syn_III HTH_11 Gram_pos_anchor DUF536 RepA_N Acetyltransf_10"
c1="Spore-coat_CotZ DUF1360 Spore-coat_CotZ" 
c2="GTP_EFTU_D3 Hexapep SecE NusG Ribosomal_L11_N Ribosomal_L1 Ribosomal_L10 Ribosomal_L12 RNA_pol_Rpb2_6 RNA_pol_Rpb1_1 EamA Ribosom_S12_S23 Ribosomal_S7 GTP_EFTU GTP_EFTU"
c3="Imm40 Cys_rich_CPCC MafB"
# d = OrderedDict({'a':1, 'b':2}),
plasmid_train = collections.OrderedDict({"1p":p1.split(" "), "2p": p2.split(" "), "3p":p3.split(" "), "4p":p4.split(" ")})
chromosome_train = collections.OrderedDict({"1c":c1.split(" "), "2c":c2.split(" "), "3c":c3.split(" ")})
#plasmid_train = collections.OrderedDict(plasmid_dataset.items()[:7000])
#chromosome_train = collections.OrderedDict(chromosome_dataset.items()[:150])


  
#plasmid_test = collections.OrderedDict(plasmid_dataset.items()[7000:])
#chromosome_test = collections.OrderedDict(chromosome_dataset.items()[150:])   


#plasmid_test = plasmid_dataset[7000:]
#chromosome_test = chromosome_dataset[150:]


#with open('plasmid_train.txt', 'wb') as outfile:
 #    for i,k in plasmid_train.items():
  #     outfile.write(i+ " " + k + "\n")    

#with open('chromosome_train.txt', 'wb') as outfile:
 #    for i,k in chromosome_train.items():
  #      outfile.write(i+ " " + k + "\n")    


#with open('plasmid_test.txt', 'wb') as outfile:
 #    for i,k in plasmid_test.items():
  #      outfile.write(i+ " " + k + "\n")    


#with open('chromosome_test.txt', 'wb') as outfile:
 #    for i,k in chromosome_test.items():
  #      outfile.write(i+ " " + k + "\n")    




print (len(plasmid_train))
print (len(chromosome_train))


#print (plasmid_train)
#print (chromosome_train)



# ok, let's assume we have two training sets. How to train classifier?
#Combine both datasets:

train=collections.OrderedDict(plasmid_train, **chromosome_train)

print(train)

print ("Initializing matrix...")


X=[len(pfam_list)*[0]]*len(train)  # initialize matirx of n rows of zeroes, len = num of pfams

print ("Populating matrix...")

counter=0
for key, value in train.items(): # take each sample  key, value in d.items():
#	print (key)
#	print (value)
	for j in value: # take each hit
	    hit_index = pfam_list.index(j)  # int(str([x for x in range(len(pfam_list)) if pfam_list[x]==j]))   #int(i for i,x in enumerate(pfam_list) if x == j) # pfam_list.index(j)  
	    X[counter][hit_index]+=1
	counter+=1

y=["Plasmid"]*len(plasmid_train)+["Chromosome"]*len(chromosome_train)

#print (X)
#print (y)

#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=42)

#print (X_train)

print ("MultinomialNB training...")

# Train classificator
clf = MultinomialNB(alpha=0.1)
clf.fit(X, y)

#print (clf)
# Test

#t = "Terminase_2 Resolvase HTH_17 Phage_AlpA VirE DUF3924 HTH_39"
#[[-3.66445918 -0.02595188]]
#[[-2.49097886 -0.08646118]]

#t= "DUF536 Acetyltransf_10 RepA_N DUF536 Gram_pos_anchor HTH_11 ACP_syn_III AAA_22 rve Resolvase Penicillinase_R Peptidase_M56 Beta-lactamase2 Resolvase DUF4969 DUF1700 PadR ABC_membrane Pkinase"
#[[-3.66512806 -0.0259343 ]]
#[ -1.19801233e+01  -6.26758063e-06]]
#input_list=(t.split(" "))

#input_list=['Biotin_lipoyl', 'CPSase_L_D2', 'UvrD_C', 'NAD_binding_1', 'PdxJ', 'ParBc', 'Pribosyltran', 'BATS', 'Patatin', 'UbiA', 'Fapy_DNA_glyco', 'Cytochrom_C_asm', 'SRP54', 'Aminotran_4', 'Fimbrial', 'Fimbrial', 'Fimbrial', 'Usher', 'PapD_N', 'FimA', 'CinA', 'Sugar_tr', 'GTP-bdg_N', 'Hfq', 'IPPT', 'Pterin_bind', 'Peptidase_M41', 'FtsJ', 'CRS1_YhbY', 'DUF498', 'ILVD_EDD', 'zf-CHCC', 'Glyco_transf_4', 'Wzy_C', 'Lip_A_acyltrans', 'Glycos_transf_N', 'Alpha_L_fucos', 'Ribosomal_S16', 'RimM', 'tRNA_m1G_MT', 'Ribosomal_L19', 'Peptidase_M24', 'THDPS_N_2', 'ArsC', 'Peptidase_M20', 'Asn_synthase', 'Beta_helix', 'NMO', 'LexA_DNA_bind', 'RecA', 'tRNA-synt_2c', 'CsrA', 'Peptidase_M3', 'PALP', 'Cu-oxidase_3', 'CopB', 'tRNA-synt_1', 'DNA_pol3_chi', 'Peptidase_M17', 'YjgP_YjgQ', 'YjgP_YjgQ', 'RDD', 'SIS', 'Acyl_transf_3', 'Spermine_synth', 'Orn_Arg_deC_N', 'adh_short', 'SPOR', 'DALR_1', 'RadC', 'DFP', 'dUTPase', 'PGM_PMM_I', 'Pribosyltran', 'DUF3289', 'Peptidase_C1', 'DUF2326', 'Exo_endo_phos', 'MFS_1', 'AnmK', 'Peptidase_M23', 'tRNA-synt_1b', 'Fapy_DNA_glyco', 'Cytochrom_C_asm', 'DUF885', 'adh_short_C2', 'GTP_EFTU', 'HlyIII', 'Glycos_transf_2', 'TatD_DNase', 'Rick_17kDa_Anti', 'GCV_H', 'GCV_T', 'NfeD', 'Band_7', 'MazG', 'Inositol_P', 'NUDIX', 'Aminotran_3', 'Methyltrans_RNA', 'DEAD', 'PTPS', 'DUF2059', 'Acyltransferase', 'DUF502', 'Methyltransf_23', 'UPF0093', 'ACCA', 'DNA_pol3_alpha', 'SAICAR_synt', 'DnaJ', 'Ribul_P_3_epim', 'Chorismate_bind', 'GATase', 'Glycos_transf_3', 'IGPS', 'HAD', 'MarR', 'Peptidase_M24', 'DUF885', 'DUF885', 'TGT', 'YajC', 'SecD-TM1', 'SecD_SecF', 'PolyA_pol', 'HPPK', 'Pantoate_transf', 'Pantoate_ligase', 'Asp_decarbox', 'PGI', 'DUF150', 'NusA_N', 'GTP_EFTU', 'RBFA', 'TruB_N', 'Ribosomal_S15', 'RNase_PH', 'HxlR', 'NmrA', 'ACR_tran', 'Biotin_lipoyl_2', 'EamA', 'DEAD', 'ETF_alpha', 'ETF', 'GDP_Man_Dehyd', 'NTP_transferase', 'dTDP_sugar_isom', 'RmlD_sub_bind', 'MannoseP_isomer', 'PGM_PMM_I', 'Peptidase_S8', 'MFS_4', 'DUF1439', 'PFK', 'ADK', 'Mur_ligase_M', 'SNARE_assoc', 'APH', 'Peptidase_M17', 'AI-2E_transport', 'Trypsin_2', 'Epimerase', 'AMP-binding', 'Aconitase', 'MazE_antitoxin', 'PIN', 'Aconitase_2_N', 'HSDR_N', 'Methylase_S', 'N6_Mtase', 'HsdM_N', 'Biotin_lipoyl_2', 'TIM', 'SecG', 'Oxidored_q4', 'Oxidored_q6', 'Complex1_30kDa', 'Complex1_49kDa', '2Fe-2S_thioredx', 'Complex1_51K', 'Molybdopterin', 'NADHdh', 'Fer4', 'Oxidored_q3', 'Oxidored_q2', 'Proton_antipo_M', 'Proton_antipo_M', 'Proton_antipo_M', 'adh_short_C2', 'CitMHS', 'Porin_O_P', 'Response_reg', '2CSK_N', 'SBP_bac_11', 'Response_reg', 'PAS_7', 'Porin_O_P', 'SDF', 'PTA_PTB', 'HSP90', 'Cons_hypoth95', 'CTP_transf_like', 'Fer4', 'G_glu_transpept', 'GGDEF', 'Peptidase_M23', 'Amidohydro_1', 'Haemolytic', 'zf-dskA_traR', 'Hydrolase_4', 'MFS_3', 'SufE', 'tRNA-synt_1e', 'OTCace_N', 'Arginosuc_synth', 'Peptidase_M20', 'NAT', 'Semialdhyde_dh', 'Lyase_1', 'AA_kinase', 'Aldedh', 'HCBP_related', 'dCMP_cyt_deam_1', 'Fe_dep_repress', 'Nramp', 'ATE_C', '4HBT_3', 'Peptidase_S8', 'Peptidase_S15', 'Acyltransferase', 'TonB_dep_Rec', 'AdoHcyase', 'AMP-binding', 'RNase_HII', 'LpxB', 'Acetyltransf_11', 'FabA', 'Hexapep', 'Bac_surface_Ag', 'Peptidase_M50', 'DXP_reductoisom', 'CTP_transf_1', 'Prenyltransf', 'RRF', 'Toluene_X', 'Rhomboid', 'DUF1820', 'AA_kinase', 'Aldolase', 'ILVD_EDD', 'Glucosamine_iso', 'Glucokinase', 'G6PD_C', 'ABC_tran', 'GCV_T_C', 'DUF1674', 'Sdh_cyt', 'Sdh_cyt', 'FAD_binding_2', 'Fer2_3', 'Sdh5', 'MacB_PCD', 'ABC_tran', 'Competence', 'MotA_ExbB', 'ExbD', 'ABC_membrane', 'LpxK', 'Phage_GPD', 'Tail_P2_I', 'Baseplate_J', 'GPW_gp25', 'HigB-like_toxin', 'HTH_3', 'MqsR_toxin', 'MqsA_antitoxin', 'Phage_base_V', 'Minor_tail_Z', 'Phage_capsid', 'Phage_portal_2', 'VRR_NUC', 'SNF2_N', 'DUF4224', 'Phage_integrase', 'CbiA', 'His_Phos_1', 'YceI', 'PLDc_2', 'HTH_31', 'Asparaginase', 'FMN_red', 'NAPRTase', 'PAP2', 'DUF2782', 'DNA_pol_A', 'DapB_C', 'CPSase_sm_chain', 'CPSase_L_D2', 'GreA_GreB_N', 'DHHA1', 'PCRF', 'tRNA-synt_2', 'Response_reg', 'Response_reg', 'ECH_1', 'AA_kinase', 'Mur_ligase_M', 'TMP-TENI', 'MTHFR', 'SIMPL', 'Maf', 'RNase_E_G', 'AsmA_2', 'PmbA_TldD', 'DUF615', 'PmbA_TldD', 'LysR_substrate', 'FMN_red', 'ADH_N', 'DJ-1_PfpI', 'NTP_transf_3', 'CM_2', 'ATP-synt_DE_N', 'ATP-synt_ab', 'ATP-synt', 'ATP-synt_ab', 'OSCP', 'ATP-synt_B', 'ATP-synt_C', 'ATP-synt_A', 'Ribosomal_S10', 'Ribosomal_L3', 'Ribosomal_L4', 'Ribosomal_L23', 'Ribosomal_L2_C', 'Ribosomal_S19', 'Ribosomal_L22', 'Ribosomal_S3_C', 'Ribosomal_L16', 'Ribosomal_L29', 'Ribosomal_S17', 'Ribosomal_L14', 'ribosomal_L24', 'Ribosomal_L5_C', 'Ribosomal_S14', 'Ribosomal_S8', 'Ribosomal_L6', 'Ribosomal_L18p', 'Ribosomal_S5_C', 'Ribosomal_L30', 'Ribosomal_L27A', 'SecY', 'Ribosomal_S13', 'Ribosomal_S11', 'Ribosomal_S4', 'RNA_pol_A_bac', 'Ribosomal_L17', 'BMFP', 'Mg_chelatase', 'Abhydrolase_1', 'Lipase_chap', 'UDPG_MGDP_dh_N', 'Trigger_C', 'CLP_protease', 'AAA_2', 'Lon_C', 'Bac_DNA_binding', 'Rotamase_3', 'DUF475', 'Ribonuc_red_lgC', 'Ribonuc_red_sm', 'Flavodoxin_1', 'Thioredoxin', 'PseudoU_synth_2', 'HAD_2', 'DUF1415', 'Ribosomal_L28', 'Ribosomal_L33', 'PLDc_2', 'GST_N', 'Ldh_1_N', 'Pro_isomerase', 'GTP_EFTU', 'DUF2127', 'HlyD_3', 'ABC_membrane', 'ABC2_membrane', 'ABC_tran', 'Neisseria_PilC', 'HrpB_C', 'GAF_2', 'Bac_surface_Ag', 'TamB', 'PEP_mutase', 'Terminase_6', 'Aconitase', 'Aconitase_C', 'PrpF', 'Arginase', 'A2M_N', 'Abhydrolase_3', 'Cupin_6', 'RNase_T', 'MS_channel', 'PPDK_N', 'Kinase-PPPase', 'DUF1249', 'NUDIX', 'Autotransporter', 'Glyco_hydro_6', 'Radical_SAM', 'BPL_LplA_LipB', 'DUF493', 'Yop-YscD_cpl', 'Phasin', 'DUF445', 'DUF21', 'Peptidase_S10', 'DNA_gyraseB', 'CTP_synth_N', 'DAHP_synth_1', 'Enolase_C', 'DivIC', 'IspD', 'YgbB', 'Smr', 'SGL', 'ETF_QO', '2OG-FeII_Oxy_2', 'ABC_trans_aux', 'MlaD', 'ABC_tran', 'MlaE', 'ThrE', 'MreB_Mbl', 'MreC', 'MreD', 'Transpeptidase', 'FTSW_RODA_SPOVE', 'Queuosine_synth', 'RelA_SpoT', 'MinE', 'CbiA', 'MinC_C', 'Acetyltransf_1', 'HhH-GPD', 'YceI', 'Ni_hydr_CYTB', 'YceI', 'HATPase_c', 'DNA_methylase', 'DUF262', 'URO-D', 'DHQ_synthase', 'SKI', 'Putative_PNPOx', 'tRNA-synt_1c', 'CutC', 'SBP_bac_11', 'BPD_transp_1', 'BPD_transp_1', 'ABC_tran', 'Lip_A_acyltrans', 'Sigma70_ner', 'DNA_topoisoIV', 'MarR', 'Abhydrolase_6', 'Aminotran_1_2', 'COX15-CtaA', 'UbiA', 'Lactamase_B', 'Lactamase_B', 'PolyA_pol', 'SLT', 'PS_Dcarbxylase', 'SCO1-SenC', 'MTS', 'Chorismate_synt', 'Semialdhyde_dhC', 'PseudoU_synth_1', 'PRAI', 'PALP', 'Trp_syntA', 'ATPase']

#input_list = (chromosome_train["NC_004556.1"])

#a=create_vector_pfams([input_list])

print(clf.predict(X))
print(clf.predict_log_proba(X))
#print(clf.get_params())

print (y)

#import pickle
# save the classifier
#with open('my_dumped_classifier.pkl', 'wb') as fid:
 #   pickle.dump(clf, fid)    
