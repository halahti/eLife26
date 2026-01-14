#This is the replication code for the data preparation for genetic variables (PGIs, PCs and IBD estimates)
#Lahtinen H., Ganna A., Kaprio J., Korhonen K., Lombardi S., Silventoinen K., Martikainen P. Heterogenous associations of polygenic indices of 35 traits with mortality. Published in eLife.
#To be published in  eLife
#Software used: PLINK 1.9-2.0, R
#see Stata commands for phenotype data prepartaion and analysis in another code file of the github repository https://github.com/halahti/eLife26/
#Hannu Lahtinen 12 January 2026



#pruning the data for principal components and relatedness estimation

#long-range ld.set (--exclude) from https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt

grun.py -n tmp_prune -q highmem.q -c " /plink-1.9/plink --bfile /combined/fr_fh_h2000_tw_autosome --maf 0.01 \
--indep-pairwise 50 5 0.2 --exclude /finriski/high_ld_fr.set  --out /combined/fr_ft_h2000_tw_pruned"



#Principal components
grun.py -n tmp_pca -q hugemem.q -c " /plink-1.9/plink --bfile /combined/fr_fh_h2000_tw_autosome  \
--pca 20   --extract /combined/fr_ft_h2000_tw_pruned.prune.in --out /combined/fr_fh_h2000_tw_pca"


#Relatedness estimates
grun.py -n tmp_pihat  -q hugemem.q -c " /plink-1.9/plink --bfile  /combined/fr_fh_h2000_tw_autosome  \
--genome --min 0.05   --extract /combined/fr_ft_h2000_tw_pruned.prune.in --out /combined/fr_fh_h2000_tw_0.05"


########
#PGIs  #
########   


for i in ACTIVITY1 ADHD1 ADVENTURE1 AFB2 ASTECZRHI1 ASTHMA1 AUDIT1 BMI2 CANNABIS2 CP2 CPD2 DEP1 DPW2  EVERSMOKE2 EXTRA2 FAMSAT1 FRIENDSAT1 HAYFEVER1 HEIGHT3 HIGHMATH1 LEFTOUT1 MENARCHE1 MIGRAINE1 MORNING1 NARCIS1 NEARSIGHTED1 NEBwomen2 NEURO2 OPEN2 READING1 RELIGATT1 RISK2 SELFHEALTH1 SELFMATH1  
do
wc -l $i-single_weights_LDpred_p1.0000e+00.txt
head $i-single_weights_LDpred_p1.0000e+00.txt
done

#####################
#a quick visit to R##
#####################

R
info<-readRDS("/homes/hlahtine/ldpred2/25503788")
info$sid<-paste0(info$chr,":",info$pos)
info2<-data.frame(info$a0,info$a1,info$sid,info$rsid)
write.table(info2,file = "/hg19ids.txt", quote = F, col.names=F,row.names=F)
q("no")

###############################
#Back in the Unix command line#
###############################

for i in ACTIVITY1 ADHD1 ADVENTURE1 AFB2 ASTECZRHI1 ASTHMA1 AUDIT1 BMI2 CANNABIS2 CP2 CPD2 DEP1 DPW2  EVERSMOKE2 EXTRA2 FAMSAT1 FRIENDSAT1 HAYFEVER1 HEIGHT3 HIGHMATH1 LEFTOUT1 MENARCHE1 MIGRAINE1 MORNING1 NARCIS1 NEARSIGHTED1 NEBwomen2 NEURO2 OPEN2 READING1 RELIGATT1 RISK2 SELFHEALTH1 SELFMATH1  
do
join -j 3 -o 2.4,1.3,1.4,1.5,1.6,1.7  <(sort -k3 $i-single_weights_LDpred_p1.0000e+00.txt) <(sort -k3 hg19ids.txt)  > $i.txt
done


cd /gwasres/PGIrepo

#ei vaan looppi suostunut tallentamaan suoraan --out-komentoon. kierretään mv:llä?

for i in ACTIVITY1 ADHD1 ADVENTURE1 AFB2 ASTECZRHI1 ASTHMA1 AUDIT1 BMI2 CANNABIS2 CP2 CPD2 DEP1 DPW2  EVERSMOKE2 EXTRA2 FAMSAT1 FRIENDSAT1 HAYFEVER1 HEIGHT3 HIGHMATH1 LEFTOUT1 MENARCHE1 MIGRAINE1 MORNING1 NARCIS1 NEARSIGHTED1 NEBwomen2 NEURO2 OPEN2 READING1 RELIGATT1 RISK2 SELFHEALTH1 SELFMATH1  
do
/plink-2-080621/plink2 --bfile /combined/hm3_autosome_R8rs \
 --double-id --rm-dup force-first --score  /$i.txt 1 3 6 --hwe 1e-8 --out pg
mv pg.sscore $i

/plink-2-080621/plink2 --bfile /twins/twin_hm3 \
 --double-id --rm-dup force-first --score  /$i.txt 1 3 6 --hwe 1e-8 --out /gwasres/PGIrepo/pgtw
mv pgtw.sscore tw$i
done



#########################################################
#there are two exceptions, 				#
#1) EA6 needs to be dited for hapmap variants		#
#2) RISK2 has direct rs-codes				#
# these hve to be formed separately			#
#########################################################




grun.py -n ru_ea6  -q highmem.q -c " /plink-2-080621/plink2 --bfile /combined/full_autosome_R8ids \
 --double-id --rm-dup force-first --score  /EA6-single_weights_SBayesR.txt 2 5 8 --hwe 1e-8 --out /gwasres/PGIrepo/ea6_frfhh2000"


/plink-2-080621/plink2 --bfile /combined/hm3_autosome_R8rs \
 --double-id --rm-dup force-first --score  /RISK2-single_weights_LDpred_p1.0000e+00.txt 3 4 7 --hwe 1e-8 --out pg
mv pg.sscore RISK2

/plink-2-080621/plink2 --bfile /twins/twin_hm3 \
 --double-id --rm-dup force-first --score  /RISK2-single_weights_LDpred_p1.0000e+00.txt 3 4 7 --hwe 1e-8 --out /gwasres/PGIrepo/pgtw
mv pgtw.sscore twRISK2


