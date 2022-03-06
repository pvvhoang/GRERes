#===============================
# For 50K samples, g+e
# 1) Sampling distributions of model parameters for the genome-residual partitioning model.
#===============================

# Input data
# .fam						fam file
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# Run on statgen server
export DATA=/data/alh-admvhp1/GRERes/Data
export OUTPUT=/data/alh-admvhp1/GRERes/Data

cd $OUTPUT

# 1.1) fam file
# scp ./../../GRE/Data/qced_rdm_005_ukbb3.fam ./
# qced_rdm_005_ukbb3.fam

# 1.2) Select 50K samples
fam_file="./qced_rdm_005_ukbb3.fam"
Rscript ../selectSamples.R --fam_file=${fam_file} --number=50000 --outFile1="selectedFam.fam"

# 1.3) Genomic ralationship matrix
../Tool/plink --bfile qced_rdm_005_ukbb3 --keep selectedFam.fam --make-grm-gz --out grm_qced_rdm_005_ukbb3_selected
gunzip grm_qced_rdm_005_ukbb3_selected.grm.gz
# remove the third column
awk '{print $1,$2,$4}' grm_qced_rdm_005_ukbb3_selected.grm > grm_qced_rdm_005_ukbb3_selected_processed.grm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.4) (**) Matrix I (identity.bmat) instead of kernel matrix K (brain.bmat)

awk '{if($1 != $2) {a = 0; print $1, $2, a;} else {a = 1; print $1, $2, a;}}' grm_qced_rdm_005_ukbb3_selected_processed.grm > identity.bmat

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1.5 (***) 
# 2.2. Compute kernel matrix Q

# i. ensure positive definite kernel matrices
# .	for A
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -thread 80 -bend 2

# .	for I
../Tool/mtg2 -p selectedFam.fam -g identity.bmat -thread 80 -bend 2

# ii. Cholesky decomposition
# .	for A
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend -thread 80 -chol 1

# .	for I
../Tool/mtg2 -p selectedFam.fam -g identity.bmat.bend -thread 80 -chol 1

# iii. compute Q
../Tool/mtg2 -p selectedFam.fam -mg identity_selected_grm_bmat.chol -thread 80 -matmat 2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
# E1) Continuous data
#===============================

# Input data
# .fam						fam file
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# # Run on statgen server
# export DATA=/data/alh-admvhp1/GRERes/Data
# export OUTPUT=/data/alh-admvhp1/GRERes/Data
# 
# cd $OUTPUT

# # tango
# export DATA=/home/users/phamvv/GRERes/Data
# export OUTPUT=/home/users/phamvv/GRERes/Data
# 
# cd $OUTPUT

# gadi
export DATA=/scratch/eu82/vp8928/GRERes/gadi_50K_bi_cov0
export OUTPUT=/scratch/eu82/vp8928/GRERes/gadi_50K_bi_cov0

cd $OUTPUT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E1.1) Simulate phenotype data & fit
# .dat: phenotypic data of the main trait
# cov = 0
cov=0
if [ -e res.greml.out ]
then
    rm res.greml.out
fi
if [ -e res.coregreml.out ]
then
    rm res.coregreml.out
fi
if [ -e res.coregreml.out.converted ]
then
    rm res.coregreml.out.converted
fi
for ((i=1;i<=50;i++))
    do
    
    echo "cov = 0 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g identity.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: identity.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside identity_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    # ../Tool/mtg2 -p selectedFam.fam -mg identity_greml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> res.greml.out
     
    # .	fit CORE GREML:
    # o	inside identity_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> res.coregreml.out
    #++++++++++++++++++++++++
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> res.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> res.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> res.coregreml.out.converted
    #++++++++++++++++++++++++

  done

# cov = 0.2
cov=0.2
if [ -e resp.greml.out ]
then
    rm resp.greml.out
fi
if [ -e resp.coregreml.out ]
then
    rm resp.coregreml.out
fi
if [ -e resp.coregreml.out.converted ]
then
    rm resp.coregreml.out.converted
fi
for ((i=1;i<=50;i++))
    do
    
    echo "cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g identity.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: identity.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside identity_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside identity_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # identity_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> resp.coregreml.out
    #++++++++++++++++++++++++
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> resp.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> resp.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> resp.coregreml.out.converted
    #++++++++++++++++++++++++
  done

# cov = -0.2
cov=-0.2
if [ -e resn.greml.out ]
then
    rm resn.greml.out
fi
if [ -e resn.coregreml.out ]
then
    rm resn.coregreml.out
fi
if [ -e resn.coregreml.out.converted ]
then
    rm resn.coregreml.out.converted
fi
for ((i=1;i<=50;i++))
    do
    
    echo "cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g identity.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: identity.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # model fitting
    # o	inside identity_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside identity_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # identity_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> resn.coregreml.out
    #++++++++++++++++++++++++
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> resn.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> resn.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> resn.coregreml.out.converted
    #++++++++++++++++++++++++
  done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# E2) For binary data, correct the transformation, scale
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E2.1) Simulate phenotype data & fit
# .dat: phenotypic data of the main trait
# cov = 0
cov=0
if [ -e bi.res.greml.out ]
then
    rm bi.res.greml.out
fi
if [ -e bi.res.coregreml.out ]
then
    rm bi.res.coregreml.out
fi
if [ -e bi.res.coregreml.out.converted ]
then
    rm bi.res.coregreml.out.converted
fi
if [ -e bi.res.greml.out.transformed ]
then
    rm bi.res.greml.out.transformed
fi
if [ -e bi.res.coregreml.out.transformed ]
then
    rm bi.res.coregreml.out.transformed
fi
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0 - i = " $i
    
    #i=1
    #k=10
    #seed=`expr $i \* $k`
    seed=-1
    Rscript ../random.R --cov=${cov} --seed=${seed} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g identity.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: identity.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside identity_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.res.greml.out
     
    # .	fit CORE GREML:
    # o	inside identity_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # identity_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out
    
    #++++++++++++++++++++++++
    # estimate h2 (refer to iii. estimate functions of model parameters)
    # Estimate proporitions of phenotypic variance explained by genetics and transcriptome,
    # correlation between genetic and transcriptomic effects on phenotypes,
    # and standard errors of these estimates
    
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.res.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.res.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.res.coregreml.out.converted
    #++++++++++++++++++++++++

  done

#-----------------------
# cov = 0, transformation

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[4]} ${line[5]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.res.greml.out.transformed
  done < bi.res.greml.out
 
while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.coregreml.out.transformed
    
    # for residual
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.res.coregreml.out.transformed
    
    # for cov
    ../Tool/mtg2 -trf_h2 ${line[10]} ${line[11]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.res.coregreml.out.transformed
  done < bi.res.coregreml.out.converted
  
#-----------------------

# cov = 0.2
cov=0.2
if [ -e bi.resp.greml.out ]
then
    rm bi.resp.greml.out
fi
if [ -e bi.resp.coregreml.out ]
then
    rm bi.resp.coregreml.out
fi
if [ -e bi.resp.coregreml.out.converted ]
then
    rm bi.resp.coregreml.out.converted
fi
if [ -e bi.resp.greml.out.transformed ]
then
    rm bi.resp.greml.out.transformed
fi
if [ -e bi.resp.coregreml.out.transformed ]
then
    rm bi.resp.coregreml.out.transformed
fi
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = 0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g identity.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: identity.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside identity_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resp.greml.out
     
    # .	fit CORE GREML:
    # o	inside identity_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # identity_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out
    
    #++++++++++++++++++++++++
    # estimate h2 (refer to iii. estimate functions of model parameters)
    # Estimate proporitions of phenotypic variance explained by genetics and transcriptome,
    # correlation between genetic and transcriptomic effects on phenotypes,
    # and standard errors of these estimates
    
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resp.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resp.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resp.coregreml.out.converted
    #++++++++++++++++++++++++

  done
  
#-----------------------
# cov = 0.2, transformation

# Transform file bi.resp.greml.out

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[4]} ${line[5]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resp.greml.out.transformed
  done < bi.resp.greml.out

# Transform file bi.resp.coregreml.out.converted

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.coregreml.out.transformed
    
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resp.coregreml.out.transformed
    
    # for cov
    ../Tool/mtg2 -trf_h2 ${line[10]} ${line[11]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resp.coregreml.out.transformed
  done < bi.resp.coregreml.out.converted
#-----------------------

# cov = -0.2
cov=-0.2
if [ -e bi.resn.greml.out ]
then
    rm bi.resn.greml.out
fi
if [ -e bi.resn.coregreml.out ]
then
    rm bi.resn.coregreml.out
fi
if [ -e bi.resn.coregreml.out.converted ]
then
    rm bi.resn.coregreml.out.converted
fi
if [ -e bi.resn.greml.out.transformed ]
then
    rm bi.resn.greml.out.transformed
fi
if [ -e bi.resn.coregreml.out.transformed ]
then
    rm bi.resn.coregreml.out.transformed
fi
for ((i=1;i<=50;i++))
    do
    
    echo "binary - cov = -0.2 - i = " $i
    
    #i=1
    Rscript ../random.R --cov=${cov} --seed=${i} --number=50000
    
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol -matvec rnd.v -thread 30
    # Out: grm_qced_rdm_005_ukbb3_selected_processed.grm.bend.chol.matvec
    ../Tool/mtg2 -p selectedFam.fam -g identity.bmat.bend.chol -matvec rnd.v2 -thread 30
    # Out: identity.bmat.bend.chol.matvec

    R CMD BATCH --no-save ../stdz_gxe2_2.r     #permutation adjusted y
    R CMD BATCH --no-save ../combine_ge_cor_2.r     #combine
    
    awk '{print $1,$2,$3}' sample.dat > sample_pheno.dat
    
    # Convert the phenotype data to binary
    Rscript ../toBinary.R --phe="sample_pheno.dat"
    
    # model fitting
    # o	inside identity_greml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # .	fit GREML:
    # -d sample_pheno.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d sample_pheno.dat -mod 1 -thread 100 -out brain_greml.out > brain_greml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_greml.out >> bi.resn.greml.out
     
    # .	fit CORE GREML:
    # o	inside identity_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity.bmat
    # identity_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d sample_pheno.dat -mod 1 -thread 100 -out brain_coregreml.out > brain_coregreml.log
    awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out
    
    #++++++++++++++++++++++++
    # estimate h2 (refer to iii. estimate functions of model parameters)
    # Estimate proporitions of phenotypic variance explained by genetics and transcriptome,
    # correlation between genetic and transcriptomic effects on phenotypes,
    # and standard errors of these estimates
    
    # o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
    # These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
    grep -vwE '(LKH|h2)' brain_coregreml.out > brain_coregreml.out2
    
    # o	Prepare a do file for MTG2 to construct functions of model parameters.
    Rscript ../createTemplate.R --inFile="brain_coregreml.out2" --outFile="brain_coregreml.do"
    
    # .	Estimation:
    ../Tool/mtg2 -delta2 brain_coregreml.do > brain_coregreml.out3 
    
    # Write h2 and se to file
    awk 'n!=1 && $1!="h2" {printf $2" "$3" "}; $1=="h2" {n=1}' brain_coregreml.out >> bi.resn.coregreml.out.converted
    awk '$1=="Ratio:" {printf $2" "$4" "}; $1=="Cor" {printf $3" "$5" "}' brain_coregreml.out3 >> bi.resn.coregreml.out.converted
    awk '$1=="LKH" {printf $2"\n"}' brain_coregreml.out >> bi.resn.coregreml.out.converted
    #++++++++++++++++++++++++

  done
  
#-----------------------
# cov = -0.2, transformation

# Transform file bi.resn.greml.out

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[4]} ${line[5]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resn.greml.out.transformed
  done < bi.resn.greml.out

# Transform file bi.resn.coregreml.out.converted

while IFS=" " read -a line;
  do
    # for genotype
    ../Tool/mtg2 -trf_h2 ${line[6]} ${line[7]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.coregreml.out.transformed
    
    ../Tool/mtg2 -trf_h2 ${line[8]} ${line[9]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8" "}' temp.out >> bi.resn.coregreml.out.transformed
    
    # for cov
    ../Tool/mtg2 -trf_h2 ${line[10]} ${line[11]} 0.1 0.1 obs > temp.out
    awk '$2 == "h2" {printf $8" "}; $2 == "SE" {printf $8"\n"}' temp.out >> bi.resn.coregreml.out.transformed
  done < bi.resn.coregreml.out.converted
#-----------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# E3) Check the phenotype
#===============================

# Run on statgen server
export DATA=/data/alh-admvhp1/GRERes/Data
export OUTPUT=/data/alh-admvhp1/GRERes/Data

cd $OUTPUT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E3.1) Check the phenotype

# phenotype file: edu_corrected

# raw data
# ln -s /data/alh-admhl/sfile_ukbb_general/info/ukb6247.csv ./ukb6247.csv

# get header
# awk 'NR==1' ukb6247.csv > ukb6247_header.csv

# R ------------------
# get only columns related to qualifications, 6138
R

cols <- c("6138-0.0", "6138-0.1", "6138-0.2", "6138-0.3", "6138-0.4", "6138-0.5", 
          "6138-1.0", "6138-1.1", "6138-1.2", "6138-1.3", "6138-1.4", "6138-1.5",
          "6138-2.0", "6138-2.1", "6138-2.2", "6138-2.3", "6138-2.4", "6138-2.5")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
# > trait_id
#  [1] 1583 1584 1585 1586 1587 1588 1589 1590 1591 1592 1593 1594 1595 1596 1597
# [16] 1598 1599 1600
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
# > head(dat)
#        V1       V2       V3       V4       V5       V6       V7       V8
# 1     eid 6138-0.0 6138-0.1 6138-0.2 6138-0.3 6138-0.4 6138-0.5 6138-1.0
# 2 2080762        3        4
# 3 3165676        1        2        3
# 4 1605581        3
# 5 5336523        1
# 6 2211334        2        3        6
#         V9      V10      V11      V12      V13      V14      V15      V16
# 1 6138-1.1 6138-1.2 6138-1.3 6138-1.4 6138-1.5 6138-2.0 6138-2.1 6138-2.2
# 2
# 3
# 4
# 5
# 6
#        V17      V18      V19
# 1 6138-2.3 6138-2.4 6138-2.5
# 2
# 3
# 4
# 5
# 6

writeLines("Getting traits ...")
sink("pheno_qualification_raw.txt")
write.table(dat[,], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

awk '{print $0} NR==6{exit}' edu_corrected
# [alh-admvhp1@hscpl-statgen Data]$ awk '{print $0} NR==6{exit}' edu_corrected
# "id" "edu"
# "1" 2080762 10
# "2" 3165676 20
# "3" 1605581 10
# "4" 5336523 20
# "5" 2211334 15

# Categories of Educational Qualification (data-field, 6138)	Code in UKBB	ISCED mapping1	Converted educational level (years)
# College or University degree	1	ISCED5	20
# A levels/AS levels or equivalent	2	ISCED3	13
# O levels/GCSEs or equivalent	3	ISCED2	10
# CSEs or equivalent	4	ISCED2	10
# NVQ or HND or HNC or equivalent	5	ISCED5	19
# Other professional qualifications e.g.: nursing, teaching	6	ISCED4	15
# None of the above	-7	ISCED1	7
# Prefer not to answer	-3	NA	NA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# E4) Qualification
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E4.1) Qualification & fit, Continuous data

# Input data
# .fam						fam file
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# Run on statgen server
export DATA=/data/alh-admvhp1/GRERes/Data
export OUTPUT=/data/alh-admvhp1/GRERes/Data

cd $OUTPUT

# # tango
# export DATA=/home/users/phamvv/GRERes/Data
# export OUTPUT=/home/users/phamvv/GRERes/Data
# 
# cd $OUTPUT

# # gadi
# export DATA=/scratch/eu82/vp8928/GRERes/gadi_50K_bi_cov0
# export OUTPUT=/scratch/eu82/vp8928/GRERes/gadi_50K_bi_cov0
# 
# cd $OUTPUT

    #sample_pheno.dat: phenotype
    
    # [alh-admvhp1@hscpl-statgen Data]$ wc -l selectedFam.fam
    # 50000 selectedFam.fam
    
    # extract phenotype
    pheno_file="edu_corrected"
    selected_id="selectedFam.fam"
    out_file="qualification_random_50K.dat"
    Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --selected_id=${selected_id} --out_file=${out_file}

    # model fitting
    # .	fit GREML:
    # -d qualification_random_50K.dat: phenotype
    ../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d qualification_random_50K.dat -mod 1 -thread 100 -out qualification_greml.out > qualification_greml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' qualification_greml.out >> res.greml.out
     
    # .	fit CORE GREML:
    # o	inside identity_coregreml.matlist:
    # grm_qced_rdm_005_ukbb3_selected_processed.grm
    # identity_selected_grm_bmat.chol.matmat2
    ../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d qualification_random_50K.dat -mod 1 -thread 100 -out qualification_coregreml.out > qualification_coregreml.log
    # awk 'n!=1 && $1!="LKH" {printf $2" "$3" "}; $1=="LKH" {n=1;printf $2"\n"}' qualification_coregreml.out >> res.coregreml.out
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E4.2) Qualification & fit, Continuous data + fixed effects

# Input data
# .fam						fam file
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# Run on statgen server
export DATA=/data/alh-admvhp1/GRERes/Data
export OUTPUT=/data/alh-admvhp1/GRERes/Data

cd $OUTPUT

# # tango
# export DATA=/home/users/phamvv/GRERes/Data
# export OUTPUT=/home/users/phamvv/GRERes/Data
# 
# cd $OUTPUT

# # gadi
# export DATA=/scratch/eu82/vp8928/GRERes/gadi_50K_bi_cov0
# export OUTPUT=/scratch/eu82/vp8928/GRERes/gadi_50K_bi_cov0
# 
# cd $OUTPUT

#sample_pheno.dat: phenotype

# [alh-admvhp1@hscpl-statgen Data]$ wc -l selectedFam.fam
# 50000 selectedFam.fam

# extract phenotype
pheno_file="edu_corrected"
other_pheno_file="ukb6247.csv" # year of birth, sex, age, etc.
covariate_file="my_sample_qcinfo.txt" # PCs
selected_id="selectedFam.fam"
out_file="qualification_random_50K_adjusted.dat"
Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --other_pheno_file=${other_pheno_file} --covariate_file=${covariate_file} --selected_id=${selected_id} --out_file=${out_file}

# update fam file
awk '{print $1,$2};' qualification_random_50K_adjusted.dat > updatedFam.fam

#------------------------------
R

fam <- read.table("selectedFam.fam", stringsAsFactors=F, header=F)
fam_id <- read.table("updatedFam.fam", stringsAsFactors=F, header=F)

filename = "updatedFam2.fam"
writeLines(paste("Writing file: ", filename, sep = ""))
sink(filename)
write.table(fam[which(fam$V2 %in% fam_id$V2),], row.names = F, quote = F, col.names = F)
sink()

q()
#------------------------------

# model fitting
# .	fit GREML:
# -d qualification_random_50K_adjusted.dat: phenotype
#../Tool/mtg2 -p updatedFam2.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d qualification_random_50K_adjusted.dat -mod 1 -thread 20 -out qualification_adjusted_greml.out > qualification_adjusted_greml.log
#../Tool/mtg2 -p updatedFam2.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d qualification_random_50K_adjusted_scale.dat -mod 1 -thread 20 -out qualification_adjusted_greml_s.out > qualification_adjusted_greml_s.log
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d trait_adjusted_f -mod 1 -thread 20 -out qualification_adjusted_greml_f.out > qualification_adjusted_greml_f.log

# .	fit CORE GREML:
# o	inside identity_coregreml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# identity_selected_grm_bmat.chol.matmat2
#../Tool/mtg2 -p updatedFam2.fam -mg identity_coregreml.matlist -d qualification_random_50K_adjusted.dat -mod 1 -thread 20 -out qualification_adjusted_coregreml.out > qualification_adjusted_coregreml.log
# ../Tool/mtg2 -p updatedFam2.fam -mg identity_coregreml.matlist -d qualification_random_50K_adjusted_scale.dat -mod 1 -thread 20 -out qualification_adjusted_coregreml_s.out > qualification_adjusted_coregreml_s.log
../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d trait_adjusted_f -mod 1 -thread 20 -out qualification_adjusted_coregreml_f.out > qualification_adjusted_coregreml_f.log

#++++++++++++++++++++++++
# o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
# These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
#grep -vwE '(LKH|h2)' qualification_adjusted_coregreml.out > qualification_adjusted_coregreml.out2
grep -vwE '(LKH|h2)' qualification_adjusted_coregreml_f.out > qualification_adjusted_coregreml_f.out2

# o	Prepare a do file for MTG2 to construct functions of model parameters.
#Rscript ../createTemplate.R --inFile="qualification_adjusted_coregreml.out2" --outFile="qualification_adjusted_coregreml.do"
Rscript ../createTemplate.R --inFile="qualification_adjusted_coregreml_f.out2" --outFile="qualification_adjusted_coregreml_f.do"

# .	Estimation:
#../Tool/mtg2 -delta2 qualification_adjusted_coregreml.do > qualification_adjusted_coregreml.out3
../Tool/mtg2 -delta2 qualification_adjusted_coregreml_f.do > qualification_adjusted_coregreml_f.out3

#++++++++++++++++++++++++
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Momin's script
#------------------------------

# model fitting
# .	fit GREML:
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d trait_adjusted -mod 1 -thread 20 -out qualification_adjusted_greml_m.out > qualification_adjusted_greml_m.log

# .	fit CORE GREML:
../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d trait_adjusted -mod 1 -thread 20 -out qualification_adjusted_coregreml_m.out > qualification_adjusted_coregreml_m.log

#++++++++++++++++++++++++

grep -vwE '(LKH|h2)' qualification_adjusted_coregreml_m.out > qualification_adjusted_coregreml_m.out2

# o	Prepare a do file for MTG2 to construct functions of model parameters.
Rscript ../createTemplate.R --inFile="qualification_adjusted_coregreml_m.out2" --outFile="qualification_adjusted_coregreml_m.do"

# .	Estimation:
../Tool/mtg2 -delta2 qualification_adjusted_coregreml_m.do > qualification_adjusted_coregreml_m.out3

#++++++++++++++++++++++++
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E4.3) Qualification & fit, Binary data + fixed effects

# Run on statgen server
export DATA=/data/alh-admvhp1/GRERes/Data
export OUTPUT=/data/alh-admvhp1/GRERes/Data

cd $OUTPUT

# extract phenotype
pheno_file="edu_corrected"
other_pheno_file="ukb6247.csv" # year of birth, sex, age, etc.
covariate_file="my_sample_qcinfo.txt" # PCs
selected_id="selectedFam.fam"
out_file="qualification_random_50K_adjusted_bi.dat"
binary=T
Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --other_pheno_file=${other_pheno_file} --covariate_file=${covariate_file} --selected_id=${selected_id} --out_file=${out_file} --binary=${binary}

# model fitting
# .	fit GREML:
# -d qualification_random_50K_adjusted_bi.dat: phenotype
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d trait_adjusted_bi_f -mod 1 -thread 20 -out qualification_adjusted_bi_greml_f.out > qualification_adjusted_bi_greml_f.log
# .	fit CORE GREML:
# o	inside identity_coregreml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# identity_selected_grm_bmat.chol.matmat2
../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d trait_adjusted_bi_f -mod 1 -thread 20 -out qualification_adjusted_bi_coregreml_f.out > qualification_adjusted_bi_coregreml_f.log

#++++++++++++++++++++++++
# o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
# These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
grep -vwE '(LKH|h2)' qualification_adjusted_bi_coregreml_f.out > qualification_adjusted_bi_coregreml_f.out2

# o	Prepare a do file for MTG2 to construct functions of model parameters.
Rscript ../createTemplate.R --inFile="qualification_adjusted_bi_coregreml_f.out2" --outFile="qualification_adjusted_bi_coregreml_f.do"

# .	Estimation:
../Tool/mtg2 -delta2 qualification_adjusted_bi_coregreml_f.do > qualification_adjusted_bi_coregreml_f.out3

#++++++++++++++++++++++++

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================

#===============================
# E6) Traits
#===============================

# Standing height 50 50-0.0	50-1.0	50-2.0 (11 12 13)
# Sitting height 20015 20015-0.0	20015-1.0	20015-2.0 (2303 2304 2305)
# *Heel bone mineral density 3148 3148-0.0 (680)
# *Heel bone mineral density (BMD) T-score, automated 78 78-0.0 (22)
# Weight 21002 21002-0.0	21002-1.0	21002-2.0 (2348 2349 2350)
# Fluid intelligence 20016 20016-0.0	20016-1.0	20016-2.0 (2306 2307 2308)
# BMI 21001 21001-0.0	21001-1.0 (2346 2347)
# Hip circumference 49 49-0.0	49-1.0	49-2.0 (8 9 10)
# Waist circumference 48 48-0.0	48-1.0	48-2.0 (5 6 7)
# Diastolic blood pressure 4079 4079-0.0	4079-0.1	4079-1.0	4079-1.1	4079-2.0	4079-2.1 (774 775 776 777 778 779)
# Body fat percentage (23099) ? (ukb29455) 23099-0.0	23099-1.0 (1437 1438)
# Overall health rating (2178) ??
# Vascular/heart problems diagnosed by doctor: High blood pressure (6150_4) 6150-0.0	6150-0.1	6150-0.2	6150-0.3	6150-1.0	6150-1.1	6150-1.2	6150-1.3	6150-2.0	6150-2.1	6150-2.2	6150-2.3 (1712 -> 1723)
# Non-cancer illness code, self-reported: hypertension (20002_1065) 87 from 20002-0.0 to 20002-2.28 (1901 -> 1987)
# Systolic blood pressure, automated reading (4080) 4080-0.0 -> 4080-2.1 (780 -> 785)
# Ever smoked (20160) 20160-0.0 (2340)
# Age when periods started (menarche) (2714) 2714-0.0	2714-1.0	2714-2.0 (618 619 620)
# Age at first live birth (2754) (2754-0.0	2754-1.0	2754-2.0) (627 628 629)
# Alcohol intake frequency (1558) ??? 1558-0.0	1558-1.0	1558-2.0 (441 442 443)
# Mood swings (1920) (1920-0.0	1920-1.0	1920-2.0) 495 496 497
# Neuroticism score (20127) ??? 20127-0.0 (2339)
# Pulse rate (4194) ? (ukb43545) 4194-0.0	4194-1.0	4194-2.0	4194-3.0 (14 15 16 17)
# Happiness (4526) ??? (4526-0.0	4526-1.0	4526-2.0) (1418 1419 1420)

# Run on statgen server
export DATA=/data/alh-admvhp1/GRERes/Data
export OUTPUT=/data/alh-admvhp1/GRERes/Data

cd $OUTPUT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.1) Get phenotype - Standing height

# raw data
# ln -s /data/alh-admhl/sfile_ukbb_general/info/ukb6247.csv ./ukb6247.csv

# get header
# awk 'NR==1' ukb6247.csv > ukb6247_header.csv

# R ------------------
# get only columns related to Standing height 50 50-0.0	50-1.0	50-2.0 (11 12 13)
R

cols <- c("50-0.0", "50-1.0", "50-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 11 12 13
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1     V2     V3     V4
# 1     eid 50-0.0 50-1.0 50-2.0
# 2 2080762  165.5
# 3 3165676    188
# 4 1605581    168
# 5 5336523    173
# 6 2211334    178

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("sh_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.2) Get phenotype - Sitting height 20015 20015-0.0	20015-1.0	20015-2.0 (2303 2304 2305)

# R ------------------
# get only related columns
R

cols <- c("20015-0.0", "20015-1.0",	"20015-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 2303 2304 2305
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2        V3        V4
# 1     eid 20015-0.0 20015-1.0 20015-2.0
# 2 2080762        84
# 3 3165676        99
# 4 1605581        99
# 5 5336523        95
# 6 2211334        98

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("sith_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.3) Get phenotype - Heel bone mineral density 3148 3148-0.0 680

# R ------------------
# get only related columns
R

cols <- c("3148-0.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 680
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2
# 1     eid 3148-0.0
# 2 2080762     0.68
# 3 3165676
# 4 1605581
# 5 5336523    0.611
# 6 2211334    0.566

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("heel_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.3*) Get phenotype - Heel bone mineral density (BMD) T-score, automated 78 78-0.0 (22)

# R ------------------
# get only related columns
R

cols <- c("78-0.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 22
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1     V2
# 1     eid 78-0.0
# 2 2080762   0.88
# 3 3165676
# 4 1605581
# 5 5336523  0.265
# 6 2211334 -0.136

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("heelauto_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.4) Get phenotype - Weight 21002 21002-0.0	21002-1.0	21002-2.0 (2348 2349 2350)

# R ------------------
# get only related columns
R

cols <- c("21002-0.0", "21002-1.0", "21002-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 2348 2349 2350
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2        V3        V4
# 1     eid 21002-0.0 21002-1.0 21002-2.0
# 2 2080762      56.1
# 3 3165676      87.3
# 4 1605581      77.8
# 5 5336523      88.5
# 6 2211334      80.4

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("wei_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.5) Get phenotype - Fluid intelligence 20016 20016-0.0	20016-1.0	20016-2.0 (2306 2307 2308)

# R ------------------
# get only related columns
R

cols <- c("20016-0.0", "20016-1.0", "20016-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 2306 2307 2308
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2        V3        V4
# 1     eid 20016-0.0 20016-1.0 20016-2.0
# 2 2080762
# 3 3165676        10
# 4 1605581         8
# 5 5336523
# 6 2211334

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("fluid_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.6) Get phenotype - BMI 21001 21001-0.0	21001-1.0 (2346 2347)

# R ------------------
# get only related columns
R

cols <- c("21001-0.0", "21001-1.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 2346 2347
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2        V3
# 1     eid 21001-0.0 21001-1.0
# 2 2080762   20.4817
# 3 3165676   24.7001
# 4 1605581   27.5652
# 5 5336523     29.57
# 6 2211334   25.3756

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("bmi_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.7) Get phenotype - Hip circumference 49 49-0.0	49-1.0	49-2.0 (8 9 10)

# R ------------------
# get only related columns
R

cols <- c("49-0.0",	"49-1.0", "49-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1]  8  9 10
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1     V2     V3     V4
# 1     eid 49-0.0 49-1.0 49-2.0
# 2 2080762     88
# 3 3165676     98
# 4 1605581    104
# 5 5336523    110
# 6 2211334    101

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("hip_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.8) Get phenotype - Waist circumference 48 48-0.0	48-1.0	48-2.0 (5 6 7)

# R ------------------
# get only related columns
R

cols <- c("48-0.0", "48-1.0", "48-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 5 6 7
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1     V2     V3     V4
# 1     eid 48-0.0 48-1.0 48-2.0
# 2 2080762     78
# 3 3165676     89
# 4 1605581     96
# 5 5336523     93
# 6 2211334     91

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("wai_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.9) Get phenotype - Diastolic blood pressure 4079 4079-0.0	4079-0.1	4079-1.0	4079-1.1	4079-2.0	4079-2.1 (774 775 776 777 778 779)

# R ------------------
# get only related columns
R

cols <- c("4079-0.0", "4079-0.1", "4079-1.0", "4079-1.1", "4079-2.0", "4079-2.1")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 774 775 776 777 778 779
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4       V5       V6       V7
# 1     eid 4079-0.0 4079-0.1 4079-1.0 4079-1.1 4079-2.0 4079-2.1
# 2 2080762       71       69
# 3 3165676       78       78
# 4 1605581       77       77
# 5 5336523       81       83
# 6 2211334       78       79

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("blood_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.10) Get phenotype - # Body fat percentage (23099) ? (ukb29455) 23099-0.0	23099-1.0 (1437 1438)

# R ------------------
# get only related columns
R

cols <- c("23099-0.0", "23099-1.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb29455.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 1437 1438
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb29455.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2        V3
# 1     eid 23099-0.0 23099-1.0
# 2 5306847      35.1
# 3 5916772      15.5
# 4 1714528      31.4
# 5 4585896      34.1
# 6 5806963      20.5

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("10body_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.11) Get phenotype - Overall health rating (2178) ??

# R ------------------
# get only related columns
R

cols <- c("4079-0.0", "4079-0.1", "4079-1.0", "4079-1.1", "4079-2.0", "4079-2.1")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 774 775 776 777 778 779
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4       V5       V6       V7
# 1     eid 4079-0.0 4079-0.1 4079-1.0 4079-1.1 4079-2.0 4079-2.1
# 2 2080762       71       69
# 3 3165676       78       78
# 4 1605581       77       77
# 5 5336523       81       83
# 6 2211334       78       79

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("blood_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.12) (Categorical) Get phenotype - # Vascular/heart problems diagnosed by doctor: High blood pressure (6150_4) 6150-0.0	6150-0.1	6150-0.2	6150-0.3	6150-1.0	6150-1.1	6150-1.2	6150-1.3	6150-2.0	6150-2.1	6150-2.2	6150-2.3 (1712 -> 1723)

# R ------------------
# get only related columns
R

cols <- c("6150-0.0", "6150-0.1", "6150-0.2", "6150-0.3", "6150-1.0", "6150-1.1", "6150-1.2", "6150-1.3", "6150-2.0", "6150-2.1", "6150-2.2", "6150-2.3")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
#  [1] 1712 1713 1714 1715 1716 1717 1718 1719 1720 1721 1722 1723
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4       V5       V6       V7       V8
# 1     eid 6150-0.0 6150-0.1 6150-0.2 6150-0.3 6150-1.0 6150-1.1 6150-1.2
# 2 2080762        4
# 3 3165676       -7
# 4 1605581        4
# 5 5336523       -7
# 6 2211334        3

# Coding	Meaning
# 1	Heart attack
# 2	Angina
# 3	Stroke
# 4	High blood pressure
# -7	None of the above
# -3	Prefer not to answer
dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-3,NA,dat[,2])

writeLines("Getting traits ...")
sink("12heart_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.13) (Categorical) 13 Non-cancer illness code, self-reported: hypertension (20002_1065) 87 from 20002-0.0 to 20002-2.28 (1901 -> 1987)

# R ------------------
# get only related columns
R

# cols <- c("4079-0.0", "4079-0.1", "4079-1.0", "4079-1.1", "4079-2.0", "4079-2.1")
# header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
# trait_id = match(cols, header_trait_file)
# trait_id
# # > trait_id
# # [1] 774 775 776 777 778 779
trait_id <- c(1901:1987)
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2        V3        V4        V5        V6        V7        V8
# 1     eid 20002-0.0 20002-0.1 20002-0.2 20002-0.3 20002-0.4 20002-0.5 20002-0.6
# 2 2080762      1065
# 3 3165676
# 4 1605581      1065      1111
# 5 5336523
# 6 2211334      1082

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
for (i in 1:nrow(dat)) {
  if(!is.na(dat[i,2])) {
    if (dat[i,2] == 1065) {
      dat[i,2] <- 1
    } else {
      dat[i,2] <- 0
    }
  }
}

writeLines("Getting traits ...")
sink("13non_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.14) # 14 Systolic blood pressure, automated reading (4080) 4080-0.0 -> 4080-2.1 (780 -> 785)

# R ------------------
# get only related columns
R

cols <- c("4080-0.0", "4080-0.1", "4080-1.0", "4080-1.1", "4080-2.0", "4080-2.1")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 780 781 782 783 784 785
trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4       V5       V6       V7
# 1     eid 4080-0.0 4080-0.1 4080-1.0 4080-1.1 4080-2.0 4080-2.1
# 2 2080762      142      137
# 3 3165676      128      127
# 4 1605581      145      129
# 5 5336523      128      124
# 6 2211334      126      135

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("14blood_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.15) (Categorical) # 15 Ever smoked (20160) 20160-0.0 (2340)

# R ------------------
# get only related columns
R

cols <- c("20160-0.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 2340

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2
# 1     eid 20160-0.0
# 2 2080762         0
# 3 3165676         1
# 4 1605581         1
# 5 5336523         0
# 6 2211334         1

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("15ever_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.16) # 16 Age when periods started (menarche) (2714) 2714-0.0	2714-1.0	2714-2.0 (618 619 620)

# R ------------------
# get only related columns
R

cols <- c("2714-0.0", "2714-1.0", "2714-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 618 619 620

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4
# 1     eid 2714-0.0 2714-1.0 2714-2.0
# 2 2080762
# 3 3165676
# 4 1605581       11
# 5 5336523
# 6 2211334

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("16age_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.17) # 17 Age at first live birth (2754) (2754-0.0	2754-1.0	2754-2.0) (627 628 629)

# R ------------------
# get only related columns
R

cols <- c("2754-0.0", "2754-1.0", "2754-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 627 628 629

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > dat[1:20,]
#         V1       V2       V3       V4
# 1      eid 2754-0.0 2754-1.0 2754-2.0
# 2  2080762
# 3  3165676
# 4  1605581
# 5  5336523
# 6  2211334
# 7  5391704
# 8  3248427       29
# 9  4646409
# 10 1041921       28

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("17age_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.18) (Categorical) # 18 Alcohol intake frequency (1558) ??? 1558-0.0	1558-1.0	1558-2.0 (441 442 443)

# Coding	Meaning
# 1	Daily or almost daily
# 2	Three or four times a week
# 3	Once or twice a week
# 4	One to three times a month
# 5	Special occasions only
# 6	Never
# -3	Prefer not to answer

# R ------------------
# get only related columns
R

cols <- c("1558-0.0", "1558-1.0", "1558-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 441 442 443

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4
# 1     eid 1558-0.0 1558-1.0 1558-2.0
# 2 2080762        2
# 3 3165676        2
# 4 1605581        2
# 5 5336523        2
# 6 2211334        1

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-3,NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==6,0,dat[,2])

writeLines("Getting traits ...")
sink("18alcohol_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.19) (Categorical) # 19 Mood swings (1920) (1920-0.0	1920-1.0	1920-2.0) 495 496 497
# 1	Yes
# 0	No
# -1	Do not know
# -3	Prefer not to answer

# R ------------------
# get only related columns
R

cols <- c("1920-0.0", "1920-1.0", "1920-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 495 496 497

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4
# 1     eid 1920-0.0 1920-1.0 1920-2.0
# 2 2080762        0
# 3 3165676        0
# 4 1605581        0
# 5 5336523        0
# 6 2211334        0

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-1,NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-3,NA,dat[,2])

writeLines("Getting traits ...")
sink("19mood_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.20) # 20 Neuroticism score (20127) ??? 20127-0.0 (2339)

# R ------------------
# get only related columns
R

cols <- c("20127-0.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 2339

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1        V2
# 1     eid 20127-0.0
# 2 2080762
# 3 3165676         1
# 4 1605581         1
# 5 5336523         1
# 6 2211334         0

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("20score_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.21) # 21 Pulse rate (4194) ? (ukb43545) 4194-0.0	4194-1.0	4194-2.0	4194-3.0 (14 15 16 17)

# R ------------------
# get only related columns
R

cols <- c("4194-0.0", "4194-1.0", "4194-2.0", "4194-3.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb43545.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
# > trait_id
# [1] 14 15 16 17

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb43545.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4       V5
# 1     eid 4194-0.0 4194-1.0 4194-2.0 4194-3.0
# 2 1000019
# 3 1000022
# 4 1000035       51
# 5 1000046
# 6 1000054       70

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])

writeLines("Getting traits ...")
sink("21rate_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-1.22) (Categorical) # 22 Happiness (4526) ??? (4526-0.0	4526-1.0	4526-2.0) (1418 1419 1420)

# Coding	Meaning
# 1	Extremely happy
# 2	Very happy
# 3	Moderately happy
# 4	Moderately unhappy
# 5	Very unhappy
# 6	Extremely unhappy
# -1	Do not know
# -3	Prefer not to answer

# R ------------------
# get only related columns
R

cols <- c("4526-0.0", "4526-1.0", "4526-2.0")
header_trait_file = as.character(read.csv(pipe(paste("head -n1 ", "ukb6247.csv", sep ="")), stringsAsFactors=F, header=F))
trait_id = match(cols, header_trait_file)
trait_id
# > trait_id
# [1] 1418 1419 1420

trait_id_text = paste("$", trait_id, sep="", collapse=",")
cmd = paste("awk -F, 'BEGIN {OFS=\",\"} {print $1,", trait_id_text, "}' ", "ukb6247.csv", sep = "")
dat = read.csv(pipe(cmd), stringsAsFactors=F, header=F) 
head(dat)
# > head(dat)
#        V1       V2       V3       V4
# 1     eid 4526-0.0 4526-1.0 4526-2.0
# 2 2080762
# 3 3165676        2
# 4 1605581        2
# 5 5336523
# 6 2211334

dat[,2] <- ifelse(dat[,2]=="",NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==4,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==5,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==6,0,dat[,2])
dat[,2] <- ifelse(dat[,2]==-1,NA,dat[,2])
dat[,2] <- ifelse(dat[,2]==-3,NA,dat[,2])

writeLines("Getting traits ...")
sink("21happiness_trait.txt")
write.table(dat[,1:2], quote=F, sep="\t", row.names=F, col.names=F)
sink()

q()
# R ------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-2) Standing height
# .dat: phenotypic data of the main trait

#sample_pheno.dat: phenotype

# [alh-admvhp1@hscpl-statgen Data]$ wc -l selectedFam.fam
# 50000 selectedFam.fam

# extract phenotype
#pre="sh" # ***update this*** Standing height
#pre="sith" # sitting height
#pre="heel" # E6-1.3) Get phenotype - Heel bone mineral density 3148 3148-0.0 680
#pre="heelauto" # E6-1.3*) Get phenotype - Heel bone mineral density (BMD) T-score, automated 78 78-0.0 (22)
#pre="wei" # E6-1.4) Get phenotype - Weight 21002 21002-0.0	21002-1.0	21002-2.0 (2348 2349 2350)
#pre="fluid" # E6-1.5) Get phenotype - Fluid intelligence 20016 20016-0.0	20016-1.0	20016-2.0 (2306 2307 2308)
# pre="bmi" # E6-1.6) Get phenotype - BMI 21001 21001-0.0	21001-1.0 (2346 2347)
#pre="hip" # E6-1.7) Get phenotype - Hip circumference 49 49-0.0	49-1.0	49-2.0 (8 9 10)
# pre="wai" # E6-1.8) Get phenotype - Waist circumference 48 48-0.0	48-1.0	48-2.0 (5 6 7)
#pre="blood" # E6-1.9) Get phenotype - Diastolic blood pressure 4079 4079-0.0	4079-0.1	4079-1.0	4079-1.1	4079-2.0	4079-2.1 (774 775 776 777 778 779)
#pre="10body" # Body fat percentage (23099) ? (ukb29455) 23099-0.0	23099-1.0 (1437 1438)
#pre="14blood" # E6-1.14) # 14 Systolic blood pressure, automated reading (4080) 4080-0.0 -> 4080-2.1 (780 -> 785)
#pre="16age" # E6-1.16) # 16 Age when periods started (menarche) (2714) 2714-0.0	2714-1.0	2714-2.0 (618 619 620) , useSex=F
#pre="17age" # E6-1.17) # 17 Age at first live birth (2754) (2754-0.0	2754-1.0	2754-2.0) (627 628 629)  , useSex=F
#pre="20score" # E6-1.20) # 20 Neuroticism score (20127) ??? 20127-0.0 (2339)
pre="21rate" # E6-1.21) # 21 Pulse rate (4194) ? (ukb43545) 4194-0.0	4194-1.0	4194-2.0	4194-3.0 (14 15 16 17)
pheno_file="${pre}_trait.txt"
# ln -s /data/alh-admhl/sfile_ukbb_general/info/ukb6247.csv /data/alh-admvhp1/GRERes/Data/ukb6247.csv
# ln -s /mnt/rdfs/WholeGenomeApproach/ukbb3/my_sample_qcinfo.txt /data/alh-admvhp1/GRERes/Data/my_sample_qcinfo.txt
other_pheno_file="ukb6247.csv" # year of birth, sex, age, etc.
covariate_file="my_sample_qcinfo.txt" # PCs
selected_id="selectedFam.fam"
out_file="${pre}_random_50K_adjusted.dat"
useSex=T
Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --other_pheno_file=${other_pheno_file} --covariate_file=${covariate_file} --selected_id=${selected_id} --out_file=${out_file} --pre=${pre} --useSex=${useSex}
# sh_random_50K_adjusted.dat

# model fitting
# .	fit GREML:
# -d qualification_random_50K_adjusted.dat: phenotype
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d ${pre}_trait_adjusted_f -mod 1 -thread 20 -out ${pre}_adjusted_greml_f.out > ${pre}_adjusted_greml_f.log

# .	fit CORE GREML:
# o	inside identity_coregreml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# identity_selected_grm_bmat.chol.matmat2
../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d ${pre}_trait_adjusted_f -mod 1 -thread 20 -out ${pre}_adjusted_coregreml_f.out > ${pre}_adjusted_coregreml_f.log

#++++++++++++++++++++++++
# o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
# These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
grep -vwE '(LKH|h2)' ${pre}_adjusted_coregreml_f.out > ${pre}_adjusted_coregreml_f.out2

# o	Prepare a do file for MTG2 to construct functions of model parameters.
Rscript ../createTemplate.R --inFile="${pre}_adjusted_coregreml_f.out2" --outFile="${pre}_adjusted_coregreml_f.do"

# .	Estimation:
../Tool/mtg2 -delta2 ${pre}_adjusted_coregreml_f.do > ${pre}_adjusted_coregreml_f.out3

#++++++++++++++++++++++++
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E6-2*) For binary
# .dat: phenotypic data of the main trait

#sample_pheno.dat: phenotype

# [alh-admvhp1@hscpl-statgen Data]$ wc -l selectedFam.fam
# 50000 selectedFam.fam

# extract phenotype
#pre="12heart" # E6-1.12) (Categorical) Get phenotype - # Vascular/heart problems diagnosed by doctor: High blood pressure (6150_4) 6150-0.0	6150-0.1	6150-0.2	6150-0.3	6150-1.0	6150-1.1	6150-1.2	6150-1.3	6150-2.0	6150-2.1	6150-2.2	6150-2.3 (1712 -> 1723)
#pre="13non" # E6-1.13) (Categorical) 13 Non-cancer illness code, self-reported: hypertension (20002_1065) 87 from 20002-0.0 to 20002-2.28 (1901 -> 1987)
#pre="15ever" # E6-1.15) (Categorical) # 15 Ever smoked (20160) 20160-0.0 (2340)
#pre="18alcohol" # E6-1.18) (Categorical) # 18 Alcohol intake frequency (1558) ??? 1558-0.0	1558-1.0	1558-2.0 (441 442 443)
#pre="19mood" # E6-1.19) (Categorical) # 19 Mood swings (1920) (1920-0.0	1920-1.0	1920-2.0) 495 496 497
pre="21happiness" # E6-1.22) (Categorical) # 22 Happiness (4526) ??? (4526-0.0	4526-1.0	4526-2.0) (1418 1419 1420)
pheno_file="${pre}_trait.txt"
# ln -s /data/alh-admhl/sfile_ukbb_general/info/ukb6247.csv /data/alh-admvhp1/GRERes/Data/ukb6247.csv
# ln -s /mnt/rdfs/WholeGenomeApproach/ukbb3/my_sample_qcinfo.txt /data/alh-admvhp1/GRERes/Data/my_sample_qcinfo.txt
other_pheno_file="ukb6247.csv" # year of birth, sex, age, etc.
covariate_file="my_sample_qcinfo.txt" # PCs
selected_id="selectedFam.fam"
out_file="${pre}_random_50K_adjusted_bi.dat"
binary=T
#thres=4
thres=1
Rscript ../selectPhenotypes.R --pheno_file=${pheno_file} --other_pheno_file=${other_pheno_file} --covariate_file=${covariate_file} --selected_id=${selected_id} --out_file=${out_file}  --pre=${pre} --binary=${binary} --thres=${thres}
# sh_random_50K_adjusted.dat

# model fitting
# .	fit GREML:
# -d qualification_random_50K_adjusted.dat: phenotype
../Tool/mtg2 -p selectedFam.fam -g grm_qced_rdm_005_ukbb3_selected_processed.grm -d ${pre}_trait_adjusted_bi_f -mod 1 -thread 20 -out ${pre}_adjusted_bi_greml_f.out > ${pre}_adjusted_bi_greml_f.log

# .	fit CORE GREML:
# o	inside identity_coregreml.matlist:
# grm_qced_rdm_005_ukbb3_selected_processed.grm
# identity_selected_grm_bmat.chol.matmat2
../Tool/mtg2 -p selectedFam.fam -mg identity_coregreml.matlist -d ${pre}_trait_adjusted_bi_f -mod 1 -thread 20 -out ${pre}_adjusted_bi_coregreml_f.out > ${pre}_adjusted_bi_coregreml_f.log

#++++++++++++++++++++++++
# o	Prepare a file that contains variance & covariance estimates & the information matrix (i.e., var-cov matrix of model parameters).
# These pieces of information are inside the CORE GREML outfile, but to be readily used by MTG2, lines with 'LKH' & 'h2' need to be removed from the outfile using the following code:
grep -vwE '(LKH|h2)' ${pre}_adjusted_bi_coregreml_f.out > ${pre}_adjusted_bi_coregreml_f.out2

# o	Prepare a do file for MTG2 to construct functions of model parameters.
Rscript ../createTemplate.R --inFile="${pre}_adjusted_bi_coregreml_f.out2" --outFile="${pre}_adjusted_bi_coregreml_f.do"

# .	Estimation:
../Tool/mtg2 -delta2 ${pre}_adjusted_bi_coregreml_f.do > ${pre}_adjusted_bi_coregreml_f.out3

#++++++++++++++++++++++++
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#===============================
#===============================
