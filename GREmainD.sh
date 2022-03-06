#===============================
# For 10K samples, g+e
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

# 1.2) Select 10K samples
fam_file="./qced_rdm_005_ukbb3.fam"
Rscript ../selectSamples.R --fam_file=${fam_file} --number=10000 --outFile1="selectedFam.fam"

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
# D1) Continuous data
#===============================

# Input data
# .fam						fam file
# .grm						genomic relationship matrix
# .dat						phenotypic data of the main trait

# Run on statgen server
export DATA=/data/alh-admvhp1/GRERes/Data
export OUTPUT=/data/alh-admvhp1/GRERes/Data

cd $OUTPUT

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# D1.1) Simulate phenotype data & fit
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
    Rscript ../random.R --cov=${cov} --seed=${i} --number=10000
    
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
    Rscript ../random.R --cov=${cov} --seed=${i} --number=10000
    
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
    Rscript ../random.R --cov=${cov} --seed=${i} --number=10000
    
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
# D2) For binary data, correct the transformation, scale
#===============================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# D2.1) Simulate phenotype data & fit
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
    Rscript ../random.R --cov=${cov} --seed=${i} --number=10000
    
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
    Rscript ../random.R --cov=${cov} --seed=${i} --number=10000
    
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
    Rscript ../random.R --cov=${cov} --seed=${i} --number=10000
    
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
