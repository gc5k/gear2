#!/bin/bash
prefix=$1
thread_num=$2
xld=$3
block_size=$4
chr=$(awk '{print $1}' ${prefix}.bim | sort -u | wc -l)

### High-resolution LD
for i in $(seq ${chr})
do
    awk -v OFS='\t' -v chr=$i -v size=$block_size '$1 == chr {count++; if(count%size==1){block++}; print $2, chr"_"block}' ${prefix}.bim > ${prefix}_chr${i}.tmp
done

cat ${prefix}_chr*.tmp > ${prefix}_block${block_size}.txt
rm ${prefix}_chr*.tmp

### 判断xld参数
if [[ $xld -eq 0 ]]
then
    ../build/gear --bfile ${prefix} --snp-tag ${prefix}_block${block_size}.txt --xld --threads ${thread_num} --out ${prefix}_${block_size} > ${prefix}_${block_size}.log 2>&1
else
    ../build/gear --bfile ${prefix} --snp-tag ${prefix}_block${block_size}.txt --xld --xld-alg 1 --iter 100 --threads ${thread_num} --out ${prefix}_${block_size} > ${prefix}_${block_size}.log 2>&1
fi


###Norm I and Norm II
for j in `seq ${chr}`
do
   plink2 --bfile ${prefix} --chr-set ${chr} --chr ${j} --silent --make-bed --out ${prefix}_${j}
   ../build/gear --bfile ${prefix}_${j} --propc --threads ${thread_num} --max-it 50 --evec 30 --out ${prefix}_${j} > ${prefix}_${j}.log 2>&1
done

if [[ $xld -eq 0 ]]
then
    ../build/gear --bfile ${prefix} --xld --threads ${thread_num} --out ${prefix}
else
    ../build/gear --bfile ${prefix} --xld --xld-alg 1 --iter 100 --threads ${thread_num} --out ${prefix}
fi


if [[ $xld -eq 0 ]]
then
    Rscript --no-save Atlas.R ${prefix}_${block_size}.xld ${prefix}.xld
else
    Rscript --no-save Atlas.R ${prefix}_${block_size}.xldr ${prefix}.xldr
fi
