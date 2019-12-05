#!/bin/bash
VCF_PATH=$1
for vcf_file in ${VCF_PATH}/*.vcf; do
  bname=$(basename $vcf_file | cut -f1 -d'_');
  # check if any variant called near the snp of interest position
  cond=$(grep -v '^#' ${vcf_file} | wc -l);
  if [ $cond -eq 0 ]; then
    echo $cond |
      awk -v OFS="\t" -v bn="$bname" '{print bn, 36993455, "G", "G", "NA", "NA"}';
  else
    grep -v '^#' ${vcf_file} |
      sed 's/[^\t]*TC=\([^;]*\);/\1\t/g; s/[^\t]*TR=\([^;]*\);.*/\1\t/g;' |
        cut -f2,4,5,8,9 |
          awk -v OFS="\t" -v bn="$bname" '{print bn, $0}';
  fi
done 
