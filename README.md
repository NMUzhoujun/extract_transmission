1、默认（输出等位基因编号）：

./extract_transmission.sh \
  -v cohort.bcf \
  -p trios.ped \
  -o transmissions.tsv.gz \
  -t 8


2、额外输出碱基（PT_base/MT_base/PU_base/MU_base）：

./extract_transmission.sh \
  -v cohort.vcf.gz \
  -p trios.ped \
  --emit-bases \
  -o transmissions_with_bases.tsv.gz


3、只跑某个染色体区间，提高吞吐/并行更方便：

./extract_transmission.sh -v cohort.bcf -p trios.ped -r chr1 \
  -o chr1.transmissions.tsv.gz -t 8
