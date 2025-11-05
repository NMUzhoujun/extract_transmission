1、默认（输出等位基因编号）：
```
./extract_transmission.sh \
  -v cohort.bcf \
  -p trios.ped \
  -o transmissions.tsv.gz \
  -t 8
```

2、额外输出碱基（PT_base/MT_base/PU_base/MU_base）：
```
./extract_transmission.sh \
  -v cohort.vcf.gz \
  -p trios.ped \
  --emit-bases \
  -o transmissions_with_bases.tsv.gz
```

3、只跑某个染色体区间，提高吞吐/并行更方便：
```
./extract_transmission.sh -v cohort.bcf -p trios.ped -r chr1 \
  -o chr1.transmissions.tsv.gz -t 8
```

##Setp2：提取指定变异及其对应基因型
```
bash ../tools/hapvcf_extract_raw_byAllele.sh \
  -v /public/home/zj2020/phasing/phasing_res/check/test2.PT.vcf.gz \
  -l var_a1.txt \
  -o test2.raw
```

2) 只导出部分样本（sample.list 每行一个样本名）：
```
bash hapvcf_extract_raw_byAllele.sh \
  -v /public/.../transmissions.chr9.MT.vcf.gz \
  -l var_a1.txt \
  -s /public/home/zj2020/phasing/sample.list \
  -o MT_chr9_customA1_subset.raw
```
##Step3：
```
~/anaconda3/envs/zj/bin/Rscript ../tools/merge_raws_and_prs.R \
  effects.tsv \
  PRS_output.tsv \
  merged_genotype.tsv \
  /public/home/zj2020/phasing/phasing_res/check
```


##转换rsid为chr:bp:ref:alt
```
rs2pos_min_whitelist.sh \
  -i /public/home/zj2020/PGS_analysis/taiwan/download_data/${file}.gz \
  -d /public/home/zj2020/StandTools/annovar/annovar/humandb/hg19_Gnomad.txt \
  -v /public/project/BC0001/BC8000/BC8000_QC/new_model/all_pass_combine.vid \
  -m /public/project/BC0001/BC8000/results/vcfnorm_CNBC/vcf_split2_anno/allvid_maf_parents.tsv \
  -o /public/home/zj2020/PGS_analysis/taiwan/output_vid/${phenotype}
```

