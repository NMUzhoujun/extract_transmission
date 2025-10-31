#!/usr/bin/env bash
# 从相位多样本 VCF/BCF + PED( Family Child Father Mother ) 一次遍历生成四个 VCF
# 依赖：bcftools>=1.14, bgzip, tabix, awk
# 用法示例见文末

set -euo pipefail

VCF=""
PED=""
OUT_PREFIX="transmissions"
THREADS=4
REGION_STR=""
REGION_FILE=""

usage() {
  cat >&2 <<EOF
Usage:
  $0 -v cohort.phased.bcf -p trios.ped [-o out_prefix] [-t 8] [-r chrN|chrN:beg-end | -R regions.bed]

Outputs (bgzip+tabix):
  <prefix>.PT.vcf.gz   (GT=PT|.)
  <prefix>.MT.vcf.gz   (GT=.|MT)
  <prefix>.PU.vcf.gz   (GT=PU|.)
  <prefix>.MU.vcf.gz   (GT=.|MU)

Notes:
  - PED 前四列：Family  Child  Father  Mother（有表头也可）
  - 子代 GT 必须是相位的“a|b”，左=父传，右=母传
  - 缺失/未相位/父母缺失 => 写 "./."
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--vcf) VCF="$2"; shift 2;;
    -p|--ped) PED="$2"; shift 2;;
    -o|--out-prefix) OUT_PREFIX="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -r|--regions) REGION_STR="$2"; shift 2;;
    -R|--regions-file) REGION_FILE="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -z "$VCF" || -z "$PED" ]] && usage
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found" >&2; exit 1; }
command -v bgzip    >/dev/null || { echo "ERROR: bgzip not found" >&2; exit 1; }
command -v tabix    >/dev/null || { echo "ERROR: tabix not found" >&2; exit 1; }

# 要求 bcftools >=1.14
if ! bcftools --version 2>/dev/null | awk '/^bcftools/{split($2,a,"."); if (a[1]<1 || (a[1]==1 && a[2]<14)) exit 1}'; then
  echo "ERROR: require bcftools >=1.14" >&2; exit 2
fi

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# 1) 读 PED：保留 Family Child Father Mother
awk -v OFS='\t' 'BEGIN{FS="[ \t]"} NR==1&&($1~/[Ff]am/||$2~/[Cc]hild/){next} {if($2!=""&&$3!=""&&$4!="") print $1,$2,$3,$4}' \
  "$PED" > "$TMP/trios.tsv"

# 2) 样本集合：子代 + 父 + 母（去重）用于 VCF 子集
awk '{print $2; print $3; print $4}' "$TMP/trios.tsv" | awk '!seen[$0]++' > "$TMP/samples.unique.list"

# 3) 子代样本顺序（列顺序）
awk '{print $2}' "$TMP/trios.tsv" | awk '!seen[$0]++' > "$TMP/children.list"

# 4) 从原始 VCF 拿 contig 头，构建四个 VCF 的 header
CONTIG_HDR="$TMP/contigs.hdr"
bcftools view -h "$VCF" | awk '/^##contig=/' > "$CONTIG_HDR" || true

make_header() {
  local label="$1"
  local header="$TMP/header.$label.vcf"
  {
    echo "##fileformat=VCFv4.2"
    cat "$CONTIG_HDR"
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype; single-haplotype view for '"$label"'">'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    while read -r s; do printf '\t%s' "$s"; done < "$TMP/children.list"
    printf '\n'
  } > "$header"
}

make_header PT
make_header MT
make_header PU
make_header MU

# 5) 构建 bcftools view 参数（区域命名不匹配时自动跳过 -r）
VIEW_ARGS=(--threads "$THREADS" -S "$TMP/samples.unique.list" -Ou "$VCF")
if [[ -n "$REGION_STR" ]]; then
  mapfile -t CONTIGS < <(bcftools index -s "$VCF" | cut -f1)
  want="$REGION_STR"; nochr="${REGION_STR#chr}"; withchr="chr${REGION_STR}"
  hit=""
  for c in "${CONTIGS[@]}"; do
    [[ "$c" == "$want" || "$c" == "$nochr" || "$c" == "$withchr" ]] && { hit="$c"; break; }
  done
  [[ -n "$hit" ]] && VIEW_ARGS=(-r "$hit" "${VIEW_ARGS[@]}") || echo "[WARN] region '$REGION_STR' not found; skip -r" >&2
fi
[[ -n "$REGION_FILE" ]] && VIEW_ARGS=(-R "$REGION_FILE" "${VIEW_ARGS[@]}")

# 6) 一次遍历：query 出每个位点的 (CHROM POS ID REF ALT) + [SAMPLE,GT] 对；AWK 同时写四个 body 文件
FMT='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%SAMPLE\t%GT]\n'
( bcftools view "${VIEW_ARGS[@]}" | bcftools query -f "$FMT" - ) \
| awk -v OFS="\t" \
      -v TRIOS="$TMP/trios.tsv" \
      -v KIDS="$TMP/children.list" \
      -v BODY_PT="$TMP/body.PT.vcf" \
      -v BODY_MT="$TMP/body.MT.vcf" \
      -v BODY_PU="$TMP/body.PU.vcf" \
      -v BODY_MU="$TMP/body.MU.vcf" '
  BEGIN{
    # 读取子代顺序
    nk=0; while ((getline s < KIDS)>0){ if(s!=""){ kid[++nk]=s; is_kid[s]=1 } } close(KIDS)
    # 读取三元组映射：child -> father, mother
    while ((getline x < TRIOS)>0){
      split(x,a,"\t"); fam=a[1]; c=a[2]; f=a[3]; m=a[4];
      dad[c]=f; mom[c]=m;
    }
    cur=""
  }
  function flush(){
    if(cur=="") return
    # 输出四条：PT/MT/PU/MU
    line=CHR OFS POS OFS ID OFS REF OFS ALT OFS ".\t.\t.\tGT"
    # 组装四个行的样本 GT 列
    PTline=line; MTline=line; PUline=line; MUline=line
    for(i=1;i<=nk;i++){
      c=kid[i]
      gt_pt=(c in GT_PT)?GT_PT[c]:"./."
      gt_mt=(c in GT_MT)?GT_MT[c]:"./."
      gt_pu=(c in GT_PU)?GT_PU[c]:"./."
      gt_mu=(c in GT_MU)?GT_MU[c]:"./."
      PTline=PTline OFS gt_pt
      MTline=MTline OFS gt_mt
      PUline=PUline OFS gt_pu
      MUline=MUline OFS gt_mu
    }
    print PTline >> BODY_PT
    print MTline >> BODY_MT
    print PUline >> BODY_PU
    print MUline >> BODY_MU
    # 清空
    delete GT; delete GT_PT; delete GT_MT; delete GT_PU; delete GT_MU
  }
  # 简单拆 GT：返回拆分数（0=缺失），a[1],a[2] 为两等位
  function splitGT(s, a,   n){ 
    if(s=="." || s=="./." || s=="." ) return 0
    n = split(s, a, /[\/|]/); 
    if(n<2) return 0; 
    return n 
  }
  {
    CHR0=$1; POS0=$2; ID0=$3; REF=$4; ALT=$5
    key=CHR0"|"POS0"|"ID0"|"REF"|"ALT
    if(key!=cur){ flush(); cur=key; CHR=CHR0; POS=POS0; ID=ID0 }
    # 采样名->GT 映射
    delete GT
    i=6
    while(i<=NF){
      s=$(i); g=$(i+1); GT[s]=g; i+=2
    }
    # 对每个子代，计算四条链的 GT
    for(i=1;i<=nk;i++){
      c=kid[i]; f=dad[c]; m=mom[c]
      cGT=(c in GT)?GT[c]:"."      # 子代
      fGT=(f in GT)?GT[f]:"."      # 父亲
      mGT=(m in GT)?GT[m]:"."      # 母亲

      # 子代相位：左=PT，右=MT
      PT="."; MT="."
      if(index(cGT,"|")>0){ split(cGT,cc,/[\|]/); PT=cc[1]; MT=cc[2] }

      # 父未传 PU：在父亲等位里找“另一个”
      PU="."
      if(PT!="."){
        if(splitGT(fGT,ff)){
          if(ff[1]==ff[2]) PU=ff[1]
          else if(PT==ff[1]) PU=ff[2]
          else if(PT==ff[2]) PU=ff[1]
        }
      }
      # 母未传 MU
      MU="."
      if(MT!="."){
        if(splitGT(mGT,mm)){
          if(mm[1]==mm[2]) MU=mm[1]
          else if(MT==mm[1]) MU=mm[2]
          else if(MT==mm[2]) MU=mm[1]
        }
      }

      # 组装四个 VCF 的 GT 列
      GT_PT[c] = (PT!=".") ? PT"|." : "./."
      GT_MT[c] = (MT!=".") ? ".|"MT : "./."
      GT_PU[c] = (PU!=".") ? PU"|." : "./."
      GT_MU[c] = (MU!=".") ? ".|"MU : "./."
    }
  }
  END{ flush() }
'

# 7) 合并 header+body，压缩索引
for L in PT MT PU MU; do
  cat "$TMP/header.$L.vcf" "$TMP/body.$L.vcf" | bgzip -c > "${OUT_PREFIX}.${L}.vcf.gz"
  tabix -p vcf "${OUT_PREFIX}.${L}.vcf.gz"
  echo "[OK] ${OUT_PREFIX}.${L}.vcf.gz"
done
