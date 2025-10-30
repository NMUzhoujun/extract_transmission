#!/usr/bin/env bash
# Extract transmitted (PT/MT) and untransmitted (PU/MU) alleles per trio
# - 输入：多样本相位VCF/BCF（子代GT为"|"相位，左=父传递，右=母传递），PED文件(每行: Family Child Father Mother)
# - 输出：每个变异 × 家系一行：CHR POS ID Family Child Father Mother PT MT PU MU（默认为等位基因编号）
# - 可选：--emit-bases 同时输出碱基（PT_base, MT_base, PU_base, MU_base）
# 依赖：bcftools (建议1.14+), bgzip/tabix(可选), awk
set -euo pipefail

usage() {
  cat >&2 <<EOF
Usage:
  $0 -v input.vcf.gz -p trios.ped [-o out.tsv.gz] [--emit-bases] [-t 8] [-r chr1:1-10M | -R regions.bed]

Required:
  -v, --vcf           输入 VCF/BCF（建议已索引 .tbi/.csi）
  -p, --ped           PED/三元组列表：FamilyID  ChildID  FatherID  MotherID

Optional:
  -o, --out           输出文件（默认 transmissions.tsv.gz）
      --emit-bases    额外输出碱基列（PT_base/MT_base/PU_base/MU_base）
  -t, --threads N     解压/读取线程（默认 4）
  -r, --regions STR   仅分析该区域（如 chr1:1-10M）
  -R, --regions-file  仅分析这些区域（BED/位置列表）
Notes:
  1) 速度：建议先将 VCF 转为 BCF 并 --threads>1；也可用 -r/-R 分块并行。
  2) 缺失或未相位时，对应 PT/MT/PU/MU 输出为 "."。
  3) 多等位位点(>1 ALT)按编号映射：0->REF, 1->ALT1, 2->ALT2, ...
EOF
  exit 1
}

VCF=""
PED=""
OUT="transmissions.tsv.gz"
THREADS=4
EMIT_BASES=0
REG_ARG=()

# ---- 解析参数 ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--vcf) VCF="$2"; shift 2;;
    -p|--ped) PED="$2"; shift 2;;
    -o|--out) OUT="$2"; shift 2;;
    --emit-bases) EMIT_BASES=1; shift;;
    -t|--threads) THREADS="$2"; shift 2;;
    -r|--regions) REG_ARG+=("-r" "$2"); shift 2;;
    -R|--regions-file) REG_ARG+=("-R" "$2"); shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -z "$VCF" || -z "$PED" ]] && usage
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found" >&2; exit 1; }

# ---- 临时工作区 ----
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# 读取 PED，生成按 [子,父,母] 顺序的样本列表（保证与后续解析的三元组一致）
# PED 列假定为：Family  Child  Father  Mother（允许有表头；父母/子为0或NA则跳过）
awk -v OFS='\t' '
  BEGIN{FS="[ \t]"}
  NR==1{
    # 尝试识别表头；简单判断：若有非标准字符则也保留
    # 不做严格校验，直接兼容
  }
  {
    fam=$1; child=$2; dad=$3; mom=$4;
    # 跳过空/0/NA
    if (child==""||child=="0"||child=="NA") next
    if (dad  ==""||dad  =="0"||dad  =="NA") next
    if (mom  ==""||mom  =="0"||mom  =="NA") next
    print fam, child, dad, mom
  }
' "$PED" > "$TMP/trios.tsv"

# 样本列表（严格按照 子,父,母 的顺序为每个家系加入3行）
awk '{print $2; print $3; print $4}' "$TMP/trios.tsv" > "$TMP/samples.list"

# ---- 主流程：一次流式读取，输出每个变异×家系 的 PT/MT/PU/MU ----
# 通过 bcftools view (支持 --threads 与区域筛选) -> bcftools query 提取 GT
# 每行：CHROM POS ID REF ALT [GT_子1 GT_父1 GT_母1 GT_子2 GT_父2 GT_母2 ...]
# 然后 awk 解析三元组，计算：
#   PT = 子GT左等位（假定已相位，左=父传递）
#   MT = 子GT右等位（母传递）
#   PU = 父另一个等位（若父杂合，则为不是PT的那个；若父同型或缺失，处理为同型或"."）
#   MU = 母另一个等位（同理）
#
# emit_bases=1 时再映射成碱基（REF/ALT）。
FMT='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'

# 为了获得更高吞吐，先子集样本再 query
# 注意：如果 VCF 很大，建议先转换： bcftools view -Ou input.vcf.gz | bcftools convert -Ob -o input.bcf; bcftools index input.bcf
( bcftools view --threads "$THREADS" -S "$TMP/samples.list" "${REG_ARG[@]}" -Ou "$VCF" \
  | bcftools query -f "$FMT" - ) \
| awk -v OFS='\t' -v emit_bases="$EMIT_BASES" '
  BEGIN{
    # 读取三元组到数组：索引从1开始，和字段组对齐
    while ((getline line < "'$TMP/trios.tsv'") > 0) {
      split(line, a, "\t");
      fam[++n] = a[1]; child[n] = a[2]; dad[n] = a[3]; mom[n] = a[4];
    }
    close("'$TMP/trios.tsv'");
    header_printed=0;
  }
  {
    CHR=$1; POS=$2; ID=$3; REF=$4; ALT=$5;

    # ALT 多等位逗号分隔，构建映射 1->ALT1, 2->ALT2, ...
    split("", AL); # 清空
    if (ALT != ".") {
      k=split(ALT, AL, ",");
    } else { k=0; }

    # 打印表头一次
    if (!header_printed) {
      if (emit_bases==1) {
        print "CHROM","POS","ID","Family","Child","Father","Mother","PT","MT","PU","MU","PT_base","MT_base","PU_base","MU_base";
      } else {
        print "CHROM","POS","ID","Family","Child","Father","Mother","PT","MT","PU","MU";
      }
      header_printed=1;
    }

    # 第一个GT字段位置
    base_idx=6;

    # 为消费性能，尽量少用正则，简单检测字符
    for (g=1; g<=n; g++) {
      ci = base_idx + (g-1)*3;
      fi = ci + 1;
      mi = ci + 2;

      cGT = $(ci);
      fGT = $(fi);
      mGT = $(mi);

      # 解析子代GT，期望 x|y
      PT="."; MT=".";
      if (index(cGT, "|")>0) {
        split(cGT, cc, "|"); PT=cc[1]; MT=cc[2];
      } else if (index(cGT, "/")>0) {
        # 未相位则无法确定父/母传递，按缺失处理
        split(cGT, cc, "/"); PT="."; MT=".";
      } else if (cGT=="." || cGT=="./.") {
        PT="."; MT=".";
      }

      # 父母等位：用[/|]都可分割
      F1="."; F2=".";
      if (index(fGT,"/")>0 || index(fGT,"|")>0) { split(fGT, ff, "[/|]"); F1=ff[1]; F2=ff[2]; }

      M1="."; M2=".";
      if (index(mGT,"/")>0 || index(mGT,"|")>0) { split(mGT, mm, "[/|]"); M1=mm[1]; M2=mm[2]; }

      # 计算父未传递 PU
      PU=".";
      if (PT!="." && F1!="." && F2!=".") {
        if (F1==F2) { PU=F1; }
        else if (PT==F1) { PU=F2; }
        else if (PT==F2) { PU=F1; }
        else { PU="."; } # 子父不相容/缺失
      }

      # 计算母未传递 MU
      MU=".";
      if (MT!="." && M1!="." && M2!=".") {
        if (M1==M2) { MU=M1; }
        else if (MT==M1) { MU=M2; }
        else if (MT==M2) { MU=M1; }
        else { MU="."; }
      }

      if (emit_bases==1) {
        # 把编号映射为碱基
        PTb = (PT=="."? ".": (PT+0==0? REF: ((PT+0)<=k? AL[PT+0]: ".")));
        MTb = (MT=="."? ".": (MT+0==0? REF: ((MT+0)<=k? AL[MT+0]: ".")));
        PUb = (PU=="."? ".": (PU+0==0? REF: ((PU+0)<=k? AL[PU+0]: ".")));
        MUb = (MU=="."? ".": (MU+0==0? REF: ((MU+0)<=k? AL[MU+0]: ".")));
        print CHR, POS, ID, fam[g], child[g], dad[g], mom[g], PT, MT, PU, MU, PTb, MTb, PUb, MUb;
      } else {
        print CHR, POS, ID, fam[g], child[g], dad[g], mom[g], PT, MT, PU, MU;
      }
    }
  }
' | bgzip -c > "$OUT"

echo "Done: $OUT"
