#!/usr/bin/env bash
set -euo pipefail

usage(){
  echo "Usage: $0 -i <weights.tsv[.gz]> -d <hg19_avsnp151.txt[.gz]> -v <vid_list> -m <vid_maf.tsv> -o <out_prefix>" >&2
  echo "  -i  输入表：第1列=MarkerID(rs)，第2列=Effect_Allele；可为 .gz" >&2
  echo "  -d  ANNOVAR avsnp库：Chr Start End Ref Alt rsID（hg19: avsnp151；hg38: avsnp155），可为 .gz" >&2
  echo "  -v  变异白名单：每行一个 VariantID（chr:pos:ref:alt），纯文本或 .gz" >&2
  echo "  -m  频率文件：两列 VariantID<TAB>freq（如 allvid_maf_parents.tsv）" >&2
  echo "  -o  输出前缀" >&2
  exit 1
}

INPUT=""; ANNODB=""; OUTPREFIX=""; VIDLIST=""; MAFFILE=""
while getopts "i:d:o:v:m:" opt; do
  case "$opt" in
    i) INPUT="$OPTARG" ;;
    d) ANNODB="$OPTARG" ;;
    o) OUTPREFIX="$OPTARG" ;;
    v) VIDLIST="$OPTARG" ;;
    m) MAFFILE="$OPTARG" ;;
    *) usage ;;
  esac
done
[[ -z "$INPUT" || -z "$ANNODB" || -z "$OUTPREFIX" || -z "$VIDLIST" || -z "$MAFFILE" ]] && usage
[[ ! -f "$INPUT"  ]] && { echo "INPUT not found: $INPUT"  >&2; exit 1; }
[[ ! -f "$ANNODB" ]] && { echo "ANNODB not found: $ANNODB" >&2; exit 1; }
[[ ! -f "$VIDLIST" ]] && { echo "VIDLIST not found: $VIDLIST" >&2; exit 1; }
[[ ! -f "$MAFFILE" ]] && { echo "MAFFILE not found: $MAFFILE" >&2; exit 1; }

FINAL="${OUTPREFIX}.final.tsv"
UNMATCHED="${OUTPREFIX}.unmatched.txt"
NOTINVID="${OUTPREFIX}.not_in_vid.txt"

# 尽量减少子进程：必要时仅做一次解压到临时文件
tmp_in="$(mktemp)"; tmp_db=""; tmp_vid=""; tmp_maf=""
trap 'rm -f "$tmp_in" ${tmp_db:-} ${tmp_vid:-} ${tmp_maf:-}' EXIT

case "$INPUT" in *.gz) gzip -cd "$INPUT" > "$tmp_in" ;; *) cp -f "$INPUT" "$tmp_in" ;; esac
if [[ "$ANNODB" =~ \.gz$ ]]; then tmp_db="$(mktemp)"; gzip -cd "$ANNODB" > "$tmp_db"; DB="$tmp_db"; else DB="$ANNODB"; fi
if [[ "$VIDLIST" =~ \.gz$ ]]; then tmp_vid="$(mktemp)"; gzip -cd "$VIDLIST" > "$tmp_vid"; VID="$tmp_vid"; else VID="$VIDLIST"; fi
if [[ "$MAFFILE" =~ \.gz$ ]]; then tmp_maf="$(mktemp)"; gzip -cd "$MAFFILE" > "$tmp_maf"; MAF="$tmp_maf"; else MAF="$MAFFILE"; fi

# 一次 AWK 完成：白名单 -> 输入 -> 频率 -> 数据库；最后一次性输出
awk -v unmatched="$UNMATCHED" -v notin="$NOTINVID" '
BEGIN{FS=OFS="\t"}
# 文件1：白名单 VariantID
ARGIND==1 {
  gsub(/^[ \t]+|[ \t]+$/, "", $0);
  if ($0!="") allow[$0]=1;
  next
}
# 文件2：输入表（你的权重表）
ARGIND==2 {
  if (FNR==1) { hdr=$0; next }
  order[++n]=$1; row[$1]=$0;           # rs -> 原行
  next
}
# 文件3：频率表 VariantID \t freq
ARGIND==3 {
  if ($1=="") next;
  maf[$1]=$2+0;                        # 频率（数值）
  next
}
# 文件4：avsnp 库 Chr Start End Ref Alt rsID
ARGIND==4 {
  if (FNR==1) { next }                 # avsnp 可能没有表头；有就跳过
  if ($0 ~ /^#/) next
  chr=$1; pos=$2; ref=toupper($4); alt=toupper($5); rs=$6
  if (rs in row) {
    split(row[rs], b, FS)              # 原表拆列
    ea=toupper(b[2])                   # Effect_Allele
    if (ea==ref || ea==alt) {
      vid = chr ":" pos ":" ref ":" alt
      if (vid in allow) {
        key = b[1] SUBSEP toupper(b[2])        # (MarkerID,Effect_Allele) 作为去重键
        f = (vid in maf ? maf[vid] : -1)       # 缺失频率记为 -1
        if (!(key in best_maf) || f > best_maf[key]) {
          best_maf[key]=f;
          best_row[key]=vid OFS row[rs];
          if (!(key in seen_key)) { seen_key[key]=1; key_order[++kN]=key; }
        }
      } else {
        print vid > notin
      }
    }
    found[rs]=1
  }
  next
}
END{
  # 输出最终表头
  print "VariantID", hdr;
  # 按首次出现键的顺序输出最佳行
  for (i=1;i<=kN;i++) { k=key_order[i]; if (k in best_row) print best_row[k]; }
  # 未匹配 rs 按输入顺序输出
  for (i=1;i<=n;i++) { r=order[i]; if (!(r in found)) print r > unmatched }
}' "$VID" "$tmp_in" "$MAF" "$DB"  > "$FINAL"

echo "[OK] FINAL:        $FINAL"
echo "[OK] UNMATCHED rs: $UNMATCHED"
echo "[OK] NOT-IN-VID:   $NOTINVID"
