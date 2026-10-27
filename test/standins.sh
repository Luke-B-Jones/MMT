#!/usr/bin/env bash
# standins.sh BINDIR PYTHON [IPDSUMMARY [MODEL]] - SMRT Link stand-ins for testing MMT where SMRT Link is not installed.
#   ccs-kinetics-bystrandify  python (pysam): one pseudo-read per strand, ip = fi / ri
#   pbmm2                     index = copy of the FASTA ; align = minimap2 map-hifi + identity filter, tags kept
#   pbindex                   no-op
#   ipdSummary                a kineticsTools ipdSummary (pip) with an explicit model, when given
#   Needs samtools + minimap2 on PATH. Never used on the HPC: there MMT_SMRT_BIN points at the real SMRT Link.
set -euo pipefail
B=$1; PY=$2; KT=${3:-}; MODEL=${4:-}
mkdir -p "$B"
cat > "$B/ccs-kinetics-bystrandify" << EOF
#!/usr/bin/env bash
[[ \${1:-} == --version ]] && { echo "ccs-kinetics-bystrandify stand-in"; exit 0; }
exec "$PY" - "\$1" "\$2" << 'X'
import sys, pysam
C = str.maketrans("ACGTN", "TGCAN")
with pysam.AlignmentFile(sys.argv[1], check_sq=False) as f, pysam.AlignmentFile(sys.argv[2], "wb", template=f) as o:
    for r in f:
        for st, tag in (("fwd", "fi"), ("rev", "ri")):
            if not r.has_tag(tag): continue
            a = pysam.AlignedSegment(o.header); a.query_name = r.query_name + "/" + st; a.flag = 4
            s = r.query_sequence; a.query_sequence = s if st == "fwd" else s.translate(C)[::-1]
            a.set_tag("ip", r.get_tag(tag)); a.set_tag("RG", r.get_tag("RG")); a.set_tag("zm", r.get_tag("zm")); o.write(a)
X
EOF
cat > "$B/pbmm2" << EOF
#!/usr/bin/env bash
[[ \${1:-} == --version ]] && { echo "pbmm2 stand-in (minimap2)"; exit 0; }
if [[ \$1 == index ]]; then cp "\${@: -2:1}" "\${@: -1}"; exit 0; fi
Y=95; a=("\$@"); for i in "\${!a[@]}"; do [[ \${a[\$i]} == -y ]] && Y=\${a[\$((i+1))]}; done
IN=\${a[-2]}; OUT=\${a[-1]}; FA=\${a[-3]}
samtools fastq -T ip,RG,zm "\$IN" 2>/dev/null | minimap2 -ax map-hifi --eqx -y -t 2 "\$FA" - 2>/dev/null | "$PY" -c "
import sys, pysam
f = pysam.AlignmentFile('-'); h = f.header.to_dict(); src = pysam.AlignmentFile('\$IN', check_sq=False).header.to_dict()
h['HD'] = src['HD']; h['RG'] = src['RG']
o = pysam.AlignmentFile('\$OUT', 'wb', header=h); y = float('\$Y') / 100
for r in f:
    if r.is_unmapped: continue
    c = {}
    for op, l in r.cigartuples: c[op] = c.get(op, 0) + l
    m = c.get(7, 0); d = m + c.get(8, 0) + c.get(1, 0) + c.get(2, 0)
    if d and m / d >= y: o.write(pysam.AlignedSegment.from_dict(r.to_dict(), o.header))
o.close()"
EOF
printf '#!/usr/bin/env bash\nexit 0\n' > "$B/pbindex"
if [[ -n $KT ]]; then
    printf '#!/usr/bin/env bash\n[[ -n ${MMT_TEST_FAIL_IPD:-} && $* == *"/$MMT_TEST_FAIL_IPD/all.bam"* ]] && { echo "injected failure" >&2; exit 3; }\nexec "%s" "$@" %s\n' \
        "$KT" "${MODEL:+--ipdModel $MODEL}" > "$B/ipdSummary"
fi
chmod +x "$B"/*
