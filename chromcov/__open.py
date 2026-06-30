from dataclasses import dataclass
from pathlib import Path
import pysam

# url = "https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/COLO829T_TEST.cram"
# cramfile = pysam.AlignmentFile(url, "rc")

# @dataclass
# class CoverageMetric:

DATA_DIR = Path().resolve() / "data"

cram_file = str(DATA_DIR / "COLO829T_TEST.cram")
cram_ref = str(DATA_DIR / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
cram_index = str(DATA_DIR / "COLO829T_TEST.cram.crai")

# pysam.faidx(cram_ref)

with pysam.AlignmentFile(cram_file, "rc", reference_filename=cram_ref, index_filename=cram_index) as aln:
	# print("Loaded!")

	# for chrom, length in zip(aln.references, aln.lengths):
	# 	bases = sum(r.query_alignment_length for r in aln.fetch(chrom) if not r.is_unmapped)
	# # 	counts = aln.count_coverage(chrom)
	# # 	print(f"{chrom}: {counts}")
	# 	print(f"{chrom}: {bases / length if length else 0}")

	header = aln.header.to_dict()
	Path("data/header.json").write_text(json.dumps(header))
	

# cramfile.count_coverage("chr1")
# cramfile.references
# chrUn_KI270338v1: an unplaced genomic scaffold within the human genome assembly GRCh38/hg38
# chr22_KI270739v1_random: unplaced genomic scaffold that originates from chromosome 22, but whose exact location or orientation within the chromosome remains unknown.
# chrEBV

# cramfile = pysam.AlignmentFile("ex1.cram", "rc")

#   docker run --platform linux/amd64 \
    # -v "$(pwd)/data:/opt/mount" \
    # quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2 \
    # mosdepth -n --fast-mode -t 4 --by 1000 \
    # --fasta /opt/mount/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
    # /opt/mount/testmos \
    # /opt/mount/COLO829T_TEST.cram

# docker run --platform linux/amd64 \
#     -v "$(pwd)/data:/opt/mount" \
#     quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2 \
#     mosdepth -n -t 4 \
#     --fasta /opt/mount/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
#     /opt/mount/testmosslow \
#     /opt/mount/COLO829T_TEST.cram
