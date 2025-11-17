import pysam
from collections import defaultdict

# load and index the CellRanger BAM file
cr_bam = pysam.AlignmentFile("cell_unaligned_fix.bam", mode="rb")
cr_idx = pysam.IndexedReads(cr_bam)
cr_idx.build()

# load and iterate through the PathSeq BAM file
pathseq_bam = pysam.AlignmentFile("pathseq_credible.bam", mode="rb")

output = []
for seg in pathseq_bam.fetch(until_eof=True):
    # returns an IteratorRowSelection object, which contains one or more AlignedSegment object
    cr_list = list(cr_idx.find(seg.query_name))
    # we assume that all records belonging to the same query name will have the same CB/UB tag
    # not all records will have the CB tag and the UB tag
    if cr_list[0].has_tag("CB") and cr_list[0].has_tag("UB"):
        CB = cr_list[0].get_tag(tag="CB")
        UB = cr_list[0].get_tag(tag="UB")
        # using set_tags removes all other tags - use set_tag instead
        seg.set_tag("CB", CB, "Z")
        seg.set_tag("UB", UB, "Z")
        # d[CB][UB].append(seg)
    # keep all PathSeq alignments
    output.append(seg)

# write all PathSeq alignments with or without tags
all_pathseq_bam = pysam.AlignmentFile("pathseq_credible_tag.bam", mode="wb", template=pathseq_bam)
for seg in output:
    all_pathseq_bam.write(seg)

all_pathseq_bam.close()
cr_bam.close()
pathseq_bam.close()
