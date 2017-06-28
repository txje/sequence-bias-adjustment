import numpy
import sys
import pysam

def bam_count(bam, ref):
  r = 0
  for read in bam.fetch(ref):
    if not read.is_unmapped:
      r += 1
  return r

# usage: python bam2npy.py <bam>

bam = pysam.AlignmentFile(sys.argv[1], 'rb')
refs = [c.split()[0] for c in open(sys.argv[2]).read().strip().split('\n')]
print "Looking for", refs
print "in", bam.references
CHROMS = [c for c in refs if c in bam.references]

tagged = None
for r in bam.fetch(CHROMS[0]):
  tagged = r.has_tag("XW")
  break

numreads = [bam_count(bam, chr) for chr in CHROMS]
for i in xrange(len(CHROMS)):
  print "%s: %i reads" % (CHROMS[i], numreads[i])
numreads = sum(numreads)

if tagged:
  ary = numpy.zeros(numreads, dtype=('u1,u4,u1,f4'))
else:
  ary = numpy.zeros(numreads, dtype=('u1,u4,u1'))

flagged_unmapped = 0
i = 0
last = 0
for chr in CHROMS:
  chr_id = CHROMS.index(chr)
  for read in bam.fetch(chr):
    if read.is_unmapped:
      flagged_unmapped += 1
      # these are obviously given a reference_id and position since we got them in the fetch
      # these are included in pysam's bam.count(ref), but - according to spec - they are no good!
      continue
    if tagged:
      ary[i] = (chr_id, read.pos, read.is_reverse, read.get_tag("XW"))
    else:
      ary[i] = (chr_id, read.pos, read.is_reverse)
    i += 1
  print "%s: %i reads written to numpy file" % (chr, i-last)
  last = i

assert i == numreads, "Number of written reads doesn't match expected"

print "%i reads marked as unmapped, but given valid reference_id and position - we did NOT use them." % flagged_unmapped
numpy.save(sys.argv[3], ary)
