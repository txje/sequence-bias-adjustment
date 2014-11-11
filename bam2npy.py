import numpy
import sys
import pysam

bam = pysam.Samfile(sys.argv[1], 'rb')
CHROMS = ['chr%i' % i for i in xrange(1,23)] + ['chrX', 'chrY']

reads = []
for chr in CHROMS:
  if not chr in bam.references:
    print "%s not found" % chr
    continue
  chr_id = CHROMS.index(chr)
  for read in bam.fetch(chr):
    reads.append((chr_id, read.pos, read.is_reverse))

ary = numpy.zeros(len(reads), dtype=('u1,u4,u1'))
for r in xrange(len(reads)):
  ary[r] = reads[r]

numpy.save(sys.argv[1], ary)
