import math
import numpy
import pysam
import argparse
import time

CHROMS = None

def main(chromfile, npyfile, bamfile, outbam, tag, mult, noy, chrom):
  global CHROMS
  CHROMS = [c.split()[0] for c in open(chromfile).read().strip().split('\n')]
  if noy and "chrY" in CHROMS:
    CHROMS.remove("chrY")

  t0 = time.time()

  print "Loading npy file."
  reads = numpy.load(npyfile)
  # format: (chromid, pos, reverse?, weight)
  t1 = time.time()
  print "%.2f seconds." % (t1 - t0)

  print "Sorting npy file."
  # sort reads by chrom, pos, rev
  # they should already be in this order...
  #reads.sort(key = lambda a: a[0]<<32 + a[1]<<1 + a[2])
  t2 = time.time()
  print "%.2f seconds." % (t2 - t1)

  # skip zeros at the beginning
  r = 0
  while reads[r][0] == 0 and reads[r][1] == 0:
    r += 1
  zeros = r

  print "Loading bam file."
  refbam = pysam.AlignmentFile(bamfile, "rb")
  if chrom is not None:
    chrom = "chr"+chrom
  numrefs = (1 if chrom else len(CHROMS))
  t2 = time.time()
  print "%.2f seconds." % (t2 - t1)
  # get usable header
  header = refbam.header
  # REORDER CHROMOSOMES TO MATCH NPY FILE (CANONICAL ORDER)
  header["SQ"].sort(key = lambda a: CHROMS.index(a["SN"]) if a["SN"] in CHROMS else 999999)
  print [s["SN"] for s in header["SQ"]]
  bam = pysam.Samfile(outbam, "wb", header=header)

  print "Writing new bam file."
  written = 0
  for c in xrange(len(CHROMS)):
    print "on to %s (%i), at read %i of %i" % (CHROMS[c], c, r, reads.shape[0])
    while r < len(reads) and reads[r][0] < c:
      #print "skipping read %i:%i:%i because chrom is %i" % (reads[r][0], reads[r][1], reads[r][2], c)
      r += 1
    if r >= len(reads):
      print "  no reads on %s (%i)" % (CHROMS[c], c)
      break
    for read in refbam.fetch(CHROMS[c]):
      if reads[r][0] > c:
        break
      read.reference_id = c # SET THE REFID TO THE CURRENT CHROM ID
      while r < len(reads) and reads[r][0] == c and (reads[r][1] < read.reference_start or (reads[r][2] == 0 and read.flag & 16 > 0)):
        r += 1
      if r < len(reads) and reads[r][0] == c and reads[r][1] == read.reference_start and (reads[r][2] * 16) == (read.flag & 16):
        name = read.query_name
        if tag:
          #read.tags.append(('XW', reads[r][3]))
          read.set_tag('XW', reads[r][3], 'f')
          bam.write(read)
          written += 1
        else: # mult must be defined
          for i in xrange(int(round(reads[r][3] * mult))):
            read.query_name = name + "_%i" % i
            bam.write(read)
            written += 1
        r += 1

  print "%i reads in npy file" % r
  print "%i reads after zero-filtering" % (r - zeros)
  print "%i bam reads written" % written

  bam.close()
  refbam.close()
  t3 = time.time()
  print "%.2f seconds." % (t3 - t2)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Simply duplicate BAM reads")
  parser.add_argument("chroms", help="Chromosome file")
  parser.add_argument("npy", help="Weighted numpy file")
  parser.add_argument("bam", help="Unweighted BAM file")
  parser.add_argument("out", help="BAM file to output")
  parser.add_argument("--noy", help="Do not use chrY", action="store_true", default=False)
  parser.add_argument("--tag", help="Use XW tag to indicate weight", action="store_true", default=False)
  parser.add_argument("--mult", help="Weight multiplier", type=int)
  parser.add_argument("--chrom", help="Convert a single chromosome")
  args = parser.parse_args()
  if args.tag and args.mult is not None:
    print "Cannot specify --tag and --mult, pick one or the other"
  else:
    main(args.chroms, args.npy, args.bam, args.out, args.tag, args.mult, args.noy, args.chrom)

# python npy2bam.py $outdir/$prefix.${k}mer_adjusted.read_weights.npy $bamfiltered $outdir/$prefix.adjusted.bam --noy --tag
