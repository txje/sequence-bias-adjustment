import numpy
import argparse
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
import pysam

CHROMS = ['chr%i' % i for i in xrange(1,23)] + ['chrX']
READ_LEN = 20
MARGIN = 5
CHROM_LENS = {"chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276, "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022, "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895, "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753, "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520, "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566}


def main(bamfile, outfile, stdout, npyfile, plotfile, dostats, read_len, chromosome, query_start, query_end, clip):
  if not stdout:
    print "Loading reads..."
  bam = pysam.AlignmentFile(bamfile, "rb")

  if not stdout:
    pileups = []
  i = 0

  has_weight = False

  for chrom in [chromosome] if chromosome != None else CHROMS:
    if query_start == None or query_end == None:
      start = 1
      end = CHROM_LENS[chrom]
    positions = end - start + 1

    pileup = numpy.zeros(positions)

    c = CHROMS.index(chrom)
    if not stdout:
      print chrom

    if not stdout:
      print "Computing pileups..."
    read_count = 0
    for read in bam.fetch(chrom):
      if read_count == 0:
        try:
          w = read.opt("XW")
          has_weight = True
        except:
          has_weight = False
        format = "\n%.2f" if has_weight else "\n%i"
      if read.reference_end + MARGIN > start and read.reference_start - MARGIN < end:
        read_count += 1
        weight = (clip if read.opt("XW") > clip else (1.0/clip if read.opt("XW") < 1.0/clip else read.opt("XW"))) if (has_weight and clip is not None) else 1
        pos = read.reference_start - start + (read_len - 1 if read.is_reverse else 0)
        pileup[max(0, pos - MARGIN) : min(positions, pos + MARGIN + 1)] += weight
      elif read_count > 0: # if this is sorted appropriately, we can quit now
        break
      i += 1
    if not stdout:
      print "%i reads" % read_count

    if not stdout:
      pileups.append(pileup)
    #else:
    #  print "fixedStep chrom=%s start=%i step=1" % (chrom, start)
    #  for p in pileup:
    #    print format % p

  total_pileup = numpy.concatenate(pileups)

  if dostats:
    print "Range: %i - %i" % (total_pileup.min(), total_pileup.max())
    print "Median: %i" % numpy.median(total_pileup)
    print "Mean: %i" % numpy.mean(total_pileup)
    print "Mode: %i" % stats.mode(total_pileup)[0][0]

  if plotfile or npyfile:
    hist = numpy.histogram(total_pileup, bins=total_pileup.max())
    if npyfile:
      numpy.save(npyfile, hist[0])
    if plotfile:
      fig = pyplot.figure()
      ax = fig.add_subplot(111)
      a = hist[0]
      a = numpy.log(a)
      a[numpy.isneginf(a)] = -1
      ax.plot(a)
      ax.xlabel("Reads/pos")
      ax.ylabel("Occurrences (log)")
      pyplot.savefig(plotfile)

  if outfile == None:
    return

  if outfile[-3:] == "csv":
    fout = open(outfile, 'w')
    fout.write("position,reads")
    for pileup in pileups:
      for p in xrange(len(pileup)):
        fout.write(format % (p + start, pileup[p]))
    fout.close()
  elif outfile[-3:] == "wig":
    fout = open(outfile, 'w')
    for c in xrange(len(pileups)):
      pileup = pileups[c]
      chrom = CHROMS[c] if chromosome == None else chromosome
      if c != 0:
        fout.write("\n")
      fout.write("fixedStep chrom=%s start=%i step=1" % (chrom, start))
      for p in xrange(len(pileup)):
        fout.write(format % pileup[p])
        #fout.write(format % (pileup[p] if clip == None else (clip if pileup[p] > clip else (1.0/clip if pileup[p] < 1.0/clip else pileup[p]))))
    fout.close()
  else:
    print "Unknown output file extension '%s'" % outfile[-4:]

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Compute read pileup, perhaps weighted")
  parser.add_argument("bam", help="BAM file, perhaps with weights")
  parser.add_argument("--out", help="Output (CSV or WIG) file")
  parser.add_argument("--stdout", help="Output WIG format to standard output", action="store_true")
  parser.add_argument("--npyout", help="Numpy array (NPY) output file")
  parser.add_argument("--plot", help="Plot pileup distribution (PNG) file")
  parser.add_argument("--stats", help="Display pileup stats?", action="store_true")
  parser.add_argument("--rlen", help="Read length", type=int, default=20)
  parser.add_argument("--chrom", help="Chromosome")
  parser.add_argument("--start", help="Start position", type=int)
  parser.add_argument("--end", help="End position", type=int)
  parser.add_argument("--clip", help="Clip read weights to 1/N-N", type=float)
  args = parser.parse_args()

  main(args.bam, args.out, args.stdout, args.npyout, args.plot, args.stats, args.rlen, args.chrom, args.start, args.end, args.clip)
