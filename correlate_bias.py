import pysam
import argparse
import matplotlib
matplotlib.use("Agg")
from matplotlib import pylab
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy
import math
import random
from string import maketrans

MARGIN = 40
MEMORY_LIMIT = 32000000000 # 32 Gb (on bigmem, hopefully)

CHROMS = None
COMPL = maketrans("ACGTacgt", "TGCATGCA")

def read_bias(f):
  data = [line.strip().split(',') for line in open(f).read().strip().split('\n')]
  bias = [{} for i in xrange(MARGIN*2+1)]
  header = data[0]
  for i in xrange(1, len(data)):
    for j in xrange(1, len(data[i])):
      bias[i-1][header[j]] = float(data[i][j])
  return len(header[1]), bias

def main(bam_npy_file, ref_file, chrom_file, bias_file, out_file, read_limit=None, tile=False, plot=False, clip=None):
  global CHROMS
  CHROMS = [c.split()[0] for c in open(chrom_file).read().strip().split('\n')]
  k, bias = read_bias(bias_file)
  bam = numpy.load(bam_npy_file)
  fa = pysam.Fastafile(ref_file)
  print "Reading FASTA file..."
  refs = [fa.fetch(ch) for ch in CHROMS]
  tile_size = (k if tile else 1)

  num_samples = int(math.ceil(float((MARGIN * 2 + 1) - k + 1) / tile_size))
  bam_reads = len(bam)
  read_count = min(MEMORY_LIMIT / 4 / num_samples, bam_reads)
  read_prob = float(read_count) / bam_reads
  stacks = numpy.zeros((num_samples, read_count), dtype="f4") # dtype == 'f4' I think

  data = []
  i = 0
  s = 0
  print "Processing reads..."
  for read in bam:
    if read[1] < MARGIN or read[1] >= len(refs[read[0]]) - MARGIN:
      continue
    seq = refs[read[0]][read[1] - MARGIN : read[1] + MARGIN + 1].upper()

    if read[2]:
      seq = seq.translate(COMPL)[::-1]

    composition = {'A':0, 'C':0, 'G':0, 'T':0}
    bad = False
    for a in seq:
      # These can be Ns or weird sequence like Kirk's (with Rs and Ms)
      if not composition.has_key(a):
        bad = True
        break
      composition[a] += 1
    if bad:
      continue

    i += 1
    if i % 1000000 == 0:
      print i
    if read_limit != None and i >= read_limit:
      break

    weight = ((clip if read[3] > clip else (1.0/clip if read[3] < 1.0/clip else read[3])) if clip is not None else read[3]) if len(read) > 3 else 1

    if s < read_count and (bam_reads == read_count or random.random() < read_prob):
      for j in xrange(num_samples):
        kmer = seq[j*tile_size:j*tile_size+k]
        a = bias[j*tile_size][kmer]
        stacks[j, s] = a * weight
      s += 1

  correlation_matrix = [[0 for j in xrange(num_samples)] for l in xrange(num_samples)]
  for j in xrange(num_samples):
    correlation_matrix[j][j] = 0

  print "Computing covariance matrices..."
  print "Taking a maximum of %i samples, to save memory" % num_samples
  print "stacks:", stacks.shape
  for j in xrange(num_samples):
    print j
    for l in xrange(j+1, num_samples):
      matrix = numpy.vstack((stacks[j, :s],stacks[l, :s]))
      cov_matrix = numpy.corrcoef(matrix) # gets a normalized covariance matrix
      cov = cov_matrix[0][1]
      # use only magnitude of correlation, not direction
      if cov < 0:
        cov = -1 * cov
      correlation_matrix[j][l] = cov
      correlation_matrix[l][j] = cov

  cov_matrix = numpy.array(correlation_matrix)
  numpy.save(out_file, cov_matrix)

  print cov_matrix.shape
  print range(0, cov_matrix.shape[0], 5)
  print range(0, cov_matrix.shape[1], 5)
  cov_matrix = cov_matrix[range(0,cov_matrix.shape[0],5)][:,range(0,cov_matrix.shape[1],5)] # fancy indexing to get nonoverlapping tiles
  print cov_matrix.shape

  if plot:
    plt.imshow(cov_matrix, interpolation="none")
    plt.jet()
    cb = plt.colorbar()
    cb.set_label("covariance")
    plt.xticks(range(cov_matrix.shape[1]), range(-8,8))
    plt.yticks(range(cov_matrix.shape[0]), range(-8,8))
    plt.savefig("%s.png" % (out_file[:out_file.rindex('.')]))


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Compute bias correlation across read-relative positions")
  parser.add_argument("bam", help="BAM.npy file")
  parser.add_argument("ref", help="Fasta file")
  parser.add_argument("chroms", help="Chromosome file")
  parser.add_argument("bias", help="Bias (CSV) file")
  parser.add_argument("out", help="Output (npy) file")
  parser.add_argument("--limit", help="Maximum reads to process", type=int)
  parser.add_argument("--tile", help="Tile based on k-mer size?", action="store_true")
  parser.add_argument("--plot", help="Plot covariance matrix?", action="store_true")
  parser.add_argument("--clip", help="Clip weights", type=int)
  args = parser.parse_args()

  main(args.bam, args.ref, args.chroms, args.bias, args.out, args.limit, args.tile, args.plot, args.clip)
