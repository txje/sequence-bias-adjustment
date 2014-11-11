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

MARGIN = 40
MEMORY_LIMIT = 32000000000 # 32 Gb (on bigmem, hopefully)

COMPLEMENT = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
CHROMS = ['chr%i' % i for i in xrange(1,23)] + ['chrX', 'chrY']

def read_bias(f):
  data = [line.strip().split(',') for line in open(f).read().strip().split('\n')]
  bias = [{} for i in xrange(MARGIN*2+1)]
  header = data[0]
  for i in xrange(1, len(data)):
    for j in xrange(1, len(data[i])):
      bias[i-1][header[j]] = float(data[i][j])
  return len(header[1]), bias

def main(bam_npy_file, ref_file, bias_file, out_file, read_limit=None, tile=False, plot=False):
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
      seq = ''.join(COMPLEMENT[a] for a in seq[::-1])

    if 'N' in seq: # kind of heavy-handed
      continue

    i += 1
    if i % 1000000 == 0:
      print i
    if read_limit != None and i >= read_limit:
      break

    weight = read[3] if len(read) > 3 else 1

    if s < read_count and (bam_reads == read_count or random.random() < read_prob):
      for j in xrange(num_samples):
        a = bias[j*tile_size][seq[j*tile_size:j*tile_size+k]]
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

  if plot:
    plt.imshow(cov_matrix, interpolation="none")
    plt.jet()
    cb = plt.colorbar()
    cb.set_label("covariance")
    plt.savefig("%s.png" % (out_file[:out_file.rindex('.')]))


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Compute bias correlation across read-relative positions")
  parser.add_argument("bam", help="BAM.npy file")
  parser.add_argument("ref", help="Fasta file")
  parser.add_argument("bias", help="Bias (CSV) file")
  parser.add_argument("out", help="Output (npy) file")
  parser.add_argument("--limit", help="Maximum reads to process", type=int)
  parser.add_argument("--tile", help="Tile based on k-mer size?", action="store_true")
  parser.add_argument("--plot", help="Plot covariance matrix?", action="store_true")
  args = parser.parse_args()

  main(args.bam, args.ref, args.bias, args.out, args.limit, args.tile, args.plot)
