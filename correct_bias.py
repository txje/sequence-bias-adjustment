# ----------------------------------------------------------------------
# v1.0
#
# Main method to correct nucleotide-specific bias in NGS
# For the full pipeline, see seqbias_pipe.sh
#
# Uses "adaptive" relative weighting
# - scales weight based on nearby weights
#
# Jeremy Wang
# ----------------------------------------------------------------------

import numpy
import string
import math
import argparse
import pysam
import time

BASELINE_MARGIN = 50
MARGIN = 40
MAX_AT_POS = 5
CHROMS = None
BIAS_THRESHOLD = 5 # fold change allowed in std of allele frequencies, relative to bias at -40
TILE_COVARIANCE_THRESHOLD = 0.15 # tiles with kmers more correlated than this will be averaged, lower will be compounded
NEIGHBOR_MARGIN = 10
COMPL = string.maketrans("ACGT", "TGCA")

def read_baseline(f):
  data = [line.strip().split(',') for line in open(f, 'r').read().strip().split('\n')]
  baseline = {}
  k = len(data[0][0])
  for i in xrange(len(data[0])):
    baseline[data[0][i]] = float(data[1][i])
  return baseline

def read_bias(f):
  data = [line.strip().split(',') for line in open(f, 'r').read().strip().split('\n')]
  k = len(data[0][1])
  bias = [{} for i in xrange(len(data)-1)]
  header = data[0][1:]
  for i in xrange(1, len(data)):
    for h in xrange(len(header)):
      bias[i-1][header[h]] = float(data[i][h+1])
  return bias

def read_fasta(fasta_file, k):
  fa = pysam.Fastafile(fasta_file)
  refs = [fa.fetch(c).upper() for c in CHROMS]
  return refs

def make_kmers(k):
  if k == 0:
    return ['']
  return ['A' + mer for mer in make_kmers(k-1)] + ['C' + mer for mer in make_kmers(k-1)] + ['G' + mer for mer in make_kmers(k-1)] + ['T' + mer for mer in make_kmers(k-1)]

def compute_groups(baseline, bias, k, cov_matrix_file):
  tiles = []
  print "Finding biased tiles."

  def allele_variance(bias):
    k = len(bias.keys()[0])
    alleles = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}
    for i in xrange(k):
      for kmer,v in bias.iteritems():
        alleles[kmer[i]][i] += v
    return sum([numpy.std(v) for v in alleles.values()]) / len(alleles.keys())

  def kmer_variance(bias):
    return numpy.std(bias.values())

  baseline_variance = allele_variance(bias[0])
  print "Threshold (std of allele frequencies):", (baseline_variance * BIAS_THRESHOLD)
  for i in xrange(0, len(bias) - k + 1):
    var = allele_variance(bias[i])
    print ("Tile %i allele freq std:" % i), var
    if var > baseline_variance * BIAS_THRESHOLD and (len(tiles) == 0 or i >= tiles[-1] + k):
      tiles.append(i)

  # compute covariance between tiles
  tile_covariance_matrix = numpy.load(cov_matrix_file)
  tile_same_matrix = [[0 for t in tiles] for u in tiles]
  for i in xrange(len(tiles)):
    t0 = tiles[i]
    for j in xrange(i+1, len(tiles)):
      t1 = tiles[j]
      print ("tiles %i and %i:" % (t0, t1)), tile_covariance_matrix[t0][t1], "same if >", TILE_COVARIANCE_THRESHOLD
      if tile_covariance_matrix[t0][t1] > TILE_COVARIANCE_THRESHOLD:
        tile_same_matrix[i][j] = True
        tile_same_matrix[j][i] = True

  print "Tile similarity matrix:"
  for row in tile_same_matrix:
    print row
  groups = []
  while len(tiles) > 0:
    my_group = [0]
    for i in xrange(1, len(tile_same_matrix[0])):
      if tile_same_matrix[0][i]:
        my_group.append(i)
    groups.append([tiles.pop(t) for t in my_group[::-1]])
    # remove in reverse order
    for i in my_group[::-1]:
      tile_same_matrix.pop(i)
      for j in xrange(len(tile_same_matrix)):
        tile_same_matrix[j].pop(i)

  print "Correction groups:", groups
  return groups

def get_weight(chr, pos, strand, k, length, ref, groups, bias, baseline):
  if strand:
    read_start = pos + length - 1
    seq = ref[chr][read_start - MARGIN:read_start + MARGIN + 1].translate(COMPL)[::-1]
  else:
    seq = ref[chr][pos - MARGIN:pos + MARGIN + 1]

  # just don't reweight reads with N (or R, or M) at all
  composition = {'A':0, 'C':0, 'G':0, 'T':0}
  bad = False
  for a in seq:
    if not composition.has_key(a):
      bad = True
      break
    composition[a] += 1
  if bad:
    return seq, -1

  if len(seq) < MARGIN * 2 + 1:
    return seq, -2

  weight = 1
  for group in groups:
    ratios = [baseline[seq[i:i+k]] / bias[i][seq[i:i+k]] for i in group]
    nonzero = [factor for factor in ratios if factor > 0] # sometimes baseline is zero
    if len(nonzero) > 0:
      weight *= sum(nonzero) / len(nonzero)

  return seq, weight

def main(bam_npy_file, fasta_file, chrom_file, baseline_file, bias_file, output_file, adjusted_file, cov_matrix_file, read_limit=None, read_len=20):
  global CHROMS
  CHROMS = [c.split()[0] for c in open(chrom_file).read().strip().split('\n')]
  baseline = read_baseline(baseline_file)
  bias = read_bias(bias_file)
  # autodetect k
  k = len(baseline.keys()[0])
  print "k: %i" % k

  if bam_npy_file[-4:] != ".npy":
    print "Input not bam.npy formatted"
    return

  read_weights = []

  eval_length = MARGIN * 2 + 1

  print "Reading FASTA seqs..."
  ref = read_fasta(fasta_file, k)

  print "Reading BAM..."
  bam = numpy.load(bam_npy_file)

  groups = compute_groups(baseline, bias, k, cov_matrix_file) # tile k-mers over regions where bias deviates significantly from baseline

  total_read_weight = 0

  '''
  # add every k-mer to frequency hashes
  frequencies = []
  for f in xrange(eval_length - k + 1):
    frequencies.append({kmer:0 for kmer in make_kmers(k)})
  '''

  read_weights = numpy.zeros(bam.size, dtype=('u1,u4,u1,f4'))
  r = 0
  num_reads = 0

  # keep track of how many reads are at each pos
  last_pos = None
  num_at_pos = 0

  ns = 0

  neighbor_cache = []
  t0 = time.time()
  for read in bam:
    '''
    chr_id = read[0]
    pos = read[1]
    reverse = read[2]
    '''

    if read_limit != None and num_reads >= read_limit:
      break
    num_reads += 1
    if num_reads % 10**6 == 0:
      t = time.time()
      print "%i reads done (%.2f%%) [%.2f reads/sec]" % (num_reads, float(num_reads)/bam.size*100, float(10**6)/(t-t0))
      t0 = t

    # keep track of how many reads are at each pos
    if read[1] == last_pos:
      num_at_pos += 1
      if num_at_pos > MAX_AT_POS:
        continue
    else:
      num_at_pos = 1
      last_pos = read[1]

    seq, weight = get_weight(read[0], read[1], read[2], k, read_len, ref, groups, bias, baseline)

    if weight == -1:
      ns += 1
      continue
    if weight == -2:
      print "Incomplete (truncated) sequence at read %i" % num_reads
      continue

    # neighborhood comparison
    tot_nearby_weight = 0
    nearby = 0
    # remove past neighbors from cache
    while len(neighbor_cache) > 0 and (neighbor_cache[0][0] != read[0] or neighbor_cache[0][1] < read[1] - NEIGHBOR_MARGIN):
      neighbor_cache.pop(0)
    for pos in xrange(read[1] - NEIGHBOR_MARGIN, read[1] + NEIGHBOR_MARGIN + 1):
      if pos == read[1]:
        continue
      for n in neighbor_cache:
        if n[1] == pos and n[2] == read[2]:
          w = n[3]
          break
      else:
        nearby_seq, w = get_weight(read[0], pos, read[2], k, read_len, ref, groups, bias, baseline)
        if w >= 0:
          neighbor_cache.append((read[0], pos, read[2], w))
      if w < 0:
        continue
      nearby += 1
      tot_nearby_weight += w
    if nearby > 0:
      avg_nearby_weight = tot_nearby_weight / nearby
      weight = weight / avg_nearby_weight
    # otherwise, we'll just use the original weight

    total_read_weight += weight

    '''
    # update_frequencies
    for i in xrange(len(frequencies)):
      frequencies[i][seq[i:i+k]] += weight
    '''

    read_weights[r] = (read[0], read[1], read[2], weight)
    r += 1

    '''
    # compute fractional kmer frequencies
    # (renormalize)
    for i in xrange(len(frequencies)):
      tot = sum(v for key,v in frequencies[i].iteritems())
      for key in frequencies[i].keys():
        frequencies[i][key] /= tot
    '''

  print "%i reads with Ns" % ns

  print "Average read weight: %.4f (weight has since been renormalized)" % (total_read_weight / r)
  print "A value which deviates significantly from 1.0 indicates multiple non-independent adjustments and is probably a cause for concern"

  print "%i reads done, of %i (the remaineder filtered)." % (r, num_reads)

  '''
  print "Writing allele frequencies..."
  print "Total %i-mer frequency: %.4f" % (k, sum(v for kmer,v in frequencies[0].iteritems()))
  fout = open(output_file, 'w')
  # need to reduce k-mer frequencies to nucleotide frequencies
  kmers = 'ACGT'
  nuc_freq = [{'A':0,'C':0,'G':0,'T':0,'N':0} for f in frequencies]
  for i in xrange(len(frequencies)):
    for kmer in frequencies[i].keys():
      nuc_freq[i][kmer[0]] += frequencies[i][kmer]
  fout.write(',' + ','.join(l for l in kmers) + "\n")
  # renormalize nucleotide frequencies for given weights (see /total_read_weight)
  fout.write('\n'.join('%i,' % (a-MARGIN) + ','.join(['%.4f' % (float(nuc_freq[a][kmer] if nuc_freq[a].has_key(kmer) else -1)) for kmer in kmers]) for a in xrange(len(nuc_freq))))
  fout.close()
  '''

  print "Writing read weights..."
  numpy.save(adjusted_file, read_weights)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Correct allele frequency bias")
  parser.add_argument("bam", help="BAM compressed .npy file")
  parser.add_argument("ref", help="Fasta file")
  parser.add_argument("chroms", help="Chromosome file")
  parser.add_argument("baseline", help="Baseline allele frequencies (CSV)")
  parser.add_argument("bias", help="Read allele frequency bias (CSV)")
  parser.add_argument("out", help="Output (CSV) allele frequencies")
  parser.add_argument("adjusted", help="Output read weights (per-read adjustment, binary)")
  parser.add_argument("covmatrix", help="Tile covariance matrix (NPY)")
  parser.add_argument("--max", help="Maximum reads to process", type=int)
  parser.add_argument("--read_len", help="Read length, default is 20bp", type=int)
  args = parser.parse_args()

  main(args.bam, args.ref, args.chroms, args.baseline, args.bias, args.out, args.adjusted, args.covmatrix, args.max, args.read_len)
