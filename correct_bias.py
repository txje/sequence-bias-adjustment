# ----------------------------------------------------------------------
# v0.1
#
# Main method to correct nucleotide-specific bias in NGS
# For the full pipeline, see seqbias_pipe.sh
#
# Jeremy Wang
# Last modified (finished): 2014-11-11
# ----------------------------------------------------------------------

import numpy
import string
import math
import argparse
import pysam

BASELINE_MARGIN = 50
MARGIN = 40
MAX_AT_POS = 5
COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
CHROMS = ['chr%i' % i for i in xrange(1,23)] + ['chrX', 'chrY']
BIAS_THRESHOLD = 5 # fold change allowed in std of allele frequencies, relative to bias at -40
TILE_COVARIANCE_THRESHOLD = 0.15 # tiles with kmers more correlated than this will be averaged, lower will be compounded

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

def get_seq(read, ref, read_len):
  if read[2]:
    read_start = read[1] + read_len - 1
    return ''.join([COMPLEMENT[a] for a in ref[read[0]][read_start - MARGIN:read_start + MARGIN + 1]][::-1])
  else:
    return ref[read[0]][read[1] - MARGIN:read[1] + MARGIN + 1]

def compute_weight(group, seq, k, baseline, bias):
  weight = 1
  for combination in group:
    total = sum([baseline[i if len(baseline) > 1 else 0][seq[i:i+k]] / bias[i][seq[i:i+k]] for i in combination])
    weight *= (total / len(combination))
  return weight

def update_frequencies(frequencies, seq, weight, k):
  for i in xrange(len(frequencies)):
    frequencies[i][seq[i:i+k]] += weight

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
  for i in xrange(0, len(bias) - k + 1, k):
    print ("Tile %i allele freq std:" % i), allele_variance(bias[i])
    if allele_variance(bias[i]) > baseline_variance * BIAS_THRESHOLD:
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
  return [groups]

def main(bam_npy_file, fasta_file, baseline_file, bias_file, output_file, adjusted_file, cov_matrix_file, read_limit=None, read_len=20):
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

  #groups = [[[20,25,30,61,66,71]], [[40],[45],[51],[56]], [[35]]] #v1
  '''
  groups = [ #v2
      [[MARGIN-20,MARGIN-15,MARGIN-12] + [MARGIN+read_len+i for i in xrange(1,12,5) if MARGIN+read_len+i <= MARGIN * 2-k]],
      [[MARGIN],[MARGIN+5,MARGIN+6],[MARGIN+read_len-9],[MARGIN+read_len-4]],
      [[MARGIN-6,MARGIN-5]]]
  groups = [
      [[MARGIN],[MARGIN+read_len-4]],
      [[MARGIN-6,MARGIN-5]]]
  '''

  groups = compute_groups(baseline, bias, k, cov_matrix_file) # tile k-mers over regions where bias deviates significantly from baseline
  prev_frequencies = None
  prev_read_weights = [] # keep a list of all past weights, so we can undo them

  for group in groups:
    total_read_weight = 0

    print "Computing bias in group", group

    # add every k-mer to frequency hashes
    frequencies = []
    for f in xrange(eval_length - k + 1):
      frequencies.append({kmer:0 for kmer in make_kmers(k)})

    read_weights = numpy.zeros(bam.size, dtype=('u1,u4,u1,f4'))
    r = 0
    num_reads = 0

    # keep track of how many reads are at each pos
    last_pos = None
    num_at_pos = 0

    ns = 0

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
        print num_reads

      # keep track of how many reads are at each pos
      if read[1] == last_pos:
        num_at_pos += 1
        if num_at_pos > MAX_AT_POS:
          continue
      else:
        num_at_pos = 1
        last_pos = read[1]

      seq = get_seq(read, ref, read_len)
      # I guess, don't reweight reads with Ns at all
      if 'N' in seq:
        ns += 1
        continue

      if len(seq) < MARGIN * 2 + 1:
        print "Incomplete (truncated) sequence at read %i" % num_reads
        continue

      weight = compute_weight(group, seq, k, ([baseline] if prev_frequencies == None else prev_frequencies), bias)
      weight *= (prev_read_weights[-1][r][3] if len(prev_read_weights) > 0 else 1)

      # if this is the last pass, undo the first weighting
      # "unbaseline"
      '''
      if group == groups[-1]:
        weight /= prev_read_weights[0][r][3]
      '''

      total_read_weight += weight

      update_frequencies(frequencies, seq, weight, k)

      read_weights[r] = (read[0], read[1], read[2], weight)
      r += 1

    '''
    # multiply current and previous frequencies
    if prev_frequencies == None:
      prev_frequencies = frequencies
    else:
      prev_frequencies = [{kmer:(frequencies[i][kmer] * prev_frequencies[i][kmer]) for kmer in frequencies[i].keys()} for i in xrange(len(frequencies))]
    '''

    # compute fractional kmer frequencies
    # (renormalize)
    for i in xrange(len(frequencies)):
      tot = sum(v for key,v in frequencies[i].iteritems())
      for key in frequencies[i].keys():
        frequencies[i][key] /= tot

    # carry frequencies and weights through
    prev_frequencies = frequencies
    prev_read_weights.append(read_weights)

    print "%i reads with Ns" % ns

  print "Average read weight: %.4f (weight has since been renormalized)" % (total_read_weight / num_reads)
  print "A value which deviates significantly from 1.0 indicates multiple non-independent adjustments and is probably a cause for concern"

  print "%i reads done." % num_reads

  print "Writing allele frequencies..."
  frequencies = prev_frequencies # prev_frequencies will be the product of all adjustments
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

  print "Writing read weights..."
  numpy.save(adjusted_file, read_weights)

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Correct allele frequency bias")
  parser.add_argument("bam", help="BAM compressed .npy file")
  parser.add_argument("ref", help="Fasta file")
  parser.add_argument("baseline", help="Baseline allele frequencies (CSV)")
  parser.add_argument("bias", help="Read allele frequency bias (CSV)")
  parser.add_argument("out", help="Output (CSV) allele frequencies")
  parser.add_argument("adjusted", help="Output read weights (per-read adjustment, binary)")
  parser.add_argument("covmatrix", help="Tile covariance matrix (NPY)")
  parser.add_argument("--max", help="Maximum reads to process", type=int)
  parser.add_argument("--read_len", help="Read length, default is 20bp", type=int)
  args = parser.parse_args()

  main(args.bam, args.ref, args.baseline, args.bias, args.out, args.adjusted, args.covmatrix, args.max, args.read_len)
