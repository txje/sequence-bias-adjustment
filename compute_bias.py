import argparse
import pysam
import numpy
from string import maketrans

MAX_AT_POS = None #5
#MAX_AT_POS = 160
#print "------------- MAX_AT_POS: 160 (usually 5..) ---------------"
CHROMS = None
COMPL = maketrans("ACGT", "TGCA")

def make_kmers(k):
  if k == 0:
    return ['']
  return ['A' + mer for mer in make_kmers(k-1)] + ['C' + mer for mer in make_kmers(k-1)] + ['G' + mer for mer in make_kmers(k-1)] + ['T' + mer for mer in make_kmers(k-1)]

def main(bam_npy_file, fasta_file, chrom_file, k, output_file, read_limit, clip, read_len, margin):

  global CHROMS
  CHROMS = [c.split()[0] for c in open(chrom_file).read().strip().split('\n')]

  if not bam_npy_file[-4:] == ".npy":
    print "BAM file must be in .bam.npy format."
    return

  print "Loading BAM file..."
  bam = numpy.load(bam_npy_file)
  fa = pysam.Fastafile(fasta_file)
  kmer_counts = [{kmer:0 for kmer in make_kmers(k)} for i in xrange(margin * 2 + 1)]
  num_reads = 0
  print "Reading ref seq..."
  refs = [fa.fetch(chr).upper() for chr in CHROMS]

  # keep track of how many reads are at each pos
  last_pos = None
  num_at_pos = 0
  over_limit = 0

  print "Processing %i reads..." % len(bam)
  for read in bam:
    if read_limit != None and num_reads >= read_limit:
      break
    num_reads += 1
    if num_reads % 10**6 == 0:
      print num_reads

    # keep track of how many reads are at each pos
    if read[1] == last_pos:
      num_at_pos += 1
      if MAX_AT_POS is not None and num_at_pos > MAX_AT_POS:
        over_limit += 1
        continue
    else:
      num_at_pos = 1
      last_pos = read[1]

    # read sequence (may be reversed)
    if read[2]:
      # on reverse strand, read starts at read.pos + read_len - 1
      seq = refs[read[0]][read[1] - margin + read_len - 1:read[1] + margin + read_len].translate(COMPL)[::-1]
    else:
      seq = refs[read[0]][read[1] - margin:read[1] + margin + 1]
    # count k-mers
    weight = (1 if len(read) < 4 else read[3])

    # -----------------------------------------------
    # clip
    # -----------------------------------------------
    if clip != None:
      if weight > clip:
        weight = clip
      if weight < (1.0/clip):
        weight = (1.0/clip)

    if len(seq) < margin * 2 + 1:
      #print read, "seq is too short:", len(seq), "should be", (margin * 2 + 1)
      pass
    if k > 1:
      for i in xrange(min(margin * 2 - k + 2, len(seq)-k+1)):
        kmer = seq[i:i+k]
        if 'N' in kmer or not kmer_counts[i].has_key(kmer):
          continue
        kmer_counts[i][kmer] += weight
    else:
      for i in xrange(min(margin * 2 + 1, len(seq))):
        if seq[i] == 'N' or not kmer_counts[i].has_key(seq[i]):
          continue
        kmer_counts[i][seq[i]] += weight

  print "%i reads done." % num_reads
  if MAX_AT_POS is not None:
    print "%i positions with >%i reads" % (over_limit, MAX_AT_POS)

  kmers = list(set([key for kmers in kmer_counts for key in kmers.keys() if 'N' not in key and len(key) == k]))
  if len(kmers[0]) == 1:
    kmers = 'ACGT'
  total_reads_used = sum([float(kmer_counts[0][k] if kmer_counts[0].has_key(k) else 0) for k in kmers])
  fout = open(output_file, 'w')
  fout.write(',' + ','.join(kmers) + "\n")
  fout.write('\n'.join('%i,' % (a-margin) + ','.join(['%.8f' % (float(kmer_counts[a][k] if kmer_counts[a].has_key(k) else 0)/total_reads_used) for k in kmers]) for a in xrange(len(kmer_counts))))
  fout.close()

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Compute k-mer frequency bias")
  parser.add_argument("bam", help="BAM.npy file")
  parser.add_argument("ref", help="Fasta file")
  parser.add_argument("chroms", help="Chromosome file")
  parser.add_argument("k", help="k-mer size", type=int)
  parser.add_argument("out", help="Output (csv) file")
  parser.add_argument("--max", help="Maximum reads to process", type=int)
  parser.add_argument("--clip", help="Clip read weights to 1/N-N", type=float)
  parser.add_argument("--read_len", help="Read length", type=int, default=20)
  parser.add_argument("--margin", help="Width around read start to capture", type=int, default=40)
  args = parser.parse_args()

  main(args.bam, args.ref, args.chroms, args.k, args.out, args.max if args.max > 0 else None, args.clip if args.clip > 0 else None, args.read_len, args.margin)
