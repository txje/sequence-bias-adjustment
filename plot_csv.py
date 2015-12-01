import math
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *

MAX_LABELS = 20

def main(csv_file, out_file, log_base=0, xmin=None, xmax=None, smooth=None, size=(800,600), stacked=False):
  data = [line.strip().split(',') for line in open(csv_file, 'r').read().strip().split('\n')]
  header = data[0][1:]
  labels = [data[i][0] for i in xrange(1, len(data))]

  data = [[float(a) for a in d[1:]] for d in data[1:]]

  # pivot data table
  data = [[data[i][j] for i in xrange(len(data))] for j in xrange(len(data[0]))]

  # smooth if necessary
  if smooth:
    smoothed_data = [[] for d in data]
    for d in xrange(len(data)):
      for i in xrange(len(data[d])):
        smoothed_data[d].append(sum(data[d][max(0, i-(smooth/2)) : min(len(data[d]), i+(smooth/2)+1)]) / float(smooth))
    data = smoothed_data

  # manually stack up data lines if necessary
  if stacked:
    for i in xrange(len(data[0])):
      for d in xrange(1, len(data)):
        data[d][i] += data[d-1][i]

  fig = figure(figsize=(s/100.0 for s in size), dpi=100)
  axis = fig.add_subplot(111)
  ticks = range(len(labels))

  # try to convert labels to integers and scatter/plot
  try:
    labels = [float(l) for l in labels]
    for col in xrange(len(header)):
      axis.plot(labels, [d if log_base == 0 else (math.log(d, log_base) if d > 0 else -1) for d in data[col]], label=header[col])
    axis.set_xlim([min(labels) if xmin == None else xmin, max(labels) if xmax == None else xmax])
  # if that doesn't work, treat them as uniformly spaced string labels
  except:
    for col in xrange(len(header)):
      axis.plot([d if log_base == 0 else (math.log(d, log_base) if d > 0 else -1) for d in data[col]], label=header[col])
    axis.set_xticks(ticks)
    axis.set_xticklabels(labels)
  handles, labels = axis.get_legend_handles_labels()
  lgd = axis.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.1,0.9))
  axis.grid('on')
  savefig(out_file, bbox_extra_artists=(lgd,), bbox_inches='tight')

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Plot any CSV file, probably...")
  parser.add_argument("csv", help="CSV file")
  parser.add_argument("--out", help="Output image file")
  parser.add_argument("--log", help="Log base for y-axis", type=int)
  parser.add_argument("--xmin", help="Minimum x value to plot", type=int)
  parser.add_argument("--xmax", help="Maximum x value to plot", type=int)
  parser.add_argument("--smooth", help="Smooth (rolling average) over window of size n", type=int)
  parser.add_argument("--size", help="Plot size as WWWxHHH")
  parser.add_argument("--stacked", help="Stacked line plot", action="store_true")
  args = parser.parse_args()

  if args.size != None:
    w,h = (int(a) for a in args.size.split('x'))
  else:
    w,h = (800, 600)

  if args.out == None:
    out = args.csv[:args.csv.rindex('.')] + ".png" # replace .csv with .png
  else:
    out = args.out

  main(args.csv, out, args.log if args.log != None else 0, args.xmin, args.xmax, args.smooth, (w,h), args.stacked)
