import argparse
import pandas as pd
import os

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def parse_arguments():
  parser = argparse.ArgumentParser()
  parser.add_argument("count_file", type=isfile, help="Count matrix file ")
  parser.add_argument("threshold", type=int, default=10,
                      help="Threshold of count (Default 10)")
  parser.add_argument("output_file", help="Output file", default="output.tsv")
  args = parser.parse_args()
  return args.count_file, args.threshold, args.output_file

count_file, threshold, output_file = parse_arguments()
count = pd.read_csv(count_file, sep="\t")
transposed_count = count.transpose()
isabovethres = count.sum(axis=0) >= threshold
transposed_filt_count = transposed_count[isabovethres]
transposed_filt_count.transpose().drop("Size",1).to_csv(output_file, sep="\t", index=False)