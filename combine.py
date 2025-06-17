#!/usr/bin/env python3
import re
import sys
import argparse
from colorama import init, Fore, Style

def parse_file(path):
    """
    Parse a file of the form:
      treatment_xxx [Pathogenic]: Causal Effect: 0.0123 (…)
    Returns a dict: { treatment_name: (label, point_estimate) }
    """
    pattern = re.compile(
        r'(?P<name>\S+)\s+\[(?P<label>Pathogenic|Benign)\]:\s+'
        r'Causal Effect:\s+(?P<val>[-\d\.]+)'
    )
    effects = {}
    with open(path, 'r') as f:
        for line in f:
            m = pattern.search(line)
            if m:
                name  = m.group('name')
                label = m.group('label')
                val   = float(m.group('val'))
                effects[name] = (label, val)
    return effects

def main():
    parser = argparse.ArgumentParser(
        description='Merge causal effects from two files and colorize output'
    )
    parser.add_argument('file1', help='First input file')
    parser.add_argument('file2', help='Second input file')
    parser.add_argument(
        '-o','--output',
        default='merged_effects.txt',
        help='Where to write the merged (plain) results'
    )
    args = parser.parse_args()

    init(autoreset=True)  # initialize colorama

    # parse both files
    eff1 = parse_file(args.file1)
    eff2 = parse_file(args.file2)

    # find treatments in common
    common = set(eff1) & set(eff2)
    if not common:
        sys.stderr.write("No matching treatments found in both files.\n")
        sys.exit(1)

    # build merged list
    merged = []
    for name in common:
        label1, v1 = eff1[name]
        label2, v2 = eff2[name]
        if label1 != label2:
            sys.stderr.write(
                f"Warning: label mismatch for {name}: {label1} vs {label2}\n"
            )
        merged.append((name, label1, v1 + v2))

    # sort by summed effect desc
    merged.sort(key=lambda x: x[2], reverse=True)

    # write plain results to output file
    with open(args.output, 'w') as fout:
        for name, label, total in merged:
            fout.write(f"{name}: {total:.4f}\n")

    # print colored to console
    for name, label, total in merged:
        color = Fore.RED   if label == 'Pathogenic' else \
                Fore.GREEN if label == 'Benign'    else Style.RESET_ALL
        print(f"{color}{name}: {total:.4f}{Style.RESET_ALL}")

if __name__ == "__main__":
    main()
