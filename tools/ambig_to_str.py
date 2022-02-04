import re
import argparse


def ambig2dic(ambig_f):
    """Read an ambig.tbl file and convert it to a dictionary."""
    ambig_regex = r"resid\s*(\d*)\s*and\s*segid\s*(\w)"
    ambig_dic = {}
    with open(ambig_f) as fh:
        for line in fh.readlines():
            matches = re.finditer(ambig_regex, line)
            for m in matches:
                resid = int(m.group(1))
                chain = m.group(2)
                if chain not in ambig_dic:
                    ambig_dic[chain] = []

                ambig_dic[chain].append(resid)
    return ambig_dic


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ambig_f")
    args = parser.parse_args()

    restraint_dic = ambig2dic(args.ambig_f)

    print(restraint_dic)
