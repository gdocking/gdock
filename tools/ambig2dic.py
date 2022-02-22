import re, sys


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
                if resid not in ambig_dic[chain]:
                    ambig_dic[chain].append(resid)
    return ambig_dic


print(ambig2dic(sys.argv[1]))
