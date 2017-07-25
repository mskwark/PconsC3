import sys


def reformat(seqf, cf):
    seq = ""
    with open(seqf) as f:
        for l in f:
            if not l.startswith(">"):
                seq += l.strip()
            else:
                if seq:
                    #TODO: handle multichain here
                    break
    clist = []
    with open(cf) as f:
        for l in f:
            l = l.strip()
            c = l.split(' ')
            clist.append(c)

    for c in clist:
        resi = int(c[0])
        resj = int(c[1])
        resni = seq[resi-1]
        resnj = seq[resj-1]
        print "s1.%s%s s1.%s%s 0 8 %s" % (resni, resi, resnj, resj, c[2])


if __name__ == "__main__":
    seqf = sys.argv[1]
    cf = sys.argv[2]
    reformat(seqf, cf)
