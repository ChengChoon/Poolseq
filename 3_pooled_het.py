import os
import sys

window_size = 10000
step_size = 5000

def sliding_window(start, end, step, size):
    intervals = []
    for i in range(start, end, step):
        if (i+size>end): # ignore last window, if too small at the end
            break
        intervals.append([i, i+size])
    return intervals

#print(sliding_window(0, 1250, 100, 400))

chrs = {}
f = open("ref.bed", "rt")
r = f.readline()
while r:
    r = r.replace("\r", "").replace("\n", "").split("\t")
    chrs[r[0]] = int(r[2])
    r = f.readline()
f.close()

f = open("ES_HT_1_F_R_400_rc", "rt")
header = f.readline()
data = {}
r = f.readline()
print("reading SNP data")
while r:
    r = r.replace("\r", "").replace("\n", "").split("\t")
    if len(r)!=13:
        r = f.readline()
        continue
    chr = r[0]
    pos = int(r[1])
    col1 = r[-4].split("/")
    col2 = r[-3].split("/")
    col3 = r[-2].split("/")
    col4 = r[-1].split("/")
    a = int(col1[0])
    b = int(col2[0])
    c = int(col3[0])
    d = int(col4[0])
    data["%s_%s" % (chr, pos)] = (a,b,c,d)
    r = f.readline()
f.close()

header = ["chr", "start", "stop", "num_SNPs", "pooled_het"]
#print("\t".join(header))
f1 = open("ES_HT_1_F_R_400_rc_het_pooled_w10k_s5k", "wt")
f1.write("\t".join(header)+"\n")

for chr, chr_size in chrs.items():
    print("processing", chr)
    windows = sliding_window(0, chr_size, step_size, window_size)
    for (start, stop) in windows:
        num_snps = 0
        wa, wb, wc, wd = 0, 0, 0, 0 # window values
        for i in range(start, stop):
            wdata = data.get("%s_%s" % (chr, i), None)
            if wdata!=None:
                ta, tb, tc, td = wdata
                wa += ta
                wb += tb
                wc += tc
                wd += wd
                num_snps += 1
        wp = wa+wb
        wq = wc+wd
        try:
            pool_het = (2*wp*wq) / ((wp+wq)*(wp+wq))
        except:
            pool_het = 0
        #print("%s\t%s\t%s\t%s\t%s" % (chr, start, stop, num_snps, pool_het))
        f1.write("%s\t%s\t%s\t%s\t%s\n" % (chr, start, stop, num_snps, pool_het))

f1.close()
