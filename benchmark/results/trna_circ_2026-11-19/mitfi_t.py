import sys, datetime as dt
ts = lambda l: dt.datetime.strptime(l[:23], "%Y-%m-%d %H:%M:%S,%f")
ml = mf = None
for l in open(sys.argv[1]):
    if "[core] tRNAs (MitoM-ML" in l: ml = ts(l)
    elif ml and "[core] tRNAs:" in l: mf = ts(l)
print(f"{(mf - ml).total_seconds():.1f}" if ml and mf else "0")
