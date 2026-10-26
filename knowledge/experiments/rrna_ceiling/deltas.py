"""Distribution of rRNA end offsets (pred - truth, in transcript direction) per tool from features.tsv."""
import sys, re, collections, csv
path, tool = sys.argv[1], sys.argv[2]
def parse(loc):
    # e.g. 1..1528+  or join? keep simple: first..last and strand
    m = re.match(r'(\d+)\.\.(\d+)([+-])', loc or '')
    if not m: return None
    return int(m.group(1)), int(m.group(2)), m.group(3)
c5 = collections.Counter(); c3 = collections.Counter(); both = collections.Counter(); n = 0
for r in csv.DictReader(open(path), delimiter='\t'):
    if r['tool'] != tool or r['kind'] != 'rRNA': continue
    t, p = parse(r['truth']), parse(r['pred'])
    if not t or not p or t[2] != p[2]: continue
    n += 1
    if t[2] == '+':
        d5, d3 = p[0] - t[0], p[1] - t[1]
    else:
        d5, d3 = t[1] - p[1], t[0] - p[0]   # 5' end is the high coordinate on minus
    k = lambda d: d if abs(d) <= 5 else ('<-5' if d < 0 else '>5')
    c5[k(d5)] += 1; c3[k(d3)] += 1; both[(k(d5), k(d3))] += 1
print(tool, 'rRNAs compared', n)
order = ['<-5', -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, '>5']
print("5' offset:", '  '.join(f'{o}:{c5[o]}' for o in order))
print("3' offset:", '  '.join(f'{o}:{c3[o]}' for o in order))
print('top (5\',3\') pairs:', both.most_common(12))
