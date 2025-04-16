import sys
import os

ifilename = sys.argv[1]
outputprefix = sys.argv[2]

nvertices = 0
nedges = 0
ngraphs = 0
vertexgraphmap = dict()
vertexlabelmap = dict()
edgeidmap = dict()
edgelabelmap = dict()
adjlists = dict()
edgeidtoedge = dict()
with open(ifilename, 'r') as cfile:
    for line in cfile:
        line = line.strip()
        ltype = line[0]
        if ltype == 't':
            vertexmap = dict()
            ngraphs += 1

        elif ltype == 'v':
            _, u, ulabel = line.split(" ")
            u = int(u)
            ulabel = int(ulabel)

            vertexmap[u] = nvertices
            vertexgraphmap[nvertices] = ngraphs - 1
            vertexlabelmap[nvertices] = ulabel
            adjlists[nvertices] = []

            nvertices += 1

        elif ltype == 'e':
            _, u, v, elabel = line.split(" ")
            u = vertexmap[int(u)]
            v = vertexmap[int(v)]
            elabel = int(elabel)
            edgeid = nedges
            nedges += 1

            adjlists[u].append((v,edgeid,elabel))
            adjlists[v].append((u,edgeid,elabel))
            

os.makedirs(outputprefix, exist_ok=True)

with open("%s/metadata" % outputprefix, 'w') as ofile:
    ofile.write("%d %d\n" % (nvertices, nedges))

print(nvertices, nedges)

fixededgeidmap = dict()
with open("%s/adjlists" % outputprefix, 'w') as ofile:
    nextedgeid = 0
    for u in range(nvertices):
        adjline = ""
        adjlist = sorted(adjlists[u])
        first = True
        for i in range(len(adjlist)):
            v,edgeid,elabel = adjlist[i]
            fixededgeid = fixededgeidmap.get(edgeid)
            if fixededgeid is None:
                fixededgeid = len(fixededgeidmap)
                fixededgeidmap[edgeid] = fixededgeid
            edgelabelmap[fixededgeid] = elabel
            if not first: ofile.write(" ")
            ofile.write("%d,%d" % (v,fixededgeid))
            first = False
        ofile.write("\n")

with open("%s/vlabels" % outputprefix, 'w') as ofile:
    for u in range(nvertices):
        ofile.write("%d\n" % vertexlabelmap[u])

with open("%s/elabels" % outputprefix, 'w') as ofile:
    for eid in range(nedges):
        ofile.write("%s\n" % edgelabelmap[eid])

with open("%s/vertexgraphmap" % outputprefix, 'w') as ofile:
    for u in range(nvertices):
        ofile.write("%d\n" % vertexgraphmap[u])