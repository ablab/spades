#!/usr/bin/python
#    Copyright 2008 Paul Harrison (pfh@logarithmic.net)

#    This file is part of Velvet.

#    Velvet is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    Velvet is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Velvet; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import sys

print 'graph G {'
print 'graph [size="12,12", overlap=true]'

length_step = 1000
min_len = 20

links = [ ]
parent = { }
def root(node):
    if node not in parent: 
        parent[node] = node
    if parent[node] == node:
        return node
    parent[node] = root(parent[node])
    return parent[node]

def merge(node1, node2):
    parent[root(node1)] = root(node2)

stats = { } #Per connected subgraph statistics node_name -> [total_length, total_cover, n_nodes, n_arcs]
stats_parent = { }
def stats_root(node):
    if node not in stats_parent: 
        stats[node] = [0,0.0,[],0]
        stats_parent[node] = node
    if stats_parent[node] == node:
        return node
    stats_parent[node] = stats_root(stats_parent[node])
    return stats_parent[node]

def stats_merge(node1, node2):
    node1 = stats_root(node1)
    node2 = stats_root(node2)
    if node1 != node2:
        stats[node1][0] += stats[node2][0]
	stats[node1][1] += stats[node2][1]
        stats[node1][2].extend(stats[node2][2])
	stats[node1][3] += stats[node2][3]
	del stats[node2]
	stats_parent[node2] = node1
    return node1

node_cover = { }
covers = [ ]
while True:
   line = sys.stdin.readline()
   if not line: break
   parts = line.rstrip().split()
   if not parts: continue
   if parts[0] == 'ARC':
       merge(parts[1], str(-int(parts[2])))
       sroot = stats_merge(parts[1].lstrip('-'),parts[2].lstrip('-'))
       stats[sroot][3] += 1
       
   if parts[0] == 'NODE':
       node_id = parts[1]
       seq = sys.stdin.readline().strip()
       cover = float(parts[2])/max(len(seq),1)
       covers.append(cover)
       
       sroot = stats_root(node_id)
       stats[sroot][0] += len(seq)
       stats[sroot][1] += float(parts[2])
       stats[sroot][2].append((len(seq), node_id))
       
       if len(seq) < min_len:
           merge('-'+node_id, node_id)
       else:
           nodes_to_link = ['-'+node_id]
           for i in xrange(0, max(1, len(seq)//length_step)):
	       string_name = node_id+'_string_%d'%i
	       node_cover[string_name] = cover
               nodes_to_link.append(string_name)
           nodes_to_link.append(node_id)
	   length = 100.0*float(len(seq)) / ((len(nodes_to_link)-1)*length_step)
           for i in xrange(len(nodes_to_link)-1):
               links.append((nodes_to_link[i], nodes_to_link[i+1], length, cover))

covers.sort()
median_cover = covers[len(covers)//2]
print >> sys.stderr, 'Minimum cover', covers[0]
print >> sys.stderr, 'Median  cover', median_cover
print >> sys.stderr, 'Maximum cover', covers[-1]

for a,b,length,cover in links:
    a = root(a)
    b = root(b)
    #if a != b:
    cval = min(255,max(0,int(cover/median_cover*128.0)))
    color = '#%02x%02x%02x' % (cval, 0, 255-cval)
    print '"%s"--"%s" [len=%f, color="%s", style="setlinewidth(2)"]' % (a,b,length,color)

for node in parent:
    if parent[node] != node: continue
    if node not in node_cover:
        color = '#000000'
    else:
        cval = min(255,max(0,int(node_cover[node]/median_cover*128.0)))
	color = '#%02x%02x%02x' % (cval, cval, 255-cval)
    #print '"%s" [style=filled color="%s" label=""]' % (node, color)
    print '"%s" [style=invis, fixedsize=true, width=0, height=0, label=""]' % (node)

print '}'


print >> sys.stderr, 'Connected subgraphs'
stats = stats.values()
stats.sort()
for length, cover, nodes, n_arcs in stats:
    print >> sys.stderr, ' * %d bases, avg cover %.1f, %d nodes, avg degree %.1f' % (length, cover/length, len(nodes), n_arcs*2.0/len(nodes))
    nodes.sort()
    nodes = [ node_id for node_len, node_id in nodes[::-1] ]
    print >> sys.stderr, '   nodes', ', '.join(nodes[:10]), (len(nodes)>10 and '...' or '')
