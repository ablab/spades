#!/usr/bin/python

import string
import math
from matplotlib import use
use('Agg')

from matplotlib import pyplot

#PREFIX = "./"
PREFIX="/smallnas/dima/algorithmic-biology/assembler/data/debruijn/ECOLI_SC_LANE_1_BH/K55/11.15_20.38.26/saves/"
MAX_DIST = 10000000
#DOWN_CUT = 20
#UP_CUT = 1000

def parse_pos(filename, down_cut, up_cut) :

	contigs = []
	file = open(PREFIX + filename,'r')
	k = 0
	while (1)  :
		line = file.readline()
		line = line.strip()
		data = string.split(line,' ')
		size = len(data)

		if ( len(line) == 0 ) :
			break

		#print line
		if (size == 2) :
			
			ref_num = int(data[1])
			if ( ref_num > 0 ) :
				for i in range(ref_num) :
					line = file.readline()
					line = line.strip()
					data_seq = string.split(line,' ')
					#print line
					edge_id = int(data[0])
					start_pos = int(data_seq[1])
					end_pos = int(data_seq[3])
					if ( end_pos - start_pos > down_cut and end_pos - start_pos < up_cut ) :
						contigs.append( (edge_id, start_pos, end_pos) )
					else :
						k += 1

		
	file.close()
	print k, "reads were discarded..."
	return contigs

def parse_cvr(filename) :

	coverage = {}
	file = open(PREFIX + filename,'r')
	while (1)  :
		line = file.readline()
		line = line.strip()
		data = string.split(line,' ')

		size = len(data)

		if ( len(line) == 0 ) :
			break

		if (size == 3) :
			
			coverage[ int(data[0]) ] = float(data[1])
	file.close()
	return coverage
 
def parse_grp(filename):

	file = open(PREFIX + filename,'r')
	edges = {}
	redges = {}
	for line in file :
		line = line.strip()
		data = string.split(line,' ')

		if (data[0] == 'Edge') :
			edges[int(data[1])] = (int(data[3]),int(data[5][:-1]))
			redges[int(data[1])] = int(data[10])

	file.close()
	#print edges
	
	return edges, redges

	
def pred_fwd( start_pos, compared, min_dist) :
 	return ( start_pos > compared and start_pos - compared < min_dist ) 
 
def pred_bwd( start_pos, compared, min_dist) :
	return ( start_pos < compared and compared - start_pos < min_dist ) 

# returns id in list that is closest to id i by pos value
def find_closest( i, start_pos, list, pred, contigs, coverage ) :

	min_dist = MAX_DIST
	min_id = i

	for j in range( len(list) ) :
		if ( i == j ) :
			continue
		if pred(start_pos, list[j], min_dist) :
			min_dist = math.fabs(start_pos - list[j])
			min_id = j

	#print i, contigs[i][1], contigs[i][2], coverage[contigs[i][0]], min_id, contigs[min_id][1], contigs[min_id][2], coverage[contigs[min_id][0]] 

	return min_id

# splits graph into component discarding reads that are potential to be unique
def get_components( edges, redges ) :

	print "splitting into components..."
	#print edges
	#import itertools

	out_count = {}
	in_count = {}
	for e in edges.keys() :
		ed = edges[e][0]
		if ed in out_count.keys() :
			out_count[ed] += 1
		else :
			out_count[ed] = 1
		ed = edges[e][1]
		if ed in in_count.keys() :
			in_count[ed] += 1
		else :
			in_count[ed] = 1
	
	components = {}
	singles = {}
	for e in edges.keys() :
		v_out = edges[e][0]
		v_in = edges[e][1]

		if (not v_out in in_count or out_count[v_out] > 1 ) and (in_count[v_in] > 1 or not v_in in out_count ) :
			if not redges[e] in singles.keys() :
				singles[e] = edges[e]
		else :
			if not redges[e] in components.keys() :
				components[e] = edges[e]

	return components, singles

# dfs traversal of components that are considered to be in the same repeat 
def traverse_components( components, singles, cov, redges ) :

	visited = set()
	for e in components :
		if not e in visited:
			ins = []
			outs = []
			visit(e, visited, components, singles, ins, outs)

			skip = False

			if len(ins) != len(outs) :
				skip = True

			if not skip:
				ins_dict = {}
				outs_dict = {}
				for i in ins :
					ins_dict[i] = cov[i]
				for i in outs :
					outs_dict[i] = cov[i]

				print dict(sorted(ins_dict.iteritems(), key=lambda (k,v): (v,k)))
				print dict(sorted(outs_dict.iteritems(), key=lambda (k,v): (v,k)))
			 
				print

# part of dfs implementation
def visit( e, visited, components, singles, ins, outs ) :

	#print "singles: ", singles
	#print "components: ", components
	if not e in visited :
		visited.add(e)

		for single_e in singles.keys() :
			#print single_e, e
			if components[e][1] == singles[single_e][0] :
				outs.append(single_e)
			if components[e][0] == singles[single_e][1] :
				ins.append(single_e)

		for adj_e in components.keys() :
			if components[e][1] == components[adj_e][0] or components[e][0] == components[adj_e][1] :
				visit (adj_e, visited, components, singles, ins, outs)

	
# simple search for "fork" components in graph, which checks one-edge repeats
def split_by_coverage( contigs, edges, long ) :

	#print "EDGES: ", edges
	#print
	long_edges = set()
	all_edges = set()
	# discards edges < long
	for c in contigs :
		if ( math.fabs(c[2] - c[1]) > long ) :
			long_edges.add( c[0] )
			#print c[0]
		all_edges.add( c[0] )
	
	in_edges = {}
	out_edges = {}

	for e in edges :
		if e in long_edges :
			for flank_e in all_edges :
				if ( edges[flank_e][1] == edges[e][0] ) :
					if e in in_edges :
						in_edges[e].append(flank_e)
					else :
						in_edges[e] = [flank_e]

				if ( edges[flank_e][0] == edges[e][1] ) :
					if e in out_edges :
						out_edges[e].append(flank_e)
					else :
						out_edges[e] = [flank_e]
	print "in:"
	for el in in_edges :
		if len(in_edges[el]) > 1 :
			print el, in_edges[el]
	print "out:"
	for el in out_edges :
		if len(out_edges[el]) > 1 :
			print el, out_edges[el]
	

	print "sorted:"
	for e in long_edges:
		ins = []
		outs = []
		if e in in_edges :
			for in_e in in_edges[e] :
				ins.append((in_e,coverage[in_e]))
			
		else :
			continue

		if e in out_edges :
			for out_e in out_edges[e] :
				outs.append((out_e,coverage[out_e]))
		else :
			continue

		if len(ins) > 1 and len(outs) > 1 :
			print e
			#print ins, outs
			ins = sorted(ins, key=lambda x: x[1])
			outs = sorted(outs, key=lambda x: x[1])
			
			print 'ins:', ins
			print 'outs:', outs
x = []
y = []
def add_pairs_to_plot( list_a, list_b ):

	l = min( len(list_a), len(list_b) )
	for i in range(l) :
		x.append(list_a[i][1])
		y.append(list_b[i][1])

def print_pairs(list_a, list_b) :
	
	l = min( len(list_a), len(list_b) )

	for i in range(l) :
		print list_a[i], list_b[i]

	for i in range(l,len(list_a)) :
		print "unused ",list_a[i]

	for i in range(l, len(list_b)) :
		print "unused",list_b[i]

def get_rank( contigs, coverage ) :

	pos_fwd = [r[1] for r in contigs]
	pos_bwd = [r[2] for r in contigs]


	incoverage = []
	outcoverage = []
	for i in range( len(contigs) ) :
		start_pos = contigs[i][1]
		c = contigs[ find_closest( i, start_pos, pos_bwd, pred_fwd, contigs, coverage ) ][0]
		if c in coverage.keys() :
			incoverage.append( (c, coverage[c]) )
		else :
			incoverage.append( (c, 0) )
		
	for i in range( len(contigs) ) :
		start_pos = contigs[i][2]
		c = contigs[ find_closest( i, start_pos, pos_fwd, pred_bwd, contigs, coverage ) ][0]
		if c in coverage.keys() :
			outcoverage.append( (c,coverage[c]) )
		else :
			outcoverage.append( (c,0) )
	
	#print incoverage
	#print outcoverage
	
	change = []
	for i in range( len(incoverage) ):
		change.append( math.fabs(incoverage[i][1] - outcoverage[i][1]) )

	rank = [0] * len(incoverage)
	for i in range( len(rank) ) :
	
		current_outcoverage = outcoverage[i][1]
		c_out = outcoverage[i][0]
		c_in = incoverage[i][0]

		for j in range( len(rank) ) :
			if i == j or c_out == contigs[j][0] or c_in == contigs[j][0] :
				continue
			if math.fabs(current_outcoverage - incoverage[j][1]) < change[i] :
				rank[i] += 1

	return rank

def flush(rank) :
	f = open('/home/ksenia/coverage/rank_result.csv','w')
	for i in range(len(rank)-1):
		f.write(str(rank[i]) + ',')
	f.write(str(rank[i]) + '\n')
	f.close()

#contigs = [(12, 2, 3), (13, 4, 11), (14, 6, 7), (15, 7, 8), (16, 9, 10), (17, 12, 13) ]
#contigs = [(12,1,2),(13,3,4),(14,5,6),(15,7,8), (16, 9,12), (17,14,18), (18,13,15) ]
#coverage = { 12:5,13:4,14:8, 15:2, 16:3, 17:4 ,18:9}

#contigs = parse_pos("distance_filling.pos",20,1000)
#coverage = parse_cvr("distance_filling.cvr")
#rank = get_rank(contigs,coverage)

TEST = False 
if TEST :
	PREFIX = "./"
	contigs = parse_pos("test.pos",0,100000000)
	coverage = parse_cvr("test.cvr")
	edges, redges = parse_grp("test.grp")
	#split_by_coverage(contigs, edges, 0)
	components, singles = get_components( edges, redges )
	traverse_components( components, singles, coverage, redges )
else :
	contigs = parse_pos("distance_filling.pos",0,10000)
	coverage = parse_cvr("distance_filling.cvr")
	edges, redges = parse_grp("distance_filling.grp")
	#split_by_coverage(contigs, edges, 100)
	components, singles = get_components( edges, redges )
	traverse_components( components, singles, coverage, redges )

#	pyplot.xlim([0,500])	
#	pyplot.ylim([0,500])	
#	pyplot.plot(x,y,marker='.',linestyle="")
#	pyplot.savefig('plot.png')
#print rank
#flush(rank)
