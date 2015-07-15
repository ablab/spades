#!/usr/bin/python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import string
import math
from matplotlib import use
use('Agg')

from matplotlib import pyplot

#PREFIX = "./"
#PREFIX="/smallnas/dima/algorithmic-biology/assembler/data/debruijn/ECOLI_SC_LANE_1_BH/K55/11.15_20.38.26/saves/"
#PREFIX="/home/ksenia/ecoli_saves/"
#PREFIX="/home/ksenia/algorithmic-biology/assembler/data/debruijn/ECOLI_IS220_QUAKE_100K/K55/01.29_17.33.47/saves/"
#PREFIX="/home/ksenia/algorithmic-biology/assembler/data/debruijn/ECOLI_SC_LANE_1_BH/K55/02.19_14.55.27/saves/"
PREFIX="/home/ksenia/algorithmic-biology/assembler/data/debruijn/ECOLI_SC_LANE_1_BH/K55/02.25_19.11.47/saves/"
PREFIX="/home/ksenia/algorithmic-biology/assembler/data/debruijn/ECOLI_SC_LANE_1_BH/K55/02.27_11.06.47/saves/"
PREFIX="/home/ksenia/algorithmic-biology/assembler/data/debruijn/ECOLI_SC_LANE_1_BH/K55/02.27_19.57.16/saves/"
#PREFIX="/home/ksenia/algorithmic-biology/assembler/data/debruijn/TOY_DATASET/K55/01.30_18.17.48/saves/"
#PREFIX="/home/antipov/algorithmic-biology/assembler/data/debruijn/SAUREUS_SC_LANE_7_BH/K55/12.17_22.23.06/saves/"
MAX_DIST = 10000000
#DOWN_CUT = 20
#UP_CUT = 1000

def parse_pos(filename) :# down_cut, up_cut) :

    contigs = []
    infile = open(PREFIX + filename,'r')
    k = 0
    while (1)  :
        line = infile.readline()
        line = line.strip()
        data = string.split(line,' ')
    #	print data
        size = len(data)

        if ( len(line) == 0 ) :
            break

        #print line
        if (size == 2) :
            ref_num = int(data[1])
            if ( ref_num > 0 ) :
                #if data[0] == '10010966' :
                #	print "!!!!!"
                #	print ref_num
                for i in range(ref_num) :
                    line = infile.readline()
		    line = line.strip()
                    data_seq = string.split(line,' ')
                    #if data[0] == '10010966' :
                    #	print data_seq
                    edge_id = int(data[0])
                    start_pos = int(data_seq[1])
                    end_pos = int(data_seq[3])
                    #if ( end_pos - start_pos > down_cut and end_pos - start_pos < up_cut ) :
                    contigs.append( (edge_id, start_pos, end_pos) )
                    #else :
                    #	k += 1

#	for e in contigs :
#		if e[0] == 10010966 :
#			print "!!!!!!!"
    infile.close()
    #print k, "reads were discarded..."
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
            edges[int(data[1])] = (int(data[3]),int(data[5][:-1]), int(data[8]))
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

#def filter_components( components, upper_bound ) :



# splits graph into component discarding reads that are potential to be unique
def get_components( edges, redges ) :
    upper_bound = 4000
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
	#if e == 17014031 :
	#	print "--------------------------",out_count[v_out], in_count[v_in], v_in in out_count
        if ((not v_out in in_count or out_count[v_out] > 1 ) and (in_count[v_in] > 1 or not v_in in out_count)):# or edges[e][2] > upper_bound:
	#	if e == 17014031 :
	#		print "singles"
            	singles[e] = edges[e]
        else :
            components[e] = edges[e]

    
    print
    
    return components, singles

def path_length( path ) :

    return sum(map( lambda x : x[1][2], path.iteritems()))

# dfs traversal of components that are considered to be in the same repeat
def traverse_components( components, singles, in_cov, out_cov, edges, redges, coverage ) :
    resolved_paths_begin = {}
    resolved_paths_end = {}
    visited = set()
    success = 0
    overall = 0
    nothing_close_number = 0
    too_short_number = 0
    cannot_resolve = 0
    long_resolved = 0
    cannot_resolve_long = 0
    upper_bound = 10000
    comp_file = open(PREFIX+"component.log", "w");
    resolved_paths = []
    failed_paths = {}
    all_paths = []
    for e in components :
        path = set()
        if not e in visited:
            ins = []
            outs = []
	    
            visit(e, visited, components, path, singles, ins, outs)
	   #if len(ins) != len(outs) :
	   # 	continue
	
	    #for e in ins :
	    #	if redges[e] in resolved_edges :
	#		continue
	 #   for e in outs :
	  #  	if redges[e] in resolved_edges :
	#		continue
#
#	    for e in path :
#	    	if redges[e] in resolved_edges :
#			continue

	    all_paths.append(list(path))
 
            for r in path:
                comp_file.write(str(r) + ' ')
            comp_file.write('\nIncoming:\n')
            for r in ins:
                comp_file.write(str(r) + ' ')
            comp_file.write('\nOutgoing:\n')
            for r in outs:
                comp_file.write(str(r) + ' ')

            comp_file.write('\n\n')
            ins_dict = {}
            outs_dict = {}
	    #print ins
            for i in ins :
                ins_dict[i] = in_cov[i]
            for i in outs :
                outs_dict[i] = out_cov[i]

            # overall counts components that do not have regions unaligned to reference
            has_empty_alignment = False
	    
            for edge in path:


                if len(get_ref_coord(edge, contigs)) == 0 :
                    has_empty_alignment = True
                    break
	    
            # do not traverse components with hanging vertices
            if len(ins) == 0 or len(outs) == 0 :
                continue
	    
            # filter loops
            if len( set(ins) & set(outs) ) > 0 :
                continue 

            #l = path_length( dict(filter ( lambda x : x[0] in path , edges.iteritems() )) )
            #print "max path length is ", len(path), '(', l, 'nt ) - ', path

	
	    ins_list = sorted(ins_dict.iteritems(), key=lambda (k,v): (-v,k))
            outs_list = sorted(outs_dict.iteritems(), key=lambda (k,v): (-v,k))

                #ins_list, outs_list =

            ps = find_closest( ins_list, outs_list)

            l = path_length( dict(filter ( lambda x : x[0] in path , edges.iteritems() )) )

            if len(ps) == 0 and not has_empty_alignment:
                    cannot_resolve += 1
                    #print "can not resolve: max path length is ", len(path), '(', l, 'nt ) - '
                    if l > 300 :
                        cannot_resolve_long += 1
                        print l, ':', path

            for e_pair in ps :
			
		    skip = False
                    coord_list_front = get_ref_coord(e_pair[0][0], contigs)
                    coord_list_back = get_ref_coord(e_pair[1][0], contigs)
		    overall += 1
		    resolved_path = resolve( e_pair, path, edges )
		    #print "resolved_path:", resolved_path	
		    print "pair:", e_pair[0], e_pair[1]
                    #if len(coord_list_front) >  0 and len(coord_list_back) > 0:
		    if True:
			#resolved_path = resolve( e_pair, path, edges )
			#print "resolved_path:", resolved_path	
			#for e in resolved_path:
			#	if redges[e] in resolved_edges :
			#		skip = True

			if not skip :
				#for e in resolved_path :
				resolved_paths.append(resolved_path)

				#print resolved_edges

				if len(resolved_path) > 2 :
                            		resolved_paths_begin[resolved_path[0]] = resolved_path
                            		resolved_paths_end[resolved_path[-1]] = resolved_path
					print "resolved_path:",resolved_path

                            		if not has_empty_alignment :
                                		ifAligned, nothing_close_number = compare_to_ref(resolved_path, contigs, ps, nothing_close_number)
                                		if  ifAligned:
							print "aligned:", resolved_path
                                    			success += 1
                                    			#print "success: max path length is ", len(path), '(', l, 'nt ) - '
                                    			if l > 300 :
                                        			long_resolved += 1
                                    			#	print ">300 -", l, path

                                		else :
							failed_paths[resolved_path[0]] = resolved_path
                                    			print "fail:  path length is ", len(resolved_path), '(', l, 'nt ) - '
				    			print "failed_path: " + str(resolved_path),
							for e in resolved_path :
								print e, ':', in_cov[e], '-', out_cov[e], ';',
							print

    #unite_all_paths(resolved_paths_begin, resolved_paths_end);
    #unite_all_paths(failed_paths, resolved_paths_end);
    
    #print len(failed_paths)
    #resolved_paths_begin = filter_reverse_complementary( resolved_paths_begin, redges ) 
    #resolved_paths_begin = filter_reverse_complementary( failed_paths, redges ) 
    #print len(failed_paths)
    #output_united_paths(resolved_paths_begin, "coverage_based.pth", edges,redges )
    
    #commented if TEST
    #output_united_paths(failed_paths, "coverage_based_failed.pth", edges,redges )
    unite_all_paths(resolved_paths_begin, resolved_paths_end);
    
    #print len(failed_paths)
    #resolved_paths_begin = filter_reverse_complementary( resolved_paths_begin, redges ) 
    #resolved_paths_begin = filter_reverse_complementary( failed_paths, redges ) 
    #print len(failed_paths)
    output_united_paths(resolved_paths_begin, "coverage_based.pth", edges,redges )
    #output_united_paths(failed_paths, "coverage_based.pth", edges,redges )

    print 'overall (not counting the paths we can not resolve):',overall

    print "all_paths: ", all_paths

    print "can not resolve", cannot_resolve
    print "can not resolve > 300", cannot_resolve_long
    print 'successfully aligned: ', success, ' of ', overall, ' resolved paths that do not have empty alignments inside'
    print 'long correctly resolved paths: ', long_resolved

    print 'resolved repeats, with edges in the path that can not be joined in reference: ', nothing_close_number
    #print 'in edge or out edge is(are) not covered by genome alignment: ', too_short_number, '(', float(too_short_number)/(overall-success), ')', ' of ', overall - success, ' unaligned repeats'
    #print 'other ', overall - success - too_short_number - nothing_close_number, ' do not have an edge in the reference alignment that would confirm in - out single edges'

def filter_reverse_complementary( resolved_paths, redges ) :

	new_resolved_paths = []
	#print resolved_paths
	#set_of_edges = set([item for sublist in resolved_paths.values() for item in sublist])
	for begining in resolved_paths.keys() :
		
		path = resolved_paths[begining]
		
		for begin_r in resolved_paths.keys() :
			if_reverse_compl = True
			path_r = resolved_paths[begin_r]
			if len(path_r) != len(path) : 
				print "different lengths"
				continue	
			
			print path
			print path_r
			path_length = len(path)
			for i in range(path_length) :
		
				print redges[path_r[path_length - 1 - i]]
				if path[i] != redges[path_r[path_length - 1 - i]] :
					if_reverse_compl = False
					break
			
			if if_reverse_compl :
				print "unique!"
				resolved_paths.pop(begining,0)
	return resolved_paths
			

def unite_all_paths( resolved_paths_begin, resolved_paths_end) :
    changed = True
    while changed:
        changed = False
        new_path = []
        for start in resolved_paths_begin:
            if start in resolved_paths_end:
                changed = True
                new_path = resolved_paths_end[start]
		new_path.extend(resolved_paths_begin[start][1:])
		#print resolved_paths_end[start]
		#print resolved_paths_begin[start]
        #       new_path = resolved_paths_end[start]
		#new_path.extend(resolved_paths_begin[start][1:])
		#print "PATH " + str(new_path)
                del resolved_paths_end[start]
                del resolved_paths_begin[start]
                break
        if changed:
            resolved_paths_begin[new_path[0]] = new_path
            resolved_paths_end[new_path[-1]] = new_path
        if changed == False:
            break;

def output_united_paths(resolved_paths_begin, outfilename, edges, redges):

    print resolved_paths_begin
    seq_file = open(PREFIX+ "distance_estimation.sqn", "r")
    sequences = {}
    while (1):
        cont_name = seq_file.readline()
        if len(cont_name) == 0:
            break
        cont = seq_file.readline().strip()
        sequences[int(cont_name.strip()[1:])] = cont
    seq_path_file =  open(PREFIX+ "distance_estimation_paths.sqn", "w")
    outFile = open (outfilename, "w")
    outFile.write(str(len(resolved_paths_begin)) + '\n')
    used = {}
    path_count = 0;
    for first in resolved_paths_begin:
        path = resolved_paths_begin[first]
        if redges[path[0]] in used:
            continue
        used[path[0]] = 1
        used[path[-1]] = 1
    #	print path
    	pl = 0
        start = True
        seq = ""
        path_count +=1
        for e in path:
            used[e] = 1
            if start:
                seq_path_file.write(">" + str(e) +"_path\n")
                seq += sequences[e]
            else:
                seq += sequences[e][55:]


            start = False

            pl += edges[e][2]
            outFile.write(str(len (path)) + " " + str(pl))
            outFile.write('\n')
            for j in range (0, len(path)):
                outFile.write(str(path[j]) + ' ')
            outFile.write('\n')
            for j in range (0, len(path)):
                outFile.write(str(edges[path[j]][2]) + ' ')
            outFile.write('\n')
            outFile.write('\n')
        seq_path_file.write(seq + '\n')
    for e in edges:
        if not e in used and not redges[e] in used:
            used[redges[e]] = 1
            used[e] = 1
            seq_path_file.write(">" + str(e) +"\n")
            seq_path_file.write( sequences[e])
            seq_path_file.write("\n")
    print str(path_count) + " is path count"

def dfs_comp( v_start, v_end, component, edges, path ) :

	#if v_start == v_end :
	#	return [8]
	for e in component:
		if v_start == edges[e][0] :
			if v_end == edges[e][1] :
				return [e]
			path = dfs_comp(edges[e][1], v_end, component, edges, path)
			if len(path) > 0:
				print path
				return [e] + path
	return []	

			

def resolve( in_out, component, edges ) :

    path = []

    edge_in = in_out[0]
    edge_out = in_out[1]

    #path.append(edge_in[0])

    v_start_repeat = edges[edge_in[0]][1]
    #print edges[edge_in[0]][0],  v_start_repeat,
    v_end_repeat = edges[edge_out[0]][0]

    path = dfs_comp(v_start_repeat, v_end_repeat, component, edges, path)
    print path
    '''print v_start_repeat, v_end_repeat
    cur_v_from = v_start_repeat
    #visited = set()
    #print component, in_out
    #ordered_component = sorted( list(component), key=lambda x: edges[x][0])
    #print component, edges[9372386], edges[9962223]
    #print ordered_component
    #print "component:", component
    #for e in component :
    #	print e, edges[e][0], edges[e][1]
    print component
    print in_out
    while cur_v_from != v_end_repeat :
    	#print cur_v_from, v_end_repeat
    	for e in component :
        #if not e in visited :
            #print 'e',e
	    #print cur_v_from
            if cur_v_from == edges[e][0] :
                path.append(e)
                cur_v_from = edges[e][1]
            if cur_v_from == v_end_repeat :
         	path.append(edge_out[0])

                break
            #print e,'-',path

    #if len(visited) == len(component) :
    #	path.append(edge_out[0])

    #print edges[edge_out[0]][1]

    #if path[-1] == path[-2] :
    #	print "path:", path'''
    return [edge_in[0]] + path + [edge_out[0]]

def get_closest( el, list_b ) :

    from math import fabs

    e_min = (0,10000000000)
    for e_b in list_b :

        if fabs(e_b[1] - el[1]) < e_min[1] :
            e_min = e_b
    return e_min

def find_closest( list_a, list_b ) :

    pairs = []

    length_bound = min(len(list_a),len(list_b))

    if len(list_a) != len(list_b) :
    	threshold = 0.93#0.97
    	threshold_diff = 0.5#0.64
    else :
    	threshold = 0.93
	threshold_diff = 0.5
    # check if resolution is possible
     	threshold = 0.97
    	threshold_diff = 0.64
    for i in range( length_bound - 1 ) :

        # minimize (dont need max min actually since sorted)
        diff_a = float( min(list_a[i][1], list_a[i + 1][1] )) / max( list_a[i][1], list_a[i + 1][1])
        diff_b = float( min(list_b[i][1], list_b[i + 1][1] )) / max( list_b[i][1], list_b[i + 1][1] )

        # maximize
        diff_ab = float ( min(list_a[i][1], list_b[i][1]) ) / max(list_a[i][1], list_b[i][1])
        #diff_ab_next = float ( min(list_a[i+1][1], list_b[i+1][1]) ) / max(list_a[i+1][1], list_b[i+1][1])

	#print diff_a, diff_b, diff_ab 

        if diff_a > threshold or diff_b > threshold or diff_ab < threshold_diff :#or diff_ab_next < threshold_diff :
            return []
        #if diff_a < threshold or diff_a > 1.0/threshold or diff_b < threshold or diff_b > 1.0/threshold or diff_a/diff_b > threshold_diff or diff_a/diff_b < 1.0/threshold_diff:
        #	return []

    i = length_bound - 1
 
    if i > -1 and float ( min(list_a[i][1], list_b[i][1]) ) / max(list_a[i][1], list_b[i][1]) < threshold_diff :
    	return []
    
    list_remained = []
    if length_bound < len(list_a) :
        list_remained = list_a

    elif length_bound < len(list_b) :
        list_remained = list_b

    if length_bound > 0 :
    	for i in range(length_bound-1, len(list_remained)-1) :
        	diff_rem = float( min(list_remained[i][1], list_remained[i + 1][1] )) / max( list_remained[i][1], list_remained[i + 1][1])
        	if diff_rem > threshold :
            		return []
    for i in range( length_bound ) :
        pairs.append((list_a[i],list_b[i]))

    # check if resolution is possible
    #for i in range(len(pairs)) :
        #if pairs[i][0]
    return pairs


'''	in_out = True
    if len(list_a) <= len(list_b) :
        l_a = list_a
        l_b = list_b
    else :
        l_a = list_b
        l_b = list_a
        in_out = False

    for e in l_a :

        e_closest = get_closest(e,l_b)
        l_b.remove(e_closest)
        if in_out:
            pairs.append((e,e_closest))
        else :
            pairs.append((e_closest,e)) '''

# part of dfs implementation
def visit( e, visited, components, path, singles, ins, outs ) :

    #print "singles: ", singles
    #print "components: ", components
    if not e in visited :
        visited.add(e)
        path.add(e)

        for single_e in singles.keys() :
            #print single_e, e
            if components[e][1] == singles[single_e][0] :
                outs.append(single_e)
            if components[e][0] == singles[single_e][1] :
                ins.append(single_e)

        for adj_e in components.keys() :
            if components[e][1] == components[adj_e][0] or components[e][0] == components[adj_e][1] :
                visit (adj_e, visited, components, path, singles, ins, outs)


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


def get_ref_coord( edge_id, reference ) :

    ref_list = []
    #print reference
    #exit()
    for e in reference :
        if e[0] == edge_id :
            ref_list.append((e[1],e[2]))

    return ref_list

def close_to(a,b) :
    from math import fabs
    epsilon = 5
    return fabs(b - a) < epsilon

def alignments_to_ref( ins_list, outs_list, path, contigs ) :

    #print ins_list
    #print outs_list
    #print path
    #print edges
    #print contigs

    ins_ref_coord = []
    outs_ref_coord = []
    component_ref_coord = []
    #print 'ins_ref:',
    for e in ins_list :
        print 	get_ref_coord( e[0], contigs ),
        #ins_ref_coord.append( get_ref_coord( e, contigs ) )
    #print

    #print 'outs_ref:',
    for e in outs_list :
        print get_ref_coord( e[0], contigs ),
        #outs_ref_coord.append( get_ref_coord( e, contigs ) )
    #print

    #print 'component:',
    for e in path:
        print get_ref_coord( e, contigs ),
        #component_ref_coord.append( get_ref_coord( e, contigs ) )
    print

def compare_to_ref( resolved_path, contigs, siblings, nothing_close_number ) :

    #for e in resolved_path :
    #	print get_ref_coord( e, contigs ),

    #return_val = True
    path_len = 0
    #if len(resolved_path) < 3  :
    #	return False, nothing_close_number, too_short_number + 1

    coord_list_front = get_ref_coord(resolved_path[0], contigs)
    #coord_list_back = get_ref_coord(e_pair[1][0], contigs)

    close = False
    for start_e in coord_list_front :
        in_pos = start_e[1]
        path_len = in_pos - start_e[0]
        for e in resolved_path[1:] :
            coord_list = get_ref_coord(e, contigs)
            close = False
            for c in coord_list :
                out_pos = c[0]
                if close_to(in_pos, out_pos) :
                    close  = True
                    #print in_pos, out_pos
                    path_len += c[1] - c[0]
                    in_pos = c[1]
                    break


            if not close :
                break
        if close :
            break
    #if match :
    #	return False, nothing_close_number, too_short_number
    if not close :
        '''print zip( resolved_path, list(coverage[x] for x in resolved_path)  )
        for pair in siblings :
            print pair[0][0], '-', get_ref_coord(pair[0][0], contigs)
            print pair[1][0], '-', get_ref_coord(pair[1][0], contigs)
        print 'siblings:', siblings
        print 'resolved path: '
        for e in resolved_path :
            print e, '-', get_ref_coord(e, contigs)
        print'''
        return False, nothing_close_number + 1

        #else :
        #	if e == resolved_path[-1] :
        #		print "lll"
        #if not match :
        #	return False, nothing_close_number, too_short_number
    last_pos = get_ref_coord(resolved_path[-1],contigs)[0]
    path_len += last_pos[1] - last_pos[0]
    #print 'length in edges:', len(resolved_path), 'length in nts:', path_len
    return True, nothing_close_number


x = []
y = []
def add_pairs_to_plot( list_a, list_b ):

    l = min( len(list_a), len(list_b) )
    for i in range(l) :
        x.append(list_a[i][1])
        y.append(list_b[i][1])

def print_pairs(list_a, list_b, edges ) :

    if len(list_a) == 0 :
        print "ins: empty",
    else :
        print "ins:",

    for i in range(len(list_a)) :
        edge = list_a[i]
        print '(', edge[0], 'cov:', edge[1], 'l:', edges[edge[0]][2], ')',
    print

    if len(list_b) == 0 :
        print "outs: empty",
    else :
        print "outs:",


    for i in range(len(list_b)) :
        edge = list_b[i]
        print '(', edge[0], 'cov:', edge[1], 'l:', edges[edge[0]][2] , ')',

    print
    print

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

TEST =False 
if TEST :
    PREFIX = "/home/ksenia/coverage/"
    contigs = parse_pos("test2.pos")
    coverage = parse_cvr("test2.cvr")
    edges, redges = parse_grp("test2.grp")
    component = set([1,2,3,5,6])
    path = []
    path = dfs_comp(1,7,component,edges,path )
    print path
    #split_by_coverage(contigs, edges, 0)
    #components, singles = get_components( edges, redges )
    #traverse_components( components, singles, coverage, coverage, edges, redges, contigs )
else :
#	upper_edges = 10000
#	print "upper bound for edges: ", upper_edges
    name = "distance_estimation"
    #name = "constructed_graph"
    print "parsing pos.."
    contigs = parse_pos(name+".pos")
#    print contigs
    print "parsing cvr.."
    #in_coverage = parse_cvr("detail_in.cvr")
    in_coverage = parse_cvr(name+".cvr")
    #out_coverage = parse_cvr("detail_out.cvr")
    out_coverage = parse_cvr(name+".cvr")
    #print in_coverage
    #print out_coverage
    print "parsing grp.."
    edges, redges = parse_grp(name+".grp")
    #print edges
    #split_by_coverage(contigs, edges, 100)
    print "getting components..."
    components, singles = get_components( edges, redges )
    print "traversing graph.."
    traverse_components( components, singles, in_coverage, out_coverage, edges, redges, contigs )

    #pyplot.xlim([0,500])
    #pyplot.ylim([0,500])
    #pyplot.plot(x,y,marker='.',linestyle="")
    #pyplot.savefig('plot.png')
#print rank
#flush(rank)
