############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys

#from rrr.py
def parser_options():
    parser = OptionParser()
    parser.add_option("-s", "", dest="saves_dir", help="Name of directory with saves")
    parser.add_option("-o", "", dest="out_dir", help="Output directory, default = out (optional)", default="out")
    parser.add_option("-g", "", dest="genome", help="File with genome (optional)")
    parser.add_option("-d", "", dest="debug_logger", help="File for debug logger (optional)", default="debug_log.txt")
    parser.add_option("-k", "", type=int, dest="k", help="k (optional)")
    parser.add_option("-D", "", type=int, dest="d", help="d (optional)")
    parser.add_option("", "--sc", dest="sc", action="store_true", help="Turn on if data is sincle-cell (optional)",
        default=False)
    return parser


def make_rectangles_from_genome(options):
    k = options.k
    ingraph = Graph()
    _, genome = fastaparser.read_fasta(options.genome).next()
    ingraph.make_graph(genome, int(k))
    ingraph.save(os.path.join(options.out_dir, "graph"))
    rs = RectangleSet(ingraph, int(options.d))
    rs.filter_without_prd()
    edges_before_loop_DG = ingraph.find_loops(10, 1000, rs)
    f_left = open(os.path.join(options.out_dir, "paired_genom_contigs_1.fasta"), "w") # TODO: what is it?
    f_right = open(os.path.join(options.out_dir, "paired_genom_contigs_2.fasta"), "w") # TODO: what is it?
    contigs_id = 0
    for key, rect in rs.rectangles.items():
        for key, diag in rect.diagonals.items():
            e1 = rect.e1.seq
            e2 = rect.e2.seq
            f_left.write(">" + str(contigs_id) + "/1\n")
            f_left.write(e1[diag.offseta:diag.offsetc])
            f_left.write("\n")
            f_right.write(">" + str(contigs_id) + "/2\n")
            f_right.write(e2[diag.offsetb:diag.offsetd])
            f_right.write("\n")
            contigs_id += 1
    bgraph = rs.bgraph_from_genome()
    bgraph.condense()
    outgraph = bgraph.project(options.out_dir, False)
    outgraph.fasta(open(os.path.join(options.out_dir, 'rectangles.fasta'), 'w'))

if __name__ == '__main__':
    ##########
    # PARAMS #
    ##########
    parser = parser_options()
    (options, args) = parser.parse_args()

    if not os.path.exists(options.out_dir):
        os.mkdir(options.out_dir)

    if options.genome and not options.saves_dir:
        if not options.k or not options.d:
            print "specify k and d"
            sys.exit(1)
        make_rectangles_from_genome(options)
        sys.exit(1)

    if not options.saves_dir:
        parser.print_help()
        sys.exit(0)

    input_dir = options.saves_dir
    outpath = options.out_dir
    reference_information_file = os.path.join(input_dir, "late_pair_info_counted_etalon_distance.txt")
    test_util = TestUtils(reference_information_file, os.path.join(outpath, options.debug_logger))
    resolve(input_dir, outpath, test_util, options.genome, options.sc)

    if test_util.has_ref_info:
        test_util.stat()

#from bigraph.py
def build_missing_rectangles(self, K, rectangles):
        return
        threshold = self.d
        self.test_utils.logger.info("treshold " + str(threshold))
        count_ovelaps = 0
        count_miss_rect = 0
        count_miss_path = 0
        true_miss_path = 0
        count_overlaps = 0
        v1s = set()
        v2s = set()
        for bv in self.vs.itervalues():
            if len(bv.inn) == 1 and len(bv.out) == 0 and len(bv.inn[0].get_seq(K, self.d)) > 3 * self.d:
                v1s.add(bv)
            if len(bv.inn) == 0 and len(bv.out) == 1 and len(bv.out[0].get_seq(K, self.d)) > 3 * self.d:
                v2s.add(bv)
        assert len(v1s) == len(v2s) # because of rev-compl
        self.test_utils.logger.info("v1s.len " + str(len(v1s)))

        all_paired_paths = []
        for v1 in v1s:
            be1 = v1.inn[0]
            diag1 = be1.diagonals[-1]
            for v2 in v2s:
                be2 = v2.out[0]
                if (be1.eid == be2.eid):
                    continue
                diag2 = be2.diagonals[0]
                paths1 = abstract_graph.find_paths(diag1.rectangle.e1.v1, diag2.rectangle.e1.v1, diag1.rectangle.e1,
                    threshold + diag1.rectangle.e1.len, DEEP_THRESHOLD)
                paths2 = abstract_graph.find_paths(diag1.rectangle.e2.v1, diag2.rectangle.e2.v1, diag1.rectangle.e2,
                    threshold + diag1.rectangle.e2.len, DEEP_THRESHOLD)
                paired_paths = find_pair_paths(paths1, paths2, diag1, diag2)
                if len(paired_paths) != 0:
                    all_paired_paths.append((paired_paths, diag1, diag2))
        self.test_utils.logger.info("all_paired_paths " + str(len(all_paired_paths)))
        can_find_one_path_more = True
        added_paths = []
        while can_find_one_path_more:
            the_best_path = None
            the_best_support = 0
            can_find_one_path_more = False

            for paired_paths in all_paired_paths:
                (best_support, best_len, best_rectangles, best_diags, best_path) = self.choose_best_path(paired_paths[0]
                    , rectangles, paired_paths[1], paired_paths[2], self.d, added_paths)
                if best_support / best_len >= 0.0 and best_support / best_len > the_best_support:
                    the_best_support = best_support / best_len
                    the_best_path = (best_support, best_len, best_rectangles, best_diags, best_path)
            if the_best_path:
                added_paths.append(the_best_path[-1])
                (best_support, best_len, best_rectangles, best_diags, best_path) = the_best_path
                can_find_one_path_more = True
                prev_diag = best_diags[0]
                true_path = True
                for diag in best_diags[1:]:
                    if prev_diag:
                        should_connect = self.test_utils.should_join(prev_diag, diag)
                        if not should_connect:
                            true_path = False
                        self.add_diagonal_and_conj(diag)
                        is_true = self.test_utils.is_true_diagonal(diag)
                        if not is_true:
                            true_path = False
                        count_miss_rect += 1
                        prev_diag = diag
                count_miss_path += 1
                if true_path:
                    true_miss_path += 1
        self.test_utils.logger.info("count_overlap " + str(count_ovelaps) + " count_miss_rect " + str(
            count_miss_rect) + " count miss path " + str(count_miss_path) + " true miss path " + str(true_miss_path))


    def choose_best_path(self, paired_paths, rectangeles_set, diag1, diag2, d, added_paths):
        best_support = 0
        best_len = 10000
        best_rectangles = []
        best_diags = []
        best_path = paired_paths[0]
        best_not_supported = 0

        for paired_path in paired_paths:
            (path1, path2, path_len) = paired_path
            if paired_path in added_paths:
                continue
            first_shift = diag1.offseta
            second_shift = diag1.offsetb
            path1.append(diag2.rectangle.e1)
            path2.append(diag2.rectangle.e2)
            rectangles = []
            diags = []
            not_supported = []
            path_support = 0
            pos_first_path = 0
            pos_second_path = 0
            first_len = first_shift
            make_less_N50 = False
            while not make_less_N50 and first_len < path_len + diag2.offseta:
                ed1 = path1[pos_first_path]
                ed2 = path2[pos_second_path]
                rectangle = Rectangle(ed1, ed2)
                rectangle.add_diagonal(d, d + first_shift - second_shift)
                rect_diag = rectangle.get_closest_diagonal(d + first_shift - second_shift)
                if (not (rect_diag.key1 == diag1.key1 and rect_diag.key2 == diag1.key2) and not(
                rect_diag.key1 == diag2.key1 and rect_diag.key2 == diag2.key2)):
                    can_use = [diag1.key1, diag1.key2, diag2.key1, diag2.key2]
                    if (rect_diag.key1 in self.vs and rect_diag.key1 not in can_use) or  (
                    rect_diag.key2 in self.vs and rect_diag.key2 not in can_use):
                        make_less_N50 = True
                        continue
                diags.append(rect_diag)
                rectangles.append(rectangle)
                rectangeles_set.use_prd_diag(rect_diag)
                #if rect_diag.prd_support < 0.00001 and (ed2.len > 10 and ed1.len > 10):
                #  make_less_N50 = True
                #  continue
                path_support += rect_diag.prd_support
                if ed2.len - second_shift < ed1.len - first_shift:
                    pos_second_path += 1
                    first_shift += ed2.len - second_shift
                    first_len += ed2.len - second_shift
                    if rect_diag.prd_support < 0.000000001:
                        not_supported.append(ed2.len - second_shift)
                    second_shift = 0
                elif ed1.len - first_shift < ed2.len - second_shift:
                    pos_first_path += 1
                    first_len += ed1.len - first_shift
                    second_shift += ed1.len - first_shift
                    if rect_diag.prd_support < 0.000000001:
                        not_supported.append(ed1.len - first_shift)
                    first_shift = 0
                else:
                    first_len += ed1.len - first_shift
                    pos_second_path += 1
                    pos_first_path += 1
                    first_shift = 0
                    second_shift = 0
            if not make_less_N50 and path_len > 1 and  path_support / path_len < 1000 and  path_support / path_len > best_support:
                best_support = path_support
                best_len = path_len
                best_rectangles = rectangles
                best_diags = diags
                best_path = (path1, path2, path_len)
                best_not_supported = not_supported
        return (best_support, best_len, best_rectangles, best_diags, best_path)


def find_pair_paths(paths1, paths2, diag1, diag2):
    paired_paths = []
    for (path1, len1) in paths1:
        for (path2, len2) in paths2:
            if path1[0] != diag1.rectangle.e1 or path2[0] != diag1.rectangle.e2:
                continue
            if len1 - diag1.offseta + diag2.offseta == len2 - diag1.offsetb + diag2.offsetb:
                paired_paths.append((path1, path2, len1))
    return paired_paths
  
  
### from graph.py
    def make_graph(self, genome, k):
        self.K = k
        kmers = self.__get_kmers_pos(genome, k)
        visit = set()
        vid = 0
        eid = 0
        edges = set()
        verts = dict()
        for key in kmers:
            if key in visit:
                continue
            body = [key[-1]]
            end_vertex = key[1:]
            while True:
                next_kmer = extend_forward(end_vertex, kmers)
                if next_kmer == None:
                    break
                body.append(next_kmer[-1])
                end_vertex = next_kmer[1:]
                visit.add(next_kmer)
                visit.add(utils.rc(next_kmer))

            begin_vertex = key[:-1]
            while True:
                next_kmer = extend_backward(begin_vertex, kmers)
                if next_kmer == None:
                    break
                body.insert(0, next_kmer[-1])
                begin_vertex = next_kmer[0:-1]
                visit.add(next_kmer)
                visit.add(utils.rc(next_kmer))

            body = begin_vertex + ''.join(body)
            if begin_vertex not in verts:
                begin_ref = self.add_vertex(vid, vid + 1)
                r_end_ref = self.add_vertex(vid + 1, vid)
                verts[begin_vertex] = begin_ref.vid
                verts[utils.rc(begin_vertex)] = r_end_ref.vid
                vid += 2
            if end_vertex not in verts:
                end_ref = self.add_vertex(vid, vid + 1)
                r_begin_ref = self.add_vertex(vid + 1, vid)
                verts[end_vertex] = end_ref.vid
                verts[utils.rc(end_vertex)] = r_begin_ref.vid
                vid += 2
            bv = verts[begin_vertex]
            ev = verts[end_vertex]
            rbv = verts[utils.rc(end_vertex)]
            rev = verts[utils.rc(begin_vertex)]
            if (bv, ev) not in edges:
                if (bv, ev) == (rbv, rev) and body == utils.rc(body):
                    self.add_edge(eid, bv, ev, len(body) - k + 1, eid)
                    edges.add((bv, ev))
                    self.add_seq(eid, body)
                    self.etalon_dist[eid] = kmers[body[:k]] + kmers[utils.rc(body)[:k]]
                    eid += 1
                else:
                    self.add_edge(eid, bv, ev, len(body) - k + 1, eid + 1)
                    self.add_edge(eid + 1, rbv, rev, len(body) - k + 1, eid)
                    edges.add((bv, ev))
                    edges.add((rbv, rev))
                    self.add_seq(eid, body)
                    self.add_seq(eid + 1, utils.rc(body))
                    self.etalon_dist[eid] = kmers[body[:k]]
                    self.etalon_dist[eid + 1] = kmers[utils.rc(body)[:k]]
                    eid += 2

alphabet = "ACGT"

def find_out_edges(vertex_body, kmer_map):
    next_kmer = [(vertex_body + base) for base in alphabet]
    return [kmer for kmer in next_kmer if kmer in kmer_map]


def find_in_edges(vertex_body, kmer_map):
    next_kmer = [(base + vertex_body) for base in alphabet]
    return [kmer for kmer in next_kmer if kmer in kmer_map]


def extend_forward(vertex_body, kmer_map):
    return extend_in_direction(vertex_body, kmer_map, True)

def extend_backward(vertex_body, kmer_map):
    return extend_in_direction(vertex_body, kmer_map, False)

def extend_in_direction(vertex_body, kmer_map, direction_forward):
    in_edge = find_in_edges(vertex_body, kmer_map)
    out_edge = find_out_edges(vertex_body, kmer_map)
    res = out_edge if direction_forward else in_edge
    if len(in_edge) == 1 and len(out_edge) == 1:
        return res[0]
    return None

    def dfs(self, e, d):
        limit1 = d - e.len
        limit2 = d
        if e.len > d:
            yield e, 0
        ls = [set() for _ in xrange(limit2)]
        ls[0].add(e.v2)
        all_dist = dict()
        if self.__from_genome():
            all_dist[(0, e.v2.vid)] = (e, self.etalon_dist[e.eid])
        for pos in xrange(limit2):
            for v in ls[pos]:
                if self.__from_genome():
                    (prev_e, dists) = all_dist[(pos, v.vid)]
                for e2 in v.out:
                    if self.__from_genome():
                        new_dists = []
                        for dist in dists:
                            if dist >= 0 and dist + prev_e.len + 1 in self.etalon_dist[e2.eid]:
                                new_dists.append(dist + prev_e.len + 1)
                            elif dist <= 0 and -(-dist + prev_e.len + 1) in self.etalon_dist[e2.eid]:
                                new_dists.append(dist - prev_e.len - 1)
                        if len(new_dists) == 0:
                            continue
                    pos2 = pos + e2.len
                    if pos2 < limit2:
                        ls[pos2].add(e2.v2)
                        if self.__from_genome():
                            if (pos2, e2.v2.vid) in all_dist:
                                new_dists += all_dist[(pos2, e2.v2.vid)][1]
                            all_dist[(pos2, e2.v2.vid)] = (e2, new_dists)
                    if pos + e2.len > limit1:
                        yield e2, pos + e.len
    def __from_genome(self):
        return len(self.etalon_dist.keys()) > 0

#from rectangle_set.py
    def filter_without_prd(self):
        self.__build_from_graph()
        self.__conjugate()
    
    def percentiles(self):
        percentiles = [utils.quantile(self.ranking, percentile).support() for percentile in xrange(101)]
        return percentiles

    def bgraph_from_genome(self):
        bg = BGraph(self.graph, self.d, self.test_utils)
        for key, rect in self.rectangles.items():
            for key, diag in rect.diagonals.items():
                bg.add_diagonal(diag)
        return bg

    def get_support(self, e1, e2):
        if (e1, e2) in self.rectangles:
            rectangle = self.rectangles[(e1, e2)]
            return max([diag.support() for diag in rectangle.diagonals.itervalues()])
        return 0

