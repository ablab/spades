#!/usr/bin/env python3
import argparse
import os
import sys

RC = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def get_rc(sequence):
    return "".join([RC[nuc] for nuc in sequence[::-1]])


def oriented_seq(sequence, orient):
    return sequence if orient == "+" else get_rc(sequence)


def parse_gfa(path):
    segments = {}
    links = []
    paths = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields[0] == "S":
                segments[fields[1]] = fields[2]
            elif fields[0] == "L":
                overlap = int(fields[5].rstrip("M"))
                links.append((fields[1], fields[2], fields[3], fields[4], overlap))
            elif fields[0] == "P":
                edge_ids = fields[2].split(",")
                edges = [(edge_id[:-1], edge_id[-1]) for edge_id in edge_ids]
                overlaps = [int(o.rstrip("M")) for o in fields[3].split(",")]
                paths.append((fields[1], edges, overlaps))
    return segments, links, paths


def reconstruct_path_sequence(segments, edges_with_orient, overlaps):
    edge_id, orient = edges_with_orient[0]
    seq = oriented_seq(segments[edge_id], orient)
    for (edge_id, orient), ov in zip(edges_with_orient[1:], overlaps):
        edge_seq = oriented_seq(segments[edge_id], orient)
        seq += edge_seq[ov:]
    return seq


def check_segments_in_paths(resolved_segments, path_sequences):
    passed = True
    for seg_id, seg_seq in resolved_segments.items():
        seg_rc = get_rc(seg_seq)
        matched = []
        for path_name, path_seq in path_sequences:
            if seg_seq in path_seq or seg_rc in path_seq:
                matched.append(path_name)
        if matched:
            print("  segment {}: found in {}".format(seg_id, ", ".join(matched)))
        else:
            print("  FAIL segment {}: not found in any path".format(seg_id))
            passed = False
    return passed


def check_link_consistency(segments, links):
    passed = True
    for from_id, from_orient, to_id, to_orient, overlap in links:
        if from_id not in segments or to_id not in segments:
            print("  SKIP link {} {} -> {} {}: missing segment"
                  .format(from_id, from_orient, to_id, to_orient))
            continue
        suffix = oriented_seq(segments[from_id], from_orient)[-overlap:]
        prefix = oriented_seq(segments[to_id], to_orient)[:overlap]
        if suffix == prefix:
            print("  link {} {} -> {} {} ({}M): OK"
                  .format(from_id, from_orient, to_id, to_orient, overlap))
        else:
            print("  FAIL link {} {} -> {} {} ({}M): overlap mismatch"
                  .format(from_id, from_orient, to_id, to_orient, overlap))
            passed = False
    return passed


def parse_edge_transform(path):
    mapping = {}
    with open(path) as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip().split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def get_merged_pairs(path):
    pairs = []
    with open(path) as handle:
        header = handle.readline().rstrip().split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields[col["Vertex result"]] in ("Completely", "Partially"):
                answer = fields[col["Answer"]]
                for token in answer.split(","):
                    in_edge, out_edge = token.split("#")
                    pairs.append((in_edge, out_edge))
    return pairs


def createparser():
    parser = argparse.ArgumentParser(
        description="Assess splitter accuracy on a simulated dataset.")
    parser.add_argument("--splitter-output", required=True,
                        help="Splitter output directory containing "
                             "resolved_graph.gfa, edge_transform.tsv, vertex_stats.tsv")
    parser.add_argument("--dataset-gfa", required=True)
    return parser


def main():
    args = createparser().parse_args()
    ok = True

    out = args.splitter_output
    resolved_gfa = os.path.join(out, "resolved_graph.gfa")
    edge_transform_path = os.path.join(out, "edge_transform.tsv")
    vertex_stats_path = os.path.join(out, "vertex_stats.tsv")

    resolved_segments, resolved_links, _ = parse_gfa(resolved_gfa)
    dataset_segments, _, dataset_paths = parse_gfa(args.dataset_gfa)

    path_sequences = []
    for name, edges, overlaps in dataset_paths:
        seq = reconstruct_path_sequence(dataset_segments, edges, overlaps)
        path_sequences.append((name, seq))

    print("=== Segment substring check ===")
    if not check_segments_in_paths(resolved_segments, path_sequences):
        ok = False

    print("\n=== Link consistency check ===")
    if not check_link_consistency(resolved_segments, resolved_links):
        ok = False

    # print("\n=== Edge transform ===")
    # edge_map = parse_edge_transform(edge_transform_path)
    # for orig, resolved in sorted(edge_map.items()):
    #     print("  {} -> {}".format(orig, resolved))

    print("\n=== Merged edge pairs ===")
    merged = get_merged_pairs(vertex_stats_path)
    for in_edge, out_edge in merged:
        print("  {} -> {}".format(in_edge, out_edge))

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
