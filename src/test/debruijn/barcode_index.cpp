//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graphio.hpp"
#include "test_utils.hpp"
#include "barcode_index/barcode_index.hpp"
#include "barcode_index/barcode_index_builder.hpp"
#include "io/binary/read_cloud.hpp"
#include "random_graph.hpp"

#include <filesystem>
#include <gtest/gtest.h>

using namespace debruijn_graph;
using namespace barcode_index;

TEST(BarcodeIndex, BasicConstruction) {
    size_t K = 55;
    Graph graph(K);
    size_t frame_size = 100;

    FrameBarcodeIndex<Graph> barcode_index(graph, frame_size);

    EXPECT_EQ(barcode_index.GetFrameSize(), frame_size);
    EXPECT_TRUE(barcode_index.empty());
    EXPECT_EQ(barcode_index.size(), 0);
}

TEST(BarcodeIndex, InitialFillMap) {
    size_t K = 55;
    Graph graph(K);

    auto v1 = graph.AddVertex();
    auto v2 = graph.AddVertex();
    auto v3 = graph.AddVertex();

    Sequence seq1 = RandomSequence(100);
    Sequence seq2 = RandomSequence(150);

    auto e1 = graph.AddEdge(v1, v2, seq1);
    auto e2 = graph.AddEdge(v2, v3, seq2);

    size_t frame_size = 50;
    FrameBarcodeIndex<Graph> barcode_index(graph, frame_size);
    barcode_index.InitialFillMap();

    EXPECT_FALSE(barcode_index.empty());
    EXPECT_EQ(barcode_index.size(), graph.e_size());
    EXPECT_EQ(barcode_index.GetBarcodeNumber(e1), 0);
    EXPECT_EQ(barcode_index.GetBarcodeNumber(e2), 0);
}

TEST(BarcodeIndex, EdgeEntryConstruction) {
    size_t K = 55;
    Graph graph(K);
    auto v1 = graph.AddVertex();
    auto v2 = graph.AddVertex();
    Sequence seq = RandomSequence(200);
    auto edge = graph.AddEdge(v1, v2, seq);

    size_t frame_size = 50;
    size_t edge_length = graph.length(edge);

    FrameEdgeEntry<Graph> entry(edge, edge_length, frame_size);

    EXPECT_EQ(entry.GetEdge(), edge);
    EXPECT_EQ(entry.GetFrameSize(), frame_size);
    EXPECT_EQ(entry.Size(), 0);
    size_t expected_frames = edge_length / frame_size + 1;
    EXPECT_EQ(entry.GetNumberOfFrames(), expected_frames);
}

TEST(BarcodeIndex, BarcodeInfoUpdate) {
    size_t num_frames = 10;
    FrameBarcodeInfo info(num_frames);

    EXPECT_EQ(info.GetCount(), 0);
    EXPECT_EQ(info.GetCovered(), 0);

    info.Update(3, 2, 4);

    EXPECT_EQ(info.GetCount(), 3);
    EXPECT_EQ(info.GetCovered(), 3);
    EXPECT_EQ(info.GetLeftMost(), 2);
    EXPECT_EQ(info.GetRightMost(), 4);

    info.Update(2, 5, 7);

    EXPECT_EQ(info.GetCount(), 5);
    EXPECT_EQ(info.GetCovered(), 6);
    EXPECT_EQ(info.GetLeftMost(), 2);
    EXPECT_EQ(info.GetRightMost(), 7);
}

TEST(BarcodeIndex, ConcurrentBufferIO) {
    size_t K = 55;
    Graph graph(K);

    auto v1 = graph.AddVertex();
    auto v2 = graph.AddVertex();
    Sequence seq = RandomSequence(200);
    auto edge = graph.AddEdge(v1, v2, seq);

    size_t frame_size = 50;
    FrameConcurrentBarcodeIndexBuffer<Graph> buffer(graph, frame_size);
    buffer.InitialFillMap();

    FrameBarcodeIndex<Graph> barcode_index(graph, frame_size);
    barcode_index.InitialFillMap();

    EXPECT_EQ(buffer.GetFrameSize(), frame_size);

    BarcodeId barcode1 = 100;
    BarcodeId barcode2 = 200;
    Range range1(10, 60);
    Range range2(80, 120);
    buffer.InsertBarcode(barcode1, edge, 5, range1);
    buffer.InsertBarcode(barcode2, edge, 3, range2);

    barcode_index.Update(buffer);

    std::filesystem::path bmap_file = std::filesystem::temp_directory_path() / "test.bmap";
    io::binary::BarcodeMapperIO<Graph> io;
    io.Save(bmap_file, barcode_index);

    FrameBarcodeIndex<Graph> loaded_index(graph, 0);
    io.Load(bmap_file, loaded_index);

    EXPECT_EQ(loaded_index.GetFrameSize(), frame_size);
    EXPECT_EQ(loaded_index.size(), barcode_index.size());
    EXPECT_EQ(loaded_index.GetBarcodeNumber(edge), 2);
}

TEST(BarcodeIndex, BarcodeFiltering) {
    size_t K = 55;
    Graph graph(K);

    auto v1 = graph.AddVertex();
    auto v2 = graph.AddVertex();
    Sequence seq = RandomSequence(200);
    auto edge = graph.AddEdge(v1, v2, seq);

    size_t frame_size = 50;
    FrameBarcodeIndex<Graph> barcode_index(graph, frame_size);
    barcode_index.InitialFillMap();

    FrameConcurrentBarcodeIndexBuffer<Graph> buffer(graph, frame_size);
    buffer.InitialFillMap();

    buffer.InsertBarcode(100, edge, 10, Range(0, 50));
    buffer.InsertBarcode(200, edge, 2, Range(60, 100));
    buffer.InsertBarcode(300, edge, 8, Range(150, 180));

    barcode_index.Update(buffer);

    EXPECT_EQ(barcode_index.GetBarcodeNumber(edge), 3);

    size_t trimming_threshold = 5;
    size_t gap_threshold = 100;
    barcode_index.Filter(trimming_threshold, gap_threshold);

    EXPECT_EQ(barcode_index.GetBarcodeNumber(edge), 1);
}

TEST(BarcodeIndex, SimpleBarcodeInfo) {
    SimpleBarcodeInfo info(5, Range(10, 50));

    EXPECT_EQ(info.GetCount(), 5);
    EXPECT_EQ(info.GetRange().start_pos, 10);
    EXPECT_EQ(info.GetRange().end_pos, 50);

    info.Update(3, Range(5, 60));

    EXPECT_EQ(info.GetCount(), 8);
    EXPECT_EQ(info.GetRange().start_pos, 5);
    EXPECT_EQ(info.GetRange().end_pos, 60);

    SimpleBarcodeInfo other(2, Range(0, 40));
    info.Update(other);

    EXPECT_EQ(info.GetCount(), 10);
    EXPECT_EQ(info.GetRange().start_pos, 0);
    EXPECT_EQ(info.GetRange().end_pos, 60);
}

