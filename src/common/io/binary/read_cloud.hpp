//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "io_base.hpp"
#include "barcode_index/barcode_index.hpp"
//fixme added temporarily to compile, remove later
#include "barcode_index/barcode_index_builder.hpp"
#include "barcode_index/barcode_info_extractor.hpp"

namespace io {

namespace binary {

template<typename Graph>
class BarcodeMapperIO : public IOSingle<barcode_index::FrameBarcodeIndex<Graph>> {
public:
    typedef barcode_index::FrameBarcodeIndex<Graph> Type;
    BarcodeMapperIO()
        : IOSingle<Type>("barcode index", ".bmap") {
    }

    void SaveImpl(BinOStream &str, const Type &mapper) override {
        str << mapper.GetFrameSize() << mapper;
    }

    void LoadImpl(BinIStream &str, Type &mapper) override {
        size_t frame_size;
        str >> frame_size;
        mapper.SetFrameSize(frame_size);
        str >> mapper;
    }
};

template<typename Graph>
struct IOTraits<barcode_index::FrameBarcodeIndex<Graph>> {
    typedef BarcodeMapperIO<Graph> Type;
};

}

}