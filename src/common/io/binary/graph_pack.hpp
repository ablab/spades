//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "basic.hpp"
#include "pipeline/graph_pack.hpp"

namespace io {

namespace binary {

/**
 * @brief  This IOer processes the graph pack including only graph-related components.
 */
class BasePackIO : public IOBase<graph_pack::GraphPack> {
public:
    using Graph = debruijn_graph::Graph;
    using Type = graph_pack::GraphPack;

    void Save(const std::string &basename, const Type &gp) override;

    bool Load(const std::string &basename, Type &gp) override;

    virtual void BinWrite(std::ostream &os, const Type &gp);

    virtual bool BinRead(std::istream &is, Type &gp);

protected:
    BasicGraphIO<Graph> graph_io_;
};

/**
 * @brief  This IOer processes all of the graph pack components.
 */
class FullPackIO : public BasePackIO {
public:
    typedef BasePackIO base;
    typedef typename graph_pack::GraphPack Type;
    void Save(const std::string &basename, const Type &gp) override;

    bool Load(const std::string &basename, Type &gp) override;

    void BinWrite(std::ostream &os, const Type &gp) override;

    bool BinRead(std::istream &is, Type &gp) override;
};

} // namespace binary

} // namespace io
