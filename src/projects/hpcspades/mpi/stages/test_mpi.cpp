//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mpi/pipeline/mpi_stage.hpp"
#include "mpi/pipeline/partask_mpi.hpp"

#include <array>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

namespace debruijn_graph {
class ArraySum {
public:
    ArraySum(const std::string &message = "") : message_{message} {};
    ArraySum(const ArraySum&) = delete;
    ArraySum(ArraySum&&) = default;

    std::string message_;
    ArraySum(std::istream &is) { std::getline(is, message_); }

    std::ostream &serialize(std::ostream &os) const { return os << message_; }

    template <typename Data>
    auto make_splitter(size_t n, const Data &data) {
        size_t N = data.size();
        auto splitter = [N, n, i = size_t(0)](std::ostream &os, size_t /*node*/) mutable -> bool {
            if (i == n) return false;
            size_t begin = i * N / n;
            size_t end = (i + 1) * N / n;
            ++i;
            os << begin << " " << end << " ";
            return true;
        };

        return splitter;
    };

    template <typename Data>
    void process(std::istream &is, std::ostream &os, const Data &data) {
        INFO("MESSAGE: " << message_);
        long long int sum = 0;
#pragma omp parallel reduction(+ : sum)
        while (true) {
            size_t begin, end;
            bool exit = false;
#pragma omp critical
            {
                if (is.peek() == EOF || !(is >> begin >> end)) {
                    exit = true;
                } else {
                    DEBUG("Extracted range: " << begin << " " << end);
                }
            }
            if (exit) break;
            for (size_t i = begin; i < end; ++i) {
                sum += data[i];
            }
        }
        INFO("Computed sum: " << sum);
        os << sum;
    }

    auto merge(const std::vector<std::istream *> &piss, ...) {
        long long int sum = 0;
        for (auto &pis : piss) {
            long long int local_sum;
            *pis >> local_sum;
            sum += local_sum;
        }

        return sum;
    };
};

class TestMPI : public spades_mpi::MPIAssemblyStage {
public:
    TestMPI() : MPIAssemblyStage("Test MPI", "test_mpi") {}

    void run(graph_pack::GraphPack& /*gp*/, const char *) override {
        INFO("TestMPI started");
        partask::TaskRegistry treg;

        const size_t N = 100000;
        std::array<int, N> data;
        std::iota(data.begin(), data.end(), 1);

        auto job = treg.add<ArraySum>(std::cref(data));
        treg.listen();

        if (treg.master()) {
            auto res = job("Message1");
            INFO("JOB RESULT: " << res);
        }

        treg.stop_listening();
    }
};
}  // namespace debruijn_graph
