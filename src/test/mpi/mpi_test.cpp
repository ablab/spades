//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <array>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>
#include <mutex>

#include "mpi/pipeline/partask_mpi.hpp"

#include "mpi/utils/logger/mpi_log_writers.hpp"

void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<mpi_console_writer>());
    attach_logger(lg);
}

class ArraySum {
public:
    ArraySum(const std::string &message = "") : message_{message} {};
    ArraySum(const ArraySum&) = delete;
    ArraySum(ArraySum&&) = default;
    ~ArraySum() noexcept {
        std::cout << "~ArraySum() node:" << partask::world_rank() << std::endl;
    }

    std::string message_;
    ArraySum(std::istream &is) {std::getline(is, message_);}

    std::ostream &serialize(std::ostream &os) const { return os << message_; }

    template <typename Data, typename... Args>
    auto make_splitter(size_t n, const Data &data, Args &&...) {
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
    }

    template <typename Data, typename... Args>
    void process(std::istream &is, std::ostream &os, const Data &data, Args &&...) {
        std::cout << "process run" << std::endl;
        std::cout << "MESSAGE: " << message_ << std::endl;
        long long int sum = 0;
#pragma omp parallel reduction(+ : sum)
        while (true) {
            size_t begin, end;
            bool exit = false;
#pragma omp critical
            {
                if (is.peek() == EOF || !(is >> begin >> end)) exit = true;
                if (!exit) std::cout << "Extracted range: " << begin << " " << end << std::endl;
            }
            if (exit) break;
            for (size_t i = begin; i < end; ++i) {
                sum += data[i];
            }
        }
        std::cout << "Computed sum: " << sum << std::endl;
        os << sum;
    }

    template <typename... Args>
    auto merge(const std::vector<std::istream *> &piss, Args &&...) {
        long long int sum = 0;
        for (auto &pis : piss) {
            long long int local_sum;
            *pis >> local_sum;
            sum += local_sum;
        }

        return sum;
    }
};

const size_t N = 100000;
std::array<int, N> data;
int main() {
    create_console_logger();
    INFO("Starting mpi test");

    {
        int a = 1, b = 2, c = 3;

        std::stringstream ss;
        auto tpl = std::make_tuple(a, b, c, std::string("21312321"));
        io::binary::BinWrite(ss, std::make_tuple(int(3), int(4), int(5), std::string("123456")));
        io::binary::BinRead(ss, tpl);
        INFO("Read: " << std::get<0>(tpl) << ":" << std::get<1>(tpl) << ":" << std::get<2>(tpl) << ":" << std::get<3>(tpl));

        io::binary::BinWrite(ss, std::tie(a, b));
        a = 42, b = 32;
        io::binary::BinRead(ss, std::tie(a, b));
        INFO("Read: " << a << ":" << b);

        io::binary::BinWrite(ss, tpl);
        io::binary::BinWrite(ss, std::vector<int>(10));
    }

    bool init = partask::init();
    INFO("MPI init: " << (init ? "done" : "failed"));
    {
        partask::ChunkedStringStream mss;
        mss << "mama,papa,pipi" << std::endl;
        std::string s;
        mss >> s;
        std::cout << "Test: " << s;
    }

    std::cout << "\n\n\n\n";


    std::unordered_map<std::string, std::vector<int>> m, m2;
    m["mama"] = {1, 2, 4};
    m["papa"] = {1, 3, 4, 2};

    std::stringstream ss;
    io::binary::BinWrite(ss, m);
    io::binary::BinRead(ss, m2);
    std::cout << m2["papa"][1] << "  <-- should be "<< m["papa"][1] << std::endl;;

    io::binary::BinWrite(std::cout, std::string("Mama"), 123, std::string("Pipi"));


    io::binary::BinWrite(std::cout, std::string("Mama"), 123, std::string("Pipi"));

    std::iota(data.begin(), data.end(), 1);
    size_t sum = std::accumulate(data.cbegin(), data.cend(), size_t(0));
    std::cout << "Actual sum: " << sum << std::endl;



    if (partask::world_rank() == 0) {
        std::string t = "0123456789";
        std::string s = "";
        for (size_t i = 0; i < 25; ++i) {
            s += t;
        }
        partask::OutputMPIStream os(1);
        // os << s;
        os << "checkX";
        os << "checkY";
        os << "checkZ";
        os << "checkA";
        os << s;
        os << "checkB";
        os << s;
        os << "checkC";
        os << "checkD";
        os << "checkE";
        os.flush();
        os.flush();
    }
    if (partask::world_rank() == 1) {
        partask::InputMPIStream is(0);
        std::string s;
        is >> s;
        std::cout << "Streams test" << std::endl;
        std::cout << s << std::endl;
        std::cout << "Streams test" << std::endl;
    }

    partask::barrier();
    std::cout << "broadcast test" << std::endl;
    partask::broadcast(m);
    std::cout << "broadcast test done" << std::endl;
    partask::barrier();

    partask::all_set_num_threads(10);

    partask::TaskRegistry reg;

    auto ptr = std::make_unique<std::string>("Mama");  // Non-copyable
    std::mutex mtx;  // Non-copyable
    auto job = reg.add<ArraySum>(std::cref(data), std::cref(ptr), std::ref(mtx));
    reg.listen();

    if (reg.master()) {
        auto res = job("Message1");
        std::cout << "JOB RESULT: " << res << std::endl;
        res = job("Message2");
        std::cout << "JOB RESULT: " << res << std::endl;
        res = job("Message3");
        std::cout << "JOB RESULT: " << res << std::endl;
        res = job();
        std::cout << "JOB RESULT: " << res << std::endl;
    }

    reg.stop_listening();
    std::cout << "Before the barrier " << __LINE__;
    partask::barrier();

    auto job2 = reg.add<ArraySum>(std::cref(data));
    reg.listen();

    if (reg.master()) {
        auto res = job();
        std::cout << "JOB RESULT: " << res << std::endl;
        res = job2();
        std::cout << "JOB RESULT: " << res << std::endl;
        res = job();
        std::cout << "JOB RESULT: " << res << std::endl;
        res = job();
        std::cout << "JOB RESULT: " << res << std::endl;
    }

    reg.stop_listening();
    reg.stop_listening();

    partask::finalize();
    return 0;
}
