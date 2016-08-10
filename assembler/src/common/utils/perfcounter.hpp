//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include <sys/time.h>
#include <string>
#include <cppformat/format.h>

struct perf_counter
{
    perf_counter()
    {
        reset();
    }

    double time() const
    {
        struct timeval now;
        gettimeofday(&now, NULL);

        return (double)(now.tv_sec - time_.tv_sec) + (double)(now.tv_usec - time_.tv_usec) * 1e-6;
    }

    double time_ms() const
    {
        return time() * 1e3;
    }

    void reset()
    {
        gettimeofday(&time_, NULL);
    }

private:
    struct timeval time_;
};


inline std::string human_readable_time(double time_in_sec)
{
//    assert(time_in_sec > 0);

    size_t msec  = size_t(time_in_sec * 1000) % 1000;
    size_t sec   = size_t(time_in_sec);
    size_t hours = sec / 3600;
    size_t mins  = (sec / 60) % 60;
    sec         %= 60;

    return fmt::format("{:3d}:{:02d}:{:02d}.{:03d}", hours, mins, sec, msec);
}

inline std::string human_readable_memory(size_t max_rss) {
    if (max_rss < 1024 * 1024) {
        return fmt::format("{:d}M", (max_rss / 1024));
    } else {
        return fmt::format("{:d}G", (max_rss / (1024 * 1024)));
    }
}

struct avg_perf_counter
{
    avg_perf_counter(/*const string& name*/)// : name_(name)
    {
        reset();
    }

//    ~avg_perf_counter() {
//        cout << "Time in counter " << name_ << ": " << human_readable_time(time()) << endl;
//    }

    int start(int ret = 0)
    {
        p_cnt_.reset();
        return ret;
    }

    int stop(int ret = 0)
    {
        counter_++;
        whole_time_ += p_cnt_.time();
        return ret;
    }
    double time() const
    {
        return whole_time_;
    }
    size_t counts()
    {
        return counter_;
    }
    double time_ms() const
    {
        return time() * 1e3;
    }

    double avg_time() const
    {
        return counter_ > 0 ? whole_time_/(double)counter_ : 0.;
    }

    double avg_time_ms() const
    {
        return avg_time() * 1e3;
    }

    void reset()
    {
        p_cnt_.reset();
        whole_time_ = 0;
        counter_ = 0;
    }

private:
  const std::string name_;
  perf_counter p_cnt_;
  double whole_time_;
  size_t counter_;

};
