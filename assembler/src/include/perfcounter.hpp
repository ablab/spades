//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once
#include <sys/time.h>

struct perf_counter
{
    perf_counter()
    {
        reset();
    }

    double time() const
    {
        timeval now;
        gettimeofday(&now, NULL);

        return (now.tv_sec - time_.tv_sec) + (now.tv_usec - time_.tv_usec) * 1e-6;
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
    timeval time_;
};


struct avg_perf_counter
{
    avg_perf_counter()
    {
        reset();
    }

    void start()
    {
        p_cnt_.reset();
    }

    {
        counter_++;
        whole_time_ += p_cnt_.time();
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
        return counter_ > 0 ? whole_time_/counter_ : 0.;
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
    perf_counter p_cnt_;
    double whole_time_;
    size_t counter_;

};


inline string human_readable_time(double time_in_sec)
{
    assert(time_in_sec > 0);

    size_t msec  = size_t(time_in_sec * 1000) % 1000;
    size_t sec   = size_t(time_in_sec);
    size_t hours = sec / 3600;
    size_t mins  = (sec / 60) % 60;
    sec         %= 60;

    return str(format("%3d:%02d:%02d.%03d") % hours % mins % sec % msec);
}
