//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "standard_vis.hpp"

namespace online_visualization {

  // write-once history
  class History
  {
    typedef HIST_ENTRY EntryT;

   public:
    void AddEntry(const string& entry) {
      add_history(entry.c_str());
      ++size_;
      VERIFY(int(size_) == history_get_history_state()->length);
    }

    const char* operator[](size_t k) const {
      VERIFY(k < size_);
      //EntryT** my_history = history_list();
      EntryT* entry = history_get(int(k + 1));
      return entry->line;
    }

    void SetEntry(size_t k, const string& entry) const {
      VERIFY(k < size_);
      //replace_history_entry(k, entry.c_str(), history_list()[k]->data);
      replace_history_entry(int(k), entry.c_str(), history_get(int(k + 1))->data);
    }

    size_t size() const {
      return size_;
    }

    const char* front() const {
      return this->operator[](0);
    }

    const char* back() const {
      return this->operator[](this->size() - 1);
    }

    static History& GetHistory() {
      static History hist;
      return hist;
    }

   private:
    size_t size_;

    History() : size_(0)
    {
    }
  };

}
