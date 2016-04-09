//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

using namespace std;

namespace dipspades {

template<typename Id>
class RedundancyMap{
    map<Id, set<Id> > red_map_;
    set<Id> all_ids_;

    void ComputeAllIDs(){
        for(auto it = red_map_.begin(); it != red_map_.end(); it++){
            all_ids_.insert(it->first);
            all_ids_.insert(it->second.begin(), it->second.end());
        }
    }

public:
    void AddNewKey(Id key){
        set<Id> empty_set;
        red_map_[key] = empty_set;

        all_ids_.insert(key);
    }

    void AddNewPair(Id key, Id value){
        red_map_[key].insert(value);

        all_ids_.insert(key);
        all_ids_.insert(value);
    }

    set<Id> GetValuesByKey(Id key){
        return red_map_[key];
    }

    void SetValuesByKey(Id key, set<Id> values){
        red_map_[key].insert(values.begin(), values.end());
        all_ids_.insert(values.begin(), values.end());
        all_ids_.insert(key);
    }

    bool ContainsKey(Id key){
        return red_map_.find(key) != red_map_.end();
    }

    size_t AllElementsNumber(){

        ComputeAllIDs();

        return all_ids_.size();
    }

    set<Id> AllElements(){

        ComputeAllIDs();

        return all_ids_;
    }

    typename map<Id, set<Id> >::iterator begin(){
        return red_map_.begin();
    }

    typename map<Id, set<Id> >::iterator end(){
        return red_map_.end();
    }

    void Clear(){
        red_map_.clear();
        all_ids_.clear();
    }
};

template<typename Id>
class RedundancyMapCondenser{

    RedundancyMap<Id> uncondensed_map_;
    RedundancyMap<Id> condensed_map_;

    size_t number_processed_;
    size_t need_processed_;

    map<Id, bool> is_processed_;

    void ProcessCondensing(){
        bool non_zero_processed = true;
        while(number_processed_ < need_processed_ && non_zero_processed){
            int num_cur_processed = 0;
            for(auto it = condensed_map_.begin(); it != condensed_map_.end(); it++){
                set<Id> cur_set = it->second;

                for(auto it_set = cur_set.begin(); it_set != cur_set.end(); it_set++){
                    if(!is_processed_[*it_set]){
                        set<Id> child_set = uncondensed_map_.GetValuesByKey(*it_set);
                        it->second.insert(child_set.begin(), child_set.end());

                        is_processed_[*it_set] = true;
                        num_cur_processed++;
                    }
                }
            }
            non_zero_processed = num_cur_processed != 0;
            number_processed_ += num_cur_processed;
            TRACE("Number of processed - " << number_processed_ << ", total number - " << need_processed_);
        }

    }

    void ClearParams(){
        number_processed_ = 0;
        need_processed_ = 0;
        is_processed_.clear();
        condensed_map_.Clear();
    }

public:
    RedundancyMap<Id> Condense(RedundancyMap<Id> uncondensed_map){
        uncondensed_map_ = uncondensed_map;
        ClearParams();

        TRACE("Start condensing");

        TRACE("Computing of main keys");
        auto all_ids_ = uncondensed_map_.AllElements();
        map<Id, bool> is_main;
        for(auto it = all_ids_.begin(); it != all_ids_.end(); it++)
            is_main[*it] = true;

        for(auto it = uncondensed_map_.begin(); it != uncondensed_map_.end(); it++){
            for(auto it_set = it->second.begin(); it_set != it->second.end(); it_set++){
                is_main[*it_set] = false;
            }
        }

        set<Id> main_keys;
        for(auto it = is_main.begin(); it != is_main.end(); it++)
            if(it->second)
                main_keys.insert(it->first);

        TRACE("Number of all keys - " << all_ids_.size());
        TRACE("Number of main keys - " << main_keys.size());
        TRACE("Condensing starts");

        need_processed_ = all_ids_.size();
        number_processed_ = 0;

        for(auto it = all_ids_.begin(); it != all_ids_.end(); it++)
            is_processed_[*it] = false;

        for(auto main_key = main_keys.begin(); main_key != main_keys.end(); main_key++){
            condensed_map_.SetValuesByKey(*main_key, uncondensed_map_.GetValuesByKey(*main_key));
            number_processed_++;
            is_processed_[*main_key] = true;
        }

        // main processing
        ProcessCondensing();

        // processing of non visiting Ids
        while(number_processed_ < need_processed_){
            size_t max_child_setsize = 0;
            Id start_id(0);
            for(auto it = is_processed_.begin(); it != is_processed_.end(); it++){
                if(!it->second && uncondensed_map_.GetValuesByKey(it->first).size() >= max_child_setsize){
                    start_id = it->first;
                    max_child_setsize = uncondensed_map_.GetValuesByKey(it->first).size();
                }
            }
            auto start_set = uncondensed_map_.GetValuesByKey(start_id);
            for(auto it = start_set.begin(); it != start_set.end(); it++)
                if(!is_processed_[*it])
                    condensed_map_.AddNewPair(start_id, *it);

            is_processed_[start_id] = true;
            number_processed_++;
            ProcessCondensing();
        }

        VERIFY(number_processed_ == need_processed_);
        return condensed_map_;
    }
};

template<typename Id>
class RedundancyMapMerger{

    bool AreMergeResultsCorrect(RedundancyMap<Id> old_map, RedundancyMap<Id> new_map){
//        cout << "Correctness - " << old_map.AllElementsNumber() << " " << new_map.AllElementsNumber() << endl;
        return old_map.AllElementsNumber() == new_map.AllElementsNumber();
    }

public:
    RedundancyMap<Id> MergeTwoMaps(RedundancyMap<Id> map1, RedundancyMap<Id> map2){

        for(auto it_old = map1.begin(); it_old != map1.end(); it_old++){
            Id old_key = it_old->first;
            auto old_set = it_old->second;

            if(map2.ContainsKey(old_key)){
                map2.SetValuesByKey(old_key, old_set);
            }
            else{
                bool is_found = false;

                for(auto it_new = map2.begin(); it_new != map2.end(); it_new++){
                    Id new_key = it_new->first;
                    auto new_set = it_new->second;
                    if(new_set.find(old_key) != new_set.end()){
                        map2.SetValuesByKey(new_key, old_set);
                        is_found = true;
                        break;
                    }
                }

                if(!is_found)
                    map2.SetValuesByKey(old_key, old_set);
            }
        }
        RedundancyMapCondenser<Id> condenser;
        map2 = condenser.Condense(map2);
        VERIFY(AreMergeResultsCorrect(map1, map2));
        return map2;
    }
};

}
