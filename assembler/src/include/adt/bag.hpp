#pragma once

template<class T, class hash = std::hash<T>>
class bag {
    typedef std::unordered_map<T, size_t, hash> Data;
    Data data_;
public:
    typedef typename Data::const_iterator const_iterator;

    void put(const T& t, size_t mult) {
        VERIFY(mult > 0);
        data_[t] += mult;
    }

    void put(const T& t) {
        put(t, 1);
    }

    bool take(const T& t, size_t mult) {
        VERIFY(mult > 0);
        /*typename map<T, size_t>::iterator*/auto it = data_.find(t);
        if (it == data_.end()) {
            return false;
        } else {
            size_t have = it->second;
            if (have < mult) {
                data_.erase(it->first);
                return false;
            } else if (have == mult) {
                data_.erase(it->first);
                return true;
            } else {
                it->second -= mult;
                return true;
            }
        }
    }

    bool take(const T& t) {
        return take(t, 1);
    }

    size_t mult(const T& t) const {
        auto it = data_.find(t);
        if (it == data_.end()) {
            return 0;
        } else {
            return it->second;
        }
    }

    void clear() {
        data_.clear();
    }

    const_iterator begin() const {
        return data_.begin();
    }

    const_iterator end() const {
        return data_.end();
    }

};
