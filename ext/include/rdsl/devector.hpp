/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * MIT License
 *
 * Copyright (c) 2022 Valasiadis Fotios
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * devector.hpp 0.0.0
 * 
 * A header-only continuous-storage double ended vector implementation for efficient insertions/removals at both the start & end of it.
 */

#ifndef DEVECTOR_RDSL_28092021
#define DEVECTOR_RDSL_28092021

#include <memory>
#include <stdexcept>
#include <iterator>

namespace rdsl{

template<class al>
using al_traits = std::allocator_traits<al>;

template<class it>
using it_traits = std::iterator_traits<it>;

template<class It>
struct is_at_least_forward{
    static constexpr bool value = false;
};

template<>
struct is_at_least_forward<std::forward_iterator_tag>{
    static constexpr bool value = true;
};

template<>
struct is_at_least_forward<std::bidirectional_iterator_tag>{
    static constexpr bool value = true;
};

template<>
struct is_at_least_forward<std::random_access_iterator_tag>{
    static constexpr bool value = true;
};

template<class It>
struct is_at_least_input{
    static constexpr bool value = is_at_least_forward<It>::value;
};

template<>
struct is_at_least_input<std::input_iterator_tag>{
    static constexpr bool value = true;
};

#if __cplusplus < 201402L
template<bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;
#else
template<bool B, class T = void>
using enable_if_t = std::enable_if_t<B, T>;
#endif

template<class It>
using is_iterator = enable_if_t<is_at_least_input<typename it_traits<It>::iterator_category>::value, int>;

struct offset_by{
    static size_t off_by(size_t free_blocks) noexcept{
        return free_blocks / 2;
    }
};

template<typename T, class Alloc = std::allocator<T>, class OffsetBy = rdsl::offset_by>
struct devector{
    using value_type = T;
    using allocator_type = Alloc;
    using offset_by_type = OffsetBy;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = typename al_traits<allocator_type>::pointer;
    using const_pointer = typename al_traits<allocator_type>::const_pointer;
    using size_type = typename al_traits<allocator_type>::size_type;
    using difference_type = typename al_traits<allocator_type>::difference_type;
    using iterator = pointer;
    using const_iterator = const_pointer;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
    struct compressed_alloc: public allocator_type{
        compressed_alloc(const allocator_type& alloc = allocator_type())
        :allocator_type(alloc)
        {}

        compressed_alloc(const compressed_alloc&) = delete;
        compressed_alloc(compressed_alloc&&) noexcept = delete;
        compressed_alloc& operator=(const compressed_alloc&) = delete;
        compressed_alloc& operator=(compressed_alloc&&) noexcept = delete;

        allocator_type& get() noexcept{ return *this; }
        const allocator_type& get() const noexcept{ return *this; }

        pointer arr;
    }alloc;

    struct compressed_offs: public offset_by_type{
        compressed_offs(const offset_by_type& offs = offset_by_type())
        :offset_by_type(offs)
        {}

        compressed_offs(const compressed_offs&) = delete;
        compressed_offs(compressed_offs&&) noexcept = delete;
        compressed_offs& operator=(const compressed_offs&) = delete;
        compressed_offs& operator=(compressed_offs&&) noexcept = delete;

        offset_by_type& get() noexcept{ return *this; }
        const offset_by_type& get() const noexcept{ return *this; }

        size_type capacity;
    }offs;

    pointer begin_;
    pointer end_;

    static constexpr float factor = 1.6f; // new_capacity = capacity * factor

    struct buffer_guard{
        pointer begin;
        pointer end;
        allocator_type& alloc;

        buffer_guard(allocator_type& alloc, pointer begin, pointer end)
        :alloc(alloc), begin(begin), end(end) {}

        buffer_guard(allocator_type& alloc, pointer start)
        :alloc(alloc), begin(start), end(start) {}

        buffer_guard(allocator_type& alloc)
        :alloc(alloc), begin(nullptr), end(nullptr) {}

        void guard(pointer begin, pointer end){
            this->begin = begin;
            this->end = end;
        }

        void guard(pointer start){
            begin = end = start;
        }

        void release(){
            begin = end;
        }

        ~buffer_guard(){
            while(begin != end){
                al_traits<allocator_type>::destroy(alloc, begin);
                ++begin;
            }
        }
    };

    struct memory_guard{
        pointer arr;
        size_type capacity;
        allocator_type& alloc;

        memory_guard(allocator_type& alloc, size_type capacity)
        :alloc(alloc), capacity(capacity), arr(alloc.allocate(capacity)) {}

        void release(){
            arr = nullptr;
        }

        ~memory_guard(){
            if(arr){
                alloc.deallocate(arr, capacity);
            }
        }
    };

public:

    size_type capacity() const noexcept{
        return offs.capacity;
    }

    size_type max_size() const{
        return al_traits<allocator_type>::max_size(alloc);
    }

    iterator begin() noexcept{
        return begin_;
    }

    const_iterator begin() const noexcept{
        return begin_;
    }

    iterator end() noexcept{
        return end_;
    }

    const_iterator end() const noexcept{
        return end_;
    }

    reverse_iterator rbegin() noexcept{
        return reverse_iterator(end_);
    }

    const_reverse_iterator rbegin() const noexcept{
        return const_reverse_iterator(end_);
    }

    reverse_iterator rend() noexcept{
        return reverse_iterator(begin_);
    }

    const_reverse_iterator rend() const noexcept{
        return const_reverse_iterator(begin_);
    }

    const_iterator cbegin() const noexcept{
        return begin();
    }

    const_iterator cend() const noexcept{
        return end();
    }

    const_reverse_iterator crbegin() const noexcept{
        return rbegin();
    }

    const_reverse_iterator crend() const noexcept{
        return rend();
    }

    size_type size() const noexcept{
        return end() - begin();
    }

    bool empty() const noexcept{
        return end() == begin();
    }

private:
    size_type free_front() const noexcept{ return begin_ - alloc.arr; }
    size_type free_back() const noexcept{ return alloc.arr + offs.capacity - end_; }
    size_type free_total() const noexcept{ return offs.capacity - size(); }

    pointer allocate_n(size_type n){
        auto ptr = alloc.allocate(n);
        offs.capacity = n;
        return ptr;
    }

    void deallocate() noexcept{
        if(offs.capacity){
            alloc.deallocate(alloc.arr, offs.capacity);
            offs.capacity = 0;
        }
    }

    void construct(size_type n, const_reference val){
        begin_ = end_ = alloc.arr + offs.off_by(offs.capacity - n);
        while(n--){
            al_traits<allocator_type>::construct(alloc, end_, val);
            ++end_;
        }
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    void construct(InputIterator first, size_type distance){
        begin_ = end_ = alloc.arr + offs.off_by(offs.capacity - distance);
        while(distance--){
            al_traits<allocator_type>::construct(alloc, end_, *first);
            ++first;
            ++end_;
        }
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    void construct_move(InputIterator first, size_type distance){
        begin_ = end_ = alloc.arr + offs.off_by(offs.capacity - distance);
        while(distance--){
            al_traits<allocator_type>::construct(alloc, end_, std::move_if_noexcept(*first));
            ++first;
            ++end_;
        }
    }

    void destroy_all() noexcept{
        while(begin_ != end_){
            al_traits<allocator_type>::destroy(alloc, begin_);
            ++begin_;
        }
    }

    void steal_ownership(devector& x) noexcept{
        offs.capacity = x.offs.capacity;
        alloc.arr = x.alloc.arr;
        begin_ = x.begin_;
        end_ = x.end_;
        
        x.alloc.arr = x.begin_ = x.end_ = nullptr;
        x.offs.capacity = 0;
    }
    
    bool in_bounds(pointer it) const noexcept{
        return it >= begin_ && it < end_;
    }

    int next_capacity() const noexcept{
        return factor * offs.capacity + 1;
    }
    /**
     * @brief Allocates a new memory chunk of *new_capacity* capacity and copies all elements into it, respecting the offset factor.
     * *new capacity* should be greater equal to size.
     */
    void reallocate(size_type new_capacity, size_type offset){
        memory_guard mem_guard(alloc, new_capacity);

        buffer_guard buf_guard(alloc, mem_guard.arr + offset);
     
        for(; begin_ != end_; ++begin_, ++buf_guard.end){
            al_traits<allocator_type>::construct(alloc, buf_guard.end, std::move_if_noexcept(*begin_));
        }
        destroy_all();
        deallocate();

        alloc.arr = mem_guard.arr;
        begin_ = buf_guard.begin;
        end_ = buf_guard.end;
        offs.capacity = new_capacity;

        buf_guard.release();
        mem_guard.release();
    }

    void reallocate(size_type new_capacity){
        reallocate(new_capacity, offs.off_by(new_capacity - size()));
    }

    template<class Pred>
    void front_shift_while(pointer& new_begin, Pred pred){
        while(!empty() && pred()){
            al_traits<allocator_type>::construct(alloc, new_begin, std::move(*begin_));
            al_traits<allocator_type>::destroy(alloc, begin_);
            ++begin_;
            ++new_begin;
        }
    }

    template<class Pred>
    void back_shift_while(pointer& new_end, Pred pred){
        while(!empty() && pred()){
            al_traits<allocator_type>::construct(alloc, new_end - 1, std::move(end_[-1]));
            al_traits<allocator_type>::destroy(alloc, end_ - 1);
            --end_;
            --new_end;
        }
    }

    /**
     * @brief Segregates the container into two parts with a gap of 'n'
     * elements starting at 'pos' while also shifting the container to 
     * the desired position according to *offset*.
     * 
     * @return pointer to the first element of the 'n' element gap. 
     */
    pointer segregate(pointer new_begin, pointer new_end, const_iterator pos, size_type n){
        pointer free_space;
        buffer_guard front_guard(alloc);
        buffer_guard back_guard(alloc);

        if(!in_bounds(new_begin)){
            front_guard.guard(new_begin);
            front_shift_while(front_guard.end, [this, pos]{ return begin_ < pos; });

            free_space = front_guard.end;

            if(!in_bounds(front_guard.end + n)){
                back_guard.guard(front_guard.end + n);
                front_shift_while(back_guard.end, []{ return true; });
            }else{
                back_guard.guard(new_end);
                back_shift_while(back_guard.begin, [this, pos]{ return end_ > pos; });
            }
        }else{
            back_guard.guard(new_end);
            back_shift_while(back_guard.begin, [this, pos]{ return end_ > pos; });

            free_space = back_guard.begin - n;

            if(!in_bounds(back_guard.begin - n - 1)){
                front_guard.guard(back_guard.begin - n);
                back_shift_while(front_guard.begin, []{ return true; });
            }else{
                front_guard.guard(new_begin);
                front_shift_while(front_guard.end, [this, pos]{ return begin_ < pos; });
            }
        }

        front_guard.release();
        back_guard.release();

        return free_space;
    }

    /**
     * @brief merges two ranges into one, destroying any elements between them.
     * Also shifts the elements in case [new_begin, new_end) isn't inside [begin, end). 
     * 
     * @return pointer 
     */
    pointer integrate(pointer new_begin, pointer new_end, const_iterator pos, size_type n){
        pointer ret = new_begin + (pos - begin_);
        buffer_guard front_guard(alloc);
        buffer_guard back_guard(alloc);
        
        if(!in_bounds(new_begin)){
            front_guard.guard(new_begin);

            while(begin_ < pos){
                al_traits<allocator_type>::construct(alloc, front_guard.end, *begin_);
                ++front_guard.end;
                pop_front();
            }

            while(n--){
                pop_front();
            }

            while(!empty()){
                al_traits<allocator_type>::construct(alloc, front_guard.end, *begin_);
                ++front_guard.end;
                pop_front();
            }
        }else if(!in_bounds(new_end)){
            back_guard.guard(new_end);

            while(end_ < pos + n){
                al_traits<allocator_type>::construct(alloc, back_guard.begin - 1, end_[-1]);
                --back_guard.begin;
                pop_back();
            }

            while(n--){
                pop_back();
            }

            while(!empty()){
                al_traits<allocator_type>::construct(alloc, back_guard.begin - 1, end_[-1]);
                --back_guard.begin;
                pop_back();
            }
        }else{
            auto const front_shift = new_begin - begin_;
            auto const back_shift = end_ - new_end;

            while(begin_ < pos){
                begin_[front_shift] = std::move(*begin_);
                pop_front();
            }

            while(begin_ < new_begin){
                pop_front();
            }

            while(end_ > pos + n){
                end_[-back_shift - 1] = std::move(end_[-1]);
                pop_back();
            }

            while(end_ > new_end){
                pop_back();
            }
        }

        begin_ = new_begin;
        end_ = new_end;
        front_guard.release();
        back_guard.release();

        return ret;
    }

    size_type capacity_to_fit(size_type n) const noexcept{
        float temp_capacity = offs.capacity ? offs.capacity : 1;
        while(temp_capacity < n){
            temp_capacity = factor * temp_capacity;
        }
        return static_cast<size_type>(temp_capacity);
    }

    template<class Insert>
    iterator insert_impl(const_iterator position, size_type n, Insert ins){
        iterator pos; // position of first newly-created element

        if(n <= free_total()){
            if(position == begin_){
                buffer_guard front_guard(alloc, begin_ - n);
                while(n--){
                    ins(front_guard.end);
                    ++front_guard.end;
                }
                begin_ = front_guard.begin;
                front_guard.release();
            }else if(position == end_){
                while(n--){
                    ins(end_);
                    ++end_;
                }
            }else{
                const pointer new_begin = alloc.arr + offs.off_by(free_total() - n);
                const pointer new_end = new_begin + n + size();

                const pointer free_space = segregate(new_begin, new_end, position, n);

                buffer_guard front_guard(alloc, new_begin, free_space);
                buffer_guard back_guard(alloc, free_space + n, new_end);

                while(n--){
                    ins(front_guard.end);
                    ++front_guard.end;
                }

                begin_ = new_begin;
                end_ = new_end;

                front_guard.release();
                back_guard.release();

                pos = free_space;
            }
        }else{
            const size_type new_size = size() + n;

            memory_guard mem_guard(alloc, capacity_to_fit(new_size));
            
            const size_type front_space = offs.off_by(mem_guard.capacity - new_size);

            buffer_guard buf_guard(alloc, mem_guard.arr + front_space);

            auto it = begin();
            for(; it < position; ++it){
                al_traits<allocator_type>::construct(alloc, buf_guard.end, std::move_if_noexcept(*it));
                ++buf_guard.end;
            }

            pos = buf_guard.end;
            while(n--){
                ins(buf_guard.end);
                ++buf_guard.end;
            }
            for(; it < end(); ++it){
                al_traits<allocator_type>::construct(alloc, buf_guard.end, std::move_if_noexcept(*it));
                ++buf_guard.end;
            }

            destroy_all();
            deallocate();

            alloc.arr = mem_guard.arr;
            offs.capacity = mem_guard.capacity;

            begin_ = buf_guard.begin;
            end_ = buf_guard.end;

            buf_guard.release();
            mem_guard.release();
        }

        return pos;
    }

public:
    
    explicit devector(const allocator_type& allocator = allocator_type(), const offset_by_type& offset_by = offset_by_type())
    :alloc(allocator), offs(offset_by)
    {
        begin_ = end_ = alloc.arr = nullptr;
        offs.capacity = 0;
    }

    explicit devector(const offset_by_type& offset_by)
    :devector(allocator_type(), offset_by)
    {}
    
    devector(
        size_type n,
        const_reference val,
        const allocator_type& allocator = allocator_type(),
        const offset_by_type& offset_by = offset_by_type()
    )
    :alloc(allocator), offs(offset_by)
    {
        alloc.arr = allocate_n(n);
        try{
            construct(n, val);
        }catch(...){
            destroy_all();
            deallocate();
            throw;
        }
    }

    devector(size_type n, const_reference val, const offset_by_type& offset_by)
    :devector(n, val, allocator_type(), offset_by)
    {}

#if __cplusplus < 201402L
    explicit devector(size_type n)
    :devector(n, value_type())
    {}
#else
    explicit devector(
        size_type n,
        const allocator_type& allocator = allocator_type(),
        const offset_by_type& offset_by = offset_by_type()
    )
    :devector(n, value_type(), allocator, offset_by)
    {}
#endif

    template<class InputIterator, is_iterator<InputIterator> = 0>
    devector(
        InputIterator first,
        InputIterator last,
        size_type distance,
        const allocator_type& allocator = allocator_type(),
        const offset_by_type& offset_by = offset_by_type()
    )
    :alloc(allocator), offs(offset_by)
    {
       alloc.arr = allocate_n(distance);
       try{
           construct(first, distance);
       }catch(...){
           destroy_all();
           deallocate();
           throw;
       }
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    devector(
        InputIterator first,
        InputIterator last,
        size_type distance,
        const offset_by_type& offset_by
    )
    :devector(first, last, distance, allocator_type(), offset_by)
    {}

    template<class InputIterator, is_iterator<InputIterator> = 0>
    devector(
        InputIterator first,
        InputIterator last,
        const allocator_type& allocator = allocator_type(),
        const offset_by_type& offset_by = offset_by_type()
    )
    :devector(allocator, offset_by)
    {
        if(is_at_least_forward<typename it_traits<InputIterator>::iterator_category>::value){
            const size_type distance = std::distance(first, last);
            alloc.arr = allocate_n(distance);
            try{
                construct(first, distance);
            }catch(...){
                destroy_all();
                deallocate();
                throw;
            }
        }else{
            while(first != last){
                push_back(*first);
                ++first;
            }
        }
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    devector(
        InputIterator first,
        InputIterator last,
        const offset_by_type& offset_by
    )
    :devector(first, last, allocator_type(), offset_by)
    {}

    devector(const devector& x, const allocator_type& allocator, const offset_by_type& offset_by)
    :devector(x.begin(), x.end(), x.size(), allocator, offset_by)
    {}

    devector(const devector& x, const allocator_type& allocator)
    :devector(x.begin(), x.end(), x.size(), allocator, x.offs)
    {}

    devector(const devector& x, const offset_by_type& offset_by)
    :devector(x.begin(), x.end(), x.size(), x.alloc, offset_by)
    {}

    devector(const devector& x)
    :devector(x, al_traits<allocator_type>::select_on_container_copy_construction(x.alloc), x.offs)
    {}

    devector(devector&& x) noexcept
    :devector(std::move(x), x.alloc, x.offs)
    {}

    devector(devector&& x, const allocator_type& allocator, const offset_by_type& offset_by) noexcept
    :offs(offset_by)
    {
        if(allocator == x.alloc){
            alloc.get() = x.alloc.get();
            steal_ownership(x);
        }else{
            alloc.get() = allocator;

            alloc.arr = allocate_n(x.size());
            construct_move(x.begin(), x.size());
        }
    }

    devector(devector&& x, const allocator_type& allocator) noexcept 
    :devector(std::move(x), allocator, x.offs)
    {}

    devector(devector&& x, const offset_by_type& offset_by) noexcept
    :devector(std::move(x), x.alloc, offset_by)
    {}

    devector(
        std::initializer_list<value_type> il,
        const allocator_type& allocator = allocator_type(),
        const offset_by_type& offset_by = offset_by_type()
    )
    :devector(il.begin(), il.end(), il.size(), allocator, offset_by)
    {}

    devector(std::initializer_list<value_type> il, const offset_by_type& offset_by)
    :devector(il.begin(), il.end(), il.size(), allocator_type(), offset_by)
    {}

    ~devector(){
        destroy_all();
        deallocate();
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    void assign (InputIterator first, size_type distance){
        destroy_all();
        if(offs.capacity < distance){
            reallocate(capacity_to_fit(distance));
        }
        construct(first, distance);
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    void assign (InputIterator first, InputIterator last){
        if(is_at_least_forward<InputIterator>::value){
            assign(first, std::distance(first,last));
        }else{
            destroy_all();
            begin_ = end_ = alloc.arr + offs.off_by(free_total());
            while(first != last){
                push_back(*first);
                ++first;
            }
        }
    }

    void assign(size_type n, const_reference val){
        destroy_all();
        if(offs.capacity < n){
            reallocate(capacity_to_fit(n));
        }
        construct(n, val);
    }

    void assign(std::initializer_list<value_type> il){
        destroy_all();
        if(offs.capacity < il.size()){
            reallocate(capacity_to_fit(il.size()));
        }
        construct(il.begin(), il.size());
    }

    devector& operator=(const devector& x){
        if(this == &x){
            return *this;
        }

        offs.get() = x.offs.get();

        if(al_traits<allocator_type>::propagate_on_container_copy_assignment::value && alloc != x.alloc){            
            destroy_all();
            deallocate();
            
            alloc.get() = x.alloc.get();

            alloc.arr = allocate_n(capacity_to_fit(x.size()));
            construct(x.begin(), x.size());

        }else{
            if(offs.capacity < x.size()){
                destroy_all();
                deallocate();
                alloc.arr = allocate_n(capacity_to_fit(x.size()));
                construct(x.begin(), x.size());
            }else{
                const pointer new_begin = alloc.arr + offs.off_by(offs.capacity - x.size());
                const pointer new_end = new_begin + x.size();

                while(!empty() && begin_ < new_begin){
                    pop_front();
                }
                
                while(!empty() && end_ > new_end){
                    pop_back();
                }

                buffer_guard guard(alloc, new_begin);

                for(auto it = x.begin_; it != x.end_; ++it, ++guard.end){
                    if(in_bounds(guard.end)){
                        *guard.end = *it;
                    }else{
                        al_traits<allocator_type>::construct(alloc, guard.end, *it);
                    }
                }

                begin_ = new_begin;
                end_ = new_end;
                guard.release();
            }
        }


        return *this;
    }

    devector& operator=(devector&& x){
        if(this == &x){
            return *this;
        }

        offs.get() = std::move(x.offs.get());

        if(!al_traits<allocator_type>::propagate_on_container_move_assignment::value){
            if(alloc != x.alloc){
                if(offs.capacity < x.size()){
                    destroy_all();
                    deallocate();
                    alloc.arr = allocate_n(capacity_to_fit(x.size()));
                    construct_move(x.begin_, x.size());
                }else{
                    const pointer new_begin = alloc.arr + offs.off_by(offs.capacity - x.size());
                    const pointer new_end = new_begin + x.size();

                    while(!empty() && begin_ < new_begin){
                        pop_front();
                    }
                    
                    while(!empty() && end_ > new_end){
                        pop_back();
                    }

                    buffer_guard guard(alloc, new_begin);

                    for(auto it = x.begin_; it != x.end_; ++it, ++guard.end){
                        if(in_bounds(guard.end)){
                            *guard.end = std::move(*it);
                        }else{
                            al_traits<allocator_type>::construct(alloc, guard.end, std::move(*it));
                        }
                    }

                    begin_ = new_begin;
                    end_ = new_end;
                    guard.release();
                }
            }else{
                destroy_all();
                deallocate();
                steal_ownership(x);
            } 
        }else{
            destroy_all();
            deallocate();
            steal_ownership(x);
            alloc.get() = std::move(x.alloc.get());
        }
        
        return *this;
    }

    devector& operator=(std::initializer_list<value_type> il){
        if(offs.capacity < il.size()){
            destroy_all();
            deallocate();
            alloc.arr = allocate_n(capacity_to_fit(il.size()));
            construct(il.begin(), il.size());
        }else{
            const pointer new_begin = alloc.arr + offs.off_by(offs.capacity - il.size());
            const pointer new_end = new_begin + il.size();

            while(!empty() && begin_ < new_begin){
                pop_front();
            }
            
            while(!empty() && end_ > new_end){
                pop_back();
            }

            buffer_guard guard(alloc, new_begin);

            for(auto it = il.begin(); it != il.end(); ++it, ++guard.end){
                if(in_bounds(guard.end)){
                    *guard.end = *it;
                }else{
                    al_traits<allocator_type>::construct(alloc, guard.end, *it);
                }
            }

            begin_ = new_begin;
            end_ = new_begin + il.size();
            guard.release();
        }

        return *this;
    }

    void resize_front(size_type n, const_reference val = value_type()){
        while(n < size()){
            pop_front();
        }
        if(n > offs.capacity){
            reallocate(capacity_to_fit(n), free_back());
        }
        while(n > size()){
            push_front(val);
        }
    }

    void resize_back(size_type n, const_reference val = value_type()){
        while(n < size()){
            pop_back();
        }
        if(n > offs.capacity){
            reallocate(capacity_to_fit(n), free_front());
        }
        while(n > size()){
            push_back(val);
        }
    }

    void resize(size_type n, const_reference val = value_type()){
        resize_back(n, val);
    }
  
    void reserve(size_type n){
        if(n > offs.capacity){
            reallocate(n);
        }
    }

    void shrink_to_fit(){
        reallocate(size());
    }

    reference operator[](size_type index){
        return begin_[index];
    }

    const_reference operator[](size_type index) const{
        return begin_[index];
    }

    reference at(size_type index){
        if(in_bounds(begin_ + index)){
            return begin_[index];
        }else{
            throw std::out_of_range("index " + std::to_string(index) + " out of range for array of size " + std::to_string(size()));
        }
    }

    const_reference at(size_type index) const{
        if(in_bounds(begin_ + index)){
            return begin_[index];
        }else{
            throw std::out_of_range("index " + std::to_string(index) + " out of range for array of size " + std::to_string(size()));
        }
    }

    reference front(){
        return *begin_;
    }

    const_reference front() const{
        return *begin_;
    }

    reference back(){
        return *(end_ - 1);
    }

    const_reference back() const{
        return *(end_ - 1);
    }

    pointer data() noexcept{
        return alloc.arr;
    }

    const_pointer data() const noexcept{
        return alloc.arr;
    }

    void push_back(const_reference val){
        if(!free_back()){
            const auto new_capacity = next_capacity();
            auto offset = offs.off_by(new_capacity - size());
            offset -= offset == (new_capacity - size());
            reallocate(new_capacity, offset);
        }

        al_traits<allocator_type>::construct(alloc, end_, val);
        ++end_;
    }

    void push_back(value_type&& val){
        if(!free_back()){
            const auto new_capacity = next_capacity();
            auto offset = offs.off_by(new_capacity - size());
            offset -= offset == (new_capacity - size());
            reallocate(new_capacity, offset);
        }

        al_traits<allocator_type>::construct(alloc, end_, std::move(val));
        ++end_;
    }

    void push_front(const_reference val){
        if(!free_front()){
            const auto new_capacity = next_capacity();
            const auto offset = offs.off_by(new_capacity - size());
            reallocate(new_capacity, offset + !offset);
        }

        al_traits<allocator_type>::construct(alloc, begin_ - 1, val);
        --begin_;
    }

    void push_front(value_type&& val){
        if(!free_front()){
            const auto new_capacity = next_capacity();
            const auto offset = offs.off_by(new_capacity - size());
            reallocate(new_capacity, offset + !offset);
        }

        al_traits<allocator_type>::construct(alloc, begin_ - 1, std::move(val));
        --begin_;
    }

    void pop_back() noexcept{
        al_traits<allocator_type>::destroy(alloc, end_ - 1);
        --end_;
    }

    void pop_front() noexcept{
        al_traits<allocator_type>::destroy(alloc, begin_);
        ++begin_;
    }

    iterator insert(const_iterator position, size_type n, const_reference val){
        return insert_impl(position, n, [&val, this](pointer p){
            al_traits<allocator_type>::construct(alloc, p, val);
        });
    }

    iterator insert(const_iterator position, const_reference val){
        return insert(position, 1, val);
    }

    iterator insert(const_iterator position, value_type&& val){
        return insert_impl(position, 1, [&val, this](pointer p) mutable{
            al_traits<allocator_type>::construct(alloc, p, std::move(val));
        });
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    iterator insert(const_iterator position, InputIterator first, size_type n){
        return insert_impl(position, n, [first, this](pointer p) mutable{
            al_traits<allocator_type>::construct(alloc, p, *first++);
        });
    }

    template<class InputIterator, is_iterator<InputIterator> = 0>
    iterator insert(const_iterator position, InputIterator first, InputIterator last){
        if(is_at_least_forward<typename it_traits<InputIterator>::iterator_category>::value){
            return insert(position, first, std::distance(first, last));
        }else{
            const auto index = position - cbegin();
            memory_guard mem_guard(alloc, cend() - position + 1);
            buffer_guard buf_guard(alloc, mem_guard.arr + mem_guard.capacity - 1);

            while(buf_guard.begin != mem_guard.arr){
                al_traits<allocator_type>::construct(alloc, buf_guard.begin - 1, std::move(back()));
                --buf_guard.begin;
                pop_back();
            }

            while(first != last){
                push_back(*first);
                ++first;
            }

            reserve(capacity() + buf_guard.end - buf_guard.begin);
            for(pointer it = buf_guard.begin; it != buf_guard.end; ++it){
                al_traits<allocator_type>::construct(alloc, end_, std::move(*it));
                ++end_;
            }
            
            return begin_ + index;
        }
    }

    iterator insert(const_iterator position, std::initializer_list<value_type> il){
        return insert(position, il.begin(), il.size());
    }

    iterator erase(const_iterator first, const_iterator last){
        iterator pos;

        if(first == begin_){
            while(begin_ < last){
                al_traits<allocator_type>::destroy(alloc, begin_);
                ++begin_;
            }
            return begin_;
        }else if(last == end_){
            while(end_ > first){
                al_traits<allocator_type>::destroy(alloc, end_ - 1);
                --end_;
            }
            return end_;
        }else{
            const size_type n = last - first;

            const pointer new_begin = alloc.arr + offs.off_by(free_total() + n);
            const pointer new_end = new_begin + size() - n;

            return integrate(new_begin, new_end, first, n);
        }
    }

    iterator erase(const_iterator position){
        return erase(position, position + 1);
    }

    void swap(devector& x){
        std::swap(alloc.arr, x.alloc.arr);
        std::swap(begin_, x.begin_);
        std::swap(end_, x.end_);
        std::swap(offs.capacity, x.offs.capacity);
        std::swap(offs, x.offs);
        if(al_traits<allocator_type>::propagate_on_container_swap){
            std::swap(alloc, x.alloc);
        }
    }

    void clear() noexcept{
        destroy_all();
    }

    template<class... Args>
    iterator emplace(const_iterator position, Args&&... args){
        return insert_impl(position, 1, [&](pointer p){
            al_traits<allocator_type>::construct(alloc, p, std::forward<Args>(args)...);
        });
    }

    template<class... Args>
    iterator emplace_back(Args&&... args){
        if(!free_back()){
            const auto new_capacity = next_capacity();
            auto offset = offs.off_by(new_capacity - size());
            offset -= offset == (new_capacity - size());
            reallocate(new_capacity, offset);
        }

        al_traits<allocator_type>::construct(alloc, end_, std::forward<Args>(args)...);
        return end_++;
    }

    template<class... Args>
    iterator emplace_front(Args&&... args){
        if(!free_front()){
            const auto new_capacity = next_capacity();
            const auto offset = offs.off_by(new_capacity - size());
            reallocate(new_capacity, offset + !offset);
        }

        al_traits<allocator_type>::construct(alloc, begin_ - 1, std::forward<Args>(args)...);
        return begin_--;
    }

    allocator_type get_allocator() const noexcept{
        return alloc;
    }

    offset_by_type get_offset_by() const noexcept{
        return offs;
    }
};

template<class T, class Alloc, class OffsetByA, class OffsetByB>
bool operator== (const devector<T, Alloc, OffsetByA>& lhs, const devector<T, Alloc, OffsetByB>& rhs){
    if(lhs.size() != rhs.size()){
        return false;
    }

    for(auto it0 = lhs.begin(), it1 = rhs.begin(); it0 != lhs.end(); ++it0, ++it1){
        if(*it0 != *it1){
            return false;
        }
    }
    
    return true;
}

template<class T, class Alloc, class OffsetByA, class OffsetByB>
bool operator!= (const devector<T, Alloc, OffsetByA>& lhs, const devector<T, Alloc, OffsetByB>& rhs){
    return !(lhs == rhs);
}

template<class T, class Alloc, class OffsetByA, class OffsetByB>
bool operator< (const devector<T, Alloc, OffsetByA>& lhs, const devector<T, Alloc, OffsetByB>& rhs){
    auto it1 = rhs.cbegin();
    for(auto it0 = lhs.cbegin(); it0 != lhs.cend(); ++it0, ++it1){
        if(it1 == rhs.cend() || *it1 < *it0){
            return false;
        }else if(*it0 < *it1){
            return true;
        }
    }

    return it1 != rhs.cend();
}

template<class T, class Alloc, class OffsetByA, class OffsetByB>
bool operator<= (const devector<T, Alloc, OffsetByA>& lhs, const devector<T, Alloc, OffsetByB>& rhs){
    return !(rhs < lhs);
}

template<class T, class Alloc, class OffsetByA, class OffsetByB>
bool operator> (const devector<T, Alloc, OffsetByA>& lhs, const devector<T, Alloc, OffsetByB>& rhs){
    return rhs < lhs;
}

template<class T, class Alloc, class OffsetByA, class OffsetByB>
bool operator>= (const devector<T, Alloc, OffsetByA>& lhs, const devector<T, Alloc, OffsetByB>& rhs){
    return !(lhs < rhs);
}

template<class T, class Alloc, class OffsetByA, class OffsetByB>
void swap(devector<T, Alloc, OffsetByA>& x, devector<T, Alloc, OffsetByB> y){
    x.swap(y);
}
} //rdsl

#endif
