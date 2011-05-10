#ifndef FILTER_ITERATOR_H_
#define FILTER_ITERATOR_H_

/**
 * Iterator with some predicate -- iterates only on elements with predicate(item) == true
 */
template<typename iterator_type, typename predicate_type>
class filter_iterator {
public:
	typedef typename iterator_type::value_type value_type;

	filter_iterator(const iterator_type& begin, const iterator_type& end, const predicate_type& pred):
		current_(begin), end_(end), pred_(pred)
	{
		while((current_ != end_) && (!pred_(*current_)))
			++current_;
	} // filter_iterator

	value_type operator*() const { return *current_; }
	value_type operator->() const { return *current_; }

	filter_iterator& operator++() { advance(); return *this; }

	bool operator==(const filter_iterator& rhs) const { return current_ == rhs.current_; }
	bool operator!=(const filter_iterator& rhs) const { return !(operator==(rhs)); }

private:
	void advance()
	{
		do
		{
			++current_;
		}
		while((current_ != end_) && (!pred_(*current_)));
	} // advance

	iterator_type current_;
	iterator_type end_;
	predicate_type pred_;
};

#endif /* FILTER_ITERATOR_H_ */
