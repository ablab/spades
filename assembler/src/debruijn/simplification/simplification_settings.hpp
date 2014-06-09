#pragma once

namespace debruijn {

namespace simplification {

class LengthThresholdFinder {
public:
    static size_t MaxTipLength(size_t read_length, size_t k, double coeff) {
        return std::max((size_t) math::round((double)std::min(k, read_length / 2) * coeff),
                        read_length);
    }

    static size_t MaxBulgeLength(size_t k, double coeff,
                                 size_t additive_coeff) {
        return std::max((size_t) math::round((double)k * coeff), k + additive_coeff);
    }

    static size_t MaxErroneousConnectionLength(size_t k, size_t param) {
        return k + param;
    }

    static size_t MaxTipOriginatedECLength(size_t read_length, size_t k,
                                           double coeff) {
        return 2 * MaxTipLength(read_length, k, coeff) - 1;
    }
};

//todo use GenomicInfo as field!
class SimplifInfoContainer {
    size_t read_length_;
    double detected_mean_coverage_;
    double detected_coverage_bound_;
    size_t iteration_count_;
    size_t iteration_;
    size_t chunk_cnt_;

public: 
    SimplifInfoContainer() : 
        read_length_(-1u),
        detected_mean_coverage_(-1.0),
        detected_coverage_bound_(-1.0),
        iteration_count_(-1u),
        iteration_(-1u),
        chunk_cnt_(-1u) {
    }

    size_t read_length() const {
        VERIFY(read_length_ != -1u);
        return read_length_;
    }

    double detected_mean_coverage() const {
        VERIFY(math::ge(detected_mean_coverage_, 0.));
        return detected_mean_coverage_;
    }

    double detected_coverage_bound() const {
        VERIFY(math::ge(detected_coverage_bound_, 0.));
        return detected_coverage_bound_;
    }

    size_t iteration_count() const {
        VERIFY(iteration_count_ != -1u);
        return iteration_count_;
    }

    size_t iteration() const {
        VERIFY(iteration_ != -1u);
        return iteration_;
    }
    
    size_t chunk_cnt() const {
        VERIFY(chunk_cnt_ != -1u);
        return chunk_cnt_;
    }

    SimplifInfoContainer& set_read_length(size_t read_length) {
        read_length_ = read_length;
        return *this;
    }

    SimplifInfoContainer& set_detected_coverage_bound(double detected_coverage_bound) {
        detected_coverage_bound_ = detected_coverage_bound;
        return *this;
    }

    SimplifInfoContainer& set_detected_mean_coverage(double detected_mean_coverage) {
    	detected_mean_coverage_ = detected_mean_coverage;
        return *this;
    }

    SimplifInfoContainer& set_iteration_count(size_t iteration_count) {
        iteration_count_ = iteration_count;
        return *this;
    }

    SimplifInfoContainer& set_iteration(size_t iteration) {
        iteration_ = iteration;
        return *this;
    }

    SimplifInfoContainer& set_chunk_cnt(size_t chunk_cnt) {
        chunk_cnt_ = chunk_cnt;
        return *this;
    }
};

template<class Graph>
class ConditionParser {
private:
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    string next_token_;
    string input_;
    const SimplifInfoContainer settings_;
    queue<string> tokenized_input_;

    size_t max_length_bound_;
    double max_coverage_bound_;

    string ReadNext() {
        if (!tokenized_input_.empty()) {
            next_token_ = tokenized_input_.front();
            tokenized_input_.pop();
        } else {
            next_token_ = "";
        }
        return next_token_;
    }

    template<typename T>
    bool RelaxMax(T& cur_max, T t) {
        if (t > cur_max) {
            cur_max = t;
            return true;
        }
        return false;
    }

    template<typename T>
    bool RelaxMin(T& cur_min, T t) {
        if (t < cur_min) {
            cur_min = t;
            return true;
        }
        return false;
    }

    double GetCoverageBound() {
        if (next_token_ == "auto") {
            return settings_.detected_coverage_bound();
        } else {
            return lexical_cast<double>(next_token_);
        }
    }

    shared_ptr<Predicate<EdgeId>> ParseCondition(size_t& min_length_bound,
                                                 double& min_coverage_bound) {
        if (next_token_ == "tc_lb") {
            double length_coeff = lexical_cast<double>(ReadNext());

            DEBUG("Creating tip length bound. Coeff " << length_coeff);
            size_t length_bound = LengthThresholdFinder::MaxTipLength(
                settings_.read_length(), g_.k(), length_coeff);

            DEBUG("Length bound" << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "to_ec_lb") {
            double length_coeff = lexical_cast<double>(ReadNext());

            DEBUG( "Creating length bound for erroneous connections originated from tip merging. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxTipOriginatedECLength(
                        settings_.read_length(), g_.k(), length_coeff);

            DEBUG("Length bound" << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "ec_lb") {
            size_t length_coeff = lexical_cast<size_t>(ReadNext());

            DEBUG("Creating ec length bound. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxErroneousConnectionLength(
                        g_.k(), length_coeff);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "lb") {
            size_t length_bound = lexical_cast<size_t>(ReadNext());

            DEBUG("Creating length bound. Value " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "cb") {
            ReadNext();
            double cov_bound = GetCoverageBound();
            DEBUG("Creating coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return make_shared<CoverageUpperBound<Graph>>(g_, cov_bound);
        } else if (next_token_ == "icb") {
            ReadNext();
            double cov_bound = GetCoverageBound();
            cov_bound = cov_bound / (double) settings_.iteration_count() * (double) (settings_.iteration() + 1);
            DEBUG("Creating iterative coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return make_shared<CoverageUpperBound<Graph>>(g_, cov_bound);
        } else if (next_token_ == "rctc") {
            ReadNext();
            DEBUG("Creating relative cov tip cond " << next_token_);
            return make_shared<RelativeCoverageTipCondition<Graph>>(
                g_, lexical_cast<double>(next_token_));
        } else if (next_token_ == "disabled") {
            DEBUG("Creating disabling condition");
            return make_shared<func::AlwaysFalse<EdgeId>>();
        } else {
            VERIFY(false);
            return make_shared<func::AlwaysTrue<EdgeId>>();
        }
    }

    shared_ptr<Predicate<EdgeId>> ParseConjunction(size_t& min_length_bound,
                                                   double& min_coverage_bound) {
        shared_ptr<Predicate<EdgeId>> answer =
                make_shared<AlwaysTrue<EdgeId>>();
        VERIFY(next_token_ == "{");
        ReadNext();
        while (next_token_ != "}") {
            answer = make_shared<AndOperator<EdgeId>>(
                answer,
                ParseCondition(min_length_bound, min_coverage_bound));
            ReadNext();
        }
        return answer;
    }

public:

    ConditionParser(const Graph& g, string input, const SimplifInfoContainer& settings)
            : g_(g),
              input_(input),
              settings_(settings),
              max_length_bound_(0),
              max_coverage_bound_(0.) {
        DEBUG("Creating parser for string " << input);
        using namespace boost;
        vector<string> tmp_tokenized_input;
        split(tmp_tokenized_input, input_, is_any_of(" ,;"), token_compress_on);
        for (auto it = tmp_tokenized_input.begin();
             it != tmp_tokenized_input.end(); ++it) {
            tokenized_input_.push(*it);
        }
        ReadNext();
    }

    shared_ptr<Predicate<EdgeId>> operator()() {
        DEBUG("Parsing");
        shared_ptr<Predicate<EdgeId>> answer = make_shared<NotOperator<EdgeId>>(
            make_shared<AlwaysTrue<EdgeId>>());
        VERIFY_MSG(next_token_ == "{", "Expected \"{\", but next token was " << next_token_);
        while (next_token_ == "{") {
            size_t min_length_bound = numeric_limits<size_t>::max();
            double min_coverage_bound = numeric_limits<double>::max();
            answer = make_shared<OrOperator<EdgeId>>(
                answer,
                ParseConjunction(min_length_bound, min_coverage_bound));
            RelaxMax(max_length_bound_, min_length_bound);
            RelaxMax(max_coverage_bound_, min_coverage_bound);
            ReadNext();
        }
        return answer;
    }

    size_t max_length_bound() const {
        return max_length_bound_;
    }

    double max_coverage_bound() const {
        return max_coverage_bound_;
    }

private:
    DECL_LOGGER("ConditionParser");
};

}

}
