#pragma once

#include "utils/logger/logger.hpp"

#include <memory>
#include <unordered_map>
#include <vector>
#include <queue>

namespace pathtree {

template <class T>
class ObjectCounter {
 public:
  static size_t object_count_max() { return count_max_; }

  static size_t object_count_current() { return count_current_; }

  template <typename... Args>
  ObjectCounter(Args &&...) noexcept {
    ++count_current_;
    count_max_ = std::max(count_max_, count_current_);
  }

  ~ObjectCounter() { --count_current_; }

 private:
  static size_t count_max_;
  static size_t count_current_;
};

template <class T>
size_t ObjectCounter<T>::count_max_ = 0;

template <class T>
size_t ObjectCounter<T>::count_current_ = 0;

template <typename T>
class Node : public ObjectCounter<Node<T>> {
  using This = Node<T>;

 public:
  std::vector<T> collect() const {
    std::vector<T> result;
    const This *p = this;
    while (p) {
      result.push_back(p->payload_);
      p = p->parent_.get();
    }

    return result;
  }

  Node(const T &payload, const std::shared_ptr<This> &parent = nullptr) : payload_{payload}, parent_{parent} {}

  static std::shared_ptr<This> child(const T &payload, const std::shared_ptr<This> &parent = nullptr) {
    return std::make_shared<This>(payload, parent);
  }

  const auto &payload() const { return payload_; }

 private:
  T payload_;
  std::shared_ptr<This> parent_;
};

template <typename T>
std::shared_ptr<Node<T>> make_child(const T &payload, const std::shared_ptr<Node<T>> &parent = nullptr) {
  return Node<T>::child(payload, parent);
}

template <typename GraphCursor>
class PathLink {
  using This = PathLink<GraphCursor>;

 public:
  double score() const {
    if (scores_.size() == 0) {
      return std::numeric_limits<double>::infinity();
    }

    return best_ancestor()->second.first;
  }

  auto best_ancestor() const {
    return std::min_element(scores_.cbegin(), scores_.cend(),
                            [](const auto &e1, const auto &e2) { return e1.second.first < e2.second.first; });
  }

  bool update(GraphCursor gp, double score, const std::shared_ptr<This> &pl) {
    auto val = std::make_pair(score, pl);
    auto it_fl = scores_.insert({gp, val});
    bool inserted = it_fl.second;
    if (!inserted) {
      const auto &it = it_fl.first;
      if (it->second.first > score) {
        it->second = std::move(val);
        return true;
      } else {
        return false;
      }
    } else {
      return true;
    }
  }

  void merge_update(const This *other, double add_fee = 0) {
    for (const auto &kv : other->scores_) {
      update(kv.first, kv.second.first + add_fee, kv.second.second);
    }
  }

  std::shared_ptr<This> merge(const This *other, double add_fee = 0) const {
    auto result = std::make_shared<This>(*this);
    result.merge_update(other, add_fee);
    return result;
  }

  std::shared_ptr<This> merge(GraphCursor gp, double score, const std::shared_ptr<This> &pl) const {
    auto result = std::make_shared<This>(*this);
    result.update(gp, score, pl);
    return result;
  }

  static std::shared_ptr<This> master_source() {
    auto result = std::make_shared<This>();
    result->scores_[GraphCursor()] = std::make_pair(0, nullptr);  // master_sourse score should be 0
    return result;
  }

  struct Comp {
    template <typename T>
    bool operator()(const T &e1, const T &e2) const {
      return std::get<double>(e1->payload()) > std::get<double>(e2->payload());
    }
  };

  std::vector<std::pair<std::vector<GraphCursor>, double>> top_k(size_t k) {
    using Node = Node<std::tuple<GraphCursor, double, This *>>;
    using SPT = std::shared_ptr<Node>;
    std::priority_queue<SPT, std::vector<SPT>, Comp> q;
    auto extract_path = [&q]() {
      SPT tail = q.top();
      q.pop();
      GraphCursor gp;
      This *p;
      double cost;
      std::tie(gp, cost, p) = tail->payload();

      std::vector<GraphCursor> path;

      for (const auto &tpl : tail->collect()) {
        path.push_back(std::get<GraphCursor>(tpl));
      }

      std::reverse(path.begin(), path.end());

      while (p) {
        auto new_tail = tail;
        p->clean_left_link_();
        auto best = p->best_ancestor();
        for (auto it = p->scores_.cbegin(); it != p->scores_.cend(); ++it) {
          double delta = it->second.first - best->second.first;
          auto spt = make_child(std::make_tuple(it->first, cost + delta, it->second.second.get()), tail);
          if (it != best) {
            q.push(spt);
          } else {
            new_tail = spt;
          }
        }
        path.push_back(best->first);
        p = best->second.second.get();
        tail = new_tail;
      }

      std::reverse(path.begin(), path.end());

      return std::make_pair(path, cost);
    };

    std::vector<std::pair<std::vector<GraphCursor>, double>> result;
    auto best = best_ancestor();
    if (best == scores_.end()) {
      return result;
      // TODO Support empty Link as a comon case and remove this workaround
    }

    double best_score = best->second.first;
    auto initial = make_child(std::make_tuple(GraphCursor(), best_score, this));
    q.push(initial);

    for (size_t i = 0; i < k && !q.empty(); ++i) {
      result.push_back(extract_path());
      TRACE((i + 1) << " top paths extracted");
      TRACE(Node::object_count_current() << " current # of best path tree nodes");
      TRACE(Node::object_count_max() << " max # of best path tree nodes");
    }

    return result;
  }

  size_t clean_left_link_() {
    // the idea is in the following:
    // if a path is a prefix or suffix of another one,
    // we do not consider them as DIFFERENT paths. We leave only one (the best) of them
    // In order to implement this we do not allow to link to be "left-terminal" (contain backref to "empty cursor") and
    // "transit" (contain backrefs to other nontrivial links).
    // Here we just remove ref to empty cursor or refs to nontrivial links dependent on what is better

    // If terminal ref is present:
    // 1) Remove all other refs worse than it
    // 2) In case of some non-terminal refs still present, remove terminal ref as well
    auto it_terminal = scores_.find(GraphCursor());
    if (it_terminal == scores_.end()) {
      return 0;
    }

    double terminal_score = it_terminal->second.first;
    std::vector<GraphCursor> to_remove;
    for (const auto &kv : scores_) {
      if (kv.first.is_empty()) {
        continue;
      }
      double score = kv.second.first;
      if (score >= terminal_score) {  // We prefer short paths
        to_remove.push_back(kv.first);
      } else {
        to_remove.push_back(GraphCursor());
      }
    }

    size_t count = 0;
    for (const auto &key : to_remove) {
      count += scores_.erase(key);
    }

    return count;
  }

  // clean_right_link
  // 1) У нас есть набор вершин, которые могут быть концами (правыми) путей.
  // 2) Изначально в этом наборе есть одна вершина --- с наименьшим штрафом (т.е. правый конец наилучшего пути)
  // 3) Вершины, которые принадлежат уже построеным оптимальным путям не могут быть началами путей,
  //    потому что для каждого пути, начинающегося из них, есть путь, который не хуже
  //    ТАК, тут нужен более мягкий критерий, выкидывать не вершины, а линки выравниваний.
  //    Если линк выравнивания пройден ранее, то не надо пытаться с него начать путь.
  // 4) Вершины, с которых начинаются пути, не могу лежать ни на каких другх путях
  //    Это тоже жестковато, представим себе скользящее окно даже на прямой.
  //    При таком агрессивном подходе мы найдем только
  //    непересекающиеся окна, а не все.
  //    Надо думать...

  std::vector<std::pair<std::string, double>> top_k_string(size_t k) {
    auto path2string = [](const auto &path) {
      std::string s;
      for (size_t i = 2; i < path.size() - 1; ++i) {
        s += path[i].letter();
      }
      return s;
    };

    std::vector<std::pair<std::string, double>> result;
    for (const auto &path_score : top_k(k)) {
      result.push_back({path2string(path_score.first), path_score.second});
    }

    return result;
  }

  std::vector<GraphCursor> best_path() const {
    std::vector<GraphCursor> path;

    const This *p = this;
    while (p) {
      auto best = p->best_ancestor();
      if (best == scores_.end()) {
        return {GraphCursor(), GraphCursor()};
        // TODO Support empty Link as a comon case and remove this workaround
      }
      path.push_back(best->first);
      p = best->second.second.get();
    }

    std::reverse(path.begin(), path.end());

    return path;
  }

  std::string best_path_string() const {
    auto path = best_path();
    std::string s;
    for (size_t i = 2; i < path.size(); ++i) {
      s += path[i].letter();
    }
    return s;
  }

  std::shared_ptr<This> clone() const { return std::make_shared<This>(*this); }

 private:
  std::unordered_map<GraphCursor, std::pair<double, std::shared_ptr<This>>> scores_;
};

}  // namespace pathtree

template<class GraphCursor>
class PathSet {
 public:
  class path_container {
   public:
    path_container(pathtree::PathLink<GraphCursor> &paths,
                   size_t k)
        : paths_(paths.top_k(k)) {}

    auto begin()  const { return paths_.begin(); }
    auto end()    const { return paths_.end();   }
    size_t size() const { return paths_.size();  }
    auto operator[](size_t n) const { return paths_[n]; }

    std::string str(const std::vector<GraphCursor> &path) const {
      std::string s;
      for (size_t i = 0; i < path.size(); ++i) {
        if (path[i].is_empty())
          continue;

        s += path[i].letter();
      }
      return s;
    }
    std::string str(size_t n) const { return str(paths_[n].first); }

 private:
    std::vector<std::pair<std::vector<GraphCursor>, double>> paths_;
  };

  pathtree::PathLink<GraphCursor> &pathlink() { return pathlink_; }
  const pathtree::PathLink<GraphCursor> &pathlink() const { return pathlink_; }

  double best_score() const { return pathlink_.score(); }
  std::vector<GraphCursor> best_path() const { return pathlink_.best_path(); }
  std::string best_path_string() const { return pathlink_.best_path_string(); }

  path_container top_k(size_t k) {
    return path_container(pathlink_, k);
  }

 private:
  pathtree::PathLink<GraphCursor> pathlink_;
};
