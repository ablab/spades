//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * devisible_tree.hpp
 *
 *  Created on: Aug 23, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef DEVISIBLE_TREE_HPP_
#define DEVISIBLE_TREE_HPP_


#include <boost/foreach.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>

#include "standard_base.hpp"


namespace omnigraph {

template <class Value>
class TreeNode {

	struct NodeComparator {
		bool operator() (TreeNode* node1, TreeNode* node2) {
			return node1->GetSize() < node2->GetSize();
		}
	};

public:
	TreeNode() : parent_(0), subtree_size_(0) { }

	void AddChild(TreeNode& node) {
		children_.push_back(&node);
		subtree_size_ += node.GetSize();
	}

	virtual size_t GetSize() const {
		return subtree_size_;
	}

	void sort() {
		std::sort(children_.begin(), children_.end(), NodeComparator());
		BOOST_FOREACH(TreeNode* node, children_) {
			node->sort();
		}
	}
//
//	virtual size_t SeparateNodes(vector<Value>& output, const size_t howMuchToSeparate) {
//		size_t separated = 0;
//
//		while (subtree_size_ > 0 && howMuchToSeparate - separated > 0) {
//			TreeNode * node = children_.back();
//
//			if (node->GetSize() <= howMuchToSeparate - separated) {
//				children_.pop_back();
//			}
//
//			// get descendants
//			size_t result = node->SeparateNodes(output, howMuchToSeparate - separated);
//			separated += result;
//			subtree_size_ -= result;
//		}
//
//		return separated;
//	}

	TreeNode* GetSubtreeWithSize(const size_t size) {
		VERIFY(size <= this->GetSize());

		if (children_.size() == 0) {
			return this;
		}

		TreeNode* biggest_child = children_.back();
		if (biggest_child->GetSize() < size) {
			return this;
		}

		TreeNode* node = biggest_child->GetSubtreeWithSize(size);
		if (node == biggest_child) {
			children_.pop_back();
		}

		subtree_size_ -= node->GetSize();
		return node;
	}

	virtual void CollectValues(vector<Value>& output) {
		BOOST_FOREACH(TreeNode* node, children_) {
			node->CollectValues(output);
		}
		subtree_size_ = 0;
		children_.clear();
	}

	TreeNode* GetRoot() {
		if (parent_ == 0) {
			return this;
		} else {
			return parent_->GetRoot();
		}
	}

	void SetParent(TreeNode& parent) {
		VERIFY(parent_ == 0);
		parent_ = &parent;
	}

	virtual ~TreeNode() { }

private:
	TreeNode * parent_;
	vector<TreeNode *> children_;
	size_t subtree_size_;
};


template <class Value>
class TreeNodeWithValue : public TreeNode<Value> {

public:
	TreeNodeWithValue(const Value& value) : value_(value) {	}

	virtual ~TreeNodeWithValue() { }
//
//	virtual size_t SeparateNodes(vector<Value>& output, const size_t howMuchToSeparate) {
//		size_t separated = TreeNode<Value>::SeparateNodes(output, howMuchToSeparate);
//
//		if (separated < howMuchToSeparate) {
//			output.push_back(value_);
//			++separated;
//		}
//
//		return separated;
//	}

	virtual void CollectValues(vector<Value>& output) {
		output.push_back(value_);
		TreeNode<Value>::CollectValues(output);
	}

	virtual size_t GetSize() const {
		return TreeNode<Value>::GetSize() + 1;
	}


private:
	Value value_;
};


template <class Graph>
class DevisibleTree {

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef TreeNodeWithValue<VertexId> Node;
	typedef TreeNode<VertexId> RootNode;


public:

	DevisibleTree(Graph & graph) {
		typedef unordered_map<VertexId, int> RankMap;
		typedef unordered_map<VertexId, VertexId> ParentMap;

		typedef boost::associative_property_map<RankMap> BoostRankMap;
		typedef boost::associative_property_map<ParentMap> BoostParentMap;

		RankMap rank_map;
		ParentMap parent_map;

		BoostRankMap boost_rank_map(rank_map);
		BoostParentMap boost_parent_map(parent_map);

		boost::disjoint_sets<BoostRankMap, BoostParentMap> dset(boost_rank_map, boost_parent_map);

		BOOST_FOREACH(const VertexId& vertex, graph) {
			dset.make_set(vertex);
		}

		BOOST_FOREACH(const VertexId& vertex, graph) {
			nodes_.push_back(Node(vertex));
			index_[vertex] = nodes_.size() - 1;
		}

// build trees
		for (auto it = graph.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			VertexId start = graph.EdgeStart(*it);
			VertexId end = graph.EdgeEnd(*it);

			VertexId start_root = dset.find_set(start);
			VertexId end_root = dset.find_set(end);

			if (start_root != end_root) {
				dset.link(start_root, end_root);

				Node& start_root_node = GetNode(start_root);
				Node& end_root_node = GetNode(end_root);

				if (start_root == dset.find_set(start_root)) {
					start_root_node.AddChild(end_root_node);
					end_root_node.SetParent(start_root_node);
				} else {
					end_root_node.AddChild(start_root_node);
					start_root_node.SetParent(end_root_node);
				}
			}
		}


// find roots of trees. Roots of disjoint dset may not be root of trees,
// but it is good point to start search.
		unordered_set<VertexId> setRoots;
		BOOST_FOREACH(const VertexId& vertex, graph) {
			setRoots.insert(dset.find_set(vertex));
		}

		BOOST_FOREACH(const VertexId& vertex, setRoots) {
			root_.AddChild(*GetNode(vertex).GetRoot());
		}

		root_.sort();
	}

	void SeparateVertices(vector<VertexId>& output, size_t size) {
//		VERIFY(GetSize() >= size);
		TreeNode<VertexId>* node = root_.GetSubtreeWithSize(min(size, GetSize()));
		output.reserve(node->GetSize());
		node->CollectValues(output);
	}

	size_t GetSize() const {
		return root_.GetSize();
	}


private:

	Node& GetNode(const VertexId& vertex) {
		auto it = index_.find(vertex);
		VERIFY(it != index_.end());

		return nodes_[it->second];
	}


private:
	unordered_map<VertexId, size_t> index_;
	RootNode root_;
	vector<Node> nodes_;
};

} // namespace omnigraph


#endif /* DEVISIBLE_TREE_HPP_ */
