//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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
#include <queue>

#include "standard_base.hpp"
#include "omni_utils.hpp"


namespace omnigraph {

template <class Value>
class TreeNode {

	struct NodeComparator {
		bool operator() (TreeNode* node1, TreeNode* node2) {
			return node1->GetSize() < node2->GetSize();
		}
	};

	struct SizePred {
		const size_t size;
		SizePred(size_t size) : size(size) { }
		bool operator()(TreeNode* node) {
			return node->GetSize() >= size;
		}
	};

public:
	TreeNode() : subtree_size_(0) { }

	void AddChild(TreeNode& node) {
		children_.push_back(&node);
	}

	void UpdateSubtreeSize() {
		BOOST_FOREACH(TreeNode* node, children_) {
			subtree_size_ += node->GetSize();
		}
	}

	virtual size_t GetSize() const {
		return subtree_size_;
	}

	TreeNode* GetSubtreeWithSize(const size_t size) {
		// find big enough child.
		auto it = std::find_if(children_.begin(), children_.end(), SizePred(size));
		if (it == children_.end()) {
			return 0;
		}

		TreeNode* node = (*it)->GetSubtreeWithSize(size);
		if (node == 0) {
			node = *it;
			children_.erase(it);
		}

		subtree_size_ -= node->GetSize();
		return node;
	}

	virtual void CollectValue(vector<Value>& output) {
	}

	virtual void CollectNodes(std::queue<TreeNode<Value>*>& nodes) {
		BOOST_FOREACH(TreeNode* node, children_) {
			nodes.push(node);
		}
		subtree_size_ = 0;
		children_.clear();
	}

	virtual ~TreeNode() { }

private:
	list<TreeNode *> children_;
	size_t subtree_size_;
};


template <class Value>
class TreeNodeWithValue : public TreeNode<Value> {

public:
	TreeNodeWithValue(const Value& value) : value_(value) {
	}

	virtual ~TreeNodeWithValue() { }

	virtual void CollectValue(vector<Value>& output) {
		output.push_back(value_);
	}

	virtual size_t GetSize() const {
		return TreeNode<Value>::GetSize() + 1;
	}

	Value GetValue() const {
		return value_;
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

	DevisibleTree(Graph & graph) : graph_(graph) {
		typedef unordered_map<VertexId, int> RankMap;
		typedef unordered_map<VertexId, VertexId> ParentMap;

		typedef boost::associative_property_map<RankMap> BoostRankMap;
		typedef boost::associative_property_map<ParentMap> BoostParentMap;

		RankMap rank_map;
		ParentMap parent_map;

		BoostRankMap boost_rank_map(rank_map);
		BoostParentMap boost_parent_map(parent_map);

		boost::disjoint_sets<BoostRankMap, BoostParentMap> dset(boost_rank_map, boost_parent_map);

		BOOST_FOREACH(const VertexId& vertex, graph_) {
			dset.make_set(vertex);
		}

		BOOST_FOREACH(const VertexId& vertex, graph_) {
			nodes_.push_back(Node(vertex));
			index_[vertex] = nodes_.size() - 1;
		}

		TRACE("Creating tree of size:" << nodes_.size());


// build trees
//		for (auto it = graph_.SmartEdgeBegin(LengthComparator<Graph>(graph_)); !it.IsEnd(); ++it) {
		BOOST_FOREACH(VertexId vertex, graph) {
			BOOST_FOREACH(EdgeId edge, graph.OutgoingEdges(vertex)) {
				VertexId start = graph_.EdgeStart(edge);
				VertexId end = graph_.EdgeEnd(edge);

				VertexId start_root = dset.find_set(start);
				VertexId end_root = dset.find_set(end);

				if (start_root != end_root) {
					dset.link(start_root, end_root);
					edges_.insert(edge);
				}
			}
		}

		TRACE("Node quantity: " << nodes_.size());
		TRACE("Edges for tree: " << edges_.size());

		unordered_set<VertexId> forest_roots;

		BOOST_FOREACH(VertexId vertex, graph_) {
			forest_roots.insert(dset.find_set(vertex));
		}

		BOOST_FOREACH(VertexId vertex, forest_roots) {
			Node& node = GetNode(vertex);
			TRACE("Adding " << vertex);
			CreateTree(node, edges_);
			root_.AddChild(node);
		}
		root_.UpdateSubtreeSize();

	}

	void SeparateVertices(vector<VertexId>& output, size_t size) {
		TreeNode<VertexId>* node = root_.GetSubtreeWithSize(min(size, GetSize()));
		if (node == 0) {
			node = &root_;
		}
		output.reserve(node->GetSize());
		CollectValues(node, output);
	}

	void CollectValues(TreeNode<VertexId>* root, vector<VertexId>& output) {
		std::queue<TreeNode<VertexId>*> nodes;
		nodes.push(root);

		while (nodes.size() > 0) {
			TreeNode<VertexId>* node = nodes.front();
			nodes.pop();
			node->CollectNodes(nodes);
			node->CollectValue(output);
		}
	}

	size_t GetSize() const {
		return root_.GetSize();
	}


private:

	const vector<EdgeId> GetEdges(VertexId vertex) const {
		vector<EdgeId> result;
		vector<EdgeId> incoming = graph_.IncomingEdges(vertex);
		vector<EdgeId> outgoing = graph_.OutgoingEdges(vertex);
		result.insert(result.end(), outgoing.begin(), outgoing.end());
		result.insert(result.end(), incoming.begin(), incoming.end());
		return result;
	}

	const vector<EdgeId> Filter(const vector<EdgeId>& vertex_edges, unordered_set<EdgeId>& edges) const {
		vector<EdgeId> result;
		BOOST_FOREACH(EdgeId edge, vertex_edges) {
			auto it = edges.find(edge);
			if (it != edges.end()) {
				TRACE("Edge " << edge << " went through the filter");
				result.push_back(edge);
				edges.erase(it);
			}
		}
		return result;
	}

	VertexId GetSecond(VertexId first, EdgeId edge) const {
		VertexId start = graph_.EdgeStart(edge);
		VertexId end = graph_.EdgeEnd(edge);
		return (first == start) ? end : start;
	}


	const vector<VertexId> GetNeighbours(VertexId vertex, unordered_set<EdgeId>& edges) const {
		const vector<EdgeId> vertex_tree_edges = Filter(GetEdges(vertex), edges);
		vector<VertexId> neighbours;
		BOOST_FOREACH(EdgeId edge, vertex_tree_edges) {
			neighbours.push_back(GetSecond(vertex, edge));
		}
		return neighbours;
	}

	enum {
	    white,
	    grey
	};
	void CreateTree(Node& root, unordered_set<EdgeId>& tree_edges) {
		typedef pair<Node*, bool> stack_elem;
		stack<stack_elem> nodes;
		nodes.push(stack_elem(&root, white));

		while(!nodes.empty()) {
			stack_elem elem = nodes.top();
			nodes.pop();
			Node* node = elem.first;

			if (elem.second == white) {
				nodes.push(stack_elem(node, grey));
				const vector<VertexId> neighbours = GetNeighbours(node->GetValue(), tree_edges);
				TRACE("Tree has " << tree_edges.size() << " edges");
				BOOST_FOREACH(VertexId neighbour, neighbours) {
					TRACE("Adding " << neighbour);
					Node& child = GetNode(neighbour);
					node->AddChild(child);
					nodes.push(stack_elem(&child, white));
				}
			} else {
				node->UpdateSubtreeSize();
			}
		}
	}

	Node& GetNode(const VertexId& vertex) {
		auto it = index_.find(vertex);
		VERIFY(it != index_.end());

		return nodes_[it->second];
	}


private:
	unordered_map<VertexId, size_t> index_;
	unordered_set<EdgeId> edges_;
	vector<Node> nodes_;
	Graph& graph_;
	RootNode root_;

private:
	DECL_LOGGER("DevisibleTree");
};

} // namespace omnigraph


#endif /* DEVISIBLE_TREE_HPP_ */
