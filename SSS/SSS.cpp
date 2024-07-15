// SSS.cpp : This file programs the main SSS framework.
// Redefinitions are in thie file. See "https://www.cnblogs.com/dan-Blog/articles/8478433.html" for solution.

#ifndef SSS_H
#define SSS_H

#include <iostream>
#include <SE3Box .h>
#include <WtFp.h>
#include <Eigen/Dense>
#include <Graph.h>
using namespace Eigen;
using namespace std;

class DeltaFeature
{
public:
	vector<Point*> Vlist;
	vector<Edge*> Elist;
	vector<Triangle*> Tlist;
	vector<Mesh*> Mlist;
	bool empty()
	{
		return Vlist.empty() && Elist.empty() && Tlist.empty() && Mlist.empty();
	}
};

class DeltaPredicate
{
	// Feature set of the root.
	DeltaFeature root_feature;
	// box_features[b->id] will be the features of box b.
	map<int, DeltaFeature> box_features;
	// Classified predicate values.
	map<int, pvalue> box_pvalues;
public:
	DeltaPredicate() {}

	// Set the feature sets.
	void set_feature(DeltaFeature F)
	{
		root_feature = F;
	}
	// Feature set of box b.
	DeltaFeature feature_of(SE3Box* b)
	{
		map<int, DeltaFeature>::iterator target = box_features.find(b->ID());
		if (target == box_features.end())
		{
			// If b is a root.
			if (b->Parent() == NULL)
			{
				box_features.insert(make_pair(b->ID(), root_feature));
				return root_feature;
			}
			// Otherwise, use its parent's feature.
			DeltaFeature phi = feature_of(b->Parent());
			box_features.insert(make_pair(b->ID(), phi));
			return phi;
		}
		else
			return target->second;
	}
	// Predicate value of box b.
	pvalue pvalue_of(SE3Box* b)
	{
		map<int, pvalue>::iterator target = box_pvalues.find(b->ID());
		if (target == box_pvalues.end())
			return UNKNOWN;
		else
			return target->second;
	}
	// Set the predicate of box b.
	void set_pvalue(SE3Box* b, pvalue Cb)
	{
		box_pvalues.insert(make_pair(b->ID(), Cb));
	}
	// Classification of b by soft pvalue.
	pvalue classify(SE3Box* b)
	{
		// Known case.
		pvalue Cb = pvalue_of(b);
		if (Cb != UNKNOWN)
			return Cb;
		DeltaFeature phi = feature_of(b);
		if (phi.empty())
		{
			set_pvalue(b, FREE);
			return FREE;
		}

		// Unknown case.
		DeltaFeature new_phi;
		DeltaWtFp* Fp;
		if (b->is_root())
			Fp = new DeltaWtFp(b->BT()->range);
		else
			Fp = new DeltaWtFp(b->BT()->range, b->BR()->range, b->BR()->WXYZ());
		for (int i = 0; i < phi.Vlist.size(); ++i)
			switch (Fp->classify(phi.Vlist[i]))
			{
			case FREE:break;
			case MIXED:new_phi.Vlist.push_back(phi.Vlist[i]); break;
			case STUCK:set_pvalue(b, STUCK); return STUCK;
			default:return UNKNOWN;
			}
		for (int i = 0; i < phi.Elist.size(); ++i)
			switch (Fp->classify(phi.Elist[i]))
			{
			case FREE:break;
			case MIXED:new_phi.Elist.push_back(phi.Elist[i]); break;
			case STUCK:set_pvalue(b, STUCK); return STUCK;
			default:return UNKNOWN;
			}
		for (int i = 0; i < phi.Tlist.size(); ++i)
			switch (Fp->classify(phi.Tlist[i]))
			{
			case FREE:break;
			case MIXED:new_phi.Tlist.push_back(phi.Tlist[i]); break;
			case STUCK:set_pvalue(b, STUCK); return STUCK;
			default:return UNKNOWN;
			}
		for (int i = 0; i < phi.Mlist.size(); ++i)
			switch (Fp->classify(phi.Mlist[i]))
			{
			case FREE:break;
			case MIXED:new_phi.Mlist.push_back(phi.Mlist[i]); break;
			case STUCK:set_pvalue(b, STUCK); return STUCK;
			default:return UNKNOWN;
			}
		if (new_phi.empty())
		{
			set_pvalue(b, FREE);
			return FREE;
		}
		else
		{
			set_pvalue(b, MIXED);
			return MIXED;
		}
	}
};

class SE3Tree
{
	// Root nodes.
	SE3Box* root;
public:
	// SE3 initialization.
	SE3Tree(MatrixId range)
	{
		R3Box* BtRoot = new R3Box(range);
		SO3Box* BrTree = new SO3Box();
		root = new SE3Box(BtRoot, BrTree);
	}
	// Find the leaf box that contains v.
	SE3Box* find(VectorXd v)
	{
		if (root->contains(v))
			return find(root, v);
		else
			return NULL;
	}
	// Find the leaf box that contains v from box B (assume that B contains v).
	SE3Box* find(SE3Box* B, VectorXd v)
	{
		if (B->is_leaf())
			return B;
		for (int i = 0; i < subsize(B); ++i)
			if (B->child(i)->contains(v))
				return find(B->child(i), v);
		return NULL;
	}
	// Root of the tree.
	SE3Box* Root()
	{
		return root;
	}
	// Number of boxes in root.
	int rootsize()
	{
		return 4;
	}
	// Number of children of b.
	int subsize(SE3Box* b)
	{
		if (b->is_root() && b->is_BRsub())
			return 4;
		else
			return 8;
	}
	// Number of directions of b.
	int dirsize(SE3Box* b)
	{
		return 7;
	}
	// Get the center point of a box.
	VectorXd center(SE3Box* b)
	{
		return b->center();
	}
	// Get the center point of the intersecting face of two neighbors.
	VectorXd cross_center(SE3Box* b0, SE3Box* b1)
	{
		// To be implemented.
		VectorXd c(7);
		return c;
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0)
	{
		root->out(os, l);
	}
};

// Config is VectorXd, Predicate is the box tree containing member UnionNode* comp, FeatureSet is vector<Feature*>, IniSize is the number of initial boxes.
template<typename Config, typename Box, typename BoxTree, typename Predicate, typename FeatureSet>
class SSS
{
	// Tree of boxes.
	BoxTree* B;
	// Soft pvalue for robot.
	Predicate* C;
	// The free graph.
	Graph<Box> G;
	// Map from box id to graph id.
	map<int, int> Gid;
	// The mixed queue.
	priority_queue<Box*> Q;

	// Procedures for mixed boxes.
	void add_mixed_node(Box* b)
	{
		Q.push(b);
	}

	// Procedures for free boxes.
	void add_free_node(Box* b)
	{
		Gid.insert(b->ID(), G.insert(b)->id);
		vector<Box*> target_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < target_neighbors.size(); ++j)
		{
			GraphNode<Box>* target_neighbor = node(target_neighbors[j]);
			if (target_neighbor != NULL)
				G.link(node(b), target_neighbor);
		}
	}

	// Find the node of a box in the free graph, return NULL if not exists.
	GraphNode<Box>* node(Box* b)
	{
		if (Gid.find(b->ID()) == Gid.end())
			return NULL;
		return G.node(Gid[b->ID()]);
	}

	// Find the class of a box in the free graph, return NULL if not exists.
	GraphNode<Box>* find(Box* b)
	{
		if (Gid.find(b->ID()) == Gid.end())
			return NULL;
		return G.find(G.node(Gid[b->ID()]));
	}

public:
	SSS(BoxTree* b, Predicate* c)
	{
		B = b;
		C = c;
	}
	// Expand(B)
	void Expand(Box* b)
	{
		// Step 1: subdivide the box.
		b->subdivide();

		// Step 2: classify soft pvalues.
		for (int i = 0; i < B->subsize(b); ++i)
			C->classify(b->child(i));

		// Step 3: maintain Q and G.
		for (int i = 0; i < B->subsize(b); ++i)
		{
			if (C->classify(b->child(i)) == MIXED)
				add_mixed_node(b->child(i));
			if (C->classify(b->child(i)) == FREE)
				add_free_node(b->child(i));
		}
	}

	// Find the channel when G.find(BoxAlpha) == G.find(BoxBeta) != NULL.
	vector<Box*> Find_Channel(Box* BoxAlpha, Box* BoxBeta)
	{
		vector<Box*> channel;
		list<GraphNode<Box>*> free_path = G.path(node(BoxAlpha), node(BoxBeta));
		for (auto it = free_path.begin(); it != free_path.end(); ++it)
			channel.push_back((*it)->content);
		return channel;
	}

	// SSS framework main process.
	vector<Config> Find_Path(Config alpha, Config beta, FeatureSet Omega, double varepsilon)
	{
		// Step 0: initializations.
		vector<Config> Path;
		Box* BoxAlpha = NULL, * BoxBeta = NULL;
		C->set_feature(Omega);

		// Step 1: classify initial boxes.
		BoxAlpha = B->find(alpha);
		BoxBeta = B->find(beta);
		if (BoxAlpha == NULL || BoxBeta == NULL)							// Initial or Target configurations not found.
			return Path;													// Empty if no-path.

		// Step 2: find the FREE box containing alpha.
		while (BoxAlpha->width() > varepsilon)
		{
			if (C->classify(BoxAlpha) == FREE)
				break;
			if (C->classify(BoxAlpha) == STUCK)
				return Path;
			Expand(BoxAlpha);
			for (int i = 0; i < B->subsize(BoxAlpha); ++i)
				if (BoxAlpha->child(i)->contains(alpha))
				{
					BoxAlpha = BoxAlpha->child(i);
					break;
				}
		}
		if (C->classify(BoxAlpha) != FREE)
			return Path;

		// Step 3: find the FREE box containing beta.
		while (BoxBeta->width() > varepsilon)
		{
			if (C->classify(BoxBeta) == FREE)
				break;
			if (C->classify(BoxBeta) == STUCK)
				return Path;
			Expand(BoxBeta);
			for (int i = 0; i < B->subsize(BoxBeta); ++i)
				if (BoxBeta->child(i)->contains(beta))
				{
					BoxBeta = BoxBeta->child(i);
					break;
				}
		}
		if (C->classify(BoxBeta) != FREE)
			return Path;

		// Step 4: expand Q.GetNext() until Box(alpha) and Box(beta) are in the same component.
		while (find(BoxAlpha) != find(BoxBeta))
		{
			if (Q.empty())
				return Path;
			if (Q.top()->width() > varepsilon)
				Expand(Q.top());
			Q.pop();
		}

		// Step 5: find the canonical path in the component of G.find(BoxAlph).
		vector<Box*> channel = Find_Channel(BoxAlpha, BoxBeta);
		for (int i = 0; i < channel.size(); ++i)
			Path.push_back(B->center(channel[i]));
		return Path;														// Empty if no-path.
	}
};

#endif