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

class Feature_Set
{
public:
	vector<Point*> Vlist;
	vector<Edge*> Elist;
	vector<Triangle*> Tlist;
	vector<Mesh*> Mlist;
};

class SE3SSS
{
	SE3Tree* T;
	vector<Point*> Vlist;
	vector<Edge*> Elist;
	vector<Triangle*> Tlist;
	vector<Mesh*> Mlist;
public:
	SE3SSS() { T = NULL; }
	SE3SSS(MatrixId B)
	{
		T = new SE3Tree(B);
	}
	void set_range(MatrixId B)
	{
		T = new SE3Tree(B);
	}

	// Set the feature sets.
	void set_feature(Feature_Set F)
	{
		Vlist = F.Vlist;
		Elist = F.Elist;
		Tlist = F.Tlist;
		Mlist = F.Mlist;
	}
	// Find the leaf box that contains v.
	SE3Box* find(VectorXd v)
	{
		return NULL;
	}
	// Number of children of b.
	int subsize(SE3Box* b)
	{
		return 8;
	}
	// Number of directions of b.
	int dirsize(SE3Box* b)
	{
		return 7;
	}
	// Classification of b by soft predicate.
	int classify(SE3Box* b)
	{
		return MIXED;
	}
	// If the box b contains a configuration v.
	bool contains(SE3Box* b, VectorXd v)
	{
		return true;
	}
	// Get the center point of a box.
	VectorXd center(SE3Box* b)
	{
		Vector3d ct = (b->BT()->range.min() + b->BT()->range.max()) / 2;
		Vector3d cr = (b->BR()->range.min() + b->BR()->range.max()) / 2;
		VectorXd c(7);
		for (int i = 0; i < 3; ++i)
			c(i) = ct(i);
		switch (b->BR()->WXYZ())
		{
		case 0: c(3) = 1; c(4) = cr(0); c(5) = cr(1); c(6) = cr(2); break;
		case 1: c(3) = cr(0); c(4) = 1; c(5) = cr(1); c(6) = cr(2); break;
		case 2: c(3) = cr(0); c(4) = cr(1); c(5) = 1; c(6) = cr(2); break;
		case 3: c(3) = cr(0); c(4) = cr(1); c(5) = cr(2); c(6) = 1; break;
		default:c(3) = 1; c(4) = 0; c(5) = 0; c(6) = 0;
		}
		return c;
	}
	// Get the center point of the intersecting face of two neighbors.
	VectorXd cross_center(SE3Box* b0, SE3Box* b1)
	{
		VectorXd c(7);
		return c;
	}
};

// Config is VectorXd, SSSTree is the box tree containing member UnionNode* comp, FeatureSet is vector<Feature*>, IniSize is the number of initial boxes.
template<typename Config, typename Box, typename SSSTree, typename FeatureSet>
class SSS
{
	// The SSS framework maintainance class.
	SSSTree* B;
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
		for (int j = 0; j < B->dirsize(b); ++j)
		{
			GraphNode<Box>* target_neighbor = node(b->neg_neighbor(j));
			if (target_neighbor != NULL)
				G.link(node(b), target_neighbor);
			target_neighbor = node(b->pos_neighbor(j));
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
	SSS(SSSTree* b)
	{
		B = b;
	}
	// Expand(B)
	void Expand(Box* b)
	{
		// Step 1: subdivide the box.
		b->subdivide();

		// Step 2: classify soft predicates.
		for (int i = 0; i < B->subsize(b); ++i)
			B->classify(b->child(i));

		// Step 3: maintain Q and G.
		for (int i = 0; i < B->subsize(b); ++i)
		{
			if (B->classify(b->child(i)) == MIXED)
				add_mixed_node(b->child(i));
			if (B->classify(b->child(i)) == FREE)
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
		vector<Config> Path;
		Box* BoxAlpha = NULL, * BoxBeta = NULL;

		// Step 1: classify initial boxes.
		B->set_feature(Omega);
		BoxAlpha = B->find(alpha);
		BoxBeta = B->find(beta);
		if (BoxAlpha == NULL || BoxBeta == NULL)							// Initial or Target configurations not found.
			return Path;													// Empty if no-path.

		// Step 2: find the FREE box containing alpha.
		while (BoxAlpha->width() > varepsilon)
		{
			if (B->classify(BoxAlpha) == FREE)
				break;
			if (B->classify(BoxAlpha) == STUCK)
				return Path;
			Expand(BoxAlpha);
			for (int i = 0; i < B->subsize(BoxAlpha); ++i)
				if (B->contains(BoxAlpha->child(i), alpha))
				{
					BoxAlpha = BoxAlpha->child(i);
					break;
				}
		}
		if (B->classify(BoxAlpha) != FREE)
			return Path;

		// Step 3: find the FREE box containing beta.
		while (BoxBeta->width() > varepsilon)
		{
			if (B->classify(BoxBeta) == FREE)
				break;
			if (B->classify(BoxBeta) == STUCK)
				return Path;
			Expand(BoxBeta);
			for (int i = 0; i < B->subsize(BoxBeta); ++i)
				if (B->contains(BoxBeta->child(i), beta))
				{
					BoxBeta = BoxBeta->child(i);
					break;
				}
		}
		if (B->classify(BoxBeta) != FREE)
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


int main()
{
	cout << "Hello World!" << endl;
    return 0;
}

#endif