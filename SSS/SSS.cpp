// SSS.cpp : This file programs the main SSS framework.

#ifndef SSS_H
#define SSS_H

#define ExpandShow 200
#define ExpandLimit 20000
#define AddFringeOnly true

#include<iostream>
#include<vector>
#include<FeatureSet.h>
#include<SE3Box.h>
#include<WtFp.h>
#include<Eigen/Dense>
#include<Graph.h>
using namespace Eigen;
using namespace std;

template<typename content>
bool operator<(std::pair<std::vector<double>, content> A, std::pair<std::vector<double>, content> B)
{
	for (int i = 0; i < min(A.first.size(), B.first.size()); ++i)
		if (A.first < B.first)
			return true;
		else if (A.first > B.first)
			return false;
		else
			continue;
	return A.first.size() < B.first.size();
}

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

	// Set the root feature sets.
	void set_feature(DeltaFeature F)
	{
		root_feature = F;
	}
	// Set the feature set for box with id.
	void set_feature(SE3Box* b, DeltaFeature new_phi)
	{
		map<int, DeltaFeature>::iterator target = box_features.find(b->ID());
		if (target == box_features.end())
			box_features.insert(make_pair(b->ID(), new_phi));
		else
			target->second = new_phi;
	}
	// Feature set of box b.
	DeltaFeature feature_of(SE3Box* b)
	{
		map<int, DeltaFeature>::iterator target = box_features.find(b->ID());
		if (target == box_features.end())
		{
			// If b is a root.
			if (b->is_root())
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
	// Minimum distance from center of b to its feature set.
	double feature_dis(SE3Box* b)
	{
		Vector3d mB = b->BT()->center();
		DeltaFeature phi = feature_of(b);
		double min_dis = b->width() + 1.0;
		double dis = 0;

		for (int i = 0; i < phi.Vlist.size(); ++i)
		{
			dis = phi.Vlist[i]->distance(mB);
			if (min_dis > dis)
				min_dis = dis;
		}
		for (int i = 0; i < phi.Elist.size(); ++i)
		{
			dis = phi.Elist[i]->distance(mB);
			if (min_dis > dis)
				min_dis = dis;
		}
		for (int i = 0; i < phi.Tlist.size(); ++i)
		{
			dis = phi.Tlist[i]->distance(mB);
			if (min_dis > dis)
				min_dis = dis;
		}
		return min_dis;
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
	void set_pvalue(SE3Box* b,pvalue Cb)
	{
		map<int, pvalue>::iterator target = box_pvalues.find(b->ID());
		if (target == box_pvalues.end())
			box_pvalues.insert(make_pair(b->ID(), Cb));
		else
			target->second = Cb;
	}
	// Construct the WtFp of box b.
	DeltaWtFp* AppFp(SE3Box* b)
	{
		DeltaWtFp* Fp = NULL;
		if (b->is_root())
			Fp = new DeltaWtFp(b->BT()->range);
		else
			Fp = new DeltaWtFp(b->BT()->range, b->BR()->range, b->BR()->WXYZ());
		return Fp;
	}
	// Construct the InFp of box b.
	DeltaInFp* AinFp(SE3Box* b)
	{
		return new DeltaInFp(b->BT()->range);
	}
	// Classification of b by soft pvalue.
	pvalue classify(SE3Box* b, bool show = false)
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

		if (show)
		{
			cout << "New classifying box: ";
			b->out();
			cout << endl;
		}

		DeltaFeature new_phi;
		DeltaFeature new_in_phi;
		DeltaWtFp* Fp = AppFp(b);
		DeltaInFp* Fq = AinFp(b);

		if (show)
		{
			Fp->out();
			Fq->out();
		}

		for (int i = 0; i < phi.Vlist.size(); ++i)
			if (Fq->classify(phi.Vlist[i]) != FREE)
			{
				new_in_phi.Vlist.push_back(phi.Vlist[i]);
				new_phi.Vlist.push_back(phi.Vlist[i]);
			}
			else if (Fp->classify(phi.Vlist[i]) != FREE)
				new_phi.Vlist.push_back(phi.Vlist[i]);
			else
				continue;
		for (int i = 0; i < phi.Elist.size(); ++i)
			if (Fq->classify(phi.Elist[i]) != FREE)
			{
				new_in_phi.Elist.push_back(phi.Elist[i]);
				new_phi.Elist.push_back(phi.Elist[i]);
			}
			else if (Fp->classify(phi.Elist[i]) != FREE)
				new_phi.Elist.push_back(phi.Elist[i]);
			else
				continue;
		for (int i = 0; i < phi.Tlist.size(); ++i)
			if (Fq->classify(phi.Tlist[i]) != FREE)
			{
				new_in_phi.Tlist.push_back(phi.Tlist[i]);
				new_phi.Tlist.push_back(phi.Tlist[i]);
			}
			else if (Fp->classify(phi.Tlist[i]) != FREE)
				new_phi.Tlist.push_back(phi.Tlist[i]);
			else
				continue;

		if (show)
		{
			new_phi.out();
			new_in_phi.out();
		}

		// Stuck check.
		if (new_in_phi.empty())
		{
			// After there is no boundary features, let's determine if the box is stuck or not.
			for (int i = 0; i < phi.Mlist.size(); ++i)
				if (Fq->classify(phi.Mlist[i]) != FREE)
				{
					set_pvalue(b, STUCK);
					set_feature(b, eset);
					if (show)
						cout << "The box is STUCK." << endl;
					return STUCK;
				}
		}

		// Free check.
		if(new_phi.empty())
		{
			set_pvalue(b, FREE);
			set_feature(b, eset);
			if (show)
				cout << "The box is FREE." << endl;
			return FREE;
		}

		for (int i = 0; i < phi.Mlist.size(); ++i)
			if (Fp->quick_classify(phi.Mlist[i]) != FREE)
				new_phi.Mlist.push_back(phi.Mlist[i]);

		set_pvalue(b, MIXED);
		set_feature(b, new_phi);
		if (show)
			cout << "The box is MIXED." << endl;
		return MIXED;

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
	SE3Box* find(VectorXd v, bool show = false)
	{
		if (root->contains(v))
			return find(root, v, show);
		else
		{
			if (show)
				cout << "(" << v.transpose() << ") not found." << endl;
			return NULL;
		}
	}
	// Find the leaf box that contains v from box B (assume that B contains v).
	SE3Box* find(SE3Box* B, VectorXd v, bool show = false)
	{
		if (B->is_leaf())
		{
			if (show)
			{
				cout << "Found box: ";
				B->out();
				cout << endl;
			}
			return B;
		}
		for (int i = 0; i < subsize(B); ++i)
			if (B->child(i)->contains(v))
				return find(B->child(i), v, show);
		if (show)
			cout << "(" << v.transpose() << ") not found." << endl;
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

enum heutype { RAND, BYID, WIDTH, TARGET, GBF, DIS, GBFDIS, WIDIS };

// Config is VectorXd, Box is SE3Box, BoxTree is SE3Tree, Predicate is DeltaPredicate, FeatureSet is DeltaFeature.
template<typename Config, typename Box, typename BoxTree, typename Predicate, typename FeatureSet>
class SSS
{
public:
	// Tree of boxes.
	BoxTree* B;
	// Soft pvalue for robot.
	Predicate* C;
	// Initial configuration.
	Config alpha;
	// Target configuration.
	Config beta;
	// Type of heuristic.
	heutype heu;
	// The free graph.
	Graph<Box> G;
	// Map from box id to graph id.
	map<int, int> Gid;
	// The mixed queue.
	priority_queue<pair<vector<double>, Box*>> Q;
	// The amount of cube that is expanded.
	int expanded;
	// The average R^3-width of cubes.
	double aver3w;
	// The average SO3-width of cubes.
	double aveso3w;
	// Amount of mixed boxes.
	int stat_mixed;
	// Amount of free boxes.
	int stat_free;
	// Amount of stuck boxes.
	int stat_stuck;
	// Average feature amount for mixed boxes.
	double avemixedfeature;
	// Average feature amount for free boxes.
	double avefreefeature;

	Vector3d config_t(Config gamma)
	{
		return Vector3d(gamma(0), gamma(1), gamma(2));
	}

	double dis_heu(Box* b)
	{
		return C->feature_dis(b);
	}

	double GBF_heu(Box* b)
	{
		return max(-log10(min((b->BT()->center() - config_t(alpha)).norm(), (b->BT()->center() - config_t(beta)).norm())), 0.1);
	}

	double width_heu(Box* b)
	{
		return max(b->BT()->width(), b->BR()->width());
	}

	double id_heu(Box* b)
	{
		return 1 / double(b->ID());
	}

	// Heuristic function. The greater, the higher priority.
	vector<double> heuristic(Box* b)
	{
		vector<double> heus;
		switch (heu)
		{
		case GBF: heus.push_back(GBF_heu(b)); break;
		case DIS: heus.push_back(dis_heu(b)); break;
		case WIDTH: heus.push_back(width_heu(b)); break;
		case BYID: heus.push_back(id_heu(b)); break;
		case GBFDIS: heus.push_back(GBF_heu(b) * dis_heu(b)); break;
		case WIDIS: heus.push_back(width_heu(b)); heus.push_back(dis_heu(b)); break;
		default: heus.push_back(rand());
		}
		return heus;
	}

	// Procedures for mixed boxes.
	void add_mixed_node(Box* b, bool show = false)
	{
		if (show)
		{
			cout << "Adding MIXED box: ";
			b->out();
			cout << endl;
		}
		if (!AddFringeOnly)
			Q.push(make_pair(heuristic(b), b));
		stat_mixed++;
		avemixedfeature += (C->feature_of(b).size() - avemixedfeature) / stat_mixed;
	}

	// Procedures for free boxes.
	void add_free_node(Box* b, bool show = false)
	{
		if (show)
		{
			cout << "Adding FREE box: ";
			b->out();
			cout << endl;
		}
		Gid.insert(make_pair(b->ID(), G.insert(b)->id));
		vector<Box*> target_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < target_neighbors.size(); ++j)
		{
			if (!AddFringeOnly)
			{
				GraphNode<Box>* target_neighbor = node(target_neighbors[j]);
				if (target_neighbor != NULL)
					G.link(node(b), target_neighbor);
			}
			else
			{
				// This change makes Q maintains purely fringe boxes.
				Box* tb = target_neighbors[j];
				switch (C->classify(tb))
				{
				case FREE: if (node(tb) != NULL) G.link(node(b), node(tb)); break;
				case MIXED: Q.push(make_pair(heuristic(tb), tb)); break;
				case STUCK: cout << "Error! FREE boxes cannot be adjacent to STUCK boxes."; break;
				case UNKNOWN: cout << "Error! Some boxes are classified as UNKNOWN!"; break;
				default: break;
				}
			}
		}
		stat_free++;
		avefreefeature += (C->feature_of(b).size() - avefreefeature) / stat_free;
	}

	// Procedures for stuck boxes.
	void add_stuck_node(Box* b, bool show = false)
	{
		if (show)
		{
			cout << "Adding STUCK box: ";
			b->out();
			cout << endl;
		}
		stat_stuck++;
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
 
	SSS(BoxTree* b, Predicate* c)
	{
		B = b;
		C = c;
		heu = RAND;
		expanded = 0;
		aver3w = 0;
		aveso3w = 0;
		stat_mixed = 0;
		stat_free = 0;
		stat_stuck = 0;
		avemixedfeature = 0;
		avefreefeature = 0;
	}

	// Show how many boxes that are expanded with hints.
	void show_expansion(string hint)
	{
		cout << expanded << " cubes have been expanded with average width " << aver3w << " * " << aveso3w << "." << endl;
		cout << "There are " << stat_mixed << " mixed boxes, " << stat_free << " free boxes and " << stat_stuck << " stuck boxes." << endl;
		cout << "The average amount of features for mixed boxes are " << avemixedfeature << "." << endl;
		//cout << "The average amount of features for free boxes are " << avefreefeature << "." << endl;
		if (hint.length() > 0)
			cout << hint << endl;
		for (int i = 0; i < 60; ++i)
			cout << "-";
		cout << endl;
	}

	// Update the statistic of expansion.
	void update_expansion(Box* b)
	{
		expanded++;
		aver3w += (b->BT()->width() - aver3w) / expanded;
		aveso3w += (b->BR()->width() - aveso3w) / expanded;
	}

	// Expand(B)
	void Expand(Box* b, bool show = false)
	{
		// Step 0: check if the box is leaf.

		if (!b->is_leaf())
		{
			if (show)
				cout << "Skipping expanding a box that is already expanded." << endl;
			return;
		}

		if (show)
		{
			cout << "Expanding box: ";
			b->out();
			cout << endl;
		}

		// Step 1: subdivide the box.
		b->subdivide();

		// Step 2: classify soft pvalues.
		for (int i = 0; i < B->subsize(b); ++i)
			C->classify(b->child(i), show);

		// Step 3: maintain Q and G.
		for (int i = 0; i < B->subsize(b); ++i)
		{
			if (C->classify(b->child(i)) == MIXED)
				add_mixed_node(b->child(i), show);
			if (C->classify(b->child(i)) == FREE)
				add_free_node(b->child(i), show);
			if (C->classify(b->child(i)) == STUCK)
				add_stuck_node(b->child(i), show);
		}

		// Step 4: show how many cubes that are expanded to give an intuition of the process.
		if (expanded % ExpandShow == 0)
			show_expansion("");
		update_expansion(b);
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
	vector<Config> Find_Path(Config _alpha, Config _beta, FeatureSet Omega, double varepsilon, heutype _h = RAND, bool show = false)
	{
		// Step 0: initializations.
		vector<Config> Path;
		Box* BoxAlpha = NULL, * BoxBeta = NULL;
		expanded = 0;
		aver3w = 0;
		aveso3w = 0;
		stat_mixed = 0;
		stat_free = 0;
		stat_stuck = 0;
		avemixedfeature = 0;
		avefreefeature = 0;
		alpha = _alpha;
		beta = _beta;
		C->set_feature(Omega);
		heu = _h;

		if (show)
			cout << "Find path initialization finished." << endl;

		// Step 1: classify initial boxes.
		if (C->classify(B->Root()) == FREE)
			add_free_node(B->Root());
		if (C->classify(B->Root()) == STUCK)
		{
			show_expansion("The environment is stuck.");
			return Path;													// Empty if no-path.
		}

		if (show)
			cout << "Classify initial boxes finished." << endl;

		// Step 2: find the FREE box containing alpha.
		if (show)
			cout << "Finding alpha: (" << alpha.transpose() << ")." << endl;
		BoxAlpha = B->find(alpha, show);
		if (BoxAlpha == NULL)
		{
			show_expansion("Configuration alpha not in environment.");
			return Path;
		}
		while (BoxAlpha->width() > varepsilon)
		{
			if (C->classify(BoxAlpha) == FREE)
				break;
			if (C->classify(BoxAlpha) == STUCK)
			{
				show_expansion("Configuration alpha is stuck.");
				return Path;
			}
			Expand(BoxAlpha,show);
			BoxAlpha = B->find(BoxAlpha, alpha, show);
		}
		if (C->classify(BoxAlpha) != FREE)
		{
			show_expansion("Configuration alpha is not free.");
			return Path;
		}

		//if (show)
		show_expansion("Found the FREE box containing alpha.");

		// Step 3: find the FREE box containing beta.
		if (show)
			cout << "Finding beta: (" << beta.transpose() << ")." << endl;
		BoxBeta = B->find(beta, show);
		if (BoxBeta == NULL)
		{
			show_expansion("Configuration beta not in environment.");
			return Path;
		}
		while (BoxBeta->width() > varepsilon)
		{
			if (C->classify(BoxBeta) == FREE)
				break;
			if (C->classify(BoxBeta) == STUCK)
			{
				show_expansion("Configuration beta is stuck.");
				return Path;
			}
			Expand(BoxBeta, show);
			BoxBeta = B->find(BoxBeta, beta, show);
		}
		if (C->classify(BoxBeta) != FREE)
		{
			show_expansion("Configuration beta is not free.");
			return Path;
		}

		//if (show)
		show_expansion("Found the FREE box containing beta.");

		// Step 4: expand Q.GetNext() until Box(alpha) and Box(beta) are in the same component.
		while (find(BoxAlpha) != find(BoxBeta))
		{
			if (Q.empty())
			{
				show_expansion("No path found.");
				return Path;
			}
			Box* Q_top = Q.top().second;
			if (Q_top->width() > varepsilon)
				Expand(Q_top, show);
			Q.pop();

			if (expanded > ExpandLimit)
			{
				show_expansion("Run over expand limit (temporary setting for experiments).");
				return Path;
			}
		}

		if (show)
			cout << "Alpha and beta have been in the same FREE component." << endl;

		// Step 5: find the canonical path in the component of G.find(BoxAlph).
		vector<Box*> channel = Find_Channel(BoxAlpha, BoxBeta);
		for (int i = 0; i < channel.size(); ++i)
			Path.push_back(B->center(channel[i]));

		show_expansion("Canonical path found.");
		return Path;														// Empty if no-path.
	}
};

void output_path(VectorXd alpha, VectorXd beta, vector<VectorXd> path)
{
	cout << "The path is the line segments connecting the following points: " << endl;
	cout << alpha.transpose() << endl;
	for (int i = 0; i < path.size(); ++i)
		cout << path[i].transpose() << endl;
	cout << beta.transpose() << endl;
}

#endif