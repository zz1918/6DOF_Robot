// SSS.cpp : This file programs the main SSS framework.

#ifndef SSS_H
#define SSS_H

#define AddFringeOnly true
#define INF 2147483647

#include<iostream>
#include<string>
#include<vector>
#include<list>
#include<queue>
#include<map>
#include<set>
#include<stack>
#include<FeatureSet.h>
#include<SE3Box.h>
#include<WtFp.h>
#include<Eigen/Dense>
#include<Graph.h>
#include<ReadSSSCommand.h>
#include<ViewerControl.h>
#include <chrono>
using namespace Eigen;
using namespace std;

int extern ExpandLimit;
int extern ExpandShow;
int extern Noisity;
SSSViewer extern viewer;
EnvironmentFeature extern env;
VectorXd extern SSSalpha, SSSbeta;
double extern envrange;
vector<Vector3d> extern SSShints;

bool is_hinted(SE3Box *b)
{
	for (int i = 0; i < SSShints.size(); ++i)
	{
		if (b->BT()->contains(SSShints[i]))
			return true;
	}
	return false;
}

//*************** Statistic variables ***************//

long long set_fringe_time = 0;
long long set_box_time = 0;
long long set_box_preparation_time = 0;
long long get_color_time = 0;
long long get_graph_node_time = 0;
long long waste_time = 0;
long long classification_time = 0;
long long classification_C1_time = 0;
long long construct_box_time = 0;
long long classify_p0_known_time = 0;
long long p0_quick_classify_time = 0;
long long p0_classify_feature_time = 0;

// The amount of cube that is expanded.
int expanded = 0;
// The average R^3-width of cubes.
double aver3w;
// The average SO3-width of cubes.
double aveso3w;
// Amount of mixed boxes.
int stat_mixed = 0;
// Amount of free boxes.
int stat_free = 0;
// Amount of stuck boxes.
int stat_stuck = 0;
// Amount of veps-small boxes.
int stat_small = 0;
// Average feature amount for mixed boxes.
double avemixedfeature;
// Average feature amount for free boxes.
double avefreefeature;
// Maximum size among all leaves.
double MaxSize;
// Minimum size among all leaves.
double MinSize;

//***************** Helper Functions ******************//

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
// Set heuristic from one number.
std::vector<double> key(double k)
{
	std::vector<double> Key;
	Key.push_back(k);
	return Key;
}

// Get the T-part of configuration.
Vector3d config_t(VectorXd gamma)
{
	return Vector3d(gamma(0), gamma(1), gamma(2));
}

// Reset the viewer as the environment.
void reset_viewer()
{
	viewer.set_env(env, envrange, SSSalpha, SSSbeta);
}

//*****************************************************//
class DeltaPredicate
{
	// Feature set of the root.
	DeltaFeature root_feature;
	// box_features[b->id] will be the features of box b.
	map<int, DeltaFeature> box_features;
	// Classified predicate values.
	map<int, pvalue> box_pvalues;
	// Classified C1 values.
	map<int, pvalue> box_C1_pvalues;
	// How many predicates?
	int size = 2;
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
	// C1 predicate value of box b.
	pvalue pvalue_C1_of(SE3Box* b)
	{
		map<int, pvalue>::iterator target = box_C1_pvalues.find(b->ID());
		if (target == box_C1_pvalues.end())
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
	// Set the C1 predicate of box b.
	void set_C1_pvalue(SE3Box* b, pvalue Cb)
	{
		map<int, pvalue>::iterator target = box_C1_pvalues.find(b->ID());
		if (target == box_C1_pvalues.end())
			box_C1_pvalues.insert(make_pair(b->ID(), Cb));
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
	// Construct the OExFp of box b.
	OExFp* OexFp(SE3Box* b)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		OExFp* newFp = new OExFp(b->BT()->range);
		auto end_time = std::chrono::high_resolution_clock::now();
		construct_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return newFp;
	}
	// Construct the CnFp of two boxes b and p.
	/*
	DeltaCnFp* ACnFp(SE3Box* b, SE3Box* p)
	{
		return new DeltaCnFp(b->BT()->range, p->BT()->range);
	}*/
	// Predicate size.
	int psize()
	{
		return size;
	}
	// Classification of b under different predicates.
	pvalue classify(SE3Box* b, int k = 0, bool show = false)
	{
		switch (k)
		{
		case 2:return FREE;
		case 1:
		{
			auto start_time = std::chrono::high_resolution_clock::now();
			pvalue C1b = C1(b, show);
			auto end_time = std::chrono::high_resolution_clock::now();
			classification_C1_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return C1b;
		}
		default:return wtC(b, show);
		}
	}
	// Classification of b by point O.
	pvalue C1(SE3Box* b, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		auto end_time = std::chrono::high_resolution_clock::now();
		// Known case.
		pvalue Cb = pvalue_C1_of(b);
		if (Cb != UNKNOWN)
		{
			end_time = std::chrono::high_resolution_clock::now();
			classify_p0_known_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return Cb;
		}
		DeltaFeature phi = feature_of(b);
		if (phi.empty())
		{
			set_C1_pvalue(b, FREE);
			end_time = std::chrono::high_resolution_clock::now();
			classify_p0_known_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return FREE;
		}

		end_time = std::chrono::high_resolution_clock::now();
		classify_p0_known_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		// Unknonw case.
		DeltaFeature new_phi;
		OExFp* Fp = OexFp(b);
		DeltaInFp* Fq = AinFp(b);
		start_time = std::chrono::high_resolution_clock::now();
		bool FpquickFree = true;
		for (int i = 0; i < phi.Mlist.size(); ++i)
			if (Fp->quick_classify(phi.Mlist[i]) == STUCK)
			{
				set_C1_pvalue(b, STUCK);
				end_time = std::chrono::high_resolution_clock::now();
				p0_quick_classify_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
				return STUCK;
			}
			else if (Fp->quick_classify(phi.Mlist[i]) != FREE)
				FpquickFree = false;
		if (FpquickFree)
		{
			set_C1_pvalue(b, FREE);
			end_time = std::chrono::high_resolution_clock::now();
			p0_quick_classify_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return FREE;
		}
		end_time = std::chrono::high_resolution_clock::now();
		p0_quick_classify_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		start_time = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < phi.Vlist.size(); ++i)
			if (Fq->classify(phi.Vlist[i]) != FREE)
				new_phi.Vlist.push_back(phi.Vlist[i]);
		for (int i = 0; i < phi.Elist.size(); ++i)
			if (Fq->classify(phi.Elist[i]) != FREE)
				new_phi.Elist.push_back(phi.Elist[i]);
		for (int i = 0; i < phi.Tlist.size(); ++i)
			if (Fq->classify(phi.Tlist[i]) != FREE)
				new_phi.Tlist.push_back(phi.Tlist[i]);
		end_time = std::chrono::high_resolution_clock::now();
		p0_classify_feature_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		// Stuck check.
		if (new_phi.empty())
		{
			// After there is no boundary features, let's determine if the box is stuck or not.
			for (int i = 0; i < phi.Mlist.size(); ++i)
				if (Fp->classify(phi.Mlist[i]) != FREE)
				{
					set_C1_pvalue(b, STUCK);
					return STUCK;
				}
			set_C1_pvalue(b, FREE);
			return FREE;
		}
		else
		{
			set_C1_pvalue(b, MIXED);
			return MIXED;
		}
	}
	// Classification of b by soft predicate.
	pvalue wtC(SE3Box* b, bool show = false)
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

		// C1 predicate.
		if (pvalue_C1_of(b) == UNKNOWN)
		{
			auto start_time = std::chrono::high_resolution_clock::now();
			C1(b, show);
			auto end_time = std::chrono::high_resolution_clock::now();
			classification_C1_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		}
		if (pvalue_C1_of(b) == STUCK)
		{
			set_pvalue(b, STUCK);
			return STUCK;
		}

		if (Noisity >= 10)
		{
			cout << "New classifying box: ";
			b->out();
			cout << endl;
		}

		DeltaFeature new_phi;
		DeltaWtFp* Fp = AppFp(b);

		// Use approximate footprint of point A to quick exclude stuckness.
		AapFp* Fa = new AapFp(Fp);
		bool Amixed = false;
		if (!Amixed)
		{
			for (int i = 0; i < phi.Vlist.size(); ++i)
				if (Fa->classify(phi.Vlist[i]) != FREE)
				{
					Amixed = true;
					break;
				}
		}
		if (!Amixed)
		{
			for (int i = 0; i < phi.Elist.size(); ++i)
				if (Fa->classify(phi.Elist[i]) != FREE)
				{
					Amixed = true;
					break;
				}
		}
		if (!Amixed)
		{
			for (int i = 0; i < phi.Tlist.size(); ++i)
				if (Fa->classify(phi.Tlist[i]) != FREE)
				{
					Amixed = true;
					break;
				}
		}
		if (!Amixed)
		{
			// After there is no boundary features, let's determine if the box is stuck or not.
			for (int i = 0; i < phi.Mlist.size(); ++i)
				if (Fa->classify(phi.Mlist[i]) != FREE)
				{
					set_pvalue(b, STUCK);
					set_feature(b, eset);
					if (Noisity >= 10)
						cout << "The box is STUCK." << endl;
					return STUCK;
				}
		}

		// Use approximate footprint of point B to quick exclude stuckness.
		BapFp* Fb = new BapFp(Fp);
		bool Bmixed = false;
		if (!Bmixed)
		{
			for (int i = 0; i < phi.Vlist.size(); ++i)
				if (Fb->classify(phi.Vlist[i]) != FREE)
				{
					Bmixed = true;
					break;
				}
		}
		if (!Bmixed)
		{
			for (int i = 0; i < phi.Elist.size(); ++i)
				if (Fb->classify(phi.Elist[i]) != FREE)
				{
					Bmixed = true;
					break;
				}
		}
		if (!Bmixed)
		{
			for (int i = 0; i < phi.Tlist.size(); ++i)
				if (Fb->classify(phi.Tlist[i]) != FREE)
				{
					Bmixed = true;
					break;
				}
		}
		if (!Bmixed)
		{
			// After there is no boundary features, let's determine if the box is stuck or not.
			for (int i = 0; i < phi.Mlist.size(); ++i)
				if (Fb->classify(phi.Mlist[i]) != FREE)
				{
					set_pvalue(b, STUCK);
					set_feature(b, eset);
					if (Noisity >= 10)
						cout << "The box is STUCK." << endl;
					return STUCK;
				}
		}

		// Both A and B are not able to assign stuck. Let's focus on the real approximate footprint.
		if (Noisity >= 15)
			Fp->out();

		for (int i = 0; i < phi.Vlist.size(); ++i)
			if (Fp->classify(phi.Vlist[i]) != FREE)
				new_phi.Vlist.push_back(phi.Vlist[i]);
		for (int i = 0; i < phi.Elist.size(); ++i)
			if (Fp->classify(phi.Elist[i]) != FREE)
				new_phi.Elist.push_back(phi.Elist[i]);
		for (int i = 0; i < phi.Tlist.size(); ++i)
			if (Fp->classify(phi.Tlist[i]) != FREE)
				new_phi.Tlist.push_back(phi.Tlist[i]);

		if (Noisity >= 15)
			new_phi.out();

		// Stuck check.
		if (new_phi.empty())
		{
			// After there is no boundary features, let's determine if the box is stuck or not.
			for (int i = 0; i < phi.Mlist.size(); ++i)
				if (Fp->classify(phi.Mlist[i]) != FREE)
				{
					set_pvalue(b, STUCK);
					set_feature(b, eset);
					if (Noisity >= 10)
						cout << "The box is STUCK." << endl;
					return STUCK;
				}
		}

		// Free check.
		if(new_phi.empty())
		{
			set_pvalue(b, FREE);
			set_feature(b, eset);
			if (Noisity >= 10)
				cout << "The box is FREE." << endl;
			return FREE;
		}

		for (int i = 0; i < phi.Mlist.size(); ++i)
			if (Fp->quick_classify(phi.Mlist[i]) != FREE)
				new_phi.Mlist.push_back(phi.Mlist[i]);

		set_pvalue(b, MIXED);
		set_feature(b, new_phi);
		if (Noisity >= 10)
			cout << "The box is MIXED." << endl;
		return MIXED;
	}
	// Classification of the connection of two boxes (if it's not stuck, then it is good).
	/*
	pvalue Cconnect(SE3Box* b, SE3Box* p, bool show = false)
	{
		DeltaFeature phi1 = feature_of(b);
		DeltaFeature phi2 = feature_of(p);
		DeltaFeature new_phi;
		DeltaCnFp* Fp = ACnFp(b, p);
		for (int i = 0; i < phi1.Vlist.size(); ++i)
			if (Fp->classify(phi1.Vlist[i]) != FREE)
				new_phi.Vlist.push_back(phi1.Vlist[i]);
		for (int i = 0; i < phi1.Elist.size(); ++i)
			if (Fp->classify(phi1.Elist[i]) != FREE)
				new_phi.Elist.push_back(phi1.Elist[i]);
		for (int i = 0; i < phi1.Tlist.size(); ++i)
			if (Fp->classify(phi1.Tlist[i]) != FREE)
				new_phi.Tlist.push_back(phi1.Tlist[i]);
		for (int i = 0; i < phi2.Vlist.size(); ++i)
			if (Fp->classify(phi2.Vlist[i]) != FREE)
				new_phi.Vlist.push_back(phi2.Vlist[i]);
		for (int i = 0; i < phi2.Elist.size(); ++i)
			if (Fp->classify(phi2.Elist[i]) != FREE)
				new_phi.Elist.push_back(phi2.Elist[i]);
		for (int i = 0; i < phi2.Tlist.size(); ++i)
			if (Fp->classify(phi2.Tlist[i]) != FREE)
				new_phi.Tlist.push_back(phi2.Tlist[i]);
		// Stuck check.
		if (new_phi.empty())
		{
			// After there is no boundary features, let's determine if the box is stuck or not.
			for (int i = 0; i < phi1.Mlist.size(); ++i)
				if (Fp->classify(phi1.Mlist[i]) != FREE)
					return STUCK;
			for (int i = 0; i < phi2.Mlist.size(); ++i)
				if (Fp->classify(phi2.Mlist[i]) != FREE)
					return STUCK;
			return FREE;
		}
		else
			return MIXED;
	}*/
};

class SE3Tree
{
	// Root nodes.
	SE3Box* root;
	// Range of the T part.
	double range;
public:
	// SE3 initialization.
	SE3Tree(MatrixId _range)
	{
		range = max(_range.max().cwiseAbs().maxCoeff(), _range.min().cwiseAbs().maxCoeff());
		R3Box* BtRoot = new R3Box(_range);
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
			if (Noisity >= 10)
				cout << "(" << v.transpose() << ") not found." << endl;
			return NULL;
		}
	}
	// Find the leaf box that contains v from box B (assume that B contains v).
	SE3Box* find(SE3Box* B, VectorXd v, bool show = false)
	{
		if (B->is_leaf())
		{
			if (Noisity >= 10)
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
		if (Noisity >= 10)
			cout << "(" << v.transpose() << ") not found." << endl;
		return NULL;
	}
	// Root of the tree.
	SE3Box* Root()
	{
		return root;
	}
	// Range of the Tree.
	double Range()
	{
		return range;
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

template<typename Box>
class GBFringe
{
	double range;
	double veps;
	// Alpha/beta boxes/fringes.
	set<int> AlphaBox, BetaBox, AlphaFringe, BetaFringe;
	// Forbidden area.
	set<int> forbid;
	// Representative alpha/beta box.
	int AlphaRepr, BetaRepr;

	// Find the global node of a box in the global graph, return NULL if not exists.
	GraphNode<Box>* gnode(Box* b, Graph<Box>& Global, map<int, int>& GlobalId)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (GlobalId.find(b->ID()) == GlobalId.end())
			return NULL;
		GraphNode<Box>* gv = Global.node(GlobalId[b->ID()]);
		auto end_time = std::chrono::high_resolution_clock::now();
		get_graph_node_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return gv;
	}
	// Get the t-color of a box b in global graph.
	Vcolor getcolor(Box* b, Graph<Box>& Global, map<int, int>& GlobalId, int t)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		Vcolor c = gnode(b, Global, GlobalId)->color(t);
		auto end_time = std::chrono::high_resolution_clock::now();
		get_color_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return c;
	}
	// Node id in global graph.
	int gid(Box* b, map<int, int>& GlobalId)
	{
		if (GlobalId.find(b->ID()) == GlobalId.end())
			return -1;
		return GlobalId[b->ID()];
	}
public:
	int set_free_num;
	int old_fringe_num;
	int new_fringe_num;
	GBFringe(double r, double ve,set<int>& forb)
	{
		range = r;
		veps = ve;
		forbid = forb;
		GBF_ini();
		set_free_num = 0;
		old_fringe_num = 0;
		new_fringe_num = 0;
		AlphaRepr = -1;
		BetaRepr = -1;
	}
	void stat_ini()
	{
		set_free_num = 0;
		old_fringe_num = 0;
		new_fringe_num = 0;
	}
	bool is_forbid(Box* b, map<int, int>& GlobalId)
	{
		return forbid.find(gid(b, GlobalId)) != forbid.end();
	}
	bool is_AlphaFringe(Box* b)
	{
		return AlphaFringe.find(b->ID()) != AlphaFringe.end();
	}
	bool is_BetaFringe(Box* b)
	{
		return BetaFringe.find(b->ID()) != BetaFringe.end();
	}
	bool is_AlphaBox(Box* b)
	{
		return AlphaBox.find(b->ID()) != AlphaBox.end();
	}
	bool is_BetaBox(Box* b)
	{
		return BetaBox.find(b->ID()) != BetaBox.end();
	}
	bool is_AlphaNeighbor(Box* b)
	{
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
			if (is_AlphaBox(b_neighbors[j]))
				return true;
		return false;
	}
	bool is_BetaNeighbor(Box* b)
	{
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
			if (is_BetaBox(b_neighbors[j]))
				return true;
		return false;
	}
	void set_AlphaFringe(Box* b)
	{
		new_fringe_num++;
		AlphaFringe.insert(b->ID());
	}
	void set_BetaFringe(Box* b)
	{
		new_fringe_num++;
		BetaFringe.insert(b->ID());
	}
	void set_AlphaBox(Box* b)
	{
		set_free_num++;
		AlphaBox.insert(b->ID());
	}
	void set_BetaBox(Box* b)
	{
		set_free_num++;
		BetaBox.insert(b->ID());
	}

	// Hint heuristic value.
	double HintHeu(Box* b)
	{
		if (is_hinted(b))
			return 1.0;
		else
			return 0.0;
	}

	// Negative exp distance to alpha fringe.
	double nExpAlphaDis(Box* b)
	{
		if (AlphaFringe.empty())
			return 0;
		double min_dis = 2 * sqrt(6) * range;
		if (AlphaRepr < 0)
		{
			int PosRepr = -1;
			for (set<int>::iterator it = AlphaFringe.begin(); it != AlphaFringe.end(); ++it)
			{
				int dis = Sep(b, to_box(*it));
				if (min_dis > dis)
				{
					min_dis = dis;
					PosRepr = *it;
				}
			}
			AlphaRepr = PosRepr;
		}
		else
			min_dis = Sep(b, to_box(AlphaRepr));
		return exp(-min_dis);
	}

	// Negative exp distance to beta fringe.
	double nExpBetaDis(Box* b)
	{
		if (BetaFringe.empty())
			return 0;
		double min_dis = 2 * sqrt(6) * range;
		if (BetaRepr < 0)
		{
			int PosRepr = -1;
			for (set<int>::iterator it = BetaFringe.begin(); it != BetaFringe.end(); ++it)
			{
				int dis = Sep(b, to_box(*it));
				if (min_dis > dis)
				{
					min_dis = dis;
					PosRepr = *it;
				}
			}
			BetaRepr = PosRepr;
		}
		else
			min_dis = Sep(b, to_box(BetaRepr));
		return exp(-min_dis);
	}

	// Update b as alpha fringe.
	void set_AlphaMix(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, map<int, int>& GlobalId, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (b->width() <= veps)
			return;
		if (is_AlphaBox(b) || is_BetaBox(b))
			cout << "Warning! Mixed box is alpha/beta box!" << endl;
		if (is_forbid(b, GlobalId))
			return;
		if (Noisity >= 9)
		{
			cout << "Pushing ";
			b->out();
			cout << " into local queue." << endl;
		}
		LQ.push(make_pair(key(nExpBetaDis(b) + HintHeu(b)), b));
		if (Noisity >= 10)
			cout << "There are " << LQ.size() << " heuristic values in queue." << endl;
		if (is_AlphaFringe(b))
			old_fringe_num++;
		else
			set_AlphaFringe(b);
		auto end_time = std::chrono::high_resolution_clock::now();
		set_fringe_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	}
	// Update b as beta fringe.
	void set_BetaMix(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, map<int, int>& GlobalId, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (b->width() <= veps)
			return;
		if (is_AlphaBox(b) || is_BetaBox(b))
			cout << "Warning! Mixed box is alpha/beta box!" << endl;
		if (is_forbid(b, GlobalId))
			return;
		if (Noisity >= 9)
		{
			cout << "Pushing ";
			b->out();
			cout << " into local queue." << endl;
		}
		LQ.push(make_pair(key(nExpAlphaDis(b) + HintHeu(b)), b));
		if (Noisity >= 10)
			cout << "There are " << LQ.size() << " heuristic values in queue." << endl;
		if (is_BetaFringe(b))
			old_fringe_num++;
		else
			set_BetaFringe(b);
		auto end_time = std::chrono::high_resolution_clock::now();
		set_fringe_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	}
	// Update b as alpha box.
	void set_AlphaFree(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, Graph<Box>& Global, map<int, int>& GlobalId, int t, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (is_AlphaBox(b) || is_BetaBox(b))
			return;

		if (is_forbid(b, GlobalId))
			return;

		set_AlphaBox(b);

		if (Noisity >= 9)
		{
			cout << "Inserting ";
			b->out();
			cout << endl;
		}

		vector<Box*> b_neighbors = b->all_adj_neighbors();
		if (Noisity >= 10)
			cout << "There are " << b_neighbors.size() << " neighbors." << endl;
		auto end_time = std::chrono::high_resolution_clock::now();
		set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		set_box_preparation_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		start_time = std::chrono::high_resolution_clock::now();
		for (int j = 0; j < b_neighbors.size(); ++j)
		{
			Box* b_neighbor = b_neighbors[j];
			switch (getcolor(b_neighbor, Global, GlobalId,t))
			{
			case GREEN:
			{
				end_time = std::chrono::high_resolution_clock::now();
				set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
				set_AlphaFree(b_neighbor, LQ, Global, GlobalId, t, show);
				start_time = std::chrono::high_resolution_clock::now();
			}
			break;
			case YELLOW:
			{
				end_time = std::chrono::high_resolution_clock::now();
				set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
				set_AlphaMix(b_neighbor, LQ, GlobalId, show);
				start_time = std::chrono::high_resolution_clock::now();
			}
			break;
			case RED: cout << "Error! FREE boxes cannot be adjacent to STUCK boxes." << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			case BLACK: cout << "Error! Some boxes are classified as UNKNOWN!" << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			default: break;
			}
		}
		end_time = std::chrono::high_resolution_clock::now();
		set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	}
	// Update b as beta box.
	void set_BetaFree(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, Graph<Box>& Global, map<int, int>& GlobalId, int t, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (is_AlphaBox(b) || is_BetaBox(b))
			return;

		if (is_forbid(b, GlobalId))
			return;

		set_BetaBox(b);

		if (Noisity >= 9)
		{
			cout << "Inserting ";
			b->out();
			cout << endl;
		}

		vector<Box*> b_neighbors = b->all_adj_neighbors();
		if (Noisity >= 10)
			cout << "There are " << b_neighbors.size() << " neighbors." << endl;
		auto end_time = std::chrono::high_resolution_clock::now();
		set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		set_box_preparation_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		for (int j = 0; j < b_neighbors.size(); ++j)
		{
			Box* b_neighbor = b_neighbors[j];
			switch (getcolor(b_neighbor, Global, GlobalId,t))
			{
			case GREEN:
			{
				end_time = std::chrono::high_resolution_clock::now();
				set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
				set_BetaFree(b_neighbor, LQ, Global, GlobalId, t, show);
				start_time = std::chrono::high_resolution_clock::now();
			}
			break;
			case YELLOW:
			{
				end_time = std::chrono::high_resolution_clock::now();
				set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
				set_BetaMix(b_neighbor, LQ, GlobalId, show);
				start_time = std::chrono::high_resolution_clock::now();
			}
			break;
			case RED: cout << "Error! FREE boxes cannot be adjacent to STUCK boxes." << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			case BLACK: cout << "Error! Some boxes are classified as UNKNOWN!" << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			default: break;
			}
		}
		end_time = std::chrono::high_resolution_clock::now();
		set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	}

	// Initialization of GBF heuristic.
	void GBF_ini()
	{
		AlphaBox.clear();
		BetaBox.clear();
		AlphaFringe.clear();
		BetaFringe.clear();
	}

};

// Quickly updating along a path by heuristics.
template<typename Box>
class GBFQuickPath
{
	double range;
	double veps;
	// Alpha/beta boxes/fringes.
	set<int> BetaBox, BetaFringe;
	// Forbidden area.
	set<int> forbid;
	// List of boxes in the box space.
	vector<Box*> BoxSpace;
	// Map boxes to the t-th box in BoxSpace that containing them.
	map<int, int> BoxParent;

	// Find the global node of a box in the global graph, return NULL if not exists.
	GraphNode<Box>* gnode(Box* b, Graph<Box>& Global, map<int, int>& GlobalId)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (GlobalId.find(b->ID()) == GlobalId.end())
			return NULL;
		GraphNode<Box>* gv = Global.node(GlobalId[b->ID()]);
		auto end_time = std::chrono::high_resolution_clock::now();
		get_graph_node_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return gv;
	}
	// Get the t-color of a box b in global graph.
	Vcolor getcolor(Box* b, Graph<Box>& Global, map<int, int>& GlobalId, int t)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		Vcolor c = gnode(b, Global, GlobalId)->color(t);
		auto end_time = std::chrono::high_resolution_clock::now();
		get_color_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return c;
	}
	// Node id in global graph.
	int gid(Box* b, map<int, int>& GlobalId)
	{
		if (GlobalId.find(b->ID()) == GlobalId.end())
			return -1;
		return GlobalId[b->ID()];
	}
public:
	int set_free_num;
	int old_fringe_num;
	int new_fringe_num;
	GBFQuickPath(double r, double ve, set<int>& forb, vector<Box*>& BSpace)
	{
		range = r;
		veps = ve;
		forbid = forb;
		BoxSpace = BSpace;
		BoxParent.clear();
		for (int i = 0; i < BoxSpace.size(); ++i)
			BoxParent.insert(make_pair(BoxSpace[i]->ID(), i));
		GBF_ini();
		set_free_num = 0;
		old_fringe_num = 0;
		new_fringe_num = 0;
	}
	// Find the t-th box in BoxSpace that containing b.
	int tParent(Box* b)
	{
		if (BoxParent.find(b->ID()) == BoxParent.end())
			return -1;
		return BoxParent[b->ID()];
	}
	// Set that the t-th box in BoxSpace is containing b.
	void setParent(Box* b, int t)
	{
		BoxParent.insert(make_pair(b->ID(), t));
	}
	void stat_ini()
	{
		set_free_num = 0;
		old_fringe_num = 0;
		new_fringe_num = 0;
	}
	bool is_forbid(Box* b, map<int, int>& GlobalId)
	{
		return forbid.find(gid(b, GlobalId)) != forbid.end();
	}
	bool is_BetaFringe(Box* b)
	{
		return BetaFringe.find(b->ID()) != BetaFringe.end();
	}
	bool is_BetaBox(Box* b)
	{
		return BetaBox.find(b->ID()) != BetaBox.end();
	}
	bool is_BetaNeighbor(Box* b)
	{
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
			if (is_BetaBox(b_neighbors[j]))
				return true;
		return false;
	}
	void set_BetaFringe(Box* b)
	{
		new_fringe_num++;
		BetaFringe.insert(b->ID());
	}
	void set_BetaBox(Box* b)
	{
		set_free_num++;
		BetaBox.insert(b->ID());
	}

	// Hint heuristic value.
	double HintHeu(Box* b)
	{
		if (is_hinted(b))
			return 1.0;
		else
			return 0.0;
	}

	// Negative exp distance to alpha fringe.
	double nExpNextDis(Box* b)
	{
		int t = tParent(b);
		if (t < 0)
			return 0;
		else if (t == 0)
			return 1;
		else if (t > BoxSpace.size())
			return 0;
		else
			return exp(-Sep(b, BoxSpace[tParent(b) - 1]));
	}

	// Update b as beta fringe.
	void set_BetaMix(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, map<int, int>& GlobalId, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (b->width() <= veps)
			return;
		if (is_BetaBox(b))
			cout << "Warning! Mixed box is alpha/beta box!" << endl;
		if (is_forbid(b, GlobalId))
			return;
		if (Noisity >= 9)
		{
			cout << "Pushing ";
			b->out();
			cout << " into local queue." << endl;
		}
		LQ.push(make_pair(key(nExpNextDis(b) + HintHeu(b) + (BoxSpace.size() - tParent(b))), b));
		if (Noisity >= 10)
			cout << "There are " << LQ.size() << " heuristic values in queue." << endl;
		if (is_BetaFringe(b))
			old_fringe_num++;
		else
			set_BetaFringe(b);
		auto end_time = std::chrono::high_resolution_clock::now();
		set_fringe_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	}
	// Update b as beta box.
	void set_BetaFree(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, Graph<Box>& Global, map<int, int>& GlobalId, int t, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (is_BetaBox(b))
			return;

		if (is_forbid(b, GlobalId))
			return;

		set_BetaBox(b);

		if (Noisity >= 9)
		{
			cout << "Inserting ";
			b->out();
			cout << endl;
		}

		vector<Box*> b_neighbors = b->all_adj_neighbors();
		if (Noisity >= 10)
			cout << "There are " << b_neighbors.size() << " neighbors." << endl;
		auto end_time = std::chrono::high_resolution_clock::now();
		set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		set_box_preparation_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		for (int j = 0; j < b_neighbors.size(); ++j)
		{
			Box* b_neighbor = b_neighbors[j];
			switch (getcolor(b_neighbor, Global, GlobalId, t))
			{
			case GREEN:
			{
				end_time = std::chrono::high_resolution_clock::now();
				set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
				set_BetaFree(b_neighbor, LQ, Global, GlobalId, t, show);
				start_time = std::chrono::high_resolution_clock::now();
			}
			break;
			case YELLOW:
			{
				end_time = std::chrono::high_resolution_clock::now();
				set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
				set_BetaMix(b_neighbor, LQ, GlobalId, show);
				start_time = std::chrono::high_resolution_clock::now();
			}
			break;
			case RED: cout << "Error! FREE boxes cannot be adjacent to STUCK boxes." << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			case BLACK: cout << "Error! Some boxes are classified as UNKNOWN!" << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			default: break;
			}
		}
		end_time = std::chrono::high_resolution_clock::now();
		set_box_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
	}

	// Initialization of GBF heuristic.
	void GBF_ini()
	{
		BetaBox.clear();
		BetaFringe.clear();
	}

};


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
	// The global graph with colors.
	Graph<Box> Global;
	// Map from box id to free graph id.
	map<int, int> Gid;
	// Map from box id to global graph id.
	map<int, int> GlobalId;
	// List of forbid area under different recursion level.
	stack<set<int> > forbids;
	// The mixed queue.
	priority_queue<pair<vector<double>, Box*> > Q;
	// Varepsilon.
	double veps;
	// Result path.
	vector<Config> Path;
	// Alpha box.
	Box* BoxAlpha = NULL;
	// Beta Box.
	Box* BoxBeta = NULL;

	// Update statistic values for boxes pushed into Q.
	void mixed_static_update(Box* b)
	{
		stat_mixed++;
		avemixedfeature += (C->feature_of(b).size() - avemixedfeature) / stat_mixed;
	}
	// Update statistic values for boxes inserted into G.
	void free_static_update(Box* b)
	{
		stat_free++;
	}
	// Update statistic values for boxes that are stuck.
	void stuck_static_update(Box* b)
	{
		stat_stuck++;
	}
	// Update statistic values for boxes that are veps-small.
	void small_static_update(Box* b)
	{
		stat_small++;
	}

	// Distance heuristic.
	double dis_heu(Box* b)
	{
		return C->feature_dis(b);
	}
	// Width heuristic.
	double width_heu(Box* b)
	{
		return max(b->BT()->width(), b->BR()->width());
	}
	// Id heuristic.
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
		case GBF: heus.push_back(0); break;
		case DIS: heus.push_back(dis_heu(b)); break;
		case WIDTH: heus.push_back(width_heu(b)); break;
		case BYID: heus.push_back(id_heu(b)); break;
		case WIDIS: heus.push_back(width_heu(b)); heus.push_back(dis_heu(b)); break;
		case RECUR:heus.push_back(width_heu(b)); break;
		default: heus.push_back(rand());
		}
		return heus;
	}

	//********************** GBF Heuristic Algorithm *********************//

	// Alpha/beta boxes/fringes.
	set<int> AlphaBox, BetaBox, AlphaFringe, BetaFringe;
	bool is_AlphaFringe(Box* b)
	{
		return AlphaFringe.find(b->ID()) != AlphaFringe.end();
	}
	bool is_BetaFringe(Box* b)
	{
		return BetaFringe.find(b->ID()) != BetaFringe.end();
	}
	bool is_AlphaBox(Box* b)
	{
		return AlphaBox.find(b->ID()) != AlphaBox.end();
	}
	bool is_BetaBox(Box* b)
	{
		return BetaBox.find(b->ID()) != BetaBox.end();
	}
	bool is_AlphaNeighbor(Box* b)
	{
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
			if (is_AlphaBox(b_neighbors[j]))
				return true;
		return false;
	}
	bool is_BetaNeighbor(Box* b)
	{
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
			if (is_BetaBox(b_neighbors[j]))
				return true;
		return false;
	}
	bool is_Alpha(Box* b)
	{
		return b->contains(alpha);
	}
	bool is_Beta(Box* b)
	{
		return b->contains(beta);
	}
	void set_AlphaFringe(Box* b)
	{
		AlphaFringe.insert(b->ID());
	}
	void set_BetaFringe(Box* b)
	{
		BetaFringe.insert(b->ID());
	}
	void set_AlphaBox(Box* b)
	{
		AlphaBox.insert(b->ID());
	}
	void set_BetaBox(Box* b)
	{
		BetaBox.insert(b->ID());
	}

	// Negative exp distance to alpha fringe.
	double nExpAlphaDis(Box* b)
	{
		if (AlphaFringe.empty())
			return 0;
		double min_dis = 2 * sqrt(6) * B->Range();
		for (set<int>::iterator it = AlphaFringe.begin(); it != AlphaFringe.end(); ++it)
			min_dis = min(min_dis, Sep(b, to_box(*it)));
		return exp(-min_dis);
	}
	// Negative exp distance to beta fringe.
	double nExpBetaDis(Box* b)
	{
		if (BetaFringe.empty())
			return 0;
		double min_dis = 2 * sqrt(6) * B->Range();
		for (set<int>::iterator it = BetaFringe.begin(); it != BetaFringe.end(); ++it)
			min_dis = min(min_dis, Sep(b, to_box(*it)));
		return exp(-min_dis);
	}

	// Update b as alpha fringe.
	void set_AlphaMix(Box* b)
	{
		if (b->width() <= veps)
			return;
		// Q can update the priority of b in this way.
		Q.push(make_pair(key(nExpBetaDis(b)), b));
		if (!is_AlphaFringe(b))
			set_AlphaFringe(b);
		if (!is_BetaFringe(b))
			mixed_static_update(b);
	}
	// Update b as beta fringe.
	void set_BetaMix(Box* b)
	{
		if (b->width() <= veps)
			return;
		// Q can update the priority of b in this way.
		Q.push(make_pair(key(nExpAlphaDis(b)), b));
		if (!is_BetaFringe(b))
			set_BetaFringe(b);
		if (!is_AlphaFringe(b))
			mixed_static_update(b);
	}
	// Update b as alpha box.
	void set_AlphaFree(Box* b)
	{
		if (is_AlphaBox(b) || is_BetaBox(b))
			return;

		set_AlphaBox(b);

		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
		{
			Box* b_neighbor = b_neighbors[j];
			switch (C->classify(b_neighbor))
			{
			case FREE: set_AlphaFree(b_neighbor); break;
			case MIXED: set_AlphaMix(b_neighbor); break;
			case STUCK: cout << "Error! FREE boxes cannot be adjacent to STUCK boxes." << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			case UNKNOWN: cout << "Error! Some boxes are classified as UNKNOWN!" << endl << "Hint: "; b_neighbor->out(); cout << endl; break;
			default: break;
			}
		}
	}
	// Update b as beta box.
	void set_BetaFree(Box* b)
	{
		if (is_AlphaBox(b) || is_BetaBox(b))
			return;

		set_BetaBox(b);

		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
		{
			Box* b_neighbor = b_neighbors[j];
			switch (C->classify(b_neighbor))
			{
			case FREE: set_BetaFree(b_neighbor); break;
			case MIXED: set_BetaMix(b_neighbor); break;
			case STUCK: cout << "Error! FREE boxes cannot be adjacent to STUCK boxes."; break;
			case UNKNOWN: cout << "Error! Some boxes are classified as UNKNOWN!"; break;
			default: break;
			}
		}
	}
	
	// Initialization of GBF heuristic.
	void GBF_ini()
	{
		AlphaBox.clear();
		BetaBox.clear();
		AlphaFringe.clear();
		BetaFringe.clear();
	}

	//********************** Update boxes by heuristics *********************//

	// This method makes Q maintains all boxes.
	void add_mixed_pure(Box* b)
	{
		Q.push(make_pair(heuristic(b), b)); 
		mixed_static_update(b);
	}
	// This method makes Q maintains purely fringe boxes.
	void add_mixed_fringe(Box* b)
	{
		vector<Box*> target_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < target_neighbors.size(); ++j)
		{
			Box* tb = target_neighbors[j];
			switch (C->classify(tb))
			{
			case FREE: add_mixed_pure(b); return;
			case MIXED: break;
			case STUCK: break;
			case UNKNOWN: cout << "Error! Some boxes are classified as UNKNOWN!"; break;
			default: break;
			}
		}
	}
	// This method makes Q maintains purely alpha and beta fringe boxes.
	void add_mixed_GBF(Box* b)
	{
		if (is_AlphaNeighbor(b))
			set_AlphaMix(b);
		if (is_BetaNeighbor(b))
			set_BetaMix(b);
	}
	// This method makes LQ maintains purely local alpha and beta fringe boxes.
	void add_mixed_recur(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, GBFQuickPath<Box>& Fringe)
	{
		if (Fringe.is_BetaNeighbor(b))
			Fringe.set_BetaMix(b, LQ, GlobalId);
	}
	// Procedures for mixed boxes.
	void add_mixed_node(Box* b, bool show = false)
	{
		if (Noisity >= 8)
		{
			cout << "Adding MIXED box: ";
			b->out();
			cout << endl;
		}
		// Never push veps-small boxes.
		if (b->width() <= veps)
			return;
		switch (heu)
		{
		case RAND:add_mixed_pure(b); break;
		case GBF: add_mixed_GBF(b); break;
		case RECUR: mixed_static_update(b); break;
		default: add_mixed_fringe(b);
		}
	}
	// Update mixed node under local priority queue LQ.
	void add_local_mixed_node(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, bool show = false)
	{
		LQ.push(make_pair(heuristic(b), b));
	}

	// This method does nothing.
	void add_free_pure(Box* b)
	{
		free_static_update(b);
	}
	// This method helps Q to maintain purely fringe boxes.
	void add_free_fringe(Box* b, vector<Box*> b_neighbors)
	{
		for (int j = 0; j < b_neighbors.size(); ++j)
		{
			// This change makes Q maintains purely fringe boxes.
			Box* b_neighbor = b_neighbors[j];
			if (C->classify(b_neighbor) == MIXED)
				Q.push(make_pair(heuristic(b_neighbor), b_neighbor));
		}
		free_static_update(b);
	}
	// This method maintains Alpha/Beta Boxes and Alpha/Beta Fringes.
	void add_free_GBF(Box* b)
	{
		if (is_Alpha(b))
			set_AlphaFree(b);
		if (is_Beta(b))
			set_BetaFree(b);
		if (is_AlphaNeighbor(b))
			set_AlphaFree(b);
		if (is_BetaNeighbor(b))
			set_BetaFree(b);
		add_free_pure(b);
	}
	// This method maintains local Alpha/Beta Boxes and Alpha/Beta Fringes.
	void add_free_recur(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, GBFQuickPath<Box>& Fringe, int t, bool show = false)
	{
		if (Noisity >= 8)
		{
			cout << "Adding FREE box: ";
			b->out();
			cout << endl;
		}

		if (is_Beta(b))
		{
			if (Noisity >= 9)
				cout << "This is the BoxBeta, inserting into beta continent." << endl;
			Fringe.set_BetaFree(b, LQ, Global, GlobalId, t, show);
		}
		if (Fringe.is_BetaNeighbor(b))
		{
			if (Noisity >= 9)
				cout << "This is a neighbor of beta box, inserting into beta continent." << endl;
			Fringe.set_BetaFree(b, LQ, Global, GlobalId, t, show);
		}
		add_free_pure(b);
		if (Noisity >= 9)
			cout << "Time: " << clock() << endl;
	}
		// Procedures for free boxes.
	void add_free_node(Box* b, bool show = false)
	{
		if (Noisity >= 8)
		{
			cout << "Adding FREE box: ";
			b->out();
			cout << endl;
		}
		// Always maintain G when classifying some node as FREE!
		Gid.insert(make_pair(b->ID(), G.insert(b)->id));
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
		{
			GraphNode<Box>* b_neighbor_node = node(b_neighbors[j]);
			if (b_neighbor_node != NULL)
				G.link(node(b), b_neighbor_node);
		}
		viewer.add_box(b);

		// Do other operations by heuristic accordingly.
		switch (heu)
		{
		case RAND:add_free_pure(b); break;
		case GBF: add_free_GBF(b); break;
		case RECUR:free_static_update(b); break;
		default: add_free_fringe(b, b_neighbors);
		}
	}

	// Procedures for stuck boxes.
	void add_stuck_node(Box* b, bool show = false)
	{
		if (Noisity >= 8)
		{
			cout << "Adding STUCK box: ";
			b->out();
			cout << endl;
		}
		stuck_static_update(b);
	}

	// Procedures for veps-small boxes.
	void add_small_node(Box* b, bool show = false)
	{
		if (Noisity >= 8)
		{
			cout << "Adding veps-small box: ";
			b->out();
			cout << endl;
		}
		small_static_update(b);
	}

	// Procedures for unknown boxes.
	void add_unknown_node(Box* b, bool show = false)
	{
		cout << "Error! Some boxes are classified as UNKNOWN!" << endl << "Hint: ";
		b->out();
		cout << endl;
	}

	//********************** Maintain Graph structures *********************//

	// Find the node of a box in the free graph, return NULL if not exists.
	GraphNode<Box>* node(Box* b)
	{
		if (Gid.find(b->ID()) == Gid.end())
			return NULL;
		return G.node(Gid[b->ID()]);
	}

	// Find the global node of a box in the global graph, return NULL if not exists.
	GraphNode<Box>* gnode(Box* b)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		if (GlobalId.find(b->ID()) == GlobalId.end())
			return NULL;
		GraphNode<Box>* gv = Global.node(GlobalId[b->ID()]);
		auto end_time = std::chrono::high_resolution_clock::now();
		get_graph_node_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return gv;
	}

	// Find the local node of a box in the local graph, return NULL if not exists.
	GraphNode<Box>* lgnode(Box* b, Graph<Box>& LG, map<int, int>& LGid)
	{
		if (LGid.find(b->ID()) == LGid.end())
			return NULL;
		return LG.node(LGid[b->ID()]);
	}

	// Insert box b into global graph.
	void ginsert(Box* b)
	{
		GlobalId.insert(make_pair(b->ID(), Global.insert(b)->id));
	}

	// Remove box b from global graph.
	void gerase(Box* b)
	{
		Global.erase(gnode(b));
	}

	// Update neighbors of box b in global graph.
	void gupdate(Box* b)
	{
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
			Global.link(gnode(b), gnode(b_neighbors[j]));
	}

	// Add a color into box b in global graph.
	void gcolor(Box* b, Vcolor c)
	{
		Global.push_color(gnode(b), c);
	}

	// Get the t-th color of box b in global graph.
	Vcolor getcolor(Box* b, int t)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		Vcolor c = gnode(b)->color(t);
		auto end_time = std::chrono::high_resolution_clock::now();
		get_color_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return c;
	}

	// Insert box b into local graph LG.
	void lginsert(Box* b, Graph<Box>& LG, map<int, int>& LGid)
	{
		LGid.insert(make_pair(b->ID(), LG.insert(b)->id));
	}

	// Update neighbors of box b in local graph LG.
	void lgupdate(Box* b, Graph<Box>& LG, map<int, int>& LGid)
	{
		vector<Box*> b_neighbors = b->all_adj_neighbors();
		for (int j = 0; j < b_neighbors.size(); ++j)
			LG.link(lgnode(b, LG, LGid), lgnode(b_neighbors[j], LG, LGid));
	}


	// Find the class of a box in the free graph, return NULL if not exists.
	GraphNode<Box>* find(Box* b)
	{
		if (Gid.find(b->ID()) == Gid.end())
			return NULL;
		return G.find(G.node(Gid[b->ID()]));
	}
 
	//********************** Constructor *********************//

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
		veps = 1.0 / 18;
		MaxSize = B->Range();
		MinSize = B->Range();
	}

	//********************** Observer functions *********************//
	
	// Show how many boxes that are expanded with hints.
	void show_expansion(string hint, bool show = false)
	{
		if (Noisity < 0)
			return;
		cout << expanded << " cubes have been expanded with average width " << aver3w << " * " << aveso3w << "." << endl;
		cout << "There are " << stat_mixed << " mixed boxes, " << stat_free << " free boxes and " << stat_stuck << " stuck boxes and " << stat_small << " veps-small boxes." << endl;
		cout << "The average amount of features for mixed boxes are " << avemixedfeature << "." << endl;
		cout << "There are " << Global.Vsize() << " boxes in the global graph." << endl;
		cout << "There are " << Global.Esize() << " edges connecting boxes in the global graph." << endl;
		cout << "Finding neighbors costs " << (find_neighbor_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Finding paths costs " << (find_path_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Setting fringes costs " << (set_fringe_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Setting continents costs " << (set_box_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Preparing for setting continents costs " << (set_box_preparation_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Getting box colors costs " << (get_color_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Getting graph node costs " << (get_graph_node_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Repeat expansions wastes " << (waste_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Classifying boxes costs " << (classification_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Classifying boxes by predicate pO costs " << (classification_C1_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Constructing Bt solid costs " << (construct_box_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Check known case of pO costs " << (classify_p0_known_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Check quick classification of pO costs " << (p0_quick_classify_time / (long double)(1000000.0)) << "s in total till now." << endl;
		cout << "Classifying of pO costs " << (p0_classify_feature_time / (long double)(1000000.0)) << "s in total till now." << endl;
		//cout << "The average amount of features for free boxes are " << avefreefeature << "." << endl;
		if (hint.length() > 0)
			cout << hint << endl;
		for (int i = 0; i < 60; ++i)
			cout << "-";
		cout << endl;
		if (Noisity >= 5)
		{
			reset_viewer();
			if (forbids.empty())
				viewer.set_graph(Global);
			else
				viewer.set_graph(Global, forbids.top());
			viewer.view();
		}
	}

	// Update the statistic of expansion.
	void update_expansion(Box* b)
	{
		expanded++;
		aver3w += (b->BT()->width() - aver3w) / expanded;
		aveso3w += (b->BR()->width() - aveso3w) / expanded;
	}

	//********************** Expand functions *********************//

	// Expand(B)
	void Expand(Box* b, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		// Step 0: check if the box is leaf.
		if (!b->is_leaf())
		{
			auto end_time = std::chrono::high_resolution_clock::now();
			waste_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return;
		}

		if (Noisity >= 7)
			cout << "Box is a leaf." << endl;

		// Step 1: subdivide the box.
		b->subdivide();
		gerase(b);

		if (Noisity >= 7)
			cout << "Box is subdivided." << endl;

		for (int i = 0; i < B->subsize(b); ++i)
			ginsert(b->child(i));

		if (Noisity >= 7)
			cout << "Children are inserted, now " << Global.Vsize() << " boxes in global graph." << endl;

		for (int i = 0; i < B->subsize(b); ++i)
			gupdate(b->child(i));

		if (Noisity >= 7)
			cout << "Box subdivision is set." << endl;

		// Step 2: classify soft pvalues and maintain global graph.
		for (int i = 0; i < B->subsize(b); ++i)
			for (int j = 0; j < C->psize(); ++j)
				if (b->child(i)->width() <= veps)
					gcolor(b->child(i), GREY);
				else
					switch (C->classify(b->child(i), j))
					{
					case FREE:gcolor(b->child(i), GREEN); break;
					case MIXED:gcolor(b->child(i), YELLOW); break;
					case STUCK:gcolor(b->child(i), RED); break;
					default:gcolor(b->child(i), BLACK);
					}

		if (Noisity >= 7)
			cout << "Children colors are set." << endl;

		// Step 3: maintain Q and G.
		for (int i = 0; i < B->subsize(b); ++i)
			switch (getcolor(b->child(i), 0))
			{
			case GREEN:add_free_node(b->child(i)); break;
			case YELLOW:add_mixed_node(b->child(i)); break;
			case RED:add_stuck_node(b->child(i)); break;
			case GREY:add_small_node(b->child(i)); break;
			default:add_unknown_node(b->child(i));
			}

		if (Noisity >= 7)
			cout << "Free graph is set." << endl;

		// Step 4: show how many cubes that are expanded to give an intuition of the process.
		if (ExpandShow == 1 || expanded % ExpandShow == 0)
			show_expansion("", show);
		update_expansion(b);
	}

	// Expand(B) under a local graph.
	void Expand_recur(Box* b, priority_queue<pair<vector<double>, Box*> >& LQ, GBFQuickPath<Box>& Fringe, int t, bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		// Step 0: check if the box is leaf avoiding forbidden area.
		if (!b->is_leaf())
		{
			auto end_time = std::chrono::high_resolution_clock::now();
			waste_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return;
		}
		if (forbids.top().find(gnode(b)->id) != forbids.top().end())
		{
			auto end_time = std::chrono::high_resolution_clock::now();
			waste_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return;
		}
		
		// No need to split this box anymore, since there is no possible R-direction in this box for the Delta robot to pass through anymore!
		/*if (t == 1 && b->BT()->width() <= 0.25)
		{
			auto end_time = std::chrono::high_resolution_clock::now();
			waste_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
			return;
		}*/

		// Step 1: subdivide the box.
		/*if (t == 1)
			b->subdivide(1);
		else*/
			b->subdivide();
		gerase(b);

		for (int i = 0; i < B->subsize(b); ++i)
			ginsert(b->child(i));
		for (int i = 0; i < B->subsize(b); ++i)
			gupdate(b->child(i));
		for (int i = 0; i < B->subsize(b); ++i)
			Fringe.setParent(b->child(i), Fringe.tParent(b));

		start_time = std::chrono::high_resolution_clock::now();
		// Step 2: classify soft pvalues and maintain global graph.
		for (int i = 0; i < B->subsize(b); ++i)
			for (int j = 0; j < C->psize(); ++j)
				if (b->child(i)->width() <= veps)
					gcolor(b->child(i), GREY);
				else if (j == 1 && b->child(i)->BT()->width() <= 0.125)
					gcolor(b->child(i), GREY);
				else
					switch (C->classify(b->child(i), j))
					{
					case FREE:gcolor(b->child(i), GREEN); break;
					case MIXED:gcolor(b->child(i), YELLOW); break;
					case STUCK:gcolor(b->child(i), RED); break;
					default:gcolor(b->child(i), BLACK);
					}
		auto end_time = std::chrono::high_resolution_clock::now();
		classification_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

		// Step 3: maintain LQ and G.
		for (int i = 0; i < B->subsize(b); ++i)
		{
			switch (getcolor(b->child(i), 0))
			{
			case GREEN:add_free_node(b->child(i)); break;
			case YELLOW:add_mixed_node(b->child(i)); break;
			case RED:add_stuck_node(b->child(i)); break;
			case GREY:add_small_node(b->child(i)); break;
			default:add_unknown_node(b->child(i));
			}
		}
		if (Noisity >= 7)
			cout << "G is set." << endl;

		for (int i = 0; i < B->subsize(b); ++i)
		{
			switch (getcolor(b->child(i), t))
			{
			case GREEN: add_free_recur(b->child(i), LQ, Fringe, t, show); break;
			case YELLOW: add_mixed_recur(b->child(i), LQ, Fringe); break;
			default:;
			}
		}

		if (Noisity >= 7)
		{
			cout << "LQ is set, " << Fringe.set_free_num << " boxes are set as alpha/beta box, ";
			cout << Fringe.old_fringe_num << " boxes are alpha/beta fringe and updated this time, ";
			cout << Fringe.new_fringe_num << " boxes are newly set as alpha/beta fringe." << endl;
		}

		// Step 4: show how many cubes that are expanded to give an intuition of the process.
		if (ExpandShow == 1 || expanded % ExpandShow == 0)
		{
			string show_hint = "";
			show_hint = show_hint + "There are " + to_string(Fringe.set_free_num) + " alpha/beta continent boxes in level t=" + to_string(t) + ".\n";
			show_hint = show_hint + "There are " + to_string(Fringe.new_fringe_num) + " newly alpha/beta fringe boxes in level t=" + to_string(t) + ".\n";
			show_hint = show_hint + "There are " + to_string(Fringe.old_fringe_num) + " updated oldly alpha/beta fringe boxes in level t=" + to_string(t) + ".\n";
			show_expansion(show_hint, show);
		}
		update_expansion(b);
	}

	//************************* General SSS framework ***********************//

	// Turn graph node path into box channel.
	vector<Box*> make_channel(list<GraphNode<Box>*>& free_path)
	{
		vector<Box*> channel;
		for (auto it = free_path.begin(); it != free_path.end(); ++it)
			channel.push_back((*it)->content);
		return channel;
	}

	// Current leaves.
	vector<Box*> current_leaves()
	{
		return Global.current_obj();
	}

	// Find the channel when G.find(BoxAlpha) == G.find(BoxBeta) != NULL.
	vector<Box*> Find_Channel()
	{
		auto lpath = G.path(node(BoxAlpha), node(BoxBeta));
		return make_channel(lpath);
	}

	// Step 0: initializations.
	bool SSS_ini(Config _alpha, Config _beta, FeatureSet Omega, double varepsilon, heutype _h = RAND, bool show = false)
	{
		Path.clear();
		BoxAlpha = NULL;
		BoxBeta = NULL;
		expanded = 0;
		aver3w = 0;
		aveso3w = 0;
		stat_mixed = 0;
		stat_free = 0;
		stat_stuck = 0;
		stat_small = 0;
		avemixedfeature = 0;
		avefreefeature = 0;
		MaxSize = B->Range();
		MinSize = B->Range();
		G.clear();
		Global.clear();
		Gid.clear();
		GlobalId.clear();
		Q = priority_queue<pair<vector<double>, Box*> >();
		alpha = _alpha;
		beta = _beta;
		C->set_feature(Omega);
		heu = _h;
		veps = varepsilon;

		if (heu == GBF)
			GBF_ini();

		if (Noisity >= 6)
			cout << "Find path initialization finished." << endl;
		return true;
	}

	// Step 1: classify initial boxes.
	bool SSS_first_classify(bool show = false)
	{
		ginsert(B->Root());
		for (int j = 0; j < C->psize(); ++j)
			switch (C->classify(B->Root(), j))
			{
			case FREE:gcolor(B->Root(), GREEN); break;
			case MIXED:gcolor(B->Root(), YELLOW); break;
			case STUCK:gcolor(B->Root(), RED); break;
			default:gcolor(B->Root(), BLACK);
			}
		if (C->classify(B->Root()) == FREE)
			add_free_node(B->Root());
		if (C->classify(B->Root()) == STUCK)
		{
			show_expansion("The environment is stuck.");
			return false;													// Empty if no-path.
		}

		if (Noisity >= 6)
			cout << "Classify initial boxes finished." << endl;
		return true;
	}

	// Step 2: find the FREE box containing alpha.
	bool SSS_find_alpha(bool show = false)
	{
		if (Noisity >= 6)
			cout << "Finding alpha: (" << alpha.transpose() << ")." << endl;
		BoxAlpha = B->find(alpha);
		if (Noisity >= 6)
			cout << "Alpha initial box is found." << endl;
		if (BoxAlpha == NULL)
		{
			show_expansion("Configuration alpha not in environment.");
			return false;
		}
		while (BoxAlpha->width() > varepsilon)
		{
			if (C->classify(BoxAlpha) == FREE)
				break;
			if (C->classify(BoxAlpha) == STUCK)
			{
				show_expansion("Configuration alpha is stuck.");
				return false;
			}
			Expand(BoxAlpha, show);
			BoxAlpha = B->find(BoxAlpha, alpha, show);
		}
		if (C->classify(BoxAlpha) != FREE)
		{
			show_expansion("Configuration alpha is not free.");
			return false;
		}

		if (show)
			show_expansion("Found the FREE box containing alpha.");
		return true;
	}

	// Step 3: find the FREE box containing beta.
	bool SSS_find_beta(bool show = false)
	{
		if (Noisity >= 6)
			cout << "Finding beta: (" << beta.transpose() << ")." << endl;
		BoxBeta = B->find(beta);
		if (BoxBeta == NULL)
		{
			show_expansion("Configuration beta not in environment.");
			return false;
		}
		while (BoxBeta->width() > varepsilon)
		{
			if (C->classify(BoxBeta) == FREE)
				break;
			if (C->classify(BoxBeta) == STUCK)
			{
				show_expansion("Configuration beta is stuck.");
				return false;
			}
			Expand(BoxBeta, show);
			BoxBeta = B->find(BoxBeta, beta, show);
		}
		if (C->classify(BoxBeta) != FREE)
		{
			show_expansion("Configuration beta is not free.");
			return false;
		}

		show_expansion("Found the FREE box containing beta.");
		return true;
	}

	// One step of expansion.
	int SSS_one_step(bool show = false)
	{
		if (find(BoxAlpha) == find(BoxBeta))
		{
			if (Noisity >= 6)
				cout << "Alpha and beta have been in the same FREE component." << endl;
			return 1;
		}
		if (Q.empty())
		{
			show_expansion("No path found.");
			return -1;
		}
		Box* Q_top = Q.top().second;
		Expand(Q_top, show);
		Q.pop();
		if (expanded > ExpandLimit)
		{
			show_expansion("Run over expand limit (temporary setting for experiments).");
			if (Noisity >= 6)
				cout << "Time Out!" << endl;
			return -1;
		}
		return 0;
	}

	// Step 4: expand Q.GetNext() until Box(alpha) and Box(beta) are in the same component.
	bool SSS_main_loop(bool show = false)
	{
		if (heu == RECUR)
			return true;
		while (true)
		{
			switch (SSS_one_step(show))
			{
			case 1:return true;
			case -1:return false;
			default:continue;
			}
		}
	}

	// Step 5: find the canonical path in the component of G.find(BoxAlph).
	bool SSS_discrete_find(bool show = false)
	{
		if (heu == RECUR)
			return Recur_discrete_find(show);
		vector<Box*> channel= Find_Channel();
		for (int i = 0; i < channel.size(); ++i)
			Path.push_back(B->center(channel[i]));
		show_expansion("Canonical path found.");
		return true;
	}

	// **************** Recursively updating SSS framework *********************** //

	// Step 1: Recursive method initialization, push C_t-mixed boxes into local queue and build local graph.
	void Recur_ini(vector<Box*> BoxSpace, priority_queue<pair<vector<double>, Box*> >& LQ, GBFQuickPath<Box>& Fringe, int t, bool show = false)
	{
		if (Noisity >= 4)
		{
			cout << "Now drawing the box space for t = " << t << "." << endl;
			reset_viewer();
			if (Noisity >= 6)
			{
				cout << "BoxAlpha: " << BoxAlpha->ID() << endl;
				cout << "BoxBeta: " << BoxBeta->ID() << endl;
				cout << "BoxSpace contains:" << endl;
			}
			for (auto it = BoxSpace.begin(); it != BoxSpace.end(); ++it)
			{
				if (Noisity >= 6)
					cout << " " << (*it)->ID();;
				if ((*it)->ID() == BoxAlpha->ID())
					continue;
				if ((*it)->ID() == BoxBeta->ID())
					continue;
				switch (getcolor(*it, t))
				{
				case GREEN: viewer.add_box(*it, Vector3d(0, 1, 0)); break;
				case YELLOW: viewer.add_box(*it, Vector3d(1, 1, 0)); break;
				case RED:viewer.add_box(*it, Vector3d(1, 0, 0)); break;
				case BLACK: viewer.add_box(*it, Vector3d(0, 0, 0)); break;
				case GREY: viewer.add_box(*it, Vector3d(0.5, 0.5, 0.5)); break;
				default:continue;
				}
			}
			if (Noisity >= 6)
				cout << endl;
			viewer.add_box(BoxAlpha, Vector3d(0, 0, 1));
			viewer.add_box(BoxBeta, Vector3d(0, 0, 1));
			viewer.view();
		}
		Fringe.set_BetaFree(BoxBeta, LQ, Global, GlobalId, t);
		if (getcolor(BoxAlpha, t) != GREEN)
			cout << "Warning! BoxAlpha is not free for t=" << t << "!" << endl;
		if (getcolor(BoxBeta, t) != GREEN)
			cout << "Warning! BoxBeta is not free for t=" << t << "!" << endl;
	}

	// One step of recursive expansion.
	pair<vector<Box*>, bool> Recur_one_step(priority_queue<pair<vector<double>, Box*> >& LQ, GBFQuickPath<Box>& Fringe, int t, bool show = false)
	{
		if (Noisity >= 6)
			cout << "Attempting to find a path for t=" << t << "." << endl;
		
		// This if is currently useless.
		if (!Global.connected(gnode(BoxAlpha), gnode(BoxBeta), forbids.top(), t + 1))
		{
			show_expansion("No possible YG-path for t=" + to_string(t) + ", return to last level of recursion.");
			return make_pair(vector<Box*>(), false);
		}

		if (t == 0)
		{
			if (G.quick_connected(node(BoxAlpha), node(BoxBeta)))
			{
				auto lpath = G.path(node(BoxAlpha), node(BoxBeta));
				return make_pair(make_channel(lpath), false);
			}
		}
		else
		{
			list<GraphNode<Box>*> try_path = Global.path(gnode(BoxAlpha), gnode(BoxBeta), forbids.top(), t);
			if (Noisity >= 6)
				cout << "Attempt find path finished." << endl;
			if (!try_path.empty())
			{
				vector<Box*> next_boxspace = make_channel(try_path);
				set<int> new_forbid = Global.Vlist;
				for (int i = 0; i < next_boxspace.size(); ++i)
					new_forbid.erase(gnode(next_boxspace[i])->id);
				forbids.push(new_forbid);
				show_expansion("Found a possible path for t=" + to_string(t) + " ,recursively calling find channel for t=" + to_string(t - 1));
				vector<Box*> new_path = RecurFindChannel(next_boxspace, t - 1, show);
				forbids.pop();
				if (!new_path.empty())
					return make_pair(new_path, false);
				else
					show_expansion("Possible path for t=" + to_string(t) + " fails, continue try for new path.");
			}
		}
		if (Noisity >= 6)
			cout << "Checked no current path." << endl;
		if (Noisity >= 6 && !LQ.empty() && !LQ.top().second->is_leaf())
			cout << "Getting rid of non-leaf boxes. There are " << LQ.size() << " boxes in local queue." << endl;
		while (!LQ.empty() && !LQ.top().second->is_leaf())
			LQ.pop();
		if (LQ.empty())
		{
			if (Noisity >= 3)
				cout << "This t = " << t << " channel has no path, return to last level to find another path." << endl;
			return make_pair(vector<Box*>(), false);
		}
		if (Noisity >= 6)
			cout << "Continue expansions." << endl;		
		Box* LQ_top = LQ.top().second;
		Expand_recur(LQ_top, LQ, Fringe, t, show);
		LQ.pop();
		if (Noisity >= 6)
			cout << "Expand finished." << endl << "-----------------------------------------------" << endl;
		if (expanded > ExpandLimit)
		{
			show_expansion("Run over expand limit (temporary setting for experiments).");
			if (Noisity >= 6)
				cout << "Time Out!" << endl;
			return make_pair(vector<Box*>(), false);
		}
		return make_pair(vector<Box*>(), true);
	}

	// Step 2: Recursive method main loop.
	vector<Box*> Recur_main_loop(priority_queue<pair<vector<double>, Box*> >& LQ, GBFQuickPath<Box>& Fringe, int t, bool show = false)
	{
		while (true)
		{
			pair<vector<Box*>, bool> One_step_result = Recur_one_step(LQ, Fringe, t, show);
			if (One_step_result.second)
				continue;
			else
				return One_step_result.first;
		}
	}

	// Recursively find channel algorithm for fixed BoxAlpha and BoxBeta in a given BoxSpace (substitution for step 4 in main algorithm).
	vector<Box*> RecurFindChannel(vector<Box*> BoxSpace, int t, bool show = false)
	{
		if (Noisity >= 6)
			cout << "Finding channel in " << BoxSpace.size() << " boxes." << endl;
		priority_queue<pair<vector<double>, Box*> > LocalQ;
		GBFQuickPath<Box> LocalFringe(B->Range(), veps, forbids.top(), BoxSpace);
		Recur_ini(BoxSpace, LocalQ, LocalFringe, t, show);
		if (Noisity >= 7)
			cout << "There are " << LocalQ.size() << " mixed boxes." << endl;
		show_expansion("Now finding channel for t = " + to_string(t));
		return Recur_main_loop(LocalQ, LocalFringe, t, show);
	}

	// If a box is a Voronoid box.
	bool isVor(Box* b)
	{
		return true;
		if (b == B->Root())
			return true;
		return C->feature_of(b).size() < C->feature_of(b->Parent()).size();
	}

	// Method to find an initial channel.
	vector<Box*> RecurIniChannel(bool show = false)
	{
		set<int> new_forbid;
		vector<Box*> leaves = current_leaves();
		if (Noisity >= 4)
			reset_viewer();
		if (Noisity >= 4)
			cout << "Drawing the initial environment." << endl;
		for (int i = 0; i < leaves.size(); ++i)
			if (!isVor(leaves[i]))
				new_forbid.insert(leaves[i]->ID());
			else if (Noisity >= 4)
				switch (getcolor(leaves[i], C->psize() - 1))
				{
				case GREEN: viewer.add_box(leaves[i], Vector3d(0, 1, 0)); break;
				case YELLOW: viewer.add_box(leaves[i], Vector3d(1, 1, 0)); break;
				case RED:viewer.add_box(leaves[i], Vector3d(1, 0, 0)); break;
				case BLACK: viewer.add_box(leaves[i], Vector3d(0, 0, 0)); break;
				case GREY: viewer.add_box(leaves[i], Vector3d(0.5, 0.5, 0.5)); break;
				default:continue;
				}
		if (Noisity >= 7)
			viewer.view();
		if (Noisity >= 4)
			cout << "Now finding an initial path." << endl;
		auto lpath = Global.path(gnode(BoxAlpha), gnode(BoxBeta), new_forbid, C->psize());
		if (Noisity >= 4 && lpath.empty())
		{
			cout << "Initial path not found." << endl;
			viewer.view();
		}
		return make_channel(lpath);
	}

	// Step 3: return the canonical path in the recursively found channel.
	bool Recur_discrete_find(bool show = false)
	{
		vector<Box*> channel;
		while (expanded <= ExpandLimit)
		{
			vector<Box*> IniChannel = RecurIniChannel(show);
			if (IniChannel.empty())
			{
				if (Noisity >= 4)
					cout << "Initial channel not found." << endl;
				break;
			}
			set<int> new_forbid = Global.Vlist;
			for (int i = 0; i < IniChannel.size(); ++i)
				new_forbid.erase(gnode(IniChannel[i])->id);
			forbids.push(new_forbid);
			channel = RecurFindChannel(IniChannel, C->psize() - 1, show);
			forbids.pop();
			if (!channel.empty())
				break;
		}
		for (int i = 0; i < channel.size(); ++i)
			Path.push_back(B->center(channel[i]));
		return true;
	}

	// Counting box size after find path.
	void BoxSizeCount()
	{
		MaxSize = 0;
		MinSize = B->Range();
		for (auto it = Global.Vlist.begin(); it != Global.Vlist.end(); ++it)
		{
			double Size = Global.node(*it)->content->width();
			if (MaxSize < Size)
				MaxSize = Size;
			if (MinSize > Size)
				MinSize = Size;
		}
	}

	// ******************************** Main function ****************************** //

	// SSS framework main process.
	vector<Config> Find_Path(Config _alpha, Config _beta, FeatureSet Omega, double varepsilon, heutype _h = RAND, bool show = false)
	{
		if (!SSS_ini(_alpha, _beta, Omega, varepsilon, _h, show))
			return Path;
		if (Noisity >= 4)
			cout << "SSS_ini finished!" << endl;
		if (!SSS_first_classify(show))
			return Path;
		if (Noisity >= 4)
			cout << "SSS_first_classify finished!" << endl;
		if (!SSS_find_alpha(show))
			return Path;
		if (Noisity >= 4)
			cout << "SSS_find_alpha finished!" << endl;
		if (!SSS_find_beta(show))
			return Path;
		if (Noisity >= 4)
			cout << "SSS_find_beta finished!" << endl;
		if (!SSS_main_loop(show))
			return Path;
		if (Noisity >= 4)
			cout << "SSS_main_loop finished!" << endl;
		if (!SSS_discrete_find(show))
			return Path;
		if (Noisity >= 4)
			cout << "SSS_discrete_find finished!" << endl;
		show_expansion("Find path algorithm finished!");
		BoxSizeCount();
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