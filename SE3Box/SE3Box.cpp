// SE3Box.cpp : This file implements the oct-tree structure for R^3 and SO(3)
// and implements the product tree structure for SE(3).

#ifndef SE3BOX_H
#define SE3BOX_H

#define MAXSHOW 2
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <interval.h>
#include <bimap.h>
#include <chrono>
using namespace std;

class R3Box;
class SO3Box;
class SE3Box;
double extern r0;
double extern varepsilon;
int extern Noisity;

vector<R3Box*> R3list;
vector<SO3Box*> SO3list;
vector<SE3Box*> SE3list;
bimap<int, int, int> SE3table;

long long find_neighbor_time = 0;

// Binary string for the tree code, 0 for \bar{1}, 1 for 1. We add a "1" at the beginning to avoid starting with "0".
class bit
{
#define EMP "E"
	int n;
public:
	// Empty string (possibly roots).
	bit()
	{
		n = 1;
	}
	// Non-empty string.
	bit(int _n)
	{
		n = _n;
	}
	// The number.
	int num()
	{
		return n;
	}
	// If this is an empty bit.
	bool is_empty()
	{
		return n < 2;
	}
	// Length of the string.
	int length()
	{
		return int(log2(n));
	}
	// The i-th position of the bit.
	int operator[](int i)
	{
		return (n >> (length() - i)) & 1;
	}
	// If this is the positive boundary bit.
	bool is_pos_boundary()
	{
		return !bool((n + 1) & n);
	}
	// If this is the negative boundary bit.
	bool is_neg_boundary()
	{
		return !bool((n - 1) & n);
	}
	// If this is the boundary bit.
	bool is_boundary()
	{
		return is_pos_boundary() || is_neg_boundary();
	}
	// If this is the boundary bit of pos.
	bool is_boundary(bool pos)
	{
		if (pos)
			return is_pos_boundary();
		else
			return is_neg_boundary();
	}
	// Positive child for this bit.
	bit pos()
	{
		return bit((n << 1) + 1);
	}
	// Negative child for this bit.
	bit neg()
	{
		return bit(n << 1);
	}
	// Flip of the bit.
	bit flip()
	{
		return bit(n ^ 1);
	}
	// Bar of the bit.
	bit bar()
	{
		return bit(n ^ ((1 << length()) - 1));
	}
	// Adj of the bit.
	bit adj()
	{
		if (n & 1)
			return bit(n + 1);
		else
			return bit(n - 1);
	}
	// Positive direction of the bit.
	bit to_pos()
	{
		return bit(n + 1);
	}
	// Negative direction of the bit.
	bit to_neg()
	{
		return bit(n - 1);
	}
	// The last bit in the positive direction.
	bit pos_most()
	{
		return bit((1 << (length() + 1)) - 1);
	}
	// The last bit in the negative direction.
	bit neg_most()
	{
		return bit(1 << length());
	}
	// If two bits are neighbor.
	bool neighbor_to(bit b)
	{
		if (length() > b.length())
			return b.neighbor_to(*this);
		int diff = b.length() - length();
		if (((n + 1) << diff) == b.n)
			return true;
		if ((n << diff) == (b.n + 1))
			return true;
		return false;
	}
	// If the bit contains another bit.
	bool contain(bit b)
	{
		if (length() > b.length())
			return false;
		int diff = b.length() - length();
		if (((n + 1) << diff) <= b.n)
			return false;
		if ((n << diff) >= (b.n + 1))
			return false;
		return true;
	}
	// If the bit is contained in the other.
	bool contained(bit b)
	{
		return b.contain(*this);
	}
	// If the bit is separated to another bit.
	bool separate(bit b)
	{
		if (length() > b.length())
			return b.separate(*this);
		int diff = b.length() - length();
		if (((n + 1) << diff) < b.n)
			return true;
		if ((n << diff) > (b.n + 1))
			return true;
		return false;
	}
	// Output the bit.
	void out(ostream& os)
	{
		if (is_empty())
			os << EMP;
		for (int i = length() - 1; i >= 0; --i)
			os << ((n >> i) & 1);
	}
#undef EMP
};
// Output a bit.
ostream& operator<<(ostream& os, bit t)
{
	t.out(os);
	return os;
}

// Get the t-th binary bit of n.
bool getBin(int n, int t)
{
	return bool((n >> t) & 1);
}

class R3Box
{
#define subsize 8
#define dim 3
	// Id number.
	int id;
	// Root node.
	R3Box* root;
	// Parent node.
	R3Box* parent;
	// Depth of the node.
	int depth;
	// Is this a leaf?
	bool leaf;
	// Children node.
	R3Box* children[subsize];
	// Neighbor node, 0 for neg, 1 for pos.
	R3Box* neighbors[dim][2];
	// Path indicator of the box.
	bit codes[dim];
public:
	// Numerical range of the box.
	MatrixId range;
	// Id number.
	int ID()
	{
		return id;
	}
	// Width of the box.
	double width()
	{
		return (range.max() - range.min()).minCoeff();
	}
	// Construct a new box.
	R3Box(MatrixId r)
	{
		range = r;
		id = R3list.size();
		R3list.push_back(this);
		root = this;
		parent = NULL;
		depth = 0;
		leaf = true;
		for (int i = 0; i < subsize; ++i)
			children[i] = NULL;
		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < 2; ++j)
				neighbors[i][j] = NULL;
		for (int i = 0; i < dim; ++i)
			codes[i] = bit();
	}
	// Root node.
	R3Box* Root()
	{
		return root;
	}
	// Set the root node.
	void set_root(R3Box* B)
	{
		root = B;
	}
	// Parent node.
	R3Box* Parent()
	{
		return parent;
	}
	// Set the parent node.
	void set_parent(R3Box* B)
	{
		parent = B;
	}
	// Depth of the node.
	int Depth()
	{
		return depth;
	}
	// Set the depth of the node.
	void set_depth(int d)
	{
		depth = d;
	}
	// Is this a leaf?
	bool Leaf()
	{
		return leaf;
	}
	// Child node.
	R3Box* child(int i)
	{
		return children[i % subsize];
	}
	// Child node by child indicator.
	R3Box* child(int t[dim])
	{
		int target = 0;
		for (int i = 0; i < dim; ++i)
			target += (t[i] % 2) << i;
		return child(target);
	}
	// Neighbor node.
	R3Box* neighbor(int i, bool pos)
	{
		if (pos)
			return neighbors[i % dim][1];
		else
			return neighbors[i % dim][0];
	}
	// Positive neighbor node.
	R3Box* pos_neighbor(int i)
	{
		return neighbors[i % dim][1];
	}
	// Negative neighbor node.
	R3Box* neg_neighbor(int i)
	{
		return neighbors[i % dim][0];
	}
	// Set neighbor node.
	void set_neighbor(R3Box* B, int i, bool pos)
	{
		if (pos)
			neighbors[i % dim][1] = B;
		else
			neighbors[i % dim][0] = B;
	}
	// Set positive neighbor node.
	void set_pos_neighbor(R3Box* B, int i)
	{
		neighbors[i % dim][1] = B;
	}
	// Set negative neighbor node.
	void set_neg_neighbor(R3Box* B, int i)
	{
		neighbors[i % dim][0] = B;
	}
	// The i-th bit of the node.
	bit code(int i)
	{
		return codes[i % dim];
	}
	// Set the i-th bit of the node.
	void set_code(bit t, int i)
	{
		codes[i % dim] = t;
	}
	// Find a target box by path indicator.
	R3Box* find(bit t[dim])
	{
		R3Box* target = root;
		while (!target->Leaf() && target->Depth() < t[0].length())
		{
			int tcode[dim];
			for (int i = 0; i < dim; ++i)
				tcode[i] = t[i][target->Depth() + 1];
			target = target->child(tcode);
		}
		return target;
	}
	// Find a target box by child indicator.
	R3Box* find(bit t[dim], R3Box* b)
	{
		R3Box* target = b;
		while (!target->Leaf() && target->Depth() < t[0].length())
		{
			int tcode[dim];
			for (int i = 0; i < dim; ++i)
				tcode[i] = t[i][target->Depth() + 1];
			target = target->child(tcode);
		}
		return target;
	}
	// Subdivision.
	void subdivide()
	{
		// Step 0: Only leaves can be subdivided.
		if (!leaf)
			return;

		// Step 1: Create subdivisions.
		// All boxes are sequenced by the lex-inequality of the range.min() coordinates.

		for (int i = 0; i < subsize; ++i)
		{
			Matrix<double, dim, 1> new_inf, new_sup;
			for (int j = 0; j < dim; ++j)
				if (getBin(i, j))
				{
					new_inf(j) = (range(j).min() + range(j).max()) / 2;
					new_sup(j) = range(j).max();
				}
				else
				{
					new_inf(j) = range(j).min();
					new_sup(j) = (range(j).min() + range(j).max()) / 2;
				}
			MatrixId new_range(new_inf, new_sup, false);
			children[i] = new R3Box(new_range);
			children[i]->set_root(root);
			children[i]->set_parent(this);
			children[i]->set_depth(depth + 1);
			for (int j = 0; j < dim; ++j)
				if (getBin(i, j))
					children[i]->set_code(codes[j].pos(), j);
				else
					children[i]->set_code(codes[j].neg(), j);
		}

		// Step 2: This is not a leaf anymore.
		leaf = false;

		// Step 3: Assign neighbors for children and re-new_child neighbors for neighbors and their children.

		for (int i = 0; i < subsize; ++i)
		{
			R3Box* new_child = children[i];
			R3Box* target = NULL;
			bit tcode[dim];
			for (int j = 0; j < dim; ++j)
			{
				// If this is not the negative boundary of j-th direction,
				// then there is a negative j-neighbor box (NULL otherwise).
				if (!new_child->code(j).is_neg_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = new_child->code(k).to_neg();
						else
							tcode[k] = new_child->code(k);
					target = find(tcode);
					new_child->set_neg_neighbor(target, j);
					if (new_child->Depth() == target->Depth())
						target->set_pos_neighbor(new_child, j);
				}
				// If this is not the positive boundary of j-th direction,
				// then there is a positive j-neighbor box (NULL otherwise).
				if (!new_child->code(j).is_pos_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = new_child->code(k).to_pos();
						else
							tcode[k] = new_child->code(k);
					target = find(tcode);
					new_child->set_pos_neighbor(target, j);
					if (new_child->Depth() == target->Depth())
						target->set_neg_neighbor(new_child, j);
				}
			}
		}
	}
	// Check if a box is a neighbor box (not necessarily principle).
	bool is_neighbor(R3Box* B)
	{
		int int_co_dim = 0;
		for (int i = 0; i < dim; ++i)
		{
			if (code(i).separate(B->code(i)))
				return false;
			if (code(i).neighbor_to(B->code(i)))
				int_co_dim += 1;
		}
		return int_co_dim == 1;
	}
	// Check if a box is containing a box.
	bool is_containing(R3Box* B)
	{
		for (int i = 0; i < dim; ++i)
			if (!code(i).contain(B->code(i)))
				return false;
		return true;
	}
	// Check if a box is intersecting to a box.
	bool is_intersect(R3Box* B)
	{
		return is_containing(B) || B->is_containing(this);
	}
	// If the box contains a configuration.
	bool contains(Vector3d v)
	{
		return range.contains(v);
	}
	// The center of the box.
	Vector3d center()
	{
		return (range.min() + range.max()) / 2;
	}
	// The range of the box.
	MatrixId Range()
	{
		return range;
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0, bool recur = true)
	{
		if (l > MAXSHOW)
			return;
		for (int i = 0; i < l; ++i)
			os << "      ";
		if (l == MAXSHOW)
		{
			os << "...";
			return;
		}
		Vector3d m = range.min();
		Vector3d M = range.max();
		os << range.transpose();
		if (leaf || !recur)
			return;
		for (int i = 0; i < subsize; ++i)
		{
			os << endl;
			children[i]->out(os, l + 1);
		}
	}
	// cout *neighbors.
	void show_neighbor(ostream& os = cout, int l = 0)
	{
		for (int i = 0; i < dim; ++i)
		{
			if (pos_neighbor(i) != NULL)
			{
				for (int j = 0; j < l; ++j)
					os << "      ";
				os << "neighbor pos " << i << ":" << endl;
				pos_neighbor(i)->show_code(os, l);
			}
			if (neg_neighbor(i) != NULL)
			{
				for (int j = 0; j < l; ++j)
					os << "      ";
				os << "neighbor neg " << i << ":" << endl;
				neg_neighbor(i)->show_code(os, l);
			}
		}
	}
	// cout *codes.
	void show_code(ostream& os = cout, int l = 0)
	{
		if (l > MAXSHOW)
			return;
		for (int j = 0; j < l; ++j)
			os << "      ";
		if (l == MAXSHOW)
		{
			os << "...";
			return;
		}
		os << code(0) << " " << code(1) << " " << code(2);
	}

#undef subsize
#undef dim
};

class SO3Box
{
#define subsize 8
#define dim 4
#define boxsize 4
#define EMP "N"
	// Id number.
	int id;
	// Root node.
	SO3Box* roots[boxsize];
	// Parent node.
	SO3Box* parent;
	// Which root it is in?
	int wxyz;
	// Depth of the node.
	int depth;
	// Is this a leaf?
	bool leaf;
	// Children node.
	SO3Box* children[subsize];
	// Neighbor node, 0 for neg, 1 for pos.
	SO3Box* neighbors[dim][2];
	// Bits of the box.
	bit codes[dim];
	// Skip the wxyz.
	int skip(int i)
	{
		return (i > wxyz) ? (i - 1) : i;
	}
public:
	// Numerical range of the box.
	MatrixId range;
	// Id number.
	int ID()
	{
		return id;
	}
	// Width of the box.
	double width()
	{
		if (wxyz < 0)
			return 2.0;
		return (range.max() - range.min()).minCoeff();
	}
	// Construct a root box.
	SO3Box()
	{
		range = MatrixId(Vector3d(-1, -1, -1), Vector3d(1, 1, 1));
		id = SO3list.size();
		SO3list.push_back(this);
		for (int i = 0; i < boxsize; ++i)
			roots[i] = new SO3Box(MatrixId(Vector3d(-1, -1, -1), Vector3d(1, 1, 1)), i);
		for (int i = 0; i < boxsize; ++i)
		{
			roots[i]->set_root(roots);
			for (int j = 0; j < dim; ++j)
				if (j == i)
					continue;
				else
				{
					roots[i]->set_pos_neighbor(roots[j], j);
					roots[i]->set_neg_neighbor(roots[j], j);
				}
		}
		parent = NULL;
		wxyz = -1;
		depth = -1;
		leaf = true;
		for (int i = 0; i < subsize; ++i)
			children[i] = NULL;
		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < 2; ++j)
				neighbors[i][j] = NULL;
		for (int i = 0; i < dim; ++i)
			codes[i] = bit();
	}
	// Construct a new box.
	SO3Box(MatrixId r, int _wxyz)
	{
		range = r;
		id = SO3list.size();
		SO3list.push_back(this);
		for (int i = 0; i < boxsize; ++i)
			roots[i] = NULL;
		parent = NULL;
		wxyz = _wxyz % boxsize;
		depth = 0;
		leaf = true;
		for (int i = 0; i < subsize; ++i)
			children[i] = NULL;
		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < 2; ++j)
				neighbors[i][j] = NULL;
		for (int i = 0; i < dim; ++i)
			codes[i] = bit();
	}
	// Root node.
	SO3Box* root(int i)
	{
		return roots[i % boxsize];
	}
	// Set the root node.
	void set_root(SO3Box* B[boxsize])
	{
		for (int i = 0; i < boxsize; ++i)
			roots[i] = B[i];
	}
	// Parent node.
	SO3Box* Parent()
	{
		return parent;
	}
	// Set the parent node.
	void set_parent(SO3Box* B)
	{
		parent = B;
	}
	// Which root it is in?
	int WXYZ()
	{
		return wxyz;
	}
	// Is this an SO3 root?
	bool is_root()
	{
		return wxyz < 0;
	}
	// Depth of the node.
	int Depth()
	{
		return depth;
	}
	// Set the depth of the node.
	void set_depth(int d)
	{
		depth = d;
	}
	// Is this a leaf?
	bool Leaf()
	{
		return leaf;
	}
	// Is this a leaf?
	bool is_leaf()
	{
		return leaf;
	}
	// Child node.
	SO3Box* child(int i)
	{
		if (wxyz < 0)
			return children[i % boxsize];
		else
			return children[i % subsize];
	}
	// Child node by code.
	SO3Box* child(int t[dim])
	{
		if (wxyz < 0)
			return NULL;
		int target = 0;
		for (int i = 0; i < dim; ++i)
			if (i < wxyz)
				target += (t[i] % 2) << i;
			else if (i > wxyz)
				target += (t[i] % 2) << (i - 1);
			else
				target += 0;
		return child(target);
	}
	// Neighbor node.
	SO3Box* neighbor(int i, bool pos)
	{
		if (wxyz < 0)
			return NULL;
		if (pos)
			return neighbors[i % dim][1];
		else
			return neighbors[i % dim][0];
	}
	// Positive neighbor node.
	SO3Box* pos_neighbor(int i)
	{
		if (wxyz < 0)
			return NULL;
		return neighbors[i % dim][1];
	}
	// Negative neighbor node.
	SO3Box* neg_neighbor(int i)
	{
		if (wxyz < 0)
			return NULL;
		return neighbors[i % dim][0];
	}
	// If the box is the negative boundary of the root to the direction i.
	bool is_neg_boundary(int i)
	{
		if (wxyz < 0)
			return true;
		return codes[i].is_neg_boundary();
	}
	// If the box is the positive boundary of the root to the direction i.
	bool is_pos_boundary(int i)
	{
		if (wxyz < 0)
			return true;
		return codes[i].is_pos_boundary();
	}
	// Set neighbor node.
	void set_neighbor(SO3Box* B, int i, bool pos)
	{
		if (pos)
			neighbors[i % dim][1] = B;
		else
			neighbors[i % dim][0] = B;
	}
	// Set positive neighbor node.
	void set_pos_neighbor(SO3Box* B, int i)
	{
		neighbors[i % dim][1] = B;
	}
	// Set negative neighbor node.
	void set_neg_neighbor(SO3Box* B, int i)
	{
		neighbors[i % dim][0] = B;
	}
	// The i-th bit of the node.
	bit code(int i)
	{
		return codes[i % dim];
	}
	// Set the i-th bit of the node.
	void set_code(bit t, int i)
	{
		codes[i % dim] = t;
	}
	// Find a target box by path indicator.
	SO3Box* find(bit t[dim], int _wxyz, bool show = false)
	{
		if (Noisity >= 12)
		{
			cout << "Finding path indicator: " << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << endl;
			cout << "In box " << _wxyz << endl;
		}
		SO3Box* target = roots[_wxyz];
		while (!target->is_leaf() && target->Depth() < t[(_wxyz + 1) % boxsize].length())
		{
			int tcode[dim];
			for (int i = 0; i < dim; ++i)
				tcode[i] = t[i][target->Depth() + 1];
			target = target->child(tcode);
		}
		if (Noisity >= 12)
		{
			cout << "Found ";
			target->show_code();
			cout << endl;
		}
		return target;
	}
	// Find a target box by child indicator.
	SO3Box* find(bit t[dim], SO3Box* b, bool show = false)
	{
		if (Noisity >= 12)
		{
			cout << "Finding path indicator: " << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << endl;
			cout << "In box " << b->WXYZ() << endl;
		}
		SO3Box* target = b;
		while (!target->is_leaf() && target->Depth() < t[(target->WXYZ() + 1) % boxsize].length())
		{
			int tcode[dim];
			for (int i = 0; i < dim; ++i)
				tcode[i] = t[i][target->Depth() + 1];
			target = target->child(tcode);
		}
		if (Noisity >= 12)
		{
			cout << "Found ";
			target->show_code();
			cout << endl;
		}
		return target;
	}
	// Subdivision.
	void subdivide(bool show = false)
	{
		// Step 0: Only leaves can be subdivided.
		if (!leaf)
			return;

		// Step 1: Create subdivisions.
		// All boxes are sequenced by the lex-inequality of the range.min() coordinates.

		if (wxyz < 0)
			for (int i = 0; i < boxsize; ++i)
			{
				children[i] = root(i);
				children[i]->set_parent(this);
				children[i]->set_depth(depth + 1);
			}
		else
			for (int i = 0; i < subsize; ++i)
			{
				Vector3d new_inf, new_sup;
				for (int j = 0; j < dim; ++j)
					if (j == wxyz)
						continue;
					else if (getBin(i, skip(j)))
					{
						new_inf(skip(j)) = (range(skip(j)).min() + range(skip(j)).max()) / 2;
						new_sup(skip(j)) = range(skip(j)).max();
					}
					else
					{
						new_inf(skip(j)) = range(skip(j)).min();
						new_sup(skip(j)) = (range(skip(j)).min() + range(skip(j)).max()) / 2;
					}
				MatrixId new_range(new_inf, new_sup, false);
				children[i] = new SO3Box(new_range, wxyz);
				children[i]->set_root(roots);
				children[i]->set_parent(this);
				children[i]->set_depth(depth + 1);
				for (int j = 0; j < dim; ++j)
					if (j == wxyz)
						children[i]->set_code(bit(), j);
					else if (getBin(i, skip(j)))
						children[i]->set_code(codes[j].pos(), j);
					else
						children[i]->set_code(codes[j].neg(), j);
			}

		// Step 2: This is not a leaf anymore.
		leaf = false;

		// Step 3: Assign neighbors for children and re-new_child neighbors for neighbors and their children.

		if (wxyz < 0)
			return;
		for (int i = 0; i < subsize; ++i)
		{
			SO3Box* new_child = children[i];
			SO3Box* target = NULL;
			bit tcode[dim];
			if (Noisity >= 14)
			{
				cout << "Assigning ";
				new_child->show_code();
				cout << endl;
			}
			for (int j = 0; j < dim; ++j)
			{
				if (Noisity >= 14)
					cout << "direction " << j << ":" << endl;
				if (j == wxyz)
					continue;
				// If this is not the negative boundary of j-th direction,
				// then there is a negative j-neighbor box (NULL otherwise).
				if (!new_child->code(j).is_neg_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = new_child->code(k).to_neg();
						else
							tcode[k] = new_child->code(k);
					target = find(tcode, wxyz, show);
					new_child->set_neg_neighbor(target, j);
					if (new_child->Depth() == target->Depth())
						target->set_pos_neighbor(new_child, j);
				}
				else
				{
					for (int k = 0; k < dim; ++k)
						if (k == wxyz)
							tcode[k] = new_child->code((k + 1) % dim).neg_most();
						else if (k == j)
							tcode[k] = bit();
						else
							tcode[k] = new_child->code(k).bar();
					target = find(tcode, j, show);
					new_child->set_neg_neighbor(target, j);
					if (new_child->Depth() == target->Depth())
						target->set_neg_neighbor(new_child, wxyz);
				}
				// If this is not the positive boundary of j-th direction,
				// then there is a positive j-neighbor box (NULL otherwise).
				if (!new_child->code(j).is_pos_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = new_child->code(k).to_pos();
						else
							tcode[k] = new_child->code(k);
					target = find(tcode, wxyz, show);
					new_child->set_pos_neighbor(target, j);
					if (new_child->Depth() == target->Depth())
						target->set_neg_neighbor(new_child, j);
				}
				else
				{
					for (int k = 0; k < dim; ++k)
						if (k == wxyz)
							tcode[k] = new_child->code((k + 1) % dim).pos_most();
						else if (k == j)
							tcode[k] = bit();
						else
							tcode[k] = new_child->code(k);
					target = find(tcode, j, show);
					new_child->set_pos_neighbor(target, j);
					if (new_child->Depth() == target->Depth())
						target->set_pos_neighbor(new_child, wxyz);
				}
			}
		}
	}
	// Check if a box is a neighbor box (not necessarily principle).
	bool is_neighbor(SO3Box* B)
	{
		if (wxyz < 0)
			return false;
		if (WXYZ() == B->WXYZ())
		{
			int int_co_dim = 0;
			for (int i = 0; i < dim; ++i)
			{
				if (i == WXYZ())
					continue;
				if (code(i).separate(B->code(i)))
					return false;
				if (code(i).neighbor_to(B->code(i)))
					int_co_dim += 1;
			}
			return int_co_dim == 1;
		}
		else
		{
			int pos_bound = -1;
			for (int i = 0; i < dim; ++i)
			{
				if (i == WXYZ())
				{
					if (!B->code(i).is_boundary())
						return false;
					if (pos_bound < 0)
					{
						if (B->code(i).is_neg_boundary())
							pos_bound = 0;
						else
							pos_bound = 1;
					}
					else
					{
						if (B->code(i).is_neg_boundary())
						{
							if (pos_bound != 0)
								return false;
						}
						else
						{
							if (pos_bound != 1)
								return false;
						}
					}
				}
				else if (i == B->WXYZ())
				{
					if (!code(i).is_boundary())
						return false;
					if (pos_bound < 0)
					{
						if (code(i).is_neg_boundary())
							pos_bound = 0;
						else
							pos_bound = 1;
					}
					else
					{
						if (code(i).is_neg_boundary())
						{
							if (pos_bound != 0)
								return false;
						}
						else
						{
							if (pos_bound != 1)
								return false;
						}
					}
				}
			}
			for (int i = 0; i < dim; ++i)
			{
				if (i != WXYZ() && i != B->WXYZ())
				{
					if (pos_bound == 1)
					{
						if ((!code(i).contain(B->code(i))) && (!code(i).contained(B->code(i))))
							return false;
					}
					else
					{
						if ((!code(i).flip().contain(B->code(i))) && (!code(i).flip().contained(B->code(i))))
							return false;
					}
				}
			}
		}
		return true;
	}
	// Check if a box is containing a box.
	bool is_containing(SO3Box* B)
	{
		if (wxyz < 0)
			return true;
		if (WXYZ() != B->WXYZ())
			return false;
		for (int i = 0; i < dim; ++i)
			if ((i != WXYZ()) && (!code(i).contain(B->code(i))))
				return false;
		return true;
	}
	// Check if a box is intersecting to a box.
	bool is_intersect(SO3Box* B)
	{
		return is_containing(B) || B->is_containing(this);
	}
	// If the box contains a configuration.
	bool contains(Vector4d v)
	{
		if (wxyz < 0)
			return true;
		else
		{
			Vector3d w;
			for (int i = 0; i < 4; ++i)
				if (i != wxyz)
					w(skip(i)) = v(i);
			return range.contains(w);
		}
	}
	// The center of the box.
	Vector4d center()
	{
		if (wxyz < 0)
			return Vector4d(1, 0, 0, 0);
		Vector3d v = (range.min() + range.max()) / 2;
		Vector4d w;
		for (int i = 0; i < 4; ++i)
			if (i != wxyz)
				w(i) = v(skip(i));
			else
				w(i) = 1;
		return w;
	}
	// The range of the box.
	MatrixId Range()
	{
		Vector3d m = range.min();
		Vector3d M = range.max();
		switch (wxyz)
		{
		case 0:return MatrixId(Vector4d(1, m(0), m(1), m(2)), Vector4d(1, M(0), M(1), M(2)), false);
		case 1:return MatrixId(Vector4d(m(0), 1, m(1), m(2)), Vector4d(M(0), 1, M(1), M(2)), false);
		case 2:return MatrixId(Vector4d(m(0), m(1), 1, m(2)), Vector4d(M(0), M(1), 1, M(2)), false);
		case 3:return MatrixId(Vector4d(m(0), m(1), m(2), 1), Vector4d(M(0), M(1), M(2), 1), false);
		default:return MatrixId(Vector4d(-1, -1, -1, -1), Vector4d(1, 1, 1, 1));
		}
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0, bool recur = true)
	{
		if (l > MAXSHOW)
			return;
		for (int i = 0; i < l; ++i)
			os << "      ";
		if (l == MAXSHOW)
		{
			os << "...";
			return;
		}
		os << Range().transpose();
		if (leaf || !recur)
			return;
		for (int i = 0; i < subsize; ++i)
		{
			os << endl;
			children[i]->out(os, l + 1);
		}
	}
	// cout *neighbors.
	void show_neighbor(ostream& os = cout, int l = 0)
	{
		for (int i = 0; i < dim; ++i)
		{
			for (int j = 0; j < l; ++j)
				os << "      ";
			os << "neighbor neg " << i << ":" << endl;
			if (neg_neighbor(i) != NULL)
				neg_neighbor(i)->show_code(os, l);
			else
				os << "NULL";
			os << endl;
			for (int j = 0; j < l; ++j)
				os << "      ";
			os << "neighbor pos " << i << ":" << endl;
			if (pos_neighbor(i) != NULL)
				pos_neighbor(i)->show_code(os, l);
			else
				os << "NULL";
			os << endl;
		}
	}
	// cout *codes.
	void show_code(ostream& os = cout, int l = 0)
	{
		if (l > MAXSHOW)
			return;
		for (int j = 0; j < l; ++j)
			os << "      ";
		if (l == MAXSHOW)
		{
			os << "...";
			return;
		}
		for (int i = 0; i < boxsize; ++i)
		{
			if (wxyz < 0 || i == wxyz)
				os << " " << EMP;
			else
				os << " " << code(i);
		}
	}
#undef EMP
#undef subsize
#undef dim
#undef boxsize
};

class SE3Box
{
#define subsize 8
#define dim 7
#define boxsize 4
#define subdim 3
protected:
	int id;
	// Translation component.
	R3Box* Bt;
	// Rotation component.
	SO3Box* Br;
	// Root node.
	SE3Box* roots[boxsize];
	// Parent node.
	SE3Box* parent;
	// Which root it is in?
	int wxyz;
	// Is this a leaf?
	bool leaf;
	// Children node.
	SE3Box* children[subsize];
	// Principle neighbor node, 0 for neg, 1 for pos.
	SE3Box* neighbors[dim][2];
	// Is this subdividing by R3?
	bool BTsub;
	// Is this subdividing by SO3?
	bool BRsub;
	// Skip the wxyz.
	int skip(int i)
	{
		return (i > wxyz) ? (i - 1) : i;
	}
public:
	// Id number.
	int ID()
	{
		return id;
	}
	// Translation component.
	R3Box* BT()
	{
		return Bt;
	}
	// Rotation component.
	SO3Box* BR()
	{
		return Br;
	}
	// Width of the box.
	double width()
	{
		return min(Bt->width(), Br->width());
	}
	// Construct a new box.
	SE3Box(R3Box* t, SO3Box* r)
	{
		id = SE3list.size();
		SE3list.push_back(this);
		SE3table.insert(t->ID(), r->ID(), id);
		Bt = t;
		Br = r;
		for (int i = 0; i < boxsize; ++i)
			roots[i] = NULL;
		parent = NULL;
		wxyz = Br->WXYZ();
		leaf = true;
		BTsub = false;
		BRsub = false;
		for (int i = 0; i < subsize; ++i)
			children[i] = NULL;
		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < 2; ++j)
				neighbors[i][j] = NULL;
	}
	// Deconstruct the box.
	~SE3Box()
	{
		for (int i = 0; i < subsize; ++i)
			if (child(i) != NULL)
				delete child(i);
		if (Bt != NULL)
			delete Bt;
		if (Br != NULL)
			delete Br;
	}
	// Root node.
	SE3Box* root(int i)
	{
		return roots[i % boxsize];
	}
	// Set the root node.
	void set_root(SE3Box* B[boxsize])
	{
		for (int i = 0; i < boxsize; ++i)
			roots[i] = B[i];
	}
	// Parent node.
	SE3Box* Parent()
	{
		return parent;
	}
	// Set the parent node.
	void set_parent(SE3Box* B)
	{
		parent = B;
	}
	// Which root it is in?
	int WXYZ()
	{
		return wxyz;
	}
	// Is this a leaf?
	bool Leaf()
	{
		return leaf;
	}
	// Is this a leaf?
	bool is_leaf()
	{
		return leaf;
	}
	// Is this an SO3 root box?
	bool is_root()
	{
		return Br->is_root();
	}
	// Is this box subdivided by BT?
	bool is_BTsub()
	{
		return BTsub;
	}
	// Is this box subdivided by BR?
	bool is_BRsub()
	{
		return BRsub;
	}
	// If the box contains a configuration.
	bool contains(VectorXd v)
	{
		Vector3d u;
		Vector4d w;
		for (int i = 0; i < subdim; ++i)
			u(i) = v(i);
		for (int i = 0; i < (dim - subdim); ++i)
			w(i) = v(i + subdim);
		return Bt->contains(u) && Br->contains(w);
	}
	// The center of the box.
	VectorXd center()
	{
		Vector3d u = Bt->center();
		Vector4d w = Br->center();
		VectorXd v(dim);
		for (int i = 0; i < dim; ++i)
			if (i < subdim)
				v(i) = u(i);
			else
				v(i) = w(i - subdim);
		return v;
	}
	// The range of the box.
	MatrixId Range()
	{
		VectorXd m(dim), M(dim);
		MatrixId Trange = BT()->Range(), Rrange = BR()->Range();
		for (int i = 0; i < dim; ++i)
		{
			if (i < subdim)
			{
				m(i) = Trange(i).min();
				M(i) = Trange(i).max();
			}
			else
			{
				m(i) = Rrange(i - subdim).min();
				M(i) = Rrange(i - subdim).max();
			}
		}
		return MatrixId(m, M);
	}
	// Child node.
	SE3Box* child(int i)
	{
		if (is_root() && BRsub)
			return children[i % boxsize];
		return children[i % subsize];
	}
	// Principle neighbor node.
	SE3Box* neighbor(int i, bool pos)
	{
		if (pos)
			return neighbors[i % dim][1];
		else
			return neighbors[i % dim][0];
	}
	// Principle neighbor for translational directions.
	SE3Box* T_neighbor(int i, bool pos)
	{
		return neighbor(i % subdim, pos);
	}
	// Principle neighbor for rotational directions.
	SE3Box* R_neighbor(int i, bool pos)
	{
		return neighbor((i % (dim - subdim)) + subdim, pos);
	}
	// Principle positive neighbor node.
	SE3Box* pos_neighbor(int i)
	{
		return neighbors[i % dim][1];
	}
	// Principle negative neighbor node.
	SE3Box* neg_neighbor(int i)
	{
		return neighbors[i % dim][0];
	}
	// Set principle neighbor node.
	void set_neighbor(SE3Box* B, int i, bool pos)
	{
		if (pos)
			neighbors[i % dim][1] = B;
		else
			neighbors[i % dim][0] = B;
	}
	// Set principle positive neighbor node.
	void set_pos_neighbor(SE3Box* B, int i)
	{
		neighbors[i % dim][1] = B;
	}
	// Set principle negative neighbor node.
	void set_neg_neighbor(SE3Box* B, int i)
	{
		neighbors[i % dim][0] = B;
	}
	// Find a target box from product table.
	SE3Box* find(int Btid, int Brid)
	{
		if (!SE3table.find(Btid, Brid))
			return NULL;
		else
			return SE3list[SE3table.coeff(Btid, Brid)];
	}
	// Determine pinciple i-neighbor for translational directions.
	SE3Box* find_T_neighbor(int i, bool pos, bool show = false)
	{
		i = i % subdim;
		R3Box* Tneighbor = BT()->neighbor(i, pos);
		if (Tneighbor == NULL)
			return NULL;
		if (Tneighbor->Depth() != BT()->Depth())
			return parent->T_neighbor(i, pos);
		if (Noisity >= 12)
		{
			cout << endl << "Finding: ";
			Tneighbor->show_code();
			cout << " *";
			BR()->show_code();
		}
		SE3Box* target = find(Tneighbor->ID(), BR()->ID());
		if (Noisity >= 12)
		{
			cout << endl << "Found: ";
			if (target == NULL)
				cout << "NULL";
			else
				target->show_code();
		}
		if (target == NULL)
			return parent->T_neighbor(i, pos);
		else
			return target;
	}
	// Determine pinciple i-neighbor for rotational directions.
	SE3Box* find_R_neighbor(int i, bool pos, bool show = false)
	{
		if (is_root())
			return NULL;
		i = i % (dim - subdim);
		if (i == wxyz)
			return NULL;
		SO3Box* Rneighbor = BR()->neighbor(i, pos);
		if (Rneighbor == NULL)
			return NULL;
		if (Rneighbor->Depth() != BR()->Depth())
			return parent->R_neighbor(i, pos);
		if (Noisity >= 12)
		{
			cout << endl << "Finding: ";
			BT()->show_code();
			cout << " *";
			Rneighbor->show_code();
		}
		SE3Box* target = find(BT()->ID(), Rneighbor->ID());
		if (Noisity >= 12)
		{
			cout << endl << "Found: ";
			if (target == NULL)
				cout << "NULL";
			else
				target->show_code();
		}
		if (target == NULL)
			return parent->R_neighbor(i, pos);
		else
			return target;
	}
	// Determine pinciple i-neighbor.
	SE3Box* find_neighbor(int i, bool pos, bool show = false)
	{
		if (i < subdim)
			return find_T_neighbor(i, pos, show);
		else
			return find_R_neighbor(i - subdim, pos, show);
	}
	// Determine if a box is aligned to this.
	bool align_to(SE3Box* B)
	{
		if (B == NULL)
			return false;
		return BT()->Depth() == B->BT()->Depth() && BR()->Depth() == B->BR()->Depth();
	}
	// When the principle neighbor is aligned, determine the direction from the neighbor to this.
	pair<int, bool> neighbor_dir(int i, bool pos)
	{
		if (i < subdim)
			return make_pair(i, !pos);
		else if (wxyz < 0)
			return make_pair(-1, false);
		else if (i == wxyz + subdim)
			return make_pair(i, pos);
		else if (!BR()->code(i - subdim).is_boundary(pos))
			return make_pair(i, !pos);
		else
			return make_pair(wxyz + subdim, pos);
	}
	// Translational partial subdivision.
	void T_Split(bool show = false)
	{
		// Step 0: Only leaves can be subdivided.
		if (!leaf)
			return;

		// Step 1: If the component is not subdivided, subdivide it first.
		if (Bt->Leaf())
			Bt->subdivide();

		// Step 2: Create subdivisions.
		for (int i = 0; i < subsize; ++i)
		{
			children[i] = new SE3Box(Bt->child(i), Br);
			children[i]->set_root(roots);
			children[i]->set_parent(this);
		}

		// Step 3: This is not a leaf anymore, it is subdivided by BT.
		leaf = false;
		BTsub = true;

		// Step 4: Assign neighbors for children and re-assign neighbors for neighbors and their children.
		for (int i = 0; i < subsize; ++i)
		{
			// The box that is going to be assigned to have a neighbor.
			SE3Box* new_child = children[i];
			// The box that is going to be assigned as the neighbor.
			SE3Box* target = NULL;
			if (Noisity >= 14)
			{
				cout << endl << "Assigning ";
				new_child->show_code();
			}
			for (int j = 0; j < dim; ++j)
			{
				if (Noisity >= 14)
					cout << endl << "direction " << j << ":";
				// Principle negative j-neighbor.
				for (int k = 0; k < 2; ++k)
				{
					bool pos = bool(k % 2);
					target = new_child->find_neighbor(j, pos, show);
					new_child->set_neighbor(target, j, pos);
					if (new_child->align_to(target))
					{
						pair<int, bool> nd = new_child->neighbor_dir(j, pos);
						target->set_neighbor(new_child, nd.first, nd.second);
					}
				}
			}
		}
	}
	// Rotational partial subdivision.
	void R_Split(bool show = false)
	{
		// Step 0: Only leaves can be subdivided.
		if (!leaf)
			return;

		// Step 1: If the component is not subdivided, subdivide it first.
		if (Br->Leaf())
			Br->subdivide();

		// Step 2: Create subdivisions.
		for (int i = 0; i < (is_root() ? boxsize : subsize); ++i)
		{
			children[i] = new SE3Box(Bt, Br->child(i));
			children[i]->set_root(roots);
			children[i]->set_parent(this);
		}

		// Step 3: This is not a leaf anymore, it is subdivided by BR.
		leaf = false;
		BRsub = true;

		// Step 4: Assign neighbors for children and re-assign neighbors for neighbors and their children.

		for (int i = 0; i < (is_root() ? boxsize : subsize); ++i)
		{
			SE3Box* new_child = children[i];
			SE3Box* target = NULL;
			if (Noisity >= 14)
			{
				cout << endl << "Assigning ";
				new_child->show_code();
			}
			for (int j = 0; j < dim; ++j)
			{
				if (Noisity >= 14)
					cout << endl << "direction " << j << ":";
				if (Noisity >= 14)
					cout << endl << "direction " << j << ":";
				// Principle negative j-neighbor.
				for (int k = 0; k < 2; ++k)
				{
					bool pos = bool(k % 2);
					target = new_child->find_neighbor(j, pos, show);
					new_child->set_neighbor(target, j, pos);
					if (new_child->align_to(target))
					{
						pair<int, bool> nd = new_child->neighbor_dir(j, pos);
						target->set_neighbor(new_child, nd.first, nd.second);
					}
				}
			}
		}
	}
	// Partial subdivision.
	void subdivide(int T = 0, bool show = false)
	{
		if (T == 1)
			T_Split(show);
		else if (T == -1)
			R_Split(show);
		else if (Bt->width() > 2 * varepsilon)
			T_Split(show);
		else
			R_Split(show);
		/*
		if (Bt->width() > 0.25)
			T_Split(show);
		else if (Br->width() * r0 > 0.25)
			R_Split(show);
		else if (Bt->width() >= Br->width() * r0)
			T_Split(show);
		else
			R_Split(show);*/
	}
	// Check if a box is a neighbor box (not necessarily principle).
	bool is_neighbor(SE3Box* B)
	{
		return(BT()->is_neighbor(B->BT()) && BR()->is_intersect(B->BR())) || (BT()->is_intersect(B->BT()) && BR()->is_neighbor(B->BR()));
	}
	// Collection of all neighbor boxes (not necessarily principle) in B.
	vector<SE3Box*> adj_neighbors(SE3Box* B)
	{
		vector<SE3Box*> neighbors;
		if (B == NULL)
			return neighbors;
		if (!is_neighbor(B))
			return neighbors;
		else if (B->Leaf())
			neighbors.push_back(B);
		else
			for (int i = 0; i < subsize; ++i)
			{
				vector<SE3Box*> subneighbors = adj_neighbors(B->child(i));
				for (int j = 0; j < subneighbors.size(); ++j)
					neighbors.push_back(subneighbors[j]);
			}
		return neighbors;
	}
	// Collection of all neighbor boxes (not necessarily principle).
	vector<SE3Box*> all_adj_neighbors(bool show = false)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		vector<SE3Box*> neighbors;
		vector<SE3Box*> subneighbors;
		for (int i = 0; i < dim; ++i)
		{
			if (Noisity >= 14)
				cout << endl << i;
			subneighbors = adj_neighbors(neg_neighbor(i));
			if (Noisity >= 14)
			{
				cout << endl;
				if (neg_neighbor(i) != NULL)
					neg_neighbor(i)->show_code();
				else
					cout << "NULL";
			}
			for (int j = 0; j < subneighbors.size(); ++j)
			{
				neighbors.push_back(subneighbors[j]);
				if (Noisity >= 14)
				{
					cout << endl;
					subneighbors[j]->show_code();
				}
			}
			subneighbors = adj_neighbors(pos_neighbor(i));
			if (Noisity >= 14)
			{
				cout << endl;
				if (pos_neighbor(i) != NULL)
					pos_neighbor(i)->show_code();
				else
					cout << "NULL";
			}
			for (int j = 0; j < subneighbors.size(); ++j)
			{
				neighbors.push_back(subneighbors[j]);
				if (Noisity >= 14)
				{
					cout << endl;
					subneighbors[j]->show_code();
				}
			}
		}
		auto end_time = std::chrono::high_resolution_clock::now();
		find_neighbor_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
		return neighbors;
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0)
	{
		if (l > MAXSHOW)
			return;
		for (int i = 0; i < l; ++i)
			os << "      ";
		if (l == MAXSHOW)
		{
			os << "...";
			return;
		}
		Bt->out(os, 0, false);
		os << " * ";
		Br->out(os, 0, false);
		if (leaf)
			return;
		for (int i = 0; i < subsize; ++i)
		{
			os << endl;
			children[i]->out(os, l + 1);
		}
	}
	// cout *neighbors.
	void show_neighbor(ostream& os = cout, int l = 0)
	{
		for (int i = 0; i < dim; ++i)
		{
			for (int j = 0; j < l; ++j)
				os << "      ";
			os << "neighbor pos " << i << ": ";
			if (pos_neighbor(i) != NULL)
			{
				pos_neighbor(i)->show_code(os, l);
				os << endl;
			}
			else
				os << "NULL" << endl;
			for (int j = 0; j < l; ++j)
				os << "      ";
			os << "neighbor neg " << i << ": ";
			if (neg_neighbor(i) != NULL)
			{
				neg_neighbor(i)->show_code(os, l);
				os << endl;
			}
			else
				os << "NULL" << endl;
		}
	}
	// cout *codes.
	void show_code(ostream& os = cout, int l = 0)
	{
		if (l > MAXSHOW)
			return;
		for (int j = 0; j < l; ++j)
			os << "      ";
		if (l == MAXSHOW)
		{
			os << "...";
			return;
		}
		Bt->show_code(os, 0);
		os << " *";
		Br->show_code(os, 0);
		os << " (Depth: " << Bt->Depth() << "*" << Br->Depth() << ")";
	}
#undef subsize
#undef dim
#undef boxsize
#undef subdim
};

// Distance between two boxes.
double Sep(SE3Box* b, SE3Box* p)
{
	return Sep(b->Range(), p->Range()).norm();
	/*
	double BTSep = Sep(b->BT()->Range(), p->BT()->Range()).norm();
	if (BTSep > 0)
		return BTSep;
	else
		return Sep(b->BR()->Range(), p->BR()->Range()).norm();*/
	//return 10 * Sep(b->BT()->Range(), p->BT()->Range()).norm() + Sep(b->BR()->Range(), p->BR()->Range()).norm();
}

SE3Box* to_box(int ID)
{
	return SE3list[ID];
}
#endif