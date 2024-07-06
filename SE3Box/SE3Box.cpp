// SE3Box.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <interval.h>
#include <bimap.h>
using namespace std;

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
	// Output the bit.
	void out(ostream &os)
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
	// Bits of the box.
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
	// Child node by code.
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
	// Find a target box by bits.
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

		// Step 3: Assign neighbors for children and re-assign neighbors for neighbors and their children.

		for (int i = 0; i < subsize; ++i)
		{
			R3Box* assign = children[i];
			R3Box* target = NULL;
			bit tcode[dim];
			for (int j = 0; j < dim; ++j)
			{
				// If this is not the negative boundary of j-th direction,
				// then there is a negative j-neighbor box (NULL otherwise).
				if (!assign->code(j).is_neg_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = assign->code(k).to_neg();
						else
							tcode[k] = assign->code(k);
					target = find(tcode);
					assign->set_neg_neighbor(target, j);
					if (assign->Depth() == target->Depth())
						target->set_pos_neighbor(assign, j);
				}
				// If this is not the positive boundary of j-th direction,
				// then there is a positive j-neighbor box (NULL otherwise).
				if (!assign->code(j).is_pos_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = assign->code(k).to_pos();
						else
							tcode[k] = assign->code(k);
					target = find(tcode);
					assign->set_pos_neighbor(target, j);
					if (assign->Depth() == target->Depth())
						target->set_neg_neighbor(assign, j);
				}
			}
		}
	}
    // cout *this.
	void out(ostream& os = cout, int l = 0, bool recur = true)
	{
		for (int i = 0; i < l; ++i)
			os << "      ";
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
	void show_code(ostream& os=cout, int l = 0)
	{
		for (int j = 0; j < l; ++j)
			os << "      ";
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
		return (range.max() - range.min()).minCoeff();
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
	SO3Box* child(int i)
	{
		return children[i % subsize];
	}
	// Child node by code.
	SO3Box* child(int t[dim])
	{
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
		if (pos)
			return neighbors[i % dim][1];
		else
			return neighbors[i % dim][0];
	}
	// Positive neighbor node.
	SO3Box* pos_neighbor(int i)
	{
		return neighbors[i % dim][1];
	}
	// Negative neighbor node.
	SO3Box* neg_neighbor(int i)
	{
		return neighbors[i % dim][0];
	}
	// If the box is the negative boundary of the root to the direction i.
	bool is_neg_boundary(int i)
	{
		return codes[i].is_neg_boundary();
	}
	// If the box is the positive boundary of the root to the direction i.
	bool is_pos_boundary(int i)
	{
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
	// Find a target box by bits.
	SO3Box* find(bit t[dim], int _wxyz, bool show = false)
	{
		if (show)
		{
			cout << "Finding code: " << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << endl;
			cout << "In box " << _wxyz << endl;
		}
		SO3Box* target = roots[_wxyz];
		while (!target->Leaf() && target->Depth() < t[(_wxyz + 1) % boxsize].length())
		{
			int tcode[dim];
			for (int i = 0; i < dim; ++i)
				tcode[i] = t[i][target->Depth() + 1];
			target = target->child(tcode);
		}
		if (show)
		{
			cout << "Found ";
			target->show_code();
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

		// Step 3: Assign neighbors for children and re-assign neighbors for neighbors and their children.

		for (int i = 0; i < subsize; ++i)
		{
			SO3Box* assign = children[i];
			SO3Box* target = NULL;
			bit tcode[dim];
			if (show)
			{
				cout << "Assigning ";
				assign->show_code();
			}
			for (int j = 0; j < dim; ++j)
			{
				if (show)
					cout << "direction " << j << ":" << endl;
				if (j == wxyz)
					continue;
				// If this is not the negative boundary of j-th direction,
				// then there is a negative j-neighbor box (NULL otherwise).
				if (!assign->code(j).is_neg_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = assign->code(k).to_neg();
						else
							tcode[k] = assign->code(k);
					target = find(tcode, wxyz, show);
					assign->set_neg_neighbor(target, j);
					if (assign->Depth() == target->Depth())
						target->set_pos_neighbor(assign, j);
				}
				else
				{
					for (int k = 0; k < dim; ++k)
						if (k == wxyz)
							tcode[k] = assign->code((k + 1) % dim).neg_most();
						else if (k == j)
							tcode[k] = bit();
						else
							tcode[k] = assign->code(k).bar();
					target = find(tcode, j, show);
					assign->set_neg_neighbor(target, j);
					if (assign->Depth() == target->Depth())
						target->set_neg_neighbor(assign, wxyz);
				}
				// If this is not the positive boundary of j-th direction,
				// then there is a positive j-neighbor box (NULL otherwise).
				if (!assign->code(j).is_pos_boundary())
				{
					for (int k = 0; k < dim; ++k)
						if (k == j)
							tcode[k] = assign->code(k).to_pos();
						else
							tcode[k] = assign->code(k);
					target = find(tcode, wxyz, show);
					assign->set_pos_neighbor(target, j);
					if (assign->Depth() == target->Depth())
						target->set_neg_neighbor(assign, j);
				}
				else
				{
					for (int k = 0; k < dim; ++k)
						if (k == wxyz)
							tcode[k] = assign->code((k + 1) % dim).pos_most();
						else if (k == j)
							tcode[k] = bit();
						else
							tcode[k] = assign->code(k);
					target = find(tcode, j, show);
					assign->set_pos_neighbor(target, j);
					if (assign->Depth() == target->Depth())
						target->set_pos_neighbor(assign, wxyz);
				}
			}
		}
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0, bool recur = true)
	{
		for (int i = 0; i < l; ++i)
			os << "      ";
		Vector3d m = range.min();
		Vector3d M = range.max();
		switch (wxyz)
		{
		case 0:os << MatrixId(Vector4d(1, m(0), m(1), m(2)), Vector4d(1, M(0), M(1), M(2)), false).transpose(); break;
		case 1:os << MatrixId(Vector4d(m(0), 1, m(1), m(2)), Vector4d(M(0), 1, M(1), M(2)), false).transpose(); break;
		case 2:os << MatrixId(Vector4d(m(0), m(1), 1, m(2)), Vector4d(M(0), M(1), 1, M(2)), false).transpose(); break;
		case 3:os << MatrixId(Vector4d(m(0), m(1), m(2), 1), Vector4d(M(0), M(1), M(2), 1), false).transpose(); break;
		default:os << range.transpose();
		}
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
			os << "neighbor pos " << i << ":" << endl;
			if (pos_neighbor(i) != NULL)
				pos_neighbor(i)->show_code(os, l);
			else
				os << "NULL" << endl;
			for (int j = 0; j < l; ++j)
				os << "      ";
			os << "neighbor neg " << i << ":" << endl;
			if (neg_neighbor(i) != NULL)
				neg_neighbor(i)->show_code(os, l);
			else
				os << "NULL" << endl;
		}
	}
	// cout *codes.
	void show_code(ostream& os = cout, int l = 0)
	{
		for (int j = 0; j < l; ++j)
			os << "      ";
		for (int i = 0; i < boxsize; ++i)
		{
			if (i == wxyz)
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

class SO3Tree
{
#define dim 4
#define boxsize 4
	// Root nodes.
	SO3Box* roots[boxsize];
public:
	// SO3 initialization.
	SO3Tree()
	{
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
	}
	// Get the root cell.
	SO3Box* root(int i)
	{
		return roots[i % boxsize];
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0)
	{
		for (int i = 0; i < boxsize; ++i)
			roots[i]->out(os, l);
	}
#undef dim
#undef boxsize
};

class SE3Box
{
#define subsize 8
#define dim 7
#define boxsize 4
#define subdim 3
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
	// Neighbor node, 0 for neg, 1 for pos.
	SE3Box* neighbors[dim][2];
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
		for (int i = 0; i < subsize; ++i)
			children[i] = NULL;
		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < 2; ++j)
				neighbors[i][j] = NULL;
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
	// Child node.
	SE3Box* child(int i)
	{
		return children[i % subsize];
	}
	// Neighbor node.
	SE3Box* neighbor(int i, bool pos)
	{
		if (pos)
			return neighbors[i % dim][1];
		else
			return neighbors[i % dim][0];
	}
	// Positive neighbor node.
	SE3Box* pos_neighbor(int i)
	{
		return neighbors[i % dim][1];
	}
	// Negative neighbor node.
	SE3Box* neg_neighbor(int i)
	{
		return neighbors[i % dim][0];
	}
	// Set neighbor node.
	void set_neighbor(SE3Box* B, int i, bool pos)
	{
		if (pos)
			neighbors[i % dim][1] = B;
		else
			neighbors[i % dim][0] = B;
	}
	// Set positive neighbor node.
	void set_pos_neighbor(SE3Box* B, int i)
	{
		neighbors[i % dim][1] = B;
	}
	// Set negative neighbor node.
	void set_neg_neighbor(SE3Box* B, int i)
	{
		neighbors[i % dim][0] = B;
	}
	// Find a target box by product table.
	SE3Box* find(int Btid, int Brid)
	{
		if (!SE3table.find(Btid, Brid))
			return NULL;
		else
			return SE3list[SE3table.coeff(Btid, Brid)];
	}
	// Translational partial subdivision.
	void R3_subdivide(bool show = false)
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

		// Step 3: This is not a leaf anymore.
		leaf = false;

		// Step 4: Assign neighbors for children and re-assign neighbors for neighbors and their children.
		for (int i = 0; i < subsize; ++i)
		{
			SE3Box* assign = children[i];
			SE3Box* target = NULL;
			if (show)
			{
				cout << "Assigning ";
				assign->show_code();
			}
			for (int j = 0; j < dim; ++j)
			{
				if (show)
					cout << "direction " << j << ":" << endl;
				if (j == wxyz + subdim)
					continue;
				// Translational direction.
				if (j < subdim)
				{
					// Negative direction.
					if (Bt->neg_neighbor(j) == NULL && (!getBin(i, j)))
						assign->set_neg_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->neg_neighbor(j)->ID(), assign->BR()->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_neg_neighbor(neg_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_neg_neighbor(target, j);
							target->set_pos_neighbor(assign, j);
						}
					}
					// Positive direction.
					if (Bt->pos_neighbor(j) == NULL && getBin(i, j))
						assign->set_pos_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->pos_neighbor(j)->ID(), assign->BR()->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_pos_neighbor(pos_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_pos_neighbor(target, j);
							target->set_neg_neighbor(assign, j);
						}
					}
				}
				// Rotational direction.
				else
				{
					// Negative direction.
					if (Br->neg_neighbor(j - subdim) == NULL)
						assign->set_neg_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->ID(), assign->BR()->neg_neighbor(j - subdim)->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_neg_neighbor(neg_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_neg_neighbor(target, j);
							if (assign->BR()->is_neg_boundary(j - subdim))
								target->set_neg_neighbor(assign, wxyz + subdim);
							else
								target->set_pos_neighbor(assign, wxyz + subdim);
						}
					}
					// Positive direction.
					if (Br->pos_neighbor(j - subdim) == NULL)
						assign->set_pos_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->ID(), assign->BR()->pos_neighbor(j - subdim)->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_pos_neighbor(pos_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_pos_neighbor(target, j);
							if (assign->BR()->is_pos_boundary(j - subdim))
								target->set_pos_neighbor(assign, wxyz + subdim);
							else
								target->set_neg_neighbor(assign, wxyz + subdim);
						}
					}
				}
			}
		}
	}
	// Rotational partial subdivision.
	void SO3_subdivide(bool show = false)
	{
		// Step 0: Only leaves can be subdivided.
		if (!leaf)
			return;

		// Step 1: If the component is not subdivided, subdivide it first.
		if (Br->Leaf())
			Br->subdivide();

		// Step 2: Create subdivisions.
		for (int i = 0; i < subsize; ++i)
		{
			children[i] = new SE3Box(Bt, Br->child(i));
			children[i]->set_root(roots);
			children[i]->set_parent(this);
		}

		// Step 3: This is not a leaf anymore.
		leaf = false;

		// Step 4: Assign neighbors for children and re-assign neighbors for neighbors and their children.
		for (int i = 0; i < subsize; ++i)
		{
			SE3Box* assign = children[i];
			SE3Box* target = NULL;
			if (show)
			{
				cout << "Assigning ";
				assign->show_code();
			}
			for (int j = 0; j < dim; ++j)
			{
				if (show)
					cout << "direction " << j << ":" << endl;
				if (j == wxyz + subdim)
					continue;
				// Translational direction.
				if (j < subdim)
				{
					// Negative direction.
					if (Bt->neg_neighbor(j) == NULL)
						assign->set_neg_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->neg_neighbor(j)->ID(), assign->BR()->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_neg_neighbor(neg_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_neg_neighbor(target, j);
							target->set_pos_neighbor(assign, j);
						}
					}
					// Positive direction.
					if (Bt->pos_neighbor(j) == NULL)
						assign->set_pos_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->pos_neighbor(j)->ID(), assign->BR()->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_pos_neighbor(pos_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_pos_neighbor(target, j);
							target->set_neg_neighbor(assign, j);
						}
					}
				}
				// Rotational direction.
				else
				{
					// Negative direction.
					if (Br->neg_neighbor(j - subdim) == NULL)
						assign->set_neg_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->ID(), assign->BR()->neg_neighbor(j - subdim)->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_neg_neighbor(neg_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_neg_neighbor(target, j);
							if (assign->BR()->is_neg_boundary(j - subdim))
								target->set_neg_neighbor(assign, wxyz + subdim);
							else
								target->set_pos_neighbor(assign, wxyz + subdim);
						}
					}
					// Positive direction.
					if (Br->pos_neighbor(j - subdim) == NULL)
						assign->set_pos_neighbor(NULL, j);
					else
					{
						target = find(assign->BT()->ID(), assign->BR()->pos_neighbor(j - subdim)->ID());
						// Not found. Then the principle neighbor is the same with parent.
						if (target == NULL)
							assign->set_pos_neighbor(pos_neighbor(j), j);
						// Found. Then "assign" and "target" is an aligned pair.
						else
						{
							assign->set_pos_neighbor(target, j);
							if (assign->BR()->is_pos_boundary(j - subdim))
								target->set_pos_neighbor(assign, wxyz + subdim);
							else
								target->set_neg_neighbor(assign, wxyz + subdim);
						}
					}
				}
			}
		}
	}
	// Partial subdivision.
	void subdivide(bool show = false)
	{
		if (Bt->width() >= Br->width())
			R3_subdivide(show);
		else
			SO3_subdivide(show);
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0)
	{
		for (int i = 0; i < l; ++i)
			os << "      ";
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
		for (int j = 0; j < l; ++j)
			os << "      ";
		Bt->show_code(os, 0);
		os << " * ";
		Br->show_code(os, 0);
	}
#undef subsize
#undef dim
#undef boxsize
#undef subdim
};

class SE3Tree
{
#define subsize 8
#define dim 7
#define boxsize 4
#define subdim 3
	// Root nodes.
	SE3Box* roots[boxsize];
public:
	// SE3 initialization.
	SE3Tree(MatrixId range)
	{
		R3Box* BtRoot = new R3Box(range);
		SO3Tree* BrTree = new SO3Tree();
		for (int i = 0; i < boxsize; ++i)
			roots[i] = new SE3Box(BtRoot, BrTree->root(i));
		for (int i = 0; i < boxsize; ++i)
		{
			roots[i]->set_root(roots);
			for (int j = 0; j < dim; ++j)
				if (j < subdim || j == i + subdim)
					continue;
				else
				{
					roots[i]->set_pos_neighbor(roots[j - subdim], j);
					roots[i]->set_neg_neighbor(roots[j - subdim], j);
				}
		}
	}
	// Get the root cell.
	SE3Box* root(int i)
	{
		return roots[i % boxsize];
	}
	// cout *this.
	void out(ostream& os = cout, int l = 0)
	{
		for (int i = 0; i < boxsize; ++i)
		{
			roots[i]->out(os, l);
			os << endl;
		}
	}
#undef subsize
#undef dim
#undef boxsize
#undef subdim
};