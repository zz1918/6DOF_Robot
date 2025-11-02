// BinStr.cpp : This file implements the binary string data structure 
// corresponding to the interval subdivisions.

#ifndef BINSTR_H
#define BINSTR_H

#include <iostream>
using namespace std;

// Binary string for the tree code, 0 for \bar{1}, 1 for 1. We add a "1" at the beginning to avoid starting with "0".
class bit
{
#define EMP "E"
#define AST "*"
	int n;								// bit number
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
		return n == 1;
	}
	// If this is an empty bit.
	bool is_emp()
	{
		return is_empty();
	}
	// If this is a special bit.
	bool is_ast()
	{
		return n == 0;
	}
	// Length of the string.
	int length()
	{
		if (n == 0)
			return -1;
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
	// If this is the boundary bit of pos ( = 0,1).
	bool is_boundary(bool pos)
	{
		if (pos)
			return is_pos_boundary();
		else
			return is_neg_boundary();
	}
	// Empty bit.
	static bit emp()
	{
		return bit();
	}
	// Special bit.
	static bit ast()
	{
		return bit(0);
	}
	// Positive child for this bit.
	bit pos()
	{
		if (is_ast())
			return bit::ast();
		return bit((n << 1) + 1);
	}
	// Negative child for this bit.
	bit neg()
	{
		if (is_ast())
			return bit::ast();
		return bit(n << 1);
	}
	// Flip of the bit.
	bit flip()
	{
		if (is_ast())
			return bit::ast();
		return bit(n ^ 1);
	}
	// Bar of the bit.
	bit bar()
	{
		if (is_ast())
			return bit::ast();
		return bit(n ^ ((1 << length()) - 1));
	}
	// Adj of the bit.
	bit adj()
	{
		if (is_ast())
			return bit::ast();
		if (n & 1)
			return bit(n + 1);
		else
			return bit(n - 1);
	}
	// Positive direction of the bit.
	bit to_pos()
	{
		if (is_ast())
			return bit::ast();
		return bit(n + 1);
	}
	// Negative direction of the bit.
	bit to_neg()
	{
		if (is_ast())
			return bit::ast();
		return bit(n - 1);
	}
	// The last bit in the positive direction.
	bit pos_most()
	{
		if (is_ast())
			return bit::ast();
		return bit((1 << (length() + 1)) - 1);
	}
	// The last bit in the negative direction.
	bit neg_most()
	{
		if (is_ast())
			return bit::ast();
		return bit(1 << length());
	}
	// If two bits are neighbor.
	bool neighbor_to(bit b)
	{
		if (is_ast())
			return false;
		if (length() > b.length())
			return b.neighbor_to(*this);
		int diff = b.length() - length();
		if (((n + 1) << diff) == b.n)
			return true;
		if ((n << diff) == (b.n + 1))
			return true;
		return false;
	}
	// If this bit contains another bit.
	bool contain(bit b)
	{
		if (is_ast())
			return false;
		if (length() > b.length())
			return false;
		int diff = b.length() - length();
		if (((n + 1) << diff) <= b.n)
			return false;
		if ((n << diff) >= (b.n + 1))
			return false;
		return true;
	}
	// If this bit is contained in another bit.
	bool contained(bit b)
	{
		if (is_ast())
			return false;
		return b.contain(*this);
	}
	// If this bit is not intersecting another bit.
	bool nintersect(bit b)
	{
		if (is_ast())
			return false;
		if (length() > b.length())
			return b.nintersect(*this);
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
		if (is_ast())
			os << AST;
		for (int i = length() - 1; i >= 0; --i)
			os << ((n >> i) & 1);
	}
#undef EMP
#undef AST
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

#endif