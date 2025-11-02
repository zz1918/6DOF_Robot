// Indicator.cpp: This files constructs all indicators.

#ifndef INDICATOR_H
#define INDICATOR_H

#include <iostream>
#include <BinStr.h>
using namespace std;

// Direction indicators.
template <int dim>
class DIR
{
public:
	int e[dim];
	DIR()
	{
		for (int i = 0; i < dim; ++i)
			e[i] = 0;
	}
	DIR(int n, bool pos)
	{
		for (int i = 0; i < dim; ++i)
			e[i] = 0;
		if (pos)
			e[n % dim] = 1;
		else
			e[n % dim] = -1;
	}
	static DIR ele(int n, bool pos)
	{
		return DIR(n, pos);
	}
	// Flip of the direction.
	DIR flip()
	{
		DIR d;
		for (int i = 0; i < dim; ++i)
			d.e[i] = -e[i];
		return d;
	}
	// If this is an elementary direction.
	bool is_ele()
	{
		int rank = 0;
		for (int i = 0; i < dim; ++i)
			rank += abs(e[dim]);
		return rank == 1;
	}
	// The corresponding position of an elementary direction.
	int elenum()
	{
		for (int i = 0; i < dim; ++i)
			if (e[i] == 1)
				return i;
			else if (e[i] == -1)
				return i + dim;
	}
};

// Elementary Directions.
template<int dim>
class ELE
{
public:
	int comp;
	bool pos;
	ELE()
	{
		comp = -1;
		pos = true;
	}
	ELE(int n, bool p)
	{
		comp = n;
		pos = p;
	}
	// Flip of the direction.
	ELE flip()
	{
		return ELE(comp, !pos);
	}
	// The corresponding position of an elementary direction.
	int elenum()
	{
		if (pos)
			return comp;
		else
			return comp + dim;
	}
};

template<int dim>
ostream& operator<<(ostream& os, DIR<dim> d)
{
	for (int i = 0; i < dim - 1; ++i)
		os << d.e[i] << " ";
	os << d.e[dim - 1];
	return os;
}

// Child indicators.
template<int dim>
class CIdct
{
public:
	int e[dim];
	CIdct()
	{
		for (int i = 0; i < dim; ++i)
			e[i] = 0;
	}
	CIdct(int n)
	{
		for (int i = 0; i < dim; ++i)
			e[i] = (((n >> i) & 1) << 1) - 1;
	}
	CIdct(int n, int skip)
	{
		for (int i = 0; i < dim; ++i)
			if (i < skip)
				e[i] = (((n >> i) & 1) << 1) - 1;
			else if (i > skip)
				e[i] = (((n >> (i - 1)) & 1) << 1) - 1;
			else
				e[i] = 2;
	}
	// Child id.
	int num()
	{
		int n = 0;
		for (int i = 0; i < dim; ++i)
			switch (e[dim - i - 1])
			{
			case -1:n = n << 1; break;
			case 1:n = (n + 1) << 1; break;
			default:;
			}
		return n >> 1;
	}
	// Flip all components.
	CIdct<dim> flip()
	{
		CIdct<dim> idt;
		for (int i = 0; i < dim; ++i)
			idt.e[i] = -e[i];
		return idt;
	}
	// Flip the j-th component.
	CIdct<dim> flip(int j)
	{
		CIdct<dim> idt;
		for (int i = 0; i < dim; ++i)
			if (i != j)
				idt.e[i] = e[i];
			else
				idt.e[i] = -e[i];
		return idt;
	}
	// Prefix of an indicator.
	template<int subdim>
	CIdct<subdim> prefix()
	{
		CIdct<subdim> idt;
		for (int i = 0; i < subdim; ++i)
			idt.e[i] = e[i];
		return idt;
	}
	// Suffix of an indicator.
	template<int subdim>
	CIdct<dim - subdim> suffix()
	{
		CIdct<dim - subdim> idt;
		for (int i = subdim; i < dim; ++i)
			idt.e[i - subdim] = e[i];
		return idt;
	}
	// The id of skip component.
	int skid()
	{
		for (int i = 0; i < dim; ++i)
			if (e[i] == 2)
				return i;
		return -1;
	}
	// Add prefix for an indicator.
	template<int supdim>
	CIdct<supdim> prefux(int skip = -1)
	{
		CIdct<supdim> idt;
		for (int i = 0; i < supdim - dim; ++i)
			if (i != skip)
				idt.e[i] = 0;
			else
				idt.e[i] = 2;
		for (int i = 0; i < dim; ++i)
			idt.e[i + supdim - dim] = e[i];
		return idt;
	}
	// Add suffix for an indicator.
	template<int supdim>
	CIdct<supdim> suffux(int skip = -1)
	{
		CIdct<supdim> idt;
		for (int i = 0; i < dim; ++i)
			idt.e[i] = e[i];
		for (int i = 0; i < supdim - dim; ++i)
			if (i != skip)
				idt.e[i + dim] = 0;
			else
				idt.e[i + dim] = 2;
		return idt;
	}
};

template<int dim>
ostream& operator<<(ostream& os, CIdct<dim> idt)
{
	for (int i = 0; i < dim - 1; ++i)
		if (idt.e[i] == 2)
			os << "*" << " ";
		else
			os << idt.e[i] << " ";
	if (idt.e[dim - 1] == 2)
		os << "*";
	else
		os << idt.e[dim - 1];
	return os;
}

// Path indicators.
template<int dim>
class PIdct
{
public:
	bit e[dim];
	PIdct(int skip = -1)
	{
		for (int i = 0; i < dim; ++i)
			e[i] = bit::emp();
		e[skip] = bit::ast();
	}
	// The path indicator of the child corresponding to a child indicator.
	PIdct<dim> child(CIdct<dim> idt)
	{
		PIdct<dim> nidt;
		for (int i = 0; i < dim; ++i)
			switch (idt.e[i])
			{
			case -1:nidt.e[i] = e[i].neg(); break;
			case 0:; nidt.e[i] = e[i]; break;
			case 1:; nidt.e[i] = e[i].pos(); break;
			default:nidt.e[i] = bit::ast();
			}
		return nidt;
	}
	// The child indicator of the path indicator at t-th level.
	CIdct<dim> cidct(int t = 0)
	{
		CIdct<dim> idt;
		for (int i = 0; i < dim; ++i)
			if (e[i].is_ast())
				idt.e[i] = 2;
			else if (e[i].is_emp())
				idt.e[i] = 0;
			else if ((e[i].num() >> t) & 1)
				idt.e[i] = 1;
			else
				idt.e[i] = -1;
		return idt;
	}
	// Prefix of an indicator.
	template<int subdim>
	PIdct<subdim> prefix()
	{
		PIdct<subdim> idt;
		for (int i = 0; i < subdim; ++i)
			idt.e[i] = e[i];
		return idt;
	}
	// Suffix of an indicator.
	template<int subdim>
	PIdct<dim - subdim> suffix()
	{
		PIdct<dim - subdim> idt;
		for (int i = subdim; i < dim; ++i)
			idt.e[i - subdim] = e[i];
		return idt;
	}
};

template<int dim>
ostream& operator<<(ostream& os, PIdct<dim> idt)
{
	for (int i = 0; i < dim - 1; ++i)
		os << idt.e[i] << " ";
	os << idt.e[dim - 1];
	return os;
}

#endif