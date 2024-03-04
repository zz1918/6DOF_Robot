// SymGroup.h : This file defines the data structure for symmetry group S_n
// Templated by n=size.

#include<iostream>
#include<iomanip>
#include"SymGroup.cpp"

int getBin(int n, int t);
int pow2(int t);
int pow2n(int n, int t);
int ne(int n);

template<int size>
class SymG;

template<int size>
std::ostream& operator<<(std::ostream& os, SymG<size> s);
