// BoxTree.h : BoxTree maintains the Box structures in R^d, where d is dim.
//

#include <iostream>
#include"BoxTree.cpp"

// Class for dim-dimensional box and its subdivision.
template<int dim>
class Box;

template<int dim>
std::ostream& operator<<(std::ostream& os, Box<dim> B);

class R3Tree;

class SO3Tree;