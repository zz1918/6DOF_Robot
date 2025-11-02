// Indicator.h

#include "Indicator.cpp"

template <int dim>
class DIR;
template<int dim>
ostream& operator<<(ostream& os, DIR<dim> d);

template <int dim>
class CIdct;
template<int dim>
ostream& operator<<(ostream& os, CIdct<dim> idt);

template <int dim>
class PIdct;
template<int dim>
ostream& operator<<(ostream& os, PIdct<dim> idt);