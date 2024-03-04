// interval.h : 
//

#include "interval.cpp"

template<typename Scarlar>
class interval;
template<typename Scarlar>
interval<Scarlar> pow(interval<Scarlar> i, int n);
template<typename Scarlar>
interval<Scarlar> pow(interval<Scarlar> i, double n);
template<typename Scarlar>
istream& operator>>(istream& is, interval<Scarlar>& i);
template<typename Scarlar>
ostream& operator<<(ostream& os, interval<Scarlar> i);
template<typename Scarlar>
class MatrixInterval;
template<typename Scarlar>
ostream& operator<<(ostream& os, MatrixInterval<Scarlar> m);