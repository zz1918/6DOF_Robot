// WtFp.h: This file implements the special Sigma_2 decomposition 
//			of the approximate footprint of the Delta robot.

#include "WtFp.cpp"

enum predicate{ MIXED, FREE, STUCK, UNKNOWN };

// *********************Help functions********************* //

// Min of 3 numbers.
double minof(double a, double b, double c);
// Max of 3 numbers.
double maxof(double a, double b, double c);

// Psuedo feature Mesh, which can only be detected by WtFp but cannot generally do parametric queries.
class Mesh;

// Approximate footprint for Delta robot.
class WtFp;
