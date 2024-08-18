// Solid.h : This file defines solids and features for collision detections.
// Features are a type of objects that constructs with their id for recording.
// Solids can detect if Sep(Solids,Feature)>t? for any feature.
// Features includes Points, Edges and Triangles which are also solids themselves.
//

#include "Solid.cpp"

// Min of 3 numbers.
double minof(double a, double b, double c);
// Max of 3 numbers.
double maxof(double a, double b, double c);

class Solid;
class Feature;
class Point;
class Edge;
class Triangle;
class Trapezoid;
class Pyramid;
class TrafficCone;
class IceCream;
