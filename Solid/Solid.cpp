// Solid.cpp : This file defines solids and features for collision detections.
// Features are a type of objects that constructs with their id for recording.
// Solids can detect if Sep(Solids,Feature)>t? for any feature.
// Features includes Points, Edges and Triangles which are also solids themselves.
//

#ifndef SOLID_H
#define SOLID_H

//#include <iostream>
#include <tuple>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;


////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Help Functions ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


// A help for quadratic equation, collecting roots of a * t^2 + 2 * b * t + c = 0 in [0,1].
vector<double> quad_root(double a, double b, double c)
{
    vector<double> t;

    // delta/4.
    double delta_4 = b * b - a * c;

    if (delta_4 < 0)
        return t;

    // Roots: t1 = c / d1 = d2 / a, t2 = c / d2 = d1 / a.
    double d1 = -b + sqrt(delta_4), d2 = -b - sqrt(delta_4);

    // Assert if 0 <= t1 <= 1, and 0 <= t2 <= 1.
    double a2 = a * a, c2 = c * c;

    // If |c|<|a|, then we use t1 = d2 / a and t2 = d1 / a.
    if (abs(c) < abs(a) || (abs(a) == abs(c) && c != 0))
    {
        // 0 <= d2 / a <= 1.
        if (d2 * a >= 0 && d2 * a <= a2)
            t.push_back(d2 / a);
        // 0 <= d1 / a <= 1.
        if (d1 * a >= 0 && d1 * a <= a2)
            t.push_back(d1 / a);
    }
    else if (abs(a) < abs(c))
    {
        // 0 <= c/d1 <= 1.
        if (d1 * c >= c2)
            t.push_back(c / d1);
        // 0 <= c/d2 <= 1.
        if (d2 * c >= c2)
            t.push_back(c / d2);
    }
    else
        t.push_back(0);
    return t;
}

// Plane (x-c).dot(n)=0 for normal vector n, line segment p+t*d for t in [0,1]. If the intersection between algebraic span is not in the segment.
bool possible_parallels(Vector3d c, Vector3d n, Vector3d p, Vector3d d)
{
    double d1 = -(p - c).dot(n), d2 = d.dot(n);
    return (d1 * d2 >= 0) && (d2 * d2 >= d1 * d2);
}

// Approximate intersection between (x-c).dot(n)=0 with p+t*d.
Vector3d PL_int(Vector3d c, Vector3d n, Vector3d p, Vector3d d)
{
    double d1 = -(p - c).dot(n), d2 = d.dot(n);
    return p + (d1 / d2) * d;
}

// Min of 3 numbers.
double minof(double a, double b, double c)
{
    return min(min(a, b), c);
}
// Max of 3 numbers.
double maxof(double a, double b, double c)
{
    return max(max(a, b), c);
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Declarations ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

class Solid;
class Feature;
class Point;
class Edge;
class Triangle;
class Trapezoid;
class Pyramid;
class TrafficCone;
class IceCream;

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Descriptions ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// A feature is a set that can be compared to.
class Feature {
public:
    int id;
    static int amount;
    Feature()
    {
        id = amount;
        amount++;
    }
};

int Feature::amount = 0;

// A solid is a set that can compare to features.
class Solid {
public:
    // f in this feature satiesfies Sep(this,f)>t.
    map<int, double> true_features;
    // f in this feature satiesfies Sep(this,f)<=t.
    map<int, double> false_features;
    // return true and record the result
    bool true_update(int F_id, double t)
    {
        map<int, double>::iterator found = true_features.find(F_id);
        if (found == true_features.end())
            true_features.insert(make_pair(F_id, t));
        else if (found->second > t)
            found->second = t;
        return true;
    }
    // return false and record the result
    bool false_update(int F_id, double t)
    {
        map<int, double>::iterator found = false_features.find(F_id);
        if (found == false_features.end())
            false_features.insert(make_pair(F_id, t));
        else if (found->second <= t)
            found->second = t;
        return false;
    }
    // return true and record the result
    bool true_update(Feature* f, double t)
    {
        return true_update(f->id, t);
    }
    // return false and record the result
    bool false_update(Feature* f, double t)
    {
        return false_update(f->id, t);
    }
    // If it is recorded false (return true if recorded).
    bool contains(int F_id, double t)
    {
        map<int, double>::iterator found = false_features.find(F_id);
        return ((found != false_features.end()) && found->second <= t);
    }
    // If it is recorded true (return true if recorded).
    bool ncontains(int F_id, double t)
    {
        map<int, double>::iterator found = true_features.find(F_id);
        return ((found != true_features.end()) && found->second >= t);
    }
    // If it is recorded false (return true if recorded).
    bool contains(Feature* f, double t)                         // If previously Sep(this,f)<=t.
    {
        return contains(f->id, t);
    }
    // If it is recorded true (return true if recorded).
    bool ncontains(Feature* f, double t)                        // If previously Sep(this,f)>t.
    {
        return ncontains(f->id, t);
    }

    // If two points are at different sides according to a vector.
    bool biside(Vector3d dir, Vector3d p, Vector3d q)
    {
        return dir.dot(p) * dir.dot(q) < 0;
    }

    // Sep(this,f)>t?
    virtual bool Sep(Point* f, double t) = 0;
    virtual bool Sep(Edge* f, double t) = 0;
    virtual bool Sep(Triangle* f, double t) = 0;
};

// Point is a feature, which is also a solid.
class Point :public Feature, public Solid {
public:
    Vector3d p;
    Point(Vector3d _p)
    {
        p = _p;
    }
    // Vector from this point to vector q.
    Vector3d direction(Vector3d q)
    {
        return q - p;
    }
    // Vector from this point to point feature f.
    Vector3d direction(Point* f)
    {
        return direction(f->p);
    }
    // Macro for direction().
    Vector3d D(Vector3d q)
    {
        return direction(q);
    }
    // Macro for direction().
    Vector3d D(Point* f)
    {
        return direction(f);
    }
    // Distance from this point to point q.
    double distance(Vector3d q)
    {
        return D(q).norm();
    }
    // Distance from this point to point feature f.
    double distance(Point* f)
    {
        return D(f).norm();
    }

    // Sep(this,f)>t?
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (distance(f) > t)
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Edge* f, double t = 0);
    // Sep(this,f)>t?
    bool Sep(Triangle* f, double t = 0);
};

// Edge is a feature, which is also a solid.
// This is actually half-edge.
class Edge :public Feature,public Solid {
public:
    Point* Boundary[2];
    Edge* inv;
    Vector3d the_direction;
    Edge(Vector3d p, Vector3d q)
    {
        Boundary[0] = new Point(p);
        Boundary[1] = new Point(q);
        the_direction = direction();
        inv = NULL;
    }
    Edge(Point* f, Point* g)
    {
        Boundary[0] = f;
        Boundary[1] = g;
        the_direction = direction();
        inv = NULL;
    }
    Edge(Point* f, Point* g,Vector3d dir)
    {
        Boundary[0] = f;
        Boundary[1] = g;
        the_direction = dir;
        inv = NULL;
    }

    // Set the inverse half-edge.
    void set_inv(Edge* f)
    {
        inv = f;
    }

    // return true and record the result
    bool true_update(int F_id, double t)
    {
        map<int, double>::iterator found = true_features.find(F_id);
        if (found == true_features.end())
            true_features.insert(make_pair(F_id, t));
        else if (found->second > t)
            found->second = t;
        if (inv == NULL)
            return true;
        // For the edge, we also need to update the inv edges.
        found = inv->true_features.find(F_id);
        if (found == inv->true_features.end())
            inv->true_features.insert(make_pair(F_id, t));
        else if (found->second > t)
            found->second = t;
        return true;
    }
    // return false and record the result
    bool false_update(int F_id, double t)
    {
        map<int, double>::iterator found = false_features.find(F_id);
        if (found == false_features.end())
            false_features.insert(make_pair(F_id, t));
        else if (found->second <= t)
            found->second = t;
        if (inv == NULL)
            return false;
        // For the edge, we also need to update the inv edges.
        found = inv->false_features.find(F_id);
        if (found == inv->false_features.end())
            inv->false_features.insert(make_pair(F_id, t));
        else if (found->second <= t)
            found->second = t;
        return false;
    }
    // return true and record the result
    bool true_update(Feature* f, double t)
    {
        return true_update(f->id, t);
    }
    // return false and record the result
    bool false_update(Feature* f, double t)
    {
        return false_update(f->id, t);
    }

    // The i-th point.
    Point* P(int i)
    {
        return Boundary[i % 2];
    }
    // Vector from P(0) to P(1).
    Vector3d direction()
    {
        return P(0)->direction(P(1));
    }
    // Macro for direction()
    Vector3d D()
    {
        return the_direction;
    }
    // Length of the edge.
    double length()
    {
        return D().norm();
    }

    //************Point Functions************//

    // Distance between point p and the line of the edge.
    double line_distance(Vector3d p)
    {
        return ((P(0)->D(p)).cross(D())).norm() / length();
    }
    // Distance between point feature f and the line of the edge.
    double line_distance(Point* f)
    {
        return line_distance(f->p);
    }

    // Signature of angle formed by p-P(i)-P(i+1) (positive if and only if it is acute angle)
    double sgn_angle(Vector3d p, int i)
    {
        return P(i)->D(p).dot(P(i)->D(P(i + 1)));
    }
    // If the projection of point p on the line of edge is in the interior of edge or not.
    bool projects_on(Vector3d p)
    {
        return (sgn_angle(p, 0) > 0) && (sgn_angle(p, 1) > 0);
    }
    // If the projection of point feature f on the line of edge is in the interior of edge or not.
    bool projects_on(Point* f)
    {
        return projects_on(f->p);
    }

    // Distance from this edge to point p.
    double distance(Vector3d p)
    {
        if (projects_on(p))
            return line_distance(p);
        else
            return min(P(0)->distance(p), P(1)->distance(p));
    }

    //************Edge Functions************//

    // Distance between line p-q and the line of the edge.
    double line_distance(Vector3d p, Vector3d q)
    {
        Vector3d Cross = D().cross(p - q);
        if (Cross.norm() == 0)
            return line_distance(p);
        else
            return abs(Cross.dot(P(0)->D(p))) / Cross.norm();
    }
    // Distance between line of point features f,g and the line of the edge.
    double line_distance(Point* f, Point* g)
    {
        return line_distance(f->p, g->p);
    }
    // Distance between line of edge feature f and the line of the edge.
    double line_distance(Edge* f)
    {
        return line_distance(f->P(0), f->P(1));
    }

    // If the projection the line of this edge on line p-q is in the edge p-q or not.
    bool projects_on(Vector3d p, Vector3d q)
    {
        Vector3d Cross = D().cross(p - q);

        // This is the vector perpendicular to the plane constructed by D() and Cross.
        Vector3d pq_dCross = D().cross(Cross);

        // p and q should be on the two sides of the plane constructed by D() and Cross.
        return biside(pq_dCross, P(0)->D(p), P(0)->D(q));
    }
    // If the projection the line of this edge on line of two point features f and g is in the edge of point features f-g or not.
    bool projects_on(Point* f, Point* g)
    {
        return projects_on(f->p, g->p);
    }
    // If the projection of the line of this edge on line of edge feature f is in the edge feature f or not.
    bool projects_on(Edge* f)
    {
        return projects_on(f->P(0), f->P(1));
    }

    //************Separation Functions************//

    // Interior criterions. Sep(this^\circ,f)>t?
    bool int_Sep(Point* f, double t)
    {
        if (line_distance(f) > t)
            return true;
        if (projects_on(f))
            return false;
        return true;
    }
    bool int_Sep(Edge* f, double t)
    {
        if (line_distance(f) > t)
            return true;
        if (projects_on(f) && f->projects_on(this))
            return false;
        return true;
    }

    // Sep(this,f)>t?
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if ((!P(0)->Sep(f, t)) || (!P(1)->Sep(f, t)))
            return false_update(f, t);
        if (int_Sep(f, t))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if ((!P(0)->Sep(f, t)) || (!P(1)->Sep(f, t)))
            return false_update(f, t);
        if ((!Sep(f->P(0), t)) || (!Sep(f->P(1), t)))
            return false_update(f, t);
        if (int_Sep(f, t))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Triangle* f, double t = 0);
};

bool Point::Sep(Edge* f, double t)
{
    return f->Sep(this, t);
}

// Triangle is a feature, which is also a solid.
class Triangle : public Feature,public Solid {
public:
    Point* Corner[3];
    Edge* Boundary[3];
    Vector3d the_normal;
    Triangle(Edge* f, Edge* g, Edge* h)
    {
        Boundary[0] = f;
        Boundary[1] = g;
        Boundary[2] = h;
        Corner[0] = f->P(0);
        Corner[1] = g->P(0);
        Corner[2] = h->P(0);
        the_normal = Normal();
    }

    // The i-th point.
    Point* P(int i)
    {
        return Corner[i % 3];
    }
    // The i-th edge.
    Edge* E(int i)
    {
        return Boundary[i % 3];
    }

    // Unit normal vector.
    Vector3d Normal()
    {
        return E(0)->D().cross(E(1)->D()).normalized();
    }
    // Macro for Normal().
    Vector3d N()
    {
        return the_normal;
    }
    // Area of the triangle.
    double area()
    {
        return E(0)->D().cross(E(1)->D()).norm() / 2;
    }

    //************Point Functions************//

    // Distance between point p and the plane of the triangle.
    double plane_distance(Vector3d p)
    {
        return abs(P(0)->D(p).dot(N()));
    }
    // Distance between point feature f and the plane of the triangle.
    double plane_distance(Point* f)
    {
        return plane_distance(f->p);
    }

    // Signature of dihedral angle formed by p-(P(i)P(i+1))-P(i+2) (positive if and only if it is acute angle)
    double sgn_angle(Vector3d p, int i)
    {
        return (P(i)->D(p).cross(E(i)->D())).dot(P(i)->D(P(i + 2)).cross(E(i)->D()));
    }
    // If the projection of point p on the plane of the triangle is in the interior of triangle or not.
    bool projects_on(Vector3d p)
    {
        return (sgn_angle(p, 0) > 0) && (sgn_angle(p, 1) > 0) && (sgn_angle(p, 2) > 0);
    }
    // If the projection of point feature f on the plane of the triangle is in the interior of triangle or not.
    bool projects_on(Point* f)
    {
        return projects_on(f->p);
    }


    // Distance from this triangle to point p.
    double distance(Vector3d p)
    {
        if (projects_on(p))
            return plane_distance(p);
        else
            return minof(E(0)->distance(p), E(1)->distance(p), E(2)->distance(p));
    }

    //************Edge Functions************//

    // The case edge f cannot reach the plane of this within the range of f.
    bool separate_to(Edge* f)
    {
        return possible_parallels(P(0)->p, N(), f->P(0)->p, f->D());
    }

    // Approximate intersection between line from p towards dir and the plane of triangle.
    Vector3d intersection(Vector3d p, Vector3d dir)
    {
        return PL_int(P(0)->p, N(), p, dir);
    }
    // Approximate intersection between line of f and the plane of triangle.
    Vector3d intersection(Edge* f)
    {
        return PL_int(P(0)->p, N(), f->P(0)->p, f->D());
    }
    // If the intersection between the plane of this triangle and the line of p-q is in segment p-q or not.
    bool projected_on(Vector3d p, Vector3d q)
    {
        return biside(N(), P(0)->D(p), P(0)->D(q));
    }
    // If the intersection between the plane of this triangle and the line of two point features is in edge formed by the two features or not.
    bool projected_on(Point* f, Point* g)
    {
        return projected_on(f->p, g->p);
    }
    // If the intersection between the plane of this triangle and the line of edge f is in f or not.
    bool projected_on(Edge* f)
    {
        return projected_on(f->P(0), f->P(1));
    }
    // If the intersection between the line of f and the plane of triangle is in triangle.
    bool intersects(Edge* f)
    {
        return projected_on(f) && projects_on(intersection(f));
    }

    //************Separation Functions************//

    // Interior criterions. Sep(this^\circ,f)>t?
    bool int_Sep(Point* f, double t)
    {
        if (plane_distance(f) > t)
            return true;
        if (projects_on(f))
            return false;
        return true;
    }
    bool int_Sep(Edge* f, double t = 0)
    {
        if (separate_to(f))
            return true;
        if (intersects(f))
            return false;
        return true;
    }

    // Sep(this,f)>t?
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        for (int i = 0; i < 3; ++i)
            if (!E(i)->Sep(f, t))
                return false_update(f, t);
        if (int_Sep(f, t))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        for (int i = 0; i < 3; ++i)
            if (!E(i)->Sep(f, t))
                return false_update(f, t);
        for (int i = 0; i < 2; ++i)
            if ((!Sep(f->P(i), t)))
                return false_update(f, t);
        if (int_Sep(f, t))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Triangle* f, double t = 0)
    {
        for (int i = 0; i < 3; ++i)
            if (!E(i)->Sep(f, t))
                return false_update(f, t);
        for (int i = 0; i < 3; ++i)
            if (!Sep(f->E(i), t))
                return false_update(f, t);
        return true_update(f, t);
    }
};

bool Point::Sep(Triangle* f, double t)
{
    return f->Sep(this, t);
}

bool Edge::Sep(Triangle* f, double t)
{
    return f->Sep(this, t);
}

// The segment solid is now merged into the edge (which is a solid as well as a feature).

// A trapezoid formed by the 4 corners where p0-p1 \\ p2-p3.
class Trapezoid :public Solid {
public:
    Point* Corner[4];
    Edge* Boundary[4];
    Vector3d the_normal;
    Trapezoid(Vector3d p, Vector3d q, Vector3d dir, double l1, double l2)
    {
        Point* P0 = new Point(p);
        Point* P1 = new Point(p + dir * l1);
        Point* Q0 = new Point(q);
        Point* Q1 = new Point(q + dir * l2);
        Corner[0] = P0;
        Corner[1] = P1;
        Corner[2] = Q1;
        Corner[3] = Q0;
        Boundary[0] = new Edge(P0, P1, dir * l1);
        Boundary[1] = new Edge(P1, Q1);
        Boundary[2] = new Edge(Q1, Q0, -dir * l2);
        Boundary[3] = new Edge(Q0, P0);
        the_normal = Normal();
    }
    Trapezoid(Point* p0, Point* p1, Point* p2, Point* p3)
    {
        Corner[0] = p0;
        Corner[1] = p1;
        Corner[2] = p2;
        Corner[3] = p3;
        Boundary[0] = new Edge(p0, p1);
        Boundary[1] = new Edge(p1, p2);
        Boundary[2] = new Edge(p2, p3);
        Boundary[3] = new Edge(p3, p0);
        the_normal = Normal();
    }

    // The i-th point.
    Point* P(int i)
    {
        return Corner[i % 4];
    }
    // The i-th edge.
    Edge* E(int i)
    {
        return Boundary[i % 4];
    }

    // Unit normal vector.
    Vector3d Normal()
    {
        return E(0)->D().cross(E(1)->D()).normalized();
    }
    // Macro for Normal().
    Vector3d N()
    {
        return the_normal;
    }
    // Area of the triangle.
    double area()
    {
        return (E(0)->D().cross(E(1)->D()).norm() + E(2)->D().cross(E(3)->D()).norm()) / 2;
    }


    //************Point Functions************//

    // Distance between point p and the plane of the triangle.
    double plane_distance(Vector3d p)
    {
        return abs(P(0)->D(p).dot(N()));
    }
    // Distance between point feature f and the plane of the triangle.
    double plane_distance(Point* f)
    {
        return plane_distance(f->p);
    }

    // Signature of dihedral angle formed by p-(P(i)P(i+1))-P(i+3) (positive if and only if it is acute angle)
    double sgn_angle(Vector3d p, int i)
    {
        return (P(i)->D(p).cross(E(i)->D())).dot(P(i)->D(P(i + 3)).cross(E(i)->D()));
    }
    // If the projection of point p on the plane of the triangle is in the interior of triangle or not.
    bool projects_on(Vector3d p)
    {
        return (sgn_angle(p, 0) > 0) && (sgn_angle(p, 1) > 0) && (sgn_angle(p, 2) > 0) && (sgn_angle(p, 3) > 0);
    }
    // If the projection of point feature f on the plane of the triangle is in the interior of triangle or not.
    bool projects_on(Point* f)
    {
        return projects_on(f->p);
    }


    //************Edge Functions************//

    // The case edge f cannot reach the plane of this within the range of f.
    bool separate_to(Edge* f)
    {
        return abs(N().dot(f->D())) < max(plane_distance(f->P(0)), plane_distance(f->P(1))) / f->length();
    }

    // Approximate intersection between line from p towards dir and the plane of triangle.
    Vector3d intersection(Vector3d p, Vector3d dir)
    {
        return p - (N().dot(P(0)->D(p)) / N().dot(dir)) * (dir);
    }
    // Approximate intersection between line of f and the plane of triangle.
    Vector3d intersection(Edge* f)
    {
        return intersection(f->P(0)->p, f->D());
    }
    // If the intersection between the plane of this triangle and the line of p-q is in segment p-q or not.
    bool projected_on(Vector3d p, Vector3d q)
    {
        return biside(N(), P(0)->D(p), P(0)->D(q));
    }
    // If the intersection between the plane of this triangle and the line of two point features is in edge formed by the two features or not.
    bool projected_on(Point* f, Point* g)
    {
        return projected_on(f->p, g->p);
    }
    // If the intersection between the plane of this triangle and the line of edge f is in f or not.
    bool projected_on(Edge* f)
    {
        return projected_on(f->P(0), f->P(1));
    }
    // If the intersection between the line of f and the plane of triangle is in triangle.
    bool intersects(Edge* f)
    {
        return projected_on(f) && projects_on(intersection(f));
    }


    //************Separation Functions************//

    // Interior criterions. Sep(this^\circ,f)>t?
    bool int_Sep(Point* f, double t)
    {
        if (plane_distance(f) > t)
            return true;
        if (projects_on(f))
            return false;
        return true;
    }
    bool int_Sep(Edge* f, double t)
    {
        if (separate_to(f))
            return true;
        if (intersects(f))
            return false;
        return true;
    }

    // Sep(this,f)>t?
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        for (int i = 0; i < 4; ++i)
            if (!E(i)->Sep(f, t))
                return false_update(f, t);
        if (int_Sep(f, t))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        for (int i = 0; i < 4; ++i)
            if (!E(i)->Sep(f, t))
                return false_update(f, t);
        for (int i = 0; i < 2; ++i)
            if ((!Sep(f->P(i), t)))
                return false_update(f, t);
        if (int_Sep(f, t))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Triangle* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        for (int i = 0; i < 4; ++i)
            if (!E(i)->Sep(f, t))
                return false_update(f, t);
        for (int i = 0; i < 3; ++i)
            if (!Sep(f->E(i), t))
                return false_update(f, t);
        return true_update(f, t);
    }
};

// A pyramid formed by two triangles p-q-r and p'-q'-r', where p' = p + dir * l1, q' = q + dir * l2, r' = r + dir * l3.
class Pyramid :public Solid {
public:
    Point* Corner[6];
    Edge* Side[9];
    Trapezoid* Vertical[3];
    Triangle* Horizontal[2];
    Pyramid(Vector3d p, Vector3d q, Vector3d r, Vector3d dir, double l1, double l2, double l3)
    {
        Point* P0 = new Point(p);
        Point* P1 = new Point(p + dir * l1);
        Point* Q0 = new Point(q);
        Point* Q1 = new Point(q + dir * l2);
        Point* R0 = new Point(r);
        Point* R1 = new Point(r + dir * l3);

        // Initialize corners
        Corner[0] = P0;
        Corner[1] = Q0;
        Corner[2] = R0;
        Corner[3] = P1;
        Corner[4] = Q1;
        Corner[5] = R1;

        // Initialize sides
        Side[0] = new Edge(P0, P1);
        Side[1] = new Edge(Q0, Q1);
        Side[2] = new Edge(R0, R1);
        Side[3] = new Edge(P0, Q0);
        Side[4] = new Edge(Q0, R0);
        Side[5] = new Edge(R0, P0);
        Side[6] = new Edge(P1, Q1);
        Side[7] = new Edge(Q1, R1);
        Side[8] = new Edge(R1, P1);

        // Initialize boundaries
        Vertical[0] = new Trapezoid(P0, P1, Q1, Q0);
        Vertical[1] = new Trapezoid(Q0, Q1, R1, R0);
        Vertical[2] = new Trapezoid(R0, R1, P1, P0);
        Horizontal[0] = new Triangle(Side[3], Side[4], Side[5]);
        Horizontal[1] = new Triangle(Side[6], Side[7], Side[8]);
    }

    Point* P(int i)
    {
        return Corner[i % 6];
    }
    Edge* E(int i)
    {
        return Side[i % 9];
    }
    Trapezoid* V(int i)
    {
        return Vertical[i % 3];
    }
    Triangle* H(int i)
    {
        return Horizontal[i % 2];
    }

    // If a point p is inside the convex hull formed by the corners or not.
    bool inside(Vector3d p)
    {
        if (H(0)->P(0)->D(p).dot(H(0)->N()) * H(1)->P(0)->D(p).dot(H(1)->N()) > 0)
            return false;
        if (V(0)->P(0)->D(p).dot(V(0)->N()) * V(1)->P(0)->D(p).dot(V(1)->N()) < 0)
            return false;
        if (V(0)->P(0)->D(p).dot(V(0)->N()) * V(2)->P(0)->D(p).dot(V(2)->N()) < 0)
            return false;
        return true;
    }
    // If a point feature may inside this pyramid.
    bool inside(Point* f)
    {
        return inside(f->p);
    }
    // If an edge feature may inside this pyramid.
    bool inside(Edge* f)
    {
        return inside(f->P(0));
    }
    // If a triangle feature may inside this pyramid.
    bool inside(Triangle* f)
    {
        return inside(f->P(0));
    }

    // Sep(this,f)>t?
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (inside(f))
            return false_update(f, t);
        for (int i = 0; i < 3; ++i)
            if (!V(i)->Sep(f, t))
                return false_update(f, t);
        for (int i = 0; i < 2; ++i)
            if (!H(i)->Sep(f, t))
                return false_update(f, t);
        return true_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (inside(f))
            return false_update(f, t);
        for (int i = 0; i < 3; ++i)
            if (!V(i)->Sep(f, t))
                return false_update(f, t);
        for (int i = 0; i < 2; ++i)
            if (!H(i)->Sep(f, t))
                return false_update(f, t);
        return true_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Triangle* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (inside(f))
            return false_update(f, t);
        for (int i = 0; i < 3; ++i)
            if (!V(i)->Sep(f, t))
                return false_update(f, t);
        for (int i = 0; i < 2; ++i)
            if (!H(i)->Sep(f, t))
                return false_update(f, t);
        return true_update(f, t);
    }
};

// A box formed by a corner and three directions to be edges.
class Cuboid: public Solid {
public:
    // Relative center.
    Vector3d rc;
    // Relative x,y,z axis.
    Vector3d rx, ry, rz;
    // 8 corners.
    Point* corner[8];
    // 6 faces.
    Trapezoid* facet[6];
    // The i-th corner.
    Point* V(int i)
    {
        return corner[i];
    }
    // The i-th face.
    Trapezoid* F(int i)
    {
        return facet[i];
    }
    // The corners of the cuboid is then rc, rc+rx, rc+ry, rc+rx+ry, rc+rz, rc+rx+rz, rc+ry+rz, rc+rx+ry+rz
    void set_VF()
    {
        corner[0] = new Point(rc);
        corner[1] = new Point(rc + rx);
        corner[2] = new Point(rc + ry);
        corner[3] = new Point(rc + rx + ry);
        corner[4] = new Point(rc + rz);
        corner[5] = new Point(rc + rx + rz);
        corner[6] = new Point(rc + ry + rz);
        corner[7] = new Point(rc + rx + ry + rz);
        facet[0] = new Trapezoid(V(0), V(1), V(3), V(2));
        facet[1] = new Trapezoid(V(1), V(0), V(4), V(5));
        facet[2] = new Trapezoid(V(4), V(0), V(2), V(6));
        facet[3] = new Trapezoid(V(5), V(4), V(6), V(7));
        facet[4] = new Trapezoid(V(3), V(1), V(5), V(7));
        facet[5] = new Trapezoid(V(2), V(3), V(7), V(6));
    }
    // Construct from relative center and relative bases.
    Cuboid(Vector3d _rc, Vector3d _rx, Vector3d _ry, Vector3d _rz)
    {
        rc = _rc;
        rx = _rx;
        ry = _ry;
        rz = _rz;
        set_VF();
    }
    // Construct from a matrix interval.
    Cuboid(Vector3d Imin, Vector3d Imax)
    {
        rc = Imin;
        rx = Vector3d(Imax(0) - Imin(0), 0, 0);
        ry = Vector3d(0, Imax(1) - Imin(1), 0);
        rz = Vector3d(0, 0, Imax(2) - Imin(2));
        set_VF();
    }
    
    // v relative to rc.
    Vector3d dir(Vector3d v)
    {
        return v - rc;
    }
    // If a point is inside the cuboid.
    bool inside(Vector3d v)
    {
        Vector3d dv = dir(v);
        if (dv.dot(rx) < 0)
            return false;
        if (dv.dot(rx) > rx.squaredNorm())
            return false;
        if (dv.dot(ry) < 0)
            return false;
        if (dv.dot(ry) > ry.squaredNorm())
            return false;
        if (dv.dot(rz) < 0)
            return false;
        if (dv.dot(rz) > rz.squaredNorm())
            return false;
        return true;
    }
    // If a point feature may inside this  cuboid.
    bool inside(Point* f)
    {
        return inside(f->p);
    }
    // If an edge feature may inside this  cuboid.
    bool inside(Edge* f)
    {
        return inside(f->P(0));
    }
    // If a triangle feature may inside this  cuboid.
    bool inside(Triangle* f)
    {
        return inside(f->P(0));
    }
    // Sep(this,f)>t?
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (inside(f))
            return false_update(f, t);
        for (int i = 0; i < 6; ++i)
            if (!F(i)->Sep(f, t))
                return false_update(f, t);
        return true_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (inside(f))
            return false_update(f, t);
        for (int i = 0; i < 6; ++i)
            if (!F(i)->Sep(f, t))
                return false_update(f, t);
        return true_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Triangle* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (inside(f))
            return false_update(f, t);
        for (int i = 0; i < 6; ++i)
            if (!F(i)->Sep(f, t))
                return false_update(f, t);
        return true_update(f, t);
    }
};

// A ball with center o and radius r.
class Ball :public Solid {
public:
    Point* origin;
    double radius;
    Ball(Vector3d o, double r)
    {
        origin = new Point(o);
        radius = r;
    }
    Ball(Point* O, double r)
    {
        origin = O;
        radius = r;
    }

    Point* O()
    {
        return origin;
    }
    double R()
    {
        return radius;
    }

    //************Separation Functions************//

    // The separation is given by the Minkowski sum lemma.
    
    // Sep(this,f)>t?
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (O()->Sep(f, t + R()))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (O()->Sep(f, t + R()))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
    // Sep(this,f)>t?
    bool Sep(Triangle* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (O()->Sep(f, t + R()))
            return true_update(f, t);
        else
            return false_update(f, t);
    }
};

// This is a 2d circular-surface of a right cone with apex v, and base disc(c,r).
class TrafficCone :public Solid {
public:
    Vector3d apex;
    Vector3d center;
    double radius;

    Edge* axis;
    Matrix3d Inv_Transition;
    double k, d;
    // k*k
    double k2;

    TrafficCone(Vector3d v, Vector3d c, double r)
    {
        apex = v;
        center = c;
        radius = r;
        axis = new Edge(v, c);
        Inv_Transition = Transition().inverse();
        d = axis->length();
        k = radius / d;
        k2 = k * k;
    }

    //************Linear Transformation Functions************//

    // SO(3) matrix that makes this into (x-x_0)^2+(y-y_0)^2=k^2(z-z_0)^2, 0<=z-z_0<=d.
    Matrix3d Transition()
    {
        Vector3d new_z = (center - apex).normalized();          // Coordinate of z under original coordinate.
        Vector3d new_x = Vector3d(new_z(2), 0, -new_z(0));
        if (new_x == Vector3d::Zero())
            new_x = Vector3d(1, 0, 0);
        else
            new_x = new_x.normalized();
        Vector3d new_y = new_z.cross(new_x);
        Matrix3d trans;
        trans.col(0) = new_x;
        trans.col(1) = new_y;
        trans.col(2) = new_z;
        return trans;
    }
    // Linear transformation that makes this into x^2+y^2=k^2z^2, 0<=z<=d.
    Vector3d Local(Vector3d p)
    {
        return Inv_Transition * (p - apex);
    }
    // Local coordinate of a point feature f.
    Vector3d Local(Point* f)
    {
        return Local(f->p);
    }
    // Local coordinates of an edge feature f.
    pair<Vector3d,Vector3d> Local(Edge* f)
    {
        return make_pair(Local(f->P(0)), Local(f->P(1)));
    }
    // Local coordinates of a triangle feature f.
    std::tuple<Vector3d, Vector3d, Vector3d> Local(Triangle* f)
    {
        return std::make_tuple(Local(f->P(0)), Local(f->P(1)), Local(f->P(2)));
    }

    // If a point p on the algebraic span of this traffic cone is in the interior of this traffic cone or not.
    bool projects_on(Vector3d Lp)
    {
        return Lp(2) <= d && Lp(2) >= 0;
    }

    //************Point Functions************//

    // If a point p is under the traffic cone.
    bool inside(Vector3d Lp)
    {
        double x = Lp(0), y = Lp(1), z = Lp(2);
        if (z < 0)
            return false;
        if (z > d)
            return false;
        if (x * x + y * y > k2 * z * z)
            return false;
        return true;
    }
    // If a point feature f is under the traffic cone.
    bool inside(Point* f)
    {
        return inside(Local(f));
    }
    // Local closest pairs between a local point Lp and the algebraic span of this traffic cone.
    vector<pair<Vector3d, Vector3d> > Lcp(Vector3d Lp)
    {
        double a = Lp(0), b = Lp(1), c = Lp(2);
        double r = sqrt(a * a + b * b);
        double t1 = (r - k * c) / (k2 * r + k * c), t2 = (r + k * c) / (k2 * r - k * c);
        Vector3d Lq1 = Vector3d(a / (1 + t1), b / (1 + t1), c / (1 - k2 * t1));
        Vector3d Lq2 = Vector3d(a / (1 + t2), b / (1 + t2), c / (1 - k2 * t2));
        vector<pair<Vector3d, Vector3d> > Lcps;
        Lcps.push_back(make_pair(Lq1, Lp));
        Lcps.push_back(make_pair(Lq2, Lp));
        return Lcps;
    }
    // Local closest pairs between a point feature and this traffic cone.
    vector<pair<Vector3d, Vector3d> > Lcp(Point* f)
    {
        return Lcp(Local(f));
    }

    //************Edge Functions************//

    // Collection of intersections between the local edge Le and this traffic cone (within the valid range 0<=t<=1).
    vector<Vector3d> Lint(pair<Vector3d, Vector3d> Le)
    {
        // l(t) := v+td. Le is 0<=t<=1.
        Vector3d v = Le.first, d = Le.second - Le.first;
        vector<Vector3d> Lints;

        // Root (t) of c_0 + 2c_1t + c_2t^2 = 0 are the intersections l(t).
        double c0 = v.dot(v) - (k2 + 1) * v(2) * v(2);
        double c1 = v.dot(d) - (k2 + 1) * v(2) * d(2);
        double c2 = d.dot(d) - (k2 + 1) * d(2) * d(2);

        // Collect possible intersections in Le.
        vector<double> t = quad_root(c2, c1, c0);
        for (int i = 0; i < t.size(); ++i)
            Lints.push_back(v + t[i] * d);
        return Lints;
    }
    // Local closest pairs between the line of Lp1-Lp2 and the algebraic span of this traffic cone.
    vector<pair<Vector3d, Vector3d> > Lcp(pair<Vector3d, Vector3d> Le)
    {
        vector<pair<Vector3d, Vector3d> > Lcps;
        Vector3d v = Le.first, d = Le.second - Le.first;
        double a = d(0), b = d(1), c = d(2), a2 = a * a, b2 = b * b, c2 = c * c;
        double kappa = k * c, kappa2 = kappa * kappa;

        // Solving quadratic equation for perpendicular conditions.
        double detA = kappa2 * (kappa2 - a2 - b2);

        if (detA > 0)
            return Lcps;

        // nx * x + nyi * y = 0 are the possible planes containing p (normal vector is (nx,nyi,0)).
        double nx = kappa2 - a2, ny1 = a * b + sqrt(-detA), ny2 = a * b - sqrt(-detA);
        Vector3d n1 = Vector3d(nx, ny1, 0).normalized(), n2 = Vector3d(nx, ny2, 0).normalized();

        // If the plane n1 detects a pair (irrelavant pairs are OK to be pushed).
        if (!possible_parallels(Vector3d::Zero(), n1, v, d))
        {
            Vector3d q = PL_int(Vector3d::Zero(), n1, v, d);
            vector<pair<Vector3d, Vector3d> > q_pair = Lcp(q);
            for (int i = 0; i < q_pair.size(); ++i)
                Lcps.push_back(q_pair[i]);
        }
        // If the plane n2 detects a pair (irrelavant pairs are OK to be pushed).
        if (!possible_parallels(Vector3d::Zero(), n2, v, d))
        {
            Vector3d q = PL_int(Vector3d::Zero(), n2, v, d);
            vector<pair<Vector3d, Vector3d> > q_pair = Lcp(q);
            for (int i = 0; i < q_pair.size(); ++i)
                Lcps.push_back(q_pair[i]);
        }
        return Lcps;
    }
    // Local closest pairs between an edge feature and this traffic cone.
    vector<pair<Vector3d, Vector3d> > Lcp(Edge* f)
    {
        pair<Vector3d, Vector3d> Le = Local(f);
        vector<pair<Vector3d, Vector3d> > CP = Lcp(Le);
        vector<Vector3d> LINT = Lint(Le);
        for (int i = 0; i < LINT.size(); ++i)
            CP.push_back(make_pair(LINT[i], LINT[i]));
        return CP;
    }


    //************Triangle Functions************//

    // There are no closest pairs between a triangle feature and this traffic cone.
    
    // If the axis intersects a triangle feature.
    bool axis_int(Triangle* f)
    {
        return !(f->int_Sep(axis, 0));
    }

    //******************Separation Functions*****************//

    // Sep(this,f)>t?
    // This separation is only given by the interior separation.
    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (inside(f))
            return false_update(f, t);
        vector<pair<Vector3d, Vector3d> > Lcps = Lcp(f);
        for (int i = 0; i < Lcps.size(); ++i)
            if (projects_on(Lcps[i].first) && (Lcps[i].first - Lcps[i].second).norm() <= t)
                return false_update(f, t);
        return true_update(f, t);
    }
    // Sep(this,f)>t?
    // This separation is only given by the interior separation.
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        for (int i = 0; i < 2; ++i)
            if (!Sep(f->P(i), t))
                return false_update(f, t);
        vector<pair<Vector3d, Vector3d> > Lcps = Lcp(f);
        for (int i = 0; i < Lcps.size(); ++i)
            if (projects_on(Lcps[i].first) && f->projects_on(Lcps[i].second) && (Lcps[i].first - Lcps[i].second).norm() <= t)
                return false_update(f, t);
        return true_update(f, t);
    }
    // Sep(this,f)>t?
    // This separation is only given by the interior separation.
    bool Sep(Triangle* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        for (int i = 0; i < 3; ++i)
            if (!Sep(f->E(i), t))
                return false_update(f, t);
        if (axis_int(f))
            return false_update(f, t);
        return true_update(f, t);
    }
};

// An ice-cream is formed by the union of a traffic cone with a ball.
// An ice-cream is defined by the apex v, center of ball o, and radius of ball r.
class IceCream :public Solid {
public:
    Point* apex;
    Ball* base;
    TrafficCone* conic;
    IceCream(Vector3d v,Vector3d o,double r)
    {
        apex = new Point(v);
        base = new Ball(o, r);
        
        // Compute the argument of the traffic cone.
        Vector3d ov = v - o;
        double h = ov.norm();
        double k = sqrt(h * h - r * r);
        Vector3d c = o + ((r * r) / (h * h)) * ov;
        conic = new TrafficCone(v, c, r * k / h);
    }
    IceCream(Point* V, Point* O, double r)
    {
        apex = V;
        base = new Ball(O, r);

        // Compute the argument of the traffic cone.
        Vector3d ov = O->D(V);
        double h = ov.norm();
        double k = sqrt(h * h - r * r);
        Vector3d c = O->p + ((r * r) / (h * h)) * ov;
        conic = new TrafficCone(V->p, c, r * k / h);
    }

    Point* V()
    {
        return apex;
    }
    Ball* B()
    {
        return base;
    }
    TrafficCone* C()
    {
        return conic;
    }


    //************Separation Functions************//

    bool Sep(Point* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (!V()->Sep(f, t))
            return false_update(f, t);
        if (!B()->Sep(f, t))
            return false_update(f, t);
        if (!C()->Sep(f, t))
            return false_update(f, t);
        return true_update(f, t);
    }
    bool Sep(Edge* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (!V()->Sep(f, t))
            return false_update(f, t);
        if (!B()->Sep(f, t))
            return false_update(f, t);
        if (!C()->Sep(f, t))
            return false_update(f, t);
        return true_update(f, t);
    }
    bool Sep(Triangle* f, double t = 0)
    {
        if (contains(f, t))
            return false;
        if (ncontains(f, t))
            return true;
        if (!V()->Sep(f, t))
            return false_update(f, t);
        if (!B()->Sep(f, t))
            return false_update(f, t);
        if (!C()->Sep(f, t))
            return false_update(f, t);
        return true_update(f, t);
    }
};

//*************************************************************************//

#endif