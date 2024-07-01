// This file defines the features and primitives.
#include<Eigen/Dense>
using namespace Eigen;

VectorXd Schmit(MatrixXd basis, VectorXd v)
{
    assert(basis.cols() == v.size());
    for (int i = 0; i < basis.rows(); ++i)
        v = v - v.dot(basis.row(i)) * basis.row(i).transpose() / basis.row(i).squarenorm();
    return v;
}

class Feature
{
public:
    virtual bool contains(Vector3d v) = 0;                             // If v is an interior point of the feature.
    virtual bool intersecting(Solid* S) = 0;
};

class Vertex :public Feature
{
    Vector3d V;
public:
    Vertex()
    {
        V = Vector3d::Zero();
    }
    Vertex(Vector3d v)
    {
        V = v;
    }
    Vector3d cor()
    {
        return V;
    }
    void set(Vector3d v)
    {
        V = v;
    }
    double distance(vector3d v)
    {
        return (v - V).norm();
    }
    double distance(Vertex v)
    {
        return distance(v.cor());
    }
    bool contains(Vector3d v)
    {
        return false;
    }
    bool intersecting(Solid* S);
};

class Edge :public Feature
{
    Matrix<double, 2, 3> V;
public:
    Edge()
    {
        V = Matrix<double, 2, 3>::Zero();
    }
    Edge(Matrix<double, 2, 3> _V)
    {
        V = _V;
    }
    Edge(Vector3d v_1, Vector3d v_2)
    {
        V.row(0) = v_1;
        V.row(1) = v_2;
    }
    Vertex vertex(int i)
    {
        return Vertex(V.row(i % 2));
    }
    Vector3d cor(int i)
    {
        return V.row(i % 2);
    }
    void set(Vector3d v, int i)
    {
        V.row(i % 2) = v;
    }
    Vector3d dir()                                                  // Direction vector.
    {
        return (V.row(1) - V.row(0)).normalized();
    }
    double length()                                                      // Length of the edge.
    {
        return (V.row(1) - V.row(0)).norm();
    }
    double line_distance(Vertex v)                                       // Line distance to v.
    {
        return dir().cross(v.cor().transpose() - cor(0)).norm();
    }
    Vertex line_projection(Vertex v)                                     // The projection of v onto this edge.
    {
        return Vertex(cor(0) + dir().dot(v.cor().transpose() - cor(0)) * dir().transpose());
    }
    double line_distance(Edge e)
    {
        if (dir().cross(e.dir()).norm() == 0)                           // Parallel means each point is an intersection.
            return Schimit(dir().transpose(), e.cor(0) - cor(0)).norm();
        Vector3d Normal = dir().cross(e.dir()).normalized();
        return abs(Normal.dot(e.cor(0) - cor(0)));
    }
    Vertex line_intersection_projection(Edge e)                         // Nearest point on the line of this edge to line of edge e.
    {
        if (dir().cross(e.dir()).norm() == 0)                           // Parallel means each point is an intersection.
            return Vertex(0);
        Vector3d Normal = dir().cross(e.dir()).normalized();
        Vector3d P = e.cor(0) - Normal.dot(e.cor(0) - cor(0)) * Normal.transpose();   // Project e.V(0) onto the plane contains this edge and normal with the Normal.
        // Find Q such that P-Q=k_1*e.dir() and V.row(0)-Q=k_2*dir().
        // This is equivalent to that P-V.row(0)=k_1*e.dir()-k_2*dir().
        // So [e.D(), D(), Normal]*[k_1,-k_2,0]^T=P-V.row(0).
        Matrix3d Local_Chart;
        Local_Chart.col(0) = e.dir();
        Local_Chart.col(1) = dir();
        Local_Chart.col(2) = Normal;
        Vector3d Coeff = Local_Chart.inverse() * (P - V.row(0).transpose());
        return Vertex(P - Coeff(0) * e.dir());
    }
    double distance(Vertex v)                                            // Segment distance to v.
    {
        if (is_on(line_projection(v)))
            return line_distance(v);
        else
            return min(vertex(0).distance(v), vertex(1).distance(v));
    }
    double distance(Edge e)
    {
        Vertex v = line_intersection_projection(e);
        if (is_on(v))
            return e.distance(v);
        else
            return min(e.distance(vertex(0)), e.distance(vertex(1)));
    }
    bool contains(Vector3d v)                                           // Point never inside the interior of edge.
    {
        return false;
    }
    bool is_on(Vertex v)                                                 // If vertex is on the edge.
    {
        return Edge(V.row(0), v.V).length() + Edge(v.V, V.row(1)).length() < length() + num_min;
    }
    bool intersecting(Solid* S);
};

class Triangle :public Feature
{
public:
    Matrix3d V;                                                     // Face set as (0,1,2) by default.
    Triangle()
    {
        V = Matrix3d::Zero();
    }
    Triangle(Matrix3d _V)
    {
        V = _V;
    }
    Triangle(Vector3d v0, Vector3d v1, Vector3d v2)
    {
        assert((v2 - v0).cross(v1 - v0).norm() > 0);
        V.row(0) = v0;
        V.row(1) = v1;
        V.row(2) = v2;
    }
    Vertex vertex(int i)
    {
        return Vertex(V.row(i % 3));
    }
    Edge edge(int i)
    {
        return Edge(V.row(i % 3), V.row((i + 1) % 3));
    }
    Vector3d Normal()
    {
        return (V.row(2) - V.row(0)).cross(V.row(1) - V.row(0)).normalized();
    }
    double area()                                                   // Length of the triangle.
    {
        return (V.row(2) - V.row(0)).cross(V.row(1) - V.row(0)).norm() / 2;
    }
    double plane_distance(Vertex v)                                 // Plane distnace to v.
    {
        return abs(Normal().dot(v.V.transpose() - V.row(0)));
    }
    Vertex plane_projection(Vertex v)                                     // Plane projection of v.
    {
        return Vertex(v.V - Normal().dot(v.V.transpose() - V.row(0)) * Normal());
    }
    double distance(Vertex v)                                       // Triangle distance to v.
    {
        if (is_on(plane_projection(v)))
            return plane_distance(v);
        else
            return min({ edge(0).distance(v),edge(1).distance(v),edge(2).distance(v) });
    }
    bool contains(Vector3d v)                                      // Point never inside the interior of triangle.
    {
        return false;
    }
    bool is_infront(Vector3d v)
    {
        return (v.transpose() - V.row(0)).dot((V.row(1) - V.row(0)).cross(V.row(2) - V.row(0))) > 0;
    }
    bool is_on(Vertex v)                                            // If vertex is on the triangle.
    {
        return Triangle(V.row(0), V.row(1), v.V).area() + Triangle(V.row(0), v.V, V.row(2)).area() + Triangle(v.V, V.row(1), V.row(2)).area() < area() + num_min;
    }
    bool intersects(Edge E)
    {
        double cos_cross_angle = Normal().dot(E.dir());
        if (cos_cross_angle == 0)
            return false;
        Vector3d intersection = E.V.row(0) - (distance(E.vertex(0)) / cos_cross_angle) * E.dir().transpose();
        if ((intersection.transpose() - E.V.row(0)).dot(E.V.row(1) - E.V.row(0)) <= 0)
            return false;
        if ((intersection.transpose() - E.V.row(1)).dot(E.V.row(0) - E.V.row(1)) <= 0)
            return false;
        Matrix3d local_chart;
        local_chart.row(0) = V.row(1) - V.row(0);
        local_chart.row(1) = V.row(2) - V.row(0);
        local_chart.row(2) = intersection.transpose() - V.row(0);
        Vector3d local_chart_solution = FullPivLU<Matrix3d>(local_chart.transpose()).kernel().col(0).normalized();
        local_chart_solution = -local_chart_solution / local_chart_solution(2);
        return (local_chart_solution(0) > 0) && (local_chart_solution(1) > 0) && (local_chart_solution(0) + local_chart_solution(1) < 1);
    }
    bool is_line_collide(Edge E)
    {
        double cos_cross_angle = Normal().dot(E.dir());
        if (cos_cross_angle == 0)
            return false;
        Vector3d intersection = E.V.row(0) - (distance(E.vertex(0)) / cos_cross_angle) * E.dir().transpose();
        Matrix3d local_chart;
        local_chart.row(0) = V.row(1) - V.row(0);
        local_chart.row(1) = V.row(2) - V.row(0);
        local_chart.row(2) = intersection.transpose() - V.row(0);
        Vector3d local_chart_solution = FullPivLU<Matrix3d>(local_chart.transpose()).kernel().col(0).normalized();
        local_chart_solution = -local_chart_solution / local_chart_solution(2);
        return (local_chart_solution(0) > 0) && (local_chart_solution(1) > 0) && (local_chart_solution(0) + local_chart_solution(1) < 1);
    }
    bool intersecting(Solid* S);
};

class Mesh :public Feature                                          // Require to be convex.
{
    Matrix<double, -1, 3> V;                                        // Vertices
    Matrix<int, -1, 3> F;                                           // Faces defined by vertices, must sequenced by dir.
    MatrixXi E;
public:
    Mesh()
    {
        E.resize(0, 0);
    }
    Mesh(MatrixXd _V, MatrixXi _F)
    {
        assert(_V.cols() == 3);
        assert(_F.cols() == 3);
        assert(_F.maxCoeff() < _V.rows());
        V = _V;
        F = _F;
        edges(F, E);
    }
    int V_size()
    {
        return V.rows();
    }
    int E_size()
    {
        return E.rows();
    }
    int F_size()
    {
        return F.rows();
    }
    Vertex vertex(int i)
    {
        return Vertex(V.row(i % V.rows()));
    } 
    Edge edge(int i)
    {
        return Edge(V.row(E(i % E.rows(), 0)), V.row(E(i % E.rows(), 1)));
    }
    Triangle Face(int i)                                            // Defines a face.
    {
        return Triangle(V.row(F(i % F.rows(), 0)), V.row(F(i % F.rows(), 1)), V.row(F(i % F.rows(), 2)));
    }
    bool contains(Vector3d v)                                      // Normal vector towards outside, require mesh to be convex.
    {
        for (Index i = 0; i < F.rows(); ++i)
            if (Face(i).is_infront(v))
                return false;
        return true;
    }
    bool intersecting(Solid* S);
};

class Circle                                                    // z'==0, x'^2+y'^2<=r^2, not assuming this as a feature now.
{
    Vector3d Center;
    Vector3d Normal;                                            // Have to be normalized.
    double Radius;
    Matrix3d Transition()                                       // Transform x'-y'-z' coordinate to x-y-z coordinate.
    {
        Vector3d new_z = Normal;
        Vector3d new_x = Vector3d(0, new_z(2), -new_z(1)).normalized();
        Vector3d new_y = new_z.cross(new_x);
        Matrix3d trans;
        trans.col(0) = new_x;
        trans.col(1) = new_y;
        trans.col(2) = new_z;
        return trans;
    }
    Matrix3d Inv_Transition;                                    // Transform x-y-z coordinate to x'-y'-z' coordinate.
    void initialization()
    {
        Inv_Transition = Transition().inverse();
    }
public:
    Circle(Vector3d _Center = Vector3d::Zero(), Vector3d _Normal = Vector3d(0, 0, 1), double _Radius = 1)
    {
        Center = _Center;
        Normal = _Normal.normalized();
        Radius = _Radius;
        initialization();
    }
    Vector3d Local(Vector3d v)                                  // x'-y'-z' coordinate.
    {
        return Inv_Transition * (v - Center);
    }
    bool intersects(Vertex Phi)
    {
        return (Phi.V - Center).dot(Normal) == 0 && (Phi.V - Center).norm() <= Radius;
    }
    bool intersects(Edge Phi)                                           // Circle collide with edge.
    {
        Edge Local_Phi(Local(Phi.V.row(0)), Local(Phi.V.row(1)));
        if (Local_Phi.V(0, 2) * Local_Phi.V(1, 2) > 0)
            return false;
        if (Local_Phi.V(0, 2) == Local_Phi.V(1, 2))
            return Phi.distance(Vertex(Vector3d::Zero())) <= Radius;
        Vector3d Local_Base = Local_Phi.V.row(0).transpose() - Local_Phi.V(0, 2) / Local_Phi.dir()(2) * Local_Phi.dir();
        return Local_Base.norm() <= Radius;
    }
    bool intersects(Triangle Phi)                                       // Circle collide with triangle.
    {
        for (Index i = 0; i < 3; ++i)
            if (intersects(Phi.edge(i)))
                return true;
        Triangle Local_Phi(Local(Phi.V.row(0)), Local(Phi.V.row(1)), Local(Phi.V.row(2)));
        if ((Local_Phi.V(0, 2) * Local_Phi.V(1, 2) > 0) && (Local_Phi.V(0, 2) * Local_Phi.V(2, 2) > 0))
            return false;
        if (Local_Phi.Normal().cross(Vector3d(0, 0, 1)).norm() == 0)
            return Local_Phi.is_on(Vertex(Vector3d::Zero()));
        vector<Vector3d> Local_Inter_Vertex;
        if (Local_Phi.V(0, 2) * Local_Phi.V(1, 2) < 0)
            Local_Inter_Vertex.push_back(Local_Phi.V.row(0).transpose() - Local_Phi.V(0, 2) / Local_Phi.edge(0).dir()(2) * Local_Phi.edge(0).dir());
        if (Local_Phi.V(1, 2) * Local_Phi.V(2, 2) < 0)
            Local_Inter_Vertex.push_back(Local_Phi.V.row(1).transpose() - Local_Phi.V(1, 2) / Local_Phi.edge(1).dir()(2) * Local_Phi.edge(1).dir());
        if (Local_Phi.V(2, 2) * Local_Phi.V(0, 2) < 0)
            Local_Inter_Vertex.push_back(Local_Phi.V.row(2).transpose() - Local_Phi.V(2, 2) / Local_Phi.edge(2).dir()(2) * Local_Phi.edge(2).dir());
        if (Local_Phi.V(0, 2) == 0)
            Local_Inter_Vertex.push_back(Local_Phi.V.row(0));
        if (Local_Phi.V(1, 2) == 0)
            Local_Inter_Vertex.push_back(Local_Phi.V.row(1));
        if (Local_Phi.V(2, 2) == 0)
            Local_Inter_Vertex.push_back(Local_Phi.V.row(2));
        if (Local_Inter_Vertex.size() != 2)                     // Triangle intersecting plane x'-O'-y' singularly.
            return false;
        return Edge(Local_Inter_Vertex[0], Local_Inter_Vertex[1]).distance(Vertex(Vector3d::Zero())) <= Radius;
    }
};

/*------------------------------Semi-Algebraic Sets Methods------------------------------*/
// CSG

class Solid                                                     // General class of closed convex semi-algebraic set.
{
public:
    virtual bool intersects(Vertex Phi, double d = 0) = 0;          // Define virtual function for collision detections with vertex feature.
    virtual bool intersects(Edge Phi, double d = 0) = 0;            // Define virtual function for collision detections with edge feature.
    virtual bool intersects(Triangle Phi, double d = 0) = 0;        // Define virtual function for collision detections with triangle feature.
    virtual bool intersects(Mesh Phi, double d = 0) = 0;            // Define virtual function for collision detections with convex mesh feature.
    virtual Vector3d interior_sample() = 0;                         // Define virtual function for FREE/STUCK detection.
};

bool Vertex::intersecting(Solid* S, double d = 0)
{
    return S->intersects(*this, d);
}
bool Edge::intersecting(Solid* S, double d = 0)
{
    return S->intersects(*this, d);
}
bool Triangle::intersecting(Solid* S, double d = 0)
{
    return S->intersects(*this, d);
}
bool Mesh::intersecting(Solid* S, double d = 0)
{
    return S->intersects(*this, d);
}

class Union_Alge :public Solid                                  // This union requires to be still convex.
{
    vector<Solid*> Primitive;
public:
    Union_Alge() {};
    void add(Solid* Alge)
    {
        Primitive.push_back(Alge);
    }
    bool intersects(Vertex Phi, double d = 0)
    {
        for (Index i = 0; i < Primitive.size(); ++i)
            if (Primitive[i]->intersects(Phi, d))
                return true;
        return false;
    }
    bool intersects(Edge Phi, double d = 0)
    {
        for (Index i = 0; i < Primitive.size(); ++i)
            if (Primitive[i]->intersects(Phi, d))
                return true;
        return false;
    }
    bool intersects(Triangle Phi, double d = 0)
    {
        for (Index i = 0; i < Primitive.size(); ++i)
            if (Primitive[i]->intersects(Phi, d))
                return true;
        return false;
    }
    bool intersects(Mesh Phi, double d = 0)
    {
        for (Index i = 0; i < Primitive.size(); ++i)
            if (Primitive[i]->intersects(Phi, d))
                return true;
        return false;
    }
    Vector3d interior_sample()
    {
        return Primitive[0]->interior_sample();
    }
};

class RPrism : public Solid                                 // This is the class for our specialized restricted-prism primitive.
{                                                           // It is constructed by connecting two triangles with paralleled edges.
    Matrix3d Centers;                                       // Center triangle, which is the mirror plane.
    double half_length[3];                                  // The half of length of the three paralleled edges.
public:
    RPrism()
    {
        Centers = Matrix3d::Zero();
        half_length[0] = 0;
        half_length[1] = 0;
        half_length[2] = 0;
    }
    RPrism(Vector3d v1, Vector3d v2, Vector3d v3, double l1, double l2, double l3)
    {
        Centers.row(0) = v1;
        Centers.row(1) = V2;
        Centers.row(2) = v3;
        half_length[0] = l1;
        half_length[1] = l2;
        half_length[2] = l3;
    }
    bool intersects(Vertex Phi, double d = 0)
    {
        return false;
    }
    bool intersects(Edge Phi, double d = 0)
    {
        return false;
    }
    bool intersects(Triangle Phi, double d = 0)
    {
        return false;
    }
    bool intersects(Mesh Phi, double d = 0)
    {
        return false;
    }
    Vector3d interior_sample()
    {
        return (Centers.row(0) + Centers.row(1) + Centers.row(2)) / 3;
    }
};

class Ball :public Solid
{
    Vector3d Center;
    double radius;
public:
    Ball(Vector3d _The_Vertex = Vector3d::Zero(), double _r = 1.0)
    {
        Center = Vector3d(0, 0, 0);
        radius = 1;
    }
    bool intersects(Vertex Phi, double d = 0)
    {
        return false;
    }
    bool intersects(Edge Phi, double d = 0)
    {
        return false;
    }
    bool intersects(Triangle Phi, double d = 0)
    {
        return false;
    }
    bool intersects(Mesh Phi, double d = 0)
    {
        return false;
    }
    Vector3d interior_sample()
    {
        return Center;
    }
};

class Cone :public Solid                                    // Cone x'^2+y'^2-k^2z'^2<=0, z'>=0, z'<=d.
{
    Vector3d The_Vertex;
    Vector3d Bottom_Center;
    double Bottom_Radius;
    Matrix3d Inv_Transition;                                    // Transform x-y-z coordinate to x'-y'-z' coordinate.
    Edge Axis;
    double k, d;
    Matrix3d Transition()                                       // Transform x'-y'-z' coordinate to x-y-z coordinate.
    {
        Vector3d new_z = (Bottom_Center - The_Vertex).normalized();     // Coordinate of z' under x-y-z coordinate.
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
    void initialization()
    {
        Axis = Edge(The_Vertex, Bottom_Center);
        Inv_Transition = Transition().inverse();
        d = Axis.length();
        k = Bottom_Radius / d;
    }
public:
    Cone(Vector3d _The_Vertex = Vector3d::Zero(), Vector3d _Bottom_Center = Vector3d(0, 0, 1), double _Bottom_Radius = 1)
    {
        assert(!_The_Vertex.isApprox(_Bottom_Center));
        assert(_Bottom_Radius > 0);
        The_Vertex = _The_Vertex;
        Bottom_Center = _Bottom_Center;
        Bottom_Radius = _Bottom_Radius;
        initialization();
    }
    Vector3d Local(Vector3d v)                                  // x'-y'-z' coordinate.
    {
        return Inv_Transition * (v - The_Vertex);
    }
    bool intersects(Vertex Phi)
    {
        Vertex Proj_Phi = Axis.line_projection(Phi);
        return (Axis.is_on(Proj_Phi)) && (Axis.line_distance(Phi) <= k * Edge(The_Vertex, Proj_Phi.V).length());
    }
    bool intersects(Edge Phi)                                   // Solve x'^2+y'^2-k^2z'^2==0 && (x'-x'_0)/A'==(y'-y'_0)/B'==(z'-z'_0)/C' and check if 0<=z'<=d and z' is between endpoints of edge.
    {
        //if (Circle(Bottom_Center, Axis.dir(), Bottom_Radius).intersects(Phi)) // If the edge collides with bottom circle.
        //    return true;
        if (intersects(Phi.vertex(0)) || intersects(Phi.vertex(1)))
            return true;
        Vector3d Local_V_0 = Local(Phi.V.row(0));
        Vector3d Local_V_1 = Local(Phi.V.row(1));
        double low_seg_bound = min(Local_V_0(2), Local_V_1(2));
        double up_seg_bound = max(Local_V_0(2), Local_V_1(2));
        Vector3d Local_Base = Local_V_0;
        Vector3d Local_dir = Edge(Local_V_0, Local_V_1).dir();
        if (Local_dir(2) == 0)                            // Edge perpendicular to axis.
            return Local_Base(2) >= 0 && Local_Base(2) <= d && Edge(Local_V_0, Local_V_1).distance(Vertex(Vector3d(0, 0, Local_Base(2)))) <= k * Local_Base(2);
        Local_dir = Local_dir / Local_dir(2); // Set C'=1.
        Local_Base = Local_Base - Local_Base(2) / Local_dir(2) * Local_dir; // Set z'_0=0.
        // x'=x'_0+A'z', y'=y'_0+B'z'.
        // x'^2+y'^2-k^2z'^2=(x'_0+A'z')^2+(y'_0+B'z')^2-k^2z'^2=(A'^2+B'^2-k^2)z'^2+2(A'x'_0+B'y'_0)z'+(x'_0^2+y'_0^2).
        // Solve (A'^2+B'^2-k^2)z'^2+2(A'x'_0+B'y'_0)z'+(x'_0^2+y'_0^2)==0.
        double Coeff_0 = Local_dir.squaredNorm() - squa(k) - 1;
        double Coeff_1 = 2 * Local_Base.dot(Local_dir);
        double Coeff_2 = Local_Base.squaredNorm();
        if (Coeff_0 == 0)
        {
            if (Coeff_1 == 0)
                return squa(Local_dir(0)) + squa(Local_dir(1)) == 0;
            double solu_z = -Coeff_2 / Coeff_1;
            return (solu_z >= 0 && solu_z <= d && solu_z >= low_seg_bound && solu_z <= up_seg_bound);
        }
        double Delta = squa(Coeff_1) - 4 * Coeff_0 * Coeff_2;
        if (Delta < 0)
            return false;
        double solu_z_1 = (-Coeff_1 + sqrt(Delta)) / (2 * Coeff_0);
        double solu_z_2 = (-Coeff_1 - sqrt(Delta)) / (2 * Coeff_0);
        // Check if exists 0<=z'<=d and z' is between endpoints of edge.
        return (solu_z_1 >= 0 && solu_z_1 <= d && solu_z_1 >= low_seg_bound && solu_z_1 <= up_seg_bound) || (solu_z_2 >= 0 && solu_z_2 <= d && solu_z_2 >= low_seg_bound && solu_z_2 <= up_seg_bound);
    }
    bool intersects(Triangle Phi)
    {
        // Collide with the boundary of Phi.
        for (Index i = 0; i < 3; ++i)
            if (intersects(Phi.edge(i)))
                return true;
        // Collide purely with the interior of Phi.
        // Since the Gaussian curvature of cone is 0, Phi cannot collide "tangently" with cone.
        // So only need to check the finite conditions.
        // Case 1: Phi crosses over the cone.
        if (Phi.intersects(Axis))
            return true;
        // Case 2: Phi crosses the bottom circle.
        if (Circle(Bottom_Center, Axis.dir(), Bottom_Radius).intersects(Phi))
            return true;
        // Otherwise, Phi cannot collide with cone.
        return false;
    }
    bool intersects(Mesh Phi)
    {
        for (Index i = 0; i < Phi.F_size(); ++i)
            if (intersects(Phi.Face(i)))
                return true;
        return false;
    }
    Vector3d interior_sample()
    {
        return (The_Vertex + Bottom_Center) / 2;
    }
};

class Cylinder :public Solid                                // Cone x'^2+y'^2-r^2<=0, z'>=-d, z'<=d.
{
    Vector3d Bottom_Center[2];
    double Bottom_Radius;
    Edge Axis;
    double d;
    Matrix3d Inv_Transition;                                    // Transform x-y-z coordinate to x'-y'-z' coordinate.
    Matrix3d Transition()                                       // Transform x'-y'-z' coordinate to x-y-z coordinate.
    {
        Vector3d new_z = (Bottom_Center[0] - Bottom_Center[1]).normalized();     // Coordinate of z' under x-y-z coordinate.
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
    void initialization()
    {
        Axis = Edge(Bottom_Center[0], Bottom_Center[1]);
        Inv_Transition = Transition().inverse();
        d = Axis.length() / 2;
    }
public:
    Cylinder(Vector3d _The_Vertex = Vector3d(0, 0, -1), Vector3d _Bottom_Center = Vector3d(0, 0, 1), double _Bottom_Radius = 1)
    {
        assert(!_The_Vertex.isApprox(_Bottom_Center));
        assert(_Bottom_Radius > 0);
        Bottom_Center[0] = _The_Vertex;
        Bottom_Center[1] = _Bottom_Center;
        Bottom_Radius = _Bottom_Radius;
        initialization();
    }
    Vector3d Local(Vector3d v)                                  // x'-y'-z' coordinate.
    {
        return Inv_Transition * (v - (Bottom_Center[0] + Bottom_Center[1]) / 2);
    }
    bool intersects(Vertex Phi)
    {
        Vector3d Local_Phi = Local(Phi.V);
        return squa(Local_Phi(0)) + squa(Local_Phi(1)) <= squa(Bottom_Radius) && Local_Phi(2) >= -d && Local_Phi(2) <= d;
    }
    bool intersects(Edge Phi)                                   // Solve x'^2+y'^2-r^2==0 && (x'-x'_0)/A'==(y'-y'_0)/B'==(z'-z'_0)/C' and check if 0<=z'<=d and z' is between endpoints of edge.
    {
        if (intersects(Phi.vertex(0)) || intersects(Phi.vertex(1)))
            return true;
        if (Circle(Bottom_Center[0], Axis.dir(), Bottom_Radius).intersects(Phi)) // If the edge collides with bottom circle (go through two bottoms).
            return true;
        Vector3d Local_V_0 = Local(Phi.V.row(0));
        Vector3d Local_V_1 = Local(Phi.V.row(1));
        double low_seg_bound = min(Local_V_0(2), Local_V_1(2));
        double up_seg_bound = max(Local_V_0(2), Local_V_1(2));
        Vector3d Local_Base = Local_V_0;
        Vector3d Local_dir = Edge(Local_V_0, Local_V_1).dir();
        if (Local_dir(2) == 0)                            // Edge perpendicular to axis.
            return Local_Base(2) >= -d && Local_Base(2) <= d && Edge(Local_V_0, Local_V_1).distance(Vertex(Vector3d(0, 0, Local_Base(2)))) <= Bottom_Radius;
        Local_dir = Local_dir / Local_dir(2); // Set C'=1.
        Local_Base = Local_Base - Local_Base(2) / Local_dir(2) * Local_dir; // Set z'_0=0.
        // x'=x'_0+A'z', y'=y'_0+B'z'.
        // x'^2+y'^2-r^2=(x'_0+A'z')^2+(y'_0+B'z')^2-r^2=(A'^2+B'^2)z'^2+2(A'x'_0+B'y'_0)z'+(x'_0^2+y'_0^2-r^2).
        // Solve (A'^2+B'^2)z'^2+2(A'x'_0+B'y'_0)z'+(x'_0^2+y'_0^2-r^2)==0.
        double Coeff_0 = Local_dir.squaredNorm() - 1;
        double Coeff_1 = 2 * Local_Base.dot(Local_dir);
        double Coeff_2 = Local_Base.squaredNorm() - squa(Bottom_Radius);
        if (Coeff_0 == 0)
        {
            if (Coeff_1 == 0)
                return squa(Local_dir(0)) + squa(Local_dir(1)) == 0;
            double solu_z = -Coeff_2 / Coeff_1;
            return (solu_z >= 0 && solu_z <= d && solu_z >= low_seg_bound && solu_z <= up_seg_bound);
        }
        double Delta = squa(Coeff_1) - 4 * Coeff_0 * Coeff_2;
        if (Delta < 0)
            return false;
        double solu_z_1 = (-Coeff_1 + sqrt(Delta)) / (2 * Coeff_0);
        double solu_z_2 = (-Coeff_1 - sqrt(Delta)) / (2 * Coeff_0);
        // Check if exists 0<=z'<=d and z' is between endpoints of edge.
        return (solu_z_1 >= -d && solu_z_1 <= d && solu_z_1 >= low_seg_bound && solu_z_1 <= up_seg_bound) || (solu_z_2 >= -d && solu_z_2 <= d && solu_z_2 >= low_seg_bound && solu_z_2 <= up_seg_bound);
    }
    bool intersects(Triangle Phi)
    {
        // Collide with the boundary of Phi.
        for (Index i = 0; i < 3; ++i)
            if (intersects(Phi.edge(i)))
                return true;
        // Collide purely with the interior of Phi.
        // Since the Gaussian curvature of cylinder is 0, Phi cannot collide "tangently" with cylinder.
        // So only need to check the finite conditions.
        // Case 1: Phi crosses over the cylinder.
        if (Phi.intersects(Axis))
            return true;
        // Case 2: Phi crosses the bottom circles.
        if (Circle(Bottom_Center[0], Axis.dir(), Bottom_Radius).intersects(Phi) || Circle(Bottom_Center[1], Axis.dir(), Bottom_Radius).intersects(Phi))
            return true;
        // Otherwise, Phi cannot collide with cylinder.
        return false;
    }
    bool intersects(Mesh Phi)
    {
        for (Index i = 0; i < Phi.F_size(); ++i)
            if (intersects(Phi.Face(i)))
                return true;
        return false;
    }
    Vector3d interior_sample()
    {
        return (Bottom_Center[0] + Bottom_Center[1]) / 2;
    }
};

class Frustum :public Solid
{

};

// Box3

/*-------------------------------These parts are abandoned-------------------------------*/

class Half_Plane :public Solid                                  // v^Tx+c <= 0 
{
public:
    Vector3d v;
    double c;
    Half_Plane(Vector3d _v = Vector3d(0, 0, 1), double _c = 0)
    {
        assert(_v.norm() > 0);
        v = _v;
        c = _c;
    }
    bool intersects(Vertex Phi)
    {
        return v.dot(Phi.V) + c <= 0;
    }
    bool intersects(Edge Phi)
    {
        for (int i = 0; i < 2; ++i)
            if (v.dot(Phi.V.row(i)) + c <= 0)
                return true;
        return false;
    }
    bool intersects(Triangle Phi)
    {
        for (int i = 0; i < 3; ++i)
            if (v.dot(Phi.V.row(i)) + c <= 0)
                return true;
        return false;
    }
    bool intersects(Mesh Phi)
    {
        for (int i = 0; i < Phi.V_size(); ++i)
            if (v.dot(Phi.vertex(i).V) + c <= 0)
                return true;
        return false;
    }
    Vector3d interior_sample()
    {
        return -((c + 1) / v.squaredNorm()) * v;
    }
};


// Algebraic operations for general convex quadratic algebraic-sets is testing now.

#define type_emptyset 0                                             // Empty set.
#define type_ellipse 1                                              // Ellipse or Point (have to contain center in order to be convex).
#define type_cone 2                                                 // Cone (have to contain axis in order to be convex).
#define type_hyperboloid_2 3                                        // Hperbloid of Two-Sheets (have not to contain center in order to be convex).
#define type_whole_space 4                                          // Whole space.
#define type_hyperboloid_1 -1                                       // Hyperbloid of One-Sheet (this will never form a convex set).
#define type_non_convex -1                                          // Non-convex cases.
#define type_impossible -1                                          // For default.

class Quadratic :public Solid                                   // x^TAx+2b^Tx+c<=0 for irreducible symmetric A.
{                                                                   // This class includes ellipses, hyperboloids and cones.
    Vector4d eigens()                                               // Coefficients for standard form ax^2+by^2+cz^2+d
    {
        Vector3cd eigen_A = A.eigenvalues();
        return Vector4d(eigen_A(0).real(), eigen_A(1).real(), eigen_A(2).real(), c + b.dot(Center()));
    }
    int type()                                                      // Check the type of the irreducible quadratic set.
    {
        int positive_inertia_index = 0;
        for (Index i = 0; i < 3; ++i)
            if (Eigens(i) > 0)
                positive_inertia_index++;
        switch (positive_inertia_index)
        {
        case 0:if (Eigens(3) <= 0) return type_whole_space; else return type_non_convex;
        case 1:return type_non_convex;
        case 2:if (Eigens(3) > 0) return type_hyperboloid_2; else if (Eigens(3) < 0)return type_hyperboloid_1; else return type_cone;
        case 3:if (Eigens(3) > 0) return type_emptyset; else return type_ellipse;
        default:return type_impossible;
        }
    }
    Vector3d axis()
    {
        switch (Type)
        {
        case type_ellipse:return FullPivLU<Matrix3d>(A - Eigens(0) * Matrix3d::Identity()).kernel().col(0).normalized();
        case type_cone:
        case type_hyperboloid_2:for (Index i = 0; i < 3; ++i)
            if (Eigens(i) < 0)
                return FullPivLU<Matrix3d>(A - Eigens(i) * Matrix3d::Identity()).kernel().col(0).normalized();
        default:return Vector3d(0, 0, 0);
        }
    }
public:
    Matrix3d A;
    Vector3d b;
    double c;
    int Type;
    Vector4d Eigens;
    Vector3d Axis;
    Quadratic(Matrix3d _A, Vector3d _b, double _c)
    {
        assert(_A.isApprox(_A.transpose()) && _A.determinant() != 0);
        A = _A;
        b = _b;
        c = _c;
        Eigens = eigens();
        Type = type();
        if (Type < 0)
            cerr << "Warning: Defining a non-convex semi-algebraic set." << endl;
        if (Type == type_emptyset)
            cerr << "Warning: Defining an empty semi-algebraic set." << endl;
        if (Type == type_whole_space)
            cerr << "Warning: Defining the whole space as a semi-algebraic set." << endl;
        Axis = axis();
    }
    Vector3d Center()                                           // Center for ellipses and hyperboloids, vertex for cones.
    {
        return -A.inverse() * b;
    }
    bool intersects(Vertex Phi)
    {
        return Phi.V.dot(A * Phi.V) + 2 * b.dot(Phi.V) + c <= 0;
    }
    bool intersects(Edge Phi)
    {
        return true;
    }
    bool intersects(Triangle Phi)
    {
        return true;
    }
    bool intersects(Mesh Phi)
    {
        return true;
    }
    Vector3d interior_sample()
    {
        switch (Type)
        {
        case type_ellipse:return Center();
        case type_cone:return Center() + Axis;
        case type_hyperboloid_2:return Center() + 2 * sqrt(Eigens(3) / max({ abs(Eigens(0)),abs(Eigens(1)),abs(Eigens(2)) })) * Axis;
        case type_whole_space:return Center();
        default:return Center();
        }
    }
};

Quadratic ball(Vector3d center, double radius)                  // (x-c)^T(x-c)-r^2<=0
{
    return Quadratic(Matrix3d::Identity(), -center, -squa(radius) + center.squaredNorm());
}

Quadratic inverse_ball(Vector3d center, double radius)           // (x-c_x)^2+(y-c_y)^2-(z-c_z)^2+r^2<=0
{
    return Quadratic(Matrix3d::Identity() - 2 * Vector3d(0, 0, 1) * Vector3d(0, 0, 1).transpose(), Vector3d(-center(0), -center(1), center(2)), squa(radius) + squa(center(0)) + squa(center(1)) - squa(center(2)));
}

class Intersect_Alge : public Solid
{

};