// This file define the class SO3 and its transformations.

#ifndef SO3_H
#define SO3_H

#define M_PI           3.14159265358979323846       /* pi */
#define doublemin 1e-8
#include<cmath>
#include<Eigen/Dense>
using namespace Eigen;

// 3d Vector Lie Bracket.
Matrix3d adv(Vector3d v)
{
    Matrix3d Adv;
    Adv << 0, -v(2), v(1),
        v(2), 0, -v(0),
        -v(1), v(0), 0;
    return Adv;
}

// Turn angle into arc.
double angle_to_arc(double theta)
{
    return theta * M_PI / 180;
}

//Turn arc into angle.
double arc_to_angle(double theta)
{
    return theta * 180 / M_PI;
}

// SO3 point represented by unit quaternion
class SO3
{
    Vector4d V;                                     // Stored by R^4 coordinate without normalization.
public:
    // For empty arguments.
    SO3()
    {
        V = Vector4d(1, 0, 0, 0);
    }
    // Construct from R^4 vector.
    SO3(double _a, double _b, double _c, double _d)
    {
        Vector4d _q(_a, _b, _c, _d);
        assert(_q.norm() > 0);
        V = _q;
    }
    // Construct from R^4 vector.
    SO3(Vector4d _q)
    {
        assert(_q.norm() > 0);
        V = _q;
    }
    // Construct from Orthogonal Matrix.
    SO3(Matrix3d r)
    {
        double a, b, c, d;
        a = sqrt(r.trace() + 1) / 2;
        if (abs(a) < doublemin)
        {
            b = sqrt(abs((r(0, 0) + 1) / 2 - a * a));
            c = sqrt(abs((r(1, 1) + 1) / 2 - a * a));
            d = sqrt(abs((r(2, 2) + 1) / 2 - a * a));
            if (r(1, 0) + r(0, 1) < 0)
                c = -c;
            if (r(2, 0) + r(0, 2) < 0)
                d = -d;
            if (r(2, 1) - r(1, 2) < 0)
                b = -b;
        }
        else
        {
            b = (r(2, 1) - r(1, 2)) / (4 * a);
            c = (r(0, 2) - r(2, 0)) / (4 * a);
            d = (r(1, 0) - r(0, 1)) / (4 * a);
        }
        V = Vector4d(a, b, c, d);
    }
    /*SO3(Matrix3d _r)
    {
        assert((_r.transpose() * _r).isApprox(Matrix3d::Identity()));   // Check if it is orthogonal.
        if (_r.isApprox(Matrix3d::Identity()))                          // If it is identity matrix.
        {
            V = Vector4d(1, 0, 0, 0);
            return;
        }
        Vector3d v = FullPivLU<Matrix3d>(_r - Matrix3d::Identity()).kernel().col(0).normalized();
        double cos_tau = 1 + (_r.trace() - 3) / 2;
        double sin_tau = -((_r * adv(v)).trace()) / 2;
        if (abs(sin_tau) < doublemin)
            V = Vector4d(0, v(0), v(1), v(2));
        else
            V = Vector4d((1 + cos_tau) / sin_tau, v(0), v(1), v(2));    // a=cot(theta/2).
    }*/
    // Construct from Axis and Angle.
    SO3(Vector3d v, double tau)
    {
        tau = angle_to_arc(tau);
        if (v.norm() < doublemin)
        {
            V = Vector4d(1, 0, 0, 0);
            return;
        }
        v = v.normalized();
        if (abs(sin(tau)) >= doublemin)
            V = Vector4d((1 + cos(tau)) / sin(tau), v(0), v(1), v(2));
        else if (abs(tau) < 1)
            V = Vector4d(1, 0, 0, 0);
        else
            V = Vector4d(0, v(0), v(1), v(2));
    }
    // Construct from pA = (cos(Aphi) , cos(Atheta) * sin(Aphi), sin(Atheta) * sin(Aphi)) and 
    // pB = ( -sin(Aphi) * cos(Btheta), cos(Atheta) * cos(Aphi) * cos(Btheta) - sin(Atheta) * sin(Btheta),	
    // sin(Atheta) * cos(Aphi) * cos(Btheta) + cos(Atheta) * sin(Btheta))
    // phi = Aphi, theta = Atheta, xi = Btheta
    SO3(double phi, double theta, double xi)
    {
        Matrix3d r;
        phi = angle_to_arc(phi);
        theta = angle_to_arc(theta);
        xi = angle_to_arc(xi);
        double cp = cos(phi), sp = sin(phi);
        double ct = cos(theta), st = sin(theta);
        double cx = cos(xi), sx = sin(xi);
        r << cp, -sp * cx, sp* sx,
            sp* ct, cp* ct* cx - st * sx, -cp * ct * sx - st * cx,
            sp* st, cp* st* cx + ct * sx, -cp * st * sx + ct * cx;
        *this = SO3(r);
    }
    /*
    SO3(double phi, double theta, double xi)
    {
        if (Vector3d(phi, theta, xi).isApprox(Vector3d::Zero()))
        {
            V = Vector4d(1, 0, 0, 0);
            return;
        }
        Vector3d v = Vector3d(-(cos(theta) * sin(xi) + sin(theta) * cos(xi)) * cos(phi), (sin(xi) - sin(theta)) * (cos(phi) - 1), (cos(xi) + cos(theta)) * (cos(phi) - 1)).normalized();
        double cos_tau = 1 + (cos(phi) + cos(theta) * cos(phi) * cos(xi) - sin(theta) * sin(xi) - sin(theta) * cos(phi) * sin(xi) + cos(theta) * cos(xi) - 3) / 2;
        double sin_tau = sqrt(1 - cos_tau * cos_tau);
        if (abs(sin_tau) < doublemin)
            V = Vector4d(0, v(0), v(1), v(2));
        else
            V = Vector4d((1 + cos_tau) / sin_tau, v(0), v(1), v(2));    // a=cot(theta/2).
    }*/
    // Quaternion representation.
    Vector4d Q()
    {
        return V.normalized();
    }
    // \wh{SO(3)} representation.
    Vector4d H()
    {
        Vector4d q = Q();
        return q / q.cwiseAbs().maxCoeff();
    }
    // Matrix representaiton.
    Matrix3d R()
    {
        Vector4d q = Q();
        Matrix3d _R;
        _R << 2 * (q(0) * q(0) + q(1) * q(1)) - 1, 2 * (q(1) * q(2) - q(0) * q(3)), 2 * (q(1) * q(3) + q(0) * q(2)),
            2 * (q(1) * q(2) + q(0) * q(3)), 2 * (q(0) * q(0) + q(2) * q(2)) - 1, 2 * (q(2) * q(3) - q(0) * q(1)),
            2 * (q(1) * q(3) - q(0) * q(2)), 2 * (q(2) * q(3) + q(0) * q(1)), 2 * (q(0) * q(0) + q(3) * q(3)) - 1;
        return _R;
    }
    // Axis - Angle representation.
    void Axis_Angle(Vector3d &v,double &theta)
    {
        v = Vector3d(V(1), V(2), V(3)).normalized();
        theta = 2 * acos(V.normalized()(0));
        if (theta > M_PI)
            theta -= 2 * M_PI;
        theta = arc_to_angle(theta);
    }
    // Axis of the rotation.
    Vector3d Axis()
    {
        Vector3d v(V(1), V(2), V(3));
        return v.normalized();
    }
    // Sine angle of quaterion.
    double Sin_theta()
    {
        Vector3d v(V(1), V(2), V(3));
        Vector4d W = V.normalized();
        Vector3d w(W(1), W(2), W(3));
        return w.norm();
    }
    // Cosine angle of quaterion.
    double Cos_theta()
    {
        Vector4d W = V.normalized();
        return W(0);
    }
    // Sine angle of the rotation.
    double Sin_tau()
    {
        return 2 * Cos_theta() * Sin_theta();
    }
    // Cosine angle of the rotation.
    double Cos_tau()
    {
        return 2 * Cos_theta() * Cos_theta() - 1;
    }
    // This SO3 element acts on a vector.
    Vector3d act_on(Vector3d v)
    {
        return R() * v;
    }
    // If two SO3 elements are the same.
    bool operator==(SO3 _Q)
    {
        Vector4d _V = _Q.Q();
        return (V(0) * _V(1) == V(1) * _V(0)) && (V(0) * _V(2) == V(2) * _V(0)) && (V(0) * _V(3) == V(3) * _V(0));
    }
};

#endif