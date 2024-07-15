// This file define the class SO3 and its transformations.

#ifndef SO3_H
#define SO3_H

#define M_PI           3.14159265358979323846       /* pi */
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

// SO3 point represented by unit quaternion
class SO3
{
    Vector4d V;                                     // Stored by R^4 coordinate without normalization.
    double r()                                      // Norm of the Storage
    {
        return V.norm();
    }
public:
    SO3(double _a = 1.0, double _b = 0.0, double _c = 0.0, double _d = 0.0)
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
    SO3(Matrix3d _r)
    {
        assert((_r.transpose() * _r).isApprox(Matrix3d::Identity()));   // Check if it is orthogonal.
        if (_r.isApprox(Matrix3d::Identity()))                          // If it is identity matrix.
        {
            V = Vector4d(1, 0, 0, 0);
            return;
        }
        Vector3d v = FullPivLU<Matrix3d>(_r - Matrix3d::Identity()).kernel().col(0).normalized();
        double cos_theta = 1 + (_r - Matrix3d::Identity())(0, 0) / (v(1) * v(1) + v(2) * v(2));
        double sin_theta = -(_r - Matrix3d::Identity() - (1 - cos_theta) * adv(v) * adv(v))(0, 1) / v(2);
        V = Vector4d((1 + cos_theta) / sin_theta, v(0), v(1), v(2));    // a=cot(theta/2).
    }
    // Construct from Axis and Angle.
    SO3(Vector3d v, double theta)
    {
        assert(v.norm() > 0 && theta <= M_PI && theta > -M_PI);
        if (theta == 0)
            V = Vector4d(1, 0, 0, 0);
        else
            V = Vector4d(1 / tan(theta / 2), v(0), v(1), v(2));
    }
    // Quaternion representation.
    Vector4d Q()
    {
        return V.normalized();
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