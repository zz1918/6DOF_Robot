// interval.cpp : This file computes the summation, substraction,
// multiplication and division of intervals and interval matrices.
// Author: zz1918@nyu.edu

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;


// All closed intervals [inf,sup].
// The Scarlar requires to be a field with 1!=0.
template<typename Scarlar>
class interval
{                           // Interval methods for totally ordered field.
    Scarlar inf;
    Scarlar sup;
public:
    interval() { inf = sup = Scarlar(0); }
    interval(Scarlar _a)
    {
        inf = sup = _a;
    }
    interval(Scarlar _a, Scarlar _b)
    {
        if (_a < _b)
        {
            inf = _a;
            sup = _b;
        }
        else
        {
            inf = _b;
            sup = _a;
        }
    }
    Scarlar min()
    {
        return inf;
    }
    Scarlar max()
    {
        return sup;
    }
    bool contains(Scarlar x)
    {
        return (inf <= x) && (x <= sup);
    }
    bool ncontains(Scarlar x)
    {
        return !contains(x);
    }
    bool is_include(interval<Scarlar> i)
    {
        return (i.min() >= inf) && (i.max() <= sup);
    }
    bool is_included(interval<Scarlar> i)
    {
        return (i.min() <= inf) && (i.max() >= sup);
    }
    bool is_intersect(interval<Scarlar> i)
    {
        return (i.max() >= inf) || (i.min() <= sup);
    }
    bool is_number()
    {
        return inf == sup;
    }
    interval<Scarlar> operator +(interval<Scarlar> i)
    {
        return interval(inf + i.min(), sup + i.max());
    }
    interval<Scarlar> operator -(interval<Scarlar> i)
    {
        return interval(inf - i.max(), sup - i.min());
    }
    interval<Scarlar> operator *(interval<Scarlar> i)
    {
        Scarlar extremes[4] = { inf * i.min(),inf * i.max(),sup * i.min(),sup * i.max() };
        sort(extremes, extremes + 4);
        return interval<Scarlar>(extremes[0], extremes[3]);
    }
    interval<Scarlar> operator /(interval<Scarlar> i)
    {
        assert(i.ncontains(Scarlar(0)));
        Scarlar extremes[4] = { inf / i.min(),inf / i.max(),sup / i.min(),sup / i.max() };
        sort(extremes, extremes + 4);
        return interval<Scarlar>(extremes[0], extremes[3]);
    }
    void operator +=(interval<Scarlar> i)
    {
        inf += i.min();
        sup += i.max();
    }
    void operator -=(interval<Scarlar> i)
    {
        inf -= i.max();
        sup -= i.min();
    }
    void operator *=(interval<Scarlar> i)
    {
        Scarlar extremes[4] = { inf * i.min(),inf * i.max(),sup * i.min(),sup * i.max() };
        sort(extremes, extremes + 4);
        inf = extremes[0];
        sup = extremes[3];
    }
    void operator /=(interval<Scarlar> i)
    {
        assert(i.ncontains(Scarlar(0)));
        Scarlar extremes[4] = { inf / i.min(),inf / i.max(),sup / i.min(),sup / i.max() };
        sort(extremes, extremes + 4);
        inf = extremes[0];
        sup = extremes[3];
    }
    interval<Scarlar> operator +(Scarlar x)
    {
        return interval(inf + x, sup + x);
    }
    interval<Scarlar> operator -(Scarlar x)
    {
        return interval(inf - x, sup - x);
    }
    interval<Scarlar> operator *(Scarlar x)
    {
        return interval(inf * x, sup * x);
    }
    interval<Scarlar> operator /(Scarlar x)
    {
        assert(x != Scarlar(0));
        return interval(inf / x, sup / x);
    }
    void operator +=(Scarlar x)
    {
        inf += x;
        sup += x;
    }
    void operator -=(Scarlar x)
    {
        inf -= x;
        sup -= x;
    }
    void operator *=(Scarlar x)
    {
        inf *= x;
        sup *= x;
    }
    void operator /=(Scarlar x)
    {
        inf /= x;
        sup /= x;
    }
    void out(ostream& os = cout)
    {
        os << "[" << inf << "," << sup << "]";
    }
};
template<typename Scarlar>
interval<Scarlar> pow(interval<Scarlar> i, unsigned int n)
{               // Interval positive-integer-pow methods for totally ordered field.
    if (n == 0)
        return interval<Scarlar>(Scarlar(1));
    bool dep_n[32];
    Scarlar pows[32][2];
    Scarlar the_pow[2];
    int length = int(log(n * 1.0) / log(2.0));
    pows[0][0] = i.min();
    pows[0][1] = i.max();
    the_pow[0] = the_pow[1] = Scarlar(1);
    for (int i = 0, k = n; i <= length; k /= 2, ++i)
        dep_n[i] = k % 2;
    for (int i = 1; i <= length; ++i)
    {
        pows[i][0] = pows[i - 1][0] * pows[i - 1][0];
        pows[i][1] = pows[i - 1][1] * pows[i - 1][1];
    }
    for (int i = 0; i <= length; ++i)
        if (dep_n[i])
        {
            the_pow[0] *= pows[i][0];
            the_pow[1] *= pows[i][1];
        }
    if (i.ncontains(Scarlar(0)))
        return interval<Scarlar>(the_pow[0], the_pow[1]);
    else
        if (n % 2)
            return interval<Scarlar>(the_pow[0], the_pow[1]);
        else
            return interval<Scarlar>(Scarlar(0), max(the_pow[0], the_pow[1]));
}
template<typename Scarlar>
interval<Scarlar> pow(interval<Scarlar> i, int n)
{
    assert(n >= 0);
    return pow(i, (unsigned int)(n));
}
template<typename Scarlar>
interval<Scarlar> pow(interval<Scarlar> i, double n)
{
    assert(i.min() >= 0);
    if (n >= 0 && n == double(int(n)))
        return pow(i, (unsigned int)(n));
    return interval<Scarlar>(pow(i.min(), n), pow(i.max(), n));
}
template<typename Scarlar>
interval<Scarlar> pow(interval<Scarlar> i, interval<Scarlar> j)
{
    assert(j.is_number());
    return pow(i, j.min());
}
#define intervali interval<int>
#define intervalf interval<float>
#define intervald interval<double>
template<typename Scarlar>
istream& operator>>(istream& is, interval<Scarlar>& i)
{
    Scarlar inf, sup;
    is >> inf >> sup;
    i = interval<Scarlar>(inf, sup);
    return is;
}
template<typename Scarlar>
ostream& operator<<(ostream& os, interval<Scarlar> i)
{
    i.out(os);
    return os;
}
template<typename Scarlar>
class MatrixInterval
{
    // In this method, the intervals are stored as two matrices,
    // the inf matrix and the sup matrix.
    // Operations are done on these two matrices.
    Matrix<Scarlar, -1, -1> inf;
    Matrix<Scarlar, -1, -1> sup;
public:
    MatrixInterval() {}
    MatrixInterval(int i, int j = 1)
    {
        inf.resize(i, j);
        sup.resize(i, j);
    }
    MatrixInterval(Matrix<Scarlar, -1, -1> X)
    {
        inf = sup = X;
    }
    MatrixInterval(Matrix<Scarlar, -1, -1> X, Matrix<Scarlar, -1, -1> Y, bool assertion = true)
    {
        assert(X.rows() == Y.rows());
        assert(X.cols() == Y.cols());
        inf.resize(X.rows(), X.cols());
        sup.resize(X.rows(), X.cols());
        if (assertion)
        {
            for (Index i = 0; i < X.rows(); ++i)
                for (Index j = 0; j < X.cols(); ++j)
                {
                    if (X(i, j) < Y(i, j))
                    {
                        inf(i, j) = X(i, j);
                        sup(i, j) = Y(i, j);
                    }
                    else
                    {
                        inf(i, j) = Y(i, j);
                        sup(i, j) = X(i, j);
                    }
                }
        }
        else
        {
            inf = X;
            sup = Y;
        }
    }
    int rows()
    {
        return inf.rows();
    }
    int cols()
    {
        return inf.cols();
    }
    int size()
    {
        return rows() * cols();
    }
    MatrixInterval<Scarlar> row(int i)
    {
        assert((i >= 0) && (i < inf.rows()));
        return MatrixInterval<Scarlar>(inf.row(i), sup.row(i));
    }
    MatrixInterval<Scarlar> col(int j)
    {
        assert((j >= 0) && (j < inf.cols()));
        return MatrixInterval<Scarlar>(inf.col(j), sup.col(j));
    }
    void set(int i, int j, interval<Scarlar> I)
    {
        inf(i, j) = I.min();
        sup(i, j) = I.max();
    }
    interval<Scarlar> operator()(int k)
    {
        assert((k >= 0) && (k < size()));
        return interval<Scarlar>(inf(k / cols(), k % cols()), sup(k / cols(), k % cols()));
    }
    interval<Scarlar> operator()(int i, int j)
    {
        assert((i >= 0) && (i < inf.rows()));
        assert((j >= 0) && (j < inf.cols()));
        return interval<Scarlar>(inf(i, j), sup(i, j));
    }
    Matrix<Scarlar,-1,-1> min()
    {
        return inf;
    }
    Matrix<Scarlar,-1,-1> max()
    {
        return sup;
    }
    bool contains(Matrix<Scarlar, 1, -1> X)
    {
        if (rows() != X.rows())
            return false;
        if (cols() != X.cols())
            return false;
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
                if ((inf(i, j) > X(i, j)) || (sup(i, j) < X(i, j)))
                    return false;
        return true;
    }
    bool ncontains(Matrix<Scarlar, 1, -1> X)
    {
        return !contains(X);
    }
    bool is_include(MatrixInterval<Scarlar> M)
    {
        if (rows() != M.rows())
            return false;
        if (cols() != M.cols())
            return false;
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
                if ((inf(i, j) > M.min()(i, j)) || (sup(i, j) < M.max()(i, j)))
                    return false;
        return true;
    }
    bool is_included(MatrixInterval<Scarlar> M)
    {
        if (rows() != M.rows())
            return false;
        if (cols() != M.cols())
            return false;
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
                if ((inf(i, j) < M.min()(i, j)) || (sup(i, j) > M.max()(i, j)))
                    return false;
        return true;
    }
    bool intersects(MatrixInterval<Scarlar> M)
    {
        if (rows() != M.rows())
            return false;
        if (cols() != M.cols())
            return false;
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
                if ((inf(i, j) <= M.max()(i, j)) || (sup(i, j) >= M.min()(i, j)))
                    return true;
        return false;
    }
    MatrixInterval<Scarlar> block(int startRow, int startCol, int blockRow, int blockCol)
    {
        return MatrixInterval<Scarlar>(inf.block(startRow, startCol, blockRow, blockCol), sup.block(startRow, startCol, blockRow, blockCol), false);
    }
    MatrixInterval<Scarlar> coblock(int _i, int _j)
    {
        assert((_i >= 0) && (_i < rows()));
        assert((_j >= 0) && (_j < cols()));
        Matrix<Scarlar, -1, -1> coblock_inf(rows() - 1, cols() - 1), coblock_sup(rows() - 1, cols() - 1);
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
                if (i != _i && j != _j)
                {
                    coblock_inf(((i < _i) ? i : (i - 1)), ((j < _j) ? j : (j - 1))) = inf(i, j);
                    coblock_sup(((i < _i) ? i : (i - 1)), ((j < _j) ? j : (j - 1))) = sup(i, j);
                }
        return MatrixInterval<Scarlar>(coblock_inf, coblock_sup, false);
    }
    MatrixInterval<Scarlar> operator+(MatrixInterval<Scarlar> M)
    {
        assert(rows() == M.rows());
        assert(cols() == M.cols());
        Matrix<Scarlar, -1, -1> sum_inf(rows(), cols());
        Matrix<Scarlar, -1, -1> sum_sup(rows(),cols());
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
            {
                sum_inf(i, j) = inf(i, j) + M.min()(i, j);
                sum_sup(i, j) = sup(i, j) + M.max()(i, j);
            }
        return MatrixInterval<Scarlar>(sum_inf, sum_sup, false);
        //return MatrixInterval<Scarlar>(min() + M.min(), max + M.max(), false);    This is not appliable due to type error.
    }
    MatrixInterval<Scarlar> operator-(MatrixInterval<Scarlar> M)
    {
        assert(rows() == M.rows());
        assert(cols() == M.cols());
        Matrix<Scarlar, -1, -1> sum_inf(rows(), cols());
        Matrix<Scarlar, -1, -1> sum_sup(rows(), cols());
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
            {
                sum_inf(i, j) = inf(i, j) - M.max()(i, j);
                sum_sup(i, j) = sup(i, j) - M.min()(i, j);
            }
        return MatrixInterval<Scarlar>(sum_inf, sum_sup, false);
        //return MatrixInterval<Scarlar>(min() - M.max(), max - M.min(), false);    This is not appliable due to type error.
    }
    MatrixInterval<Scarlar> operator*(MatrixInterval<Scarlar> M)
    {
        assert(cols() == M.rows());
        Matrix<Scarlar, -1, -1> sum_inf(rows(), cols());
        Matrix<Scarlar, -1, -1> sum_sup(rows(), cols());
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < M.cols(); ++j)
            {
                interval<Scarlar> sum(0, 0);
                for (Index k = 0; k < cols(); ++k)
                    sum = sum + (operator()(i, k) * M(k, j));
                sum_inf(i, j) = sum.min();
                sum_sup(i, j) = sum.max();
            }
        return MatrixInterval<Scarlar>(sum_inf, sum_sup, false);
    }
    MatrixInterval<Scarlar> cwiseproduct(MatrixInterval<Scarlar> M)
    {
        assert(rows() == M.rows());
        assert(cols() == M.cols());
        Matrix<Scarlar, -1, -1> sum_inf(rows(), cols());
        Matrix<Scarlar, -1, -1> sum_sup(rows(), cols());
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
            {
                Scarlar extremes[4] = { inf(i,j) * M.min()(i,j),inf(i,j) * M.max()(i,j),sup(i,j) * M.min()(i,j),sup(i,j) * M.max()(i,j) };
                sort(extremes, extremes + 4);
                sum_inf(i, j) = extremes[0];
                sum_sup(i, j) = extremes[3];
            }
        return MatrixInterval<Scarlar>(sum_inf, sum_sup, false);
    }
    MatrixInterval<Scarlar> cwisedivide(MatrixInterval<Scarlar> M)
    {
        assert(rows() == M.rows());
        assert(cols() == M.cols());
        Matrix<Scarlar, -1, -1> sum_inf(rows(), cols());
        Matrix<Scarlar, -1, -1> sum_sup(rows(), cols());
        for (Index i = 0; i < rows(); ++i)
            for (Index j = 0; j < cols(); ++j)
            {
                assert(M.min()(i, j).ncontains(Scarlar(0)));
                Scarlar extremes[4] = { inf(i,j) / M.min()(i,j),inf(i,j) / M.max()(i,j),sup(i,j) / M.min()(i,j),sup(i,j) / M.max()(i,j) };
                sort(extremes, extremes + 4);
                sum_inf(i, j) = extremes[0];
                sum_sup(i, j) = extremes[3];
            }
        return MatrixInterval<Scarlar>(sum_inf, sum_sup, false);
    }
    MatrixInterval<Scarlar> operator+(Matrix<Scarlar, -1, -1> X)
    {
        return operator+(MatrixInterval<Scarlar>(X));
    }
    MatrixInterval<Scarlar> operator-(Matrix<Scarlar, -1, -1> X)
    {
        return operator-(MatrixInterval<Scarlar>(X));
    }
    MatrixInterval<Scarlar> operator*(Matrix<Scarlar, -1, -1> X)
    {
        return operator*(MatrixInterval<Scarlar>(X));
    }
    MatrixInterval<Scarlar> transpose()
    {
        return MatrixInterval<Scarlar>(inf.transpose(), sup.transpose(), false);
    }
    interval<Scarlar> trace()
    {
        interval<Scarlar> sum(0, 0);
        for (Index i = 0; i < std::min(rows(), cols()); ++i)
            sum = sum + operator()(i, i);
        return sum;
    }
    interval<Scarlar> determinant()
    {
        assert(rows() == cols());
        if (rows() == 0)
            return interval<Scarlar>(1, 1);
        interval<Scarlar> det(0, 0);
        for (Index i = 0; i < rows(); ++i)
            if (i % 2 == 0)
                det = det + operator()(i, 0) * coblock(i, 0).determinant();
            else
                det = det - operator()(i, 0) * coblock(i, 0).determinant();
        return det;
    }
    interval<Scarlar> dot(MatrixInterval<Scarlar> M)
    {
        if (cols() == 1)
            return transpose().dot(M);
        if (M.rows() == 1)
            return dot(M.transpose());
        assert(rows() == 1);
        assert(M.cols() == 1);
        return operator*(M)(0, 0);
    }
    MatrixInterval<Scarlar> cross(MatrixInterval<Scarlar> M)
    {
        if (cols() == 1)
            return transpose().cross(M);
        if (M.rows() == 1)
            return cross(M.transpose());
        assert((rows() == 1) && (cols() == 3));
        assert((M.rows() == 3) && (M.cols() == 1));
        interval<Scarlar> x, y, z;
        Matrix<Scarlar, 3, 1> cross_inf, cross_sup;
        x = operator()(1) * M(2) - operator()(2) * M(1);
        y = operator()(2) * M(0) - operator()(0) * M(2);
        z = operator()(0) * M(1) - operator()(1) * M(0);
        cross_inf(0) = x.min();
        cross_inf(1) = y.min();
        cross_inf(2) = z.min();
        cross_sup(0) = x.max();
        cross_sup(1) = y.max();
        cross_sup(2) = z.max();
        return MatrixInterval<Scarlar>(cross_inf, cross_sup, false);
    }
    void add_row()
    {
        inf.conservativeResize(inf.rows() + 1, inf.cols());
        sup.conservativeResize(sup.rows() + 1, sup.cols());
    }
    void add_col()
    {
        inf.conservativeResize(inf.rows(), inf.cols() + 1);
        sup.conservativeResize(sup.rows(), sup.cols() + 1);
    }
    void out(ostream& os = cout)
    {
        for (Index i = 0; i < rows() - 1; ++i)
        {
            for (Index j = 0; j < cols() - 1; ++j)
            {
                operator()(i, j).out(os);
                os << " ";
            }
            operator()(i, cols() - 1).out(os);
            os << endl;
        }
        for (Index j = 0; j < cols() - 1; ++j)
        {
            operator()(rows() - 1, j).out(os);
            os << " ";
        }
        operator()(rows() - 1, cols() - 1).out(os);
    }
};
template<typename Scarlar>
ostream& operator<<(ostream& os, MatrixInterval<Scarlar> m)
{
    m.out(os);
    return os;
}
#define MatrixIi MatrixInterval<int>
#define MatrixIf MatrixInterval<float>
#define MatrixId MatrixInterval<double>