// SymGroup.cpp : This file defines the data structure for symmetry group S_n
// Templated by n=size.

#include<iostream>
#include<iomanip>

// Get the t-th binary bit of n.
int getBin(int n, int t)
{
    return (n >> t) & 1;
}
// pow(2,int t)
int pow2(int t)
{
    return 1 << t;
}
// n*pow(2,int t)
int pow2n(int n, int t)
{
    return n << t;
}
// 1-n
int ne(int n)
{
    return 1 - n;
}

template<int size>
class SymG
{
    int map[size];
    bool is_id;
public:
    // Initialize by identity.
    SymG()
    {
        for (int i = 0; i < size; ++i)
            map[i] = i;
        is_id = true;
    }
    // Consecutive single loop
    SymG(int k, int n = 1)
    {
        if (k + n < size)
            for (int i = 0; i < size; ++i)
                if ((i < k) || (i >= k + n))
                    map[i] = i;
                else if (i == k + n - 1)
                    map[i] = k;
                else
                    map[i] = i + 1;
        else
            for (int i = 0; i < size; ++i)
                if ((i < k) && (i >= k + n - size))
                    map[i] = i;
                else if (i == k + n - 1)
                    map[i] = k;
                else if (i == size - 1)
                    map[i] = 0;
                else
                    map[i] = i + 1;
        is_id = (n <= 1);
    }
    // Direct assign.
    SymG(int* m)
    {
        memcpy(map, m, size * sizeof(m[0]));
        is_id = true;
        for (int i = 0; i < size; ++i)
            if (map[i] != i)
                is_id = false;
    }
    // is_id
    bool is_ID()
    {
        return is_id;
    }
    // Entry acts on integer.
    int act(int n)
    {
        return map[n % size];
    }
    // Entry acts on bindary components.
    int binAct(int n)
    {
        if (is_id)
            return n;
        int An = n;                                         // Unchanged parts.
        int Bn = 0;                                         // Changed parts.
        for (int i = 0; i < size; ++i)
            if (map[i] != i)
            {
                An = An & ~pow2(i);
                Bn = Bn | pow2n(getBin(n, i), map[i]);
            }
        return An | Bn;
    }
    // Inverse entry.
    SymG<size> inverse()
    {
        if (is_id)
            return id();
        int m[size];
        for (int i = 0; i < size; ++i)
            m[map[i]] = i;
        return SymG<size>(m);
    }
    // Multiply
    SymG<size> operator*(SymG<size> s)
    {
        if (is_id)
            return s;
        if (s.is_ID())
            return *this;
        int m[size];
        for (int i = 0; i < size; ++i)
            m[i] = map[s.act(i)];
        return SymG<size>(m);
    }
    // Divide
    SymG<size> operator/(SymG<size> s)
    {
        return operator*(s.inverse());
    }
    // Identity.
    static SymG<size> id()
    {
        return SymG<size>();
    }
    // output
    void out(std::ostream& os = std::cout)
    {
        for (int i = 0; i < size - 1; ++i)
            os << map[i] << " ";
        os << map[size - 1] << endl;
    }
};

template<int size>
std::ostream& operator<<(std::ostream& os, SymG<size> s)
{
    s.out(os);
    return os;
}
