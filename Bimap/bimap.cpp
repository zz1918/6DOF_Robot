// bimap.cpp : Hash table for (key1,key2)|->(content).
//

#ifndef BIMAP_H
#define BIMAP_H

#include<map>
using namespace std;

template<typename key1, typename key2, typename type>
class bimap {
    map<key1, map<key2, type>> bim;
public:
    bimap() {}
    bool find(key1 k, key2 l)
    {
        auto row_key = bim.find(k);
        if (row_key == bim.end())
            return false;
        auto col_key = row_key->second.find(l);
        if (col_key == row_key->second.end())
            return false;
        return true;
    }
    bool insert(key1 k, key2 l, type c)
    {
        if (find(k, l))
            return false;
        bim[k][l] = c;
        return true;
    }
    type coeff(key1 k, key2 l)
    {
        return bim[k][l];
    }
    void clear()
    {
        bim.clear();
    }
};


#endif