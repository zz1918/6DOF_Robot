// BoxTree.cpp : BoxTree maintains the Box structures in R^d, where d is dim.
//

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <SymGroup.h>
#include <interval.h>
using namespace std;
using namespace Eigen;

// Class for dim-dimensional box and its subdivision.
template<int dim>
class Box
{
#define subsize 1<<dim
#define BDd Box<dim>
    MatrixId range;
    BDd* neighbor[2 * dim];                                       // Principle neighbor.
    BDd* children[subsize];                                       // Children.
    BDd* parent;                                                  // Parent.
    int mark;                                                     // Mark for something needed, i.e. the SO3 box.
    bool flip[2 * dim][dim];                                      // Flip array.
    SymG<dim> shift[2 * dim];                                     // Shift array.
    bool leaf;
    int depth;
    void ini_array()
    {
        for (int i = 0; i < 2 * dim; ++i)
            neighbor[i] = NULL;
        for (int i = 0; i < subsize; ++i)
            children[i] = NULL;
        memset(flip, 0, sizeof(flip));
        for (int i = 0; i < dim; ++i)
        {
            flip[i][i] = true;
            flip[i + dim][i] = true;
        }
        for (int i = 0; i < 2 * dim; ++i)
            shift[i] = SymG<dim>::id();
    }
public:
    // Box(range, mark, *parent, depth)
    Box(MatrixId r = MatrixId(dim), int m = -1, BDd* p = NULL, int d = 0)
    {
        range = r;
        parent = p;
        leaf = true;
        depth = d;
        mark = m;
        ini_array();
    }
    // range
    MatrixId Range()
    {
        return range;
    }
    // mark
    int Mark()
    {
        return mark;
    }
    // children[i]
    BDd* child(int i)
    {
        return children[i];
    }
    // *parent
    BDd* the_parent()
    {
        return parent;
    }
    // neighbor[(i + dim) % dim]
    BDd* pos_neighbor(int i)
    {
        return neighbor[(i + dim) % dim];                         // 0,1,2,3,...,-0,-1,-2,-3,...
    }
    // neighbor[((i + dim) % dim) + dim]
    BDd* neg_neighbor(int i)
    {
        return neighbor[((i + dim) % dim) + dim];                 // 0,1,2,3,...,-0,-1,-2,-3,...
    }
    // neighbor[(i + 2 * dim) % (2 * dim)]
    BDd* the_neighbor(int i)
    {
        return neighbor[(i + 2 * dim) % (2 * dim)];
    }
    // leaf
    bool is_leaf()
    {
        return leaf;
    }
    // root
    bool is_root()
    {
        return parent == NULL;
    }
    // If neighbor(i) is sibling
    bool is_sibling(int i)
    {
        if (the_neighbor(i) == NULL)
            return false;
        return the_neighbor(i)->the_parent() == parent;
    }
    // neighbor[(i + dim) % dim] = b
    void set_pos_neighbor(int i, BDd* b)
    {
        neighbor[(i + dim) % dim] = b;
    }
    // neighbor[((i + dim) % dim) + dim] = 
    void set_neg_neighbor(int i, BDd* b)
    {
        neighbor[((i + dim) % dim) + dim] = b;
    }
    // neighbor[(i + 2* dim) % (2 * dim)] = b
    void set_neighbor(int i, BDd* b)
    {
        if (i < dim)
            set_pos_neighbor(i, b);
        else
            set_neg_neighbor(i - dim, b);
    }
    // memcpy(flip[i], f, sizeof(f));
    void set_flip(int i, bool f[dim])
    {
        memcpy(flip[i], f, dim * sizeof(f[0]));
    }
    // memcpy(flip, f, sizeof(f))
    void set_flip(bool f[2 * dim][dim])
    {
        memcpy(flip, f, 2 * dim * dim * sizeof(f[0][0]));
    }
    // shift[i]=s
    void set_shift(int i, SymG<dim> s = SymG<dim>::id())
    {
        shift[i] = s;
    }
    // memcpy(shift, s, sizeof(s))
    void set_shift(SymG<dim> s[2 * dim])
    {
        memcpy(shift, s, 2 * dim * sizeof(s[0]));
    }
    // Assign "this" to as i-th neighbor of b and to all its i-pos children.
    void assign_neighbor(int i, BDd* b)
    {
        if (i > 2 * dim)
            i = i % (2 * dim);
        b->set_neighbor(i, this);
        if (b->is_leaf())
            return;
        for (int k = 0; k < subsize; ++k)
            if (getBin(k, i % dim) && (i < dim))
                assign_neighbor(i, b->child(k));
            else if (ne(getBin(k, i % dim)) && (i >= dim))
                assign_neighbor(i, b->child(k));
    }
    // Subdivision
    void subdivide()
    {
        // Step 1: Create subdivisions.
        // All boxes are sequenced by the lex-inequality of the range.min() coordinates.
        for (int i = 0; i < subsize; ++i)
        {
            Matrix<double, dim, 1> new_inf, new_sup;
            for (int j = 0; j < dim; ++j)
                if (getBin(i,j))
                {
                    new_inf(j) = (range(j).min() + range(j).max()) / 2;
                    new_sup(j) = range(j).max();
                }
                else
                {
                    new_inf(j) = range(j).min();
                    new_sup(j) = ((range(j).min() + range(j).max()) / 2);
                }
            MatrixId new_range(new_inf, new_sup, false);
            children[i] = new BDd(new_range, mark, this, depth + 1);
            bool new_flip[2 * dim][dim];
            memcpy(new_flip, flip, 2 * dim * dim * sizeof(flip[0][0]));
            for (int j = 0; j < dim; ++j)
                if (getBin(i, j))
                    new_flip[j + dim][j] = true;
                else
                    new_flip[j][j] = true;
            children[i]->set_flip(new_flip);
            for(int j=0;j<dim;++j)
                if (getBin(i, j))
                {
                    children[i]->set_shift(j, shift[j]);
                    children[i]->set_shift(j + dim);
                }
                else
                {
                    children[i]->set_shift(j);
                    children[i]->set_shift(j + dim, shift[j + dim]);
                }
        }
        // Step 2: Assign neighbors for children.
        // Step 3: Re-assign neighbors for neighbors and their children.
        for (int i = 0; i < subsize; ++i)
            for (int j = 0; j < dim; ++j)
            {
                if (getBin(i, j))
                {
                    child(i)->set_neg_neighbor(j, children[i - pow2(j)]);
                    if (pos_neighbor(j) == NULL || pos_neighbor(j)->is_leaf())
                        child(i)->set_pos_neighbor(j, pos_neighbor(j));
                    else
                    {
                        int ti = 0, tj = j + dim;
                        for (int k = 0; k < dim; ++k)
                            ti = ti | ((flip[j][k]) ? pow2n(ne(getBin(i, k)), k) : pow2n(getBin(i, k), k));
                        ti = shift[j].binAct(ti);
                        tj = shift[j].act(j) + (flip[j][j] ? dim : 0);
                        //cout << "Assigning child " << i << ", in direction positive " << j << " to the child " << ti << " of neighbor by direction " << tj << endl;
                        child(i)->set_pos_neighbor(j, pos_neighbor(j)->child(ti));                                      // For step 2.
                        child(i)->assign_neighbor(tj, pos_neighbor(j)->child(ti));                                      // For step 3.
                    }
                }
                else
                {
                    child(i)->set_pos_neighbor(j, children[i + pow2(j)]);
                    if (neg_neighbor(j) == NULL || neg_neighbor(j)->is_leaf())
                        child(i)->set_neg_neighbor(j, neg_neighbor(j));
                    else
                    {
                        int ti = 0, tj = j;
                        for (int k = 0; k < dim; ++k)
                            ti = ti | ((flip[j + dim][k]) ? pow2n(ne(getBin(i, k)), k) : pow2n(getBin(i, k), k));
                        ti = shift[j + dim].binAct(ti);
                        tj = shift[j + dim].act(j) + (flip[j + dim][j] ? 0 : dim);
                        //cout << "Assigning child " << i << ", in direction negative " << j << " to the child " << ti << " of neighbor by direction " << tj << endl;
                        child(i)->set_neg_neighbor(j, neg_neighbor(j)->child(ti));                                      // For step 2.
                        child(i)->assign_neighbor(tj, neg_neighbor(j)->child(ti));                                      // For step 3.
                    }
                }
            }
        // Step 4: This is not a leaf anymore.
        leaf = false;
    }
    // cout *this.
    void out(ostream& os = cout, int l = 0)
    {
        for (int i = 0; i < l; ++i)
            os << "      ";
        if (dim == 3)
        {
            Vector3d m = range.min();
            Vector3d M = range.max();
            switch (mark)
            {
            case 0:os << MatrixId(Vector4d(1, m(0), m(1), m(2)), Vector4d(1, M(0), M(1), M(2)), false).transpose() << endl; break;
            case 1:os << MatrixId(Vector4d(m(0), 1, m(1), m(2)), Vector4d(M(0), 1, M(1), M(2)), false).transpose() << endl; break;
            case 2:os << MatrixId(Vector4d(m(0), m(1), 1, m(2)), Vector4d(M(0), M(1), 1, M(2)), false).transpose() << endl; break;
            case 3:os << MatrixId(Vector4d(m(0), m(1), m(2), 1), Vector4d(M(0), M(1), M(2), 1), false).transpose() << endl; break;
            default:os << range.transpose() << endl;
            }
        }
        else
            os << range.transpose() << endl;
        if (leaf)
            return;
        for (int i = 0; i < subsize; ++i)
            children[i]->out(os, l + 1);
    }
    // cout *neighbors.
    void show_neighbor(ostream& os = cout, int l = 0)
    {
        for (int i = 0; i < 2 * dim; ++i)
        {
            int k = (i + dim) % (2 * dim);
            if (the_neighbor(k) != NULL)
            {
                for (int j = 0; j < l; ++j)
                    os << "      ";
                os << "neighbor " << (k < dim ? "pos" : "neg") << " " << (k < dim ? k : k - dim) << ":" << endl;
                the_neighbor(k)->out(os, l);
            }
        }
    }
#undef subsize
#undef BDd
};

template<int dim>
ostream& operator<<(ostream& os, Box<dim> B)
{
    B.out(os);
    return os;
}

class R3Tree
{
    Box<3>* root;
public:
    Box<3>* pointer;
    R3Tree(Vector3d a = Vector3d(-1, -1, -1), Vector3d b = Vector3d(1, 1, 1))
    {
        root = new Box<3>(MatrixId(a, b));
        pointer = root;
    }
};

class SO3Tree
{
    Box<3>* cells[4];
public:
    Box<3>* pointer;
    SO3Tree()
    {
        bool flip_types[4][3] = { {false,false,false},{false,true,true},{true,false,true},{true,true,false} };
        for (int i = 0; i < 4; ++i)
            cells[i] = new Box<3>(MatrixId(Vector3d(-1, -1, -1), Vector3d(1, 1, 1)), i);
        pointer = cells[0];

        // Cell[0] (C_w)
        // w 1-3-5-7 -> x 1-3-5-7
        cells[0]->set_neighbor(0, cells[1]);
        cells[0]->set_flip(0, flip_types[0]);
        cells[0]->set_shift(0, SymG<3>::id());
        // w 2-3-6-7 -> y 1-3-5-7
        cells[0]->set_neighbor(1, cells[2]);
        cells[0]->set_flip(1, flip_types[0]);
        cells[0]->set_shift(1, SymG<3>(0, 2));
        // w 4-5-6-7 -> z 1-3-5-7
        cells[0]->set_neighbor(2, cells[3]);
        cells[0]->set_flip(2, flip_types[0]);
        cells[0]->set_shift(2, SymG<3>(0, 3));
        // w 0-2-4-6 -> x 6-4-2-0
        cells[0]->set_neighbor(3, cells[1]);
        cells[0]->set_flip(3, flip_types[1]);
        cells[0]->set_shift(3, SymG<3>::id());
        // w 0-1-4-5 -> y 6-4-2-0
        cells[0]->set_neighbor(4, cells[2]);
        cells[0]->set_flip(4, flip_types[2]);
        cells[0]->set_shift(4, SymG<3>(0, 2));
        // w 0-1-2-3 -> z 6-4-2-0
        cells[0]->set_neighbor(5, cells[3]);
        cells[0]->set_flip(5, flip_types[3]);
        cells[0]->set_shift(5, SymG<3>(0, 3));

        // Cell[1] (C_x)
        // x 1-3-5-7 -> w 1-3-5-7
        cells[1]->set_neighbor(0, cells[0]);
        cells[1]->set_flip(0, flip_types[0]);
        cells[1]->set_shift(0, SymG<3>::id());
        // x 2-3-6-7 -> y 2-3-6-7
        cells[1]->set_neighbor(1, cells[2]);
        cells[1]->set_flip(1, flip_types[0]);
        cells[1]->set_shift(1, SymG<3>::id());
        // x 4-5-6-7 -> z 2-3-6-7
        cells[1]->set_neighbor(2, cells[3]);
        cells[1]->set_flip(2, flip_types[0]);
        cells[1]->set_shift(2, SymG<3>(1, 2));
        // x 0-2-4-6 -> w 6-4-2-0
        cells[1]->set_neighbor(3, cells[0]);
        cells[1]->set_flip(3, flip_types[1]);
        cells[1]->set_shift(3, SymG<3>::id());
        // x 0-1-4-5 -> y 5-4-1-0
        cells[1]->set_neighbor(4, cells[2]);
        cells[1]->set_flip(4, flip_types[2]);
        cells[1]->set_shift(4, SymG<3>::id());
        // x 0-1-2-3 -> z 5-4-1-0
        cells[1]->set_neighbor(5, cells[3]);
        cells[1]->set_flip(5, flip_types[3]);
        cells[1]->set_shift(5, SymG<3>(1, 2));

        // Cell[2] (C_y)
        // y 1-3-5-7 -> w 2-3-6-7
        cells[2]->set_neighbor(0, cells[0]);
        cells[2]->set_flip(0, flip_types[0]);
        cells[2]->set_shift(0, SymG<3>(0, 2).inverse());
        // y 2-3-6-7 -> x 2-3-6-7
        cells[2]->set_neighbor(1, cells[1]);
        cells[2]->set_flip(1, flip_types[0]);
        cells[2]->set_shift(1, SymG<3>::id());
        // y 4-5-6-7 -> z 4-5-6-7
        cells[2]->set_neighbor(2, cells[3]);
        cells[2]->set_flip(2, flip_types[0]);
        cells[2]->set_shift(2, SymG<3>::id());
        // y 0-2-4-6 -> w 5-4-1-0
        cells[2]->set_neighbor(3, cells[0]);
        cells[2]->set_flip(3, flip_types[1]);
        cells[2]->set_shift(3, SymG<3>(0, 2).inverse());
        // y 0-1-4-5 -> x 5-4-1-0
        cells[2]->set_neighbor(4, cells[1]);
        cells[2]->set_flip(4, flip_types[2]);
        cells[2]->set_shift(4, SymG<3>::id());
        // y 0-1-2-3 -> z 3-2-1-0
        cells[2]->set_neighbor(5, cells[3]);
        cells[2]->set_flip(5, flip_types[3]);
        cells[2]->set_shift(5, SymG<3>::id());

        // Cell[3] (C_z)
        // z 1-3-5-7 -> w 4-5-6-7
        cells[3]->set_neighbor(0, cells[0]);
        cells[3]->set_flip(0, flip_types[0]);
        cells[3]->set_shift(0, SymG<3>(0, 3).inverse());
        // z 2-3-6-7 -> x 4-5-6-7
        cells[3]->set_neighbor(1, cells[1]);
        cells[3]->set_flip(1, flip_types[0]);
        cells[3]->set_shift(1, SymG<3>(1, 2).inverse());
        // z 4-5-6-7 -> y 4-5-6-7
        cells[3]->set_neighbor(2, cells[2]);
        cells[3]->set_flip(2, flip_types[0]);
        cells[3]->set_shift(2, SymG<3>::id());
        // z 0-2-4-6 -> w 3-2-1-0
        cells[3]->set_neighbor(3, cells[0]);
        cells[3]->set_flip(3, flip_types[1]);
        cells[3]->set_shift(3, SymG<3>(0, 3).inverse());
        // z 0-1-4-5 -> x 3-2-1-0
        cells[3]->set_neighbor(4, cells[1]);
        cells[3]->set_flip(4, flip_types[2]);
        cells[3]->set_shift(4, SymG<3>(1, 2).inverse());
        // z 0-1-2-3 -> y 3-2-1-0
        cells[3]->set_neighbor(5, cells[2]);
        cells[3]->set_flip(5, flip_types[3]);
        cells[3]->set_shift(5, SymG<3>::id());
    }
    Box<3>* cell(int i)
    {
        return cells[i % 4];
    }
    static MatrixId SO3range(Box<3>* b)
    {
        Vector3d m = b->Range().min();
        Vector3d M = b->Range().max();
        switch (b->Mark())
        {
        case 0:return MatrixId(Vector4d(1, m(0), m(1), m(2)), Vector4d(1, M(0), M(1), M(2)), false);
        case 1:return MatrixId(Vector4d(m(0), 1, m(1), m(2)), Vector4d(M(0), 1, M(1), M(2)), false);
        case 2:return MatrixId(Vector4d(m(0), m(1), 1, m(2)), Vector4d(M(0), M(1), 1, M(2)), false);
        case 3:return MatrixId(Vector4d(m(0), m(1), m(2), 1), Vector4d(M(0), M(1), M(2), 1), false);
        default: return b->Range();
        }
    }
};

class SE3Tree
{
    R3Tree* p;
    SO3Tree* q;
};