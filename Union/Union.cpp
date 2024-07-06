// Union.cpp
// This file gives a template of union find method.
// The type in the union requires a member "UnionNode* comp" to store the component that it is in the union.
// Universal settings: if comp=NULL, then it is undefined.

class UnionNode
{
#define ROOT this
public:
    UnionNode* parent;
    int size;
    int rank;
    UnionNode()
    {
        parent = ROOT;
        size = 1;
        rank = 0;
    }
    UnionNode(UnionNode* p)
    {
        parent = p;
        size = 1;
        rank = 0;
        p->size += 1;
        p->rank = max(p->rank, 1);
    }
    bool is_root()
    {
        return parent == ROOT;
    }
    UnionNode* default_find()
    {
        if (is_root())
            return this;
        else
            return parent->default_find();
    }
    UnionNode* root_find()
    {
        if (is_root())
            return this;
        else
        {
            UnionNode* root = parent->root_find();
            if (!parent->is_root())
                parent->size -= size;
            parent = root;
            return root;
        }
    }
    UnionNode* parent_find()
    {
        if (is_root())
            return this;
        else
        {
            if (parent->is_root())
                return parent;
            else
            {
                parent->size -= size;
                parent->parent->size += size;
                parent = parent->parent;
                return parent->parent_find();
            }
        }
    }
    void link(UnionNode* p)
    {
        if (!is_root())
            parent->link(p);
        else
        {
            parent = p;
            p->size += size;
            p->rank = max(p->rank, rank + 1);
        }
    }
#undef ROOT
};

// The type in the union requires a public member "UnionNode* comp" to store the component that it is in the union.
template<typename type>
class Union
{
    int link_strategy, find_strategy;
public:
    enum ls { SIZE, RANK };
    enum fs { ROOT, PARENT };
    Union(int s = RANK, int t = ROOT)
    {
        link_strategy = s;
        find_strategy = t;
    }
    // Register the type t into the union.
    void insert(type* t)
    {
        t->comp = new UnionNode();
    }
    // Find the equivalent class of type t in the union.
    UnionNode* find(type* t)
    {
        // if comp is NULL, then t is not in this union.
        if (t->comp == NULL)
            return NULL;
        switch (find_strategy)
        {
        case ROOT:return t->comp->root_find();
        case PARENT:return t->comp->parent_find();
        default:return t->comp->default_find();
        }
    }
    // Link the equivalent class of type p with type q.
    void link(type* p, type* q)
    {
        UnionNode* fp = find(p), * fq = find(q);
        if (fp == NULL || fq == NULL || fp == fq)
            return;
        switch (link_strategy)
        {
        case SIZE: size_link(fp, fq); break;
        case RANK: rank_link(fp, fq); break;
        default: default_link(fp, fq);
        }
    }
    // Different linking strategies.
    void default_link(UnionNode* p, UnionNode* q)
    {
        p->link(q);
    }
    void size_link(UnionNode* p, UnionNode* q)
    {
        if (p->size > q->size)
            p->link(q);
        else
            q->link(p);
    }
    void rank_link(UnionNode* p, UnionNode* q)
    {
        if (p->rank > q->rank)
            p->link(q);
        else
            q->link(p);
    }
};