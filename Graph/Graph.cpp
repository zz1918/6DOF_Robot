// Graph.cpp : This file implements the graph structure with union-find algorithm and find path algorithms.
// connected(u,v) will use the union-find algorithm to check if u and v are in the same union.
// path(u,v) will output a path connecting u and v (empty if no path).

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <bimap.h>
using namespace std;

template<typename type>
class GraphNode
{
#define ROOT this
public:
    // Content of the node.
    type* content;
    // Id of the node in the graph.
    int id;
    // The parent of the node.
    GraphNode* parent;
    // Rank of the node.
    int rank;
    // Nodes connecting to this node (successors).
    vector<GraphNode<type>*> successor;
    // Constructed by content without parent.
    GraphNode(type* t, int _id = -1)
    {
        content = t;
        parent = ROOT;
        id = _id;
        rank = 0;
    }
    // Constructed by content with parent.
    GraphNode(type* t, GraphNode<type>* p, int _id = -1)
    {
        content = t;
        parent = p;
        rank = 0;
        id = _id;
        p->rank = max(p->rank, 1);
    }

    //************************* Graph Structure ************************//

    // If the node is in the graph.
    bool inserted()
    {
        return id >= 0;
    }
    // If the node is not in the graph.
    bool uninserted()
    {
        return id < 0;
    }
    // Make a node as successor.
    void make_edge(GraphNode<type>* v)
    {
        successor.push_back(v);
    }
    // Return the i-th successor.
    GraphNode<type>* succ(int i)
    {
        return successor[i];
    }

    //*************************** Union Find ***************************//
         
    // If the node is the root in the union-find structure.
    bool is_root()
    {
        return parent == ROOT;
    }
    // Find the root of the node in the union-find structure.
    GraphNode<type>* find()
    {
        if (is_root())
            return this;
        else
        {
            GraphNode<type>* root = parent->find();
            parent = root;
            return root;
        }
    }
    // Link this node to node p.
    void link_to(GraphNode<type>* v)
    {
        if (!is_root())
            parent->link_to(v);
        else
        {
            parent = v;
            v->rank = max(v->rank, rank + 1);
        }
    }
#undef ROOT
};

template<typename type>
class GraphEdge
{
public:
    // End nodes of the edge.
    GraphNode<type>* first, * second;
    // Id of the edge in the graph.
    int id;
    // Constructed by two nodes.
    GraphEdge(GraphNode<type>* u, GraphNode<type>* v, int _id = -1)
    {
        // Set end points and id.
        first = u;
        second = v;
        id = _id;

        // Make u and v as successors.
        u->make_edge(v);
        v->make_edge(u);

        // Maintain union-find structure.
        GraphNode<type>* fu = u->find(), * fv = v->find();
        if (fu != NULL && fv != NULL && fu != fv)
        {
            if (fu->rank > fv->rank)
                fv->link_to(fu);
            else
                fu->link_to(fv);
        }
    }

    //************************* Graph Structure ************************//

    // If the node is in the graph.
    bool inserted()
    {
        return id >= 0;
    }
    // If the node is not in the graph.
    bool uninserted()
    {
        return id < 0;
    }
};

template<typename type>
class Graph
{
    // Node list.
    vector<GraphNode<type>*> V;
    // Edge list.
    vector<GraphEdge<type>*> E;
    // End-points table.
    bimap<int, int, int> VE_table;
    // Visited nodes for DFS/BFS visit.
    set<int> mark;
    // Pi-table for DFS/BFS visit.
    map<int, int> pi;
    // If DFS/BFS is finished.
    bool found;
public:
    Graph()
    {
        found = false;
    }
    // Get the i-th node.
    GraphNode<type>* node(int i)
    {
        if (i >= V.size())
            return NULL;
        return V[i];
    }
    // Get the i-th edge.
    GraphEdge<type>* edge(int i)
    {
        if (i >= E.size())
            return NULL;
        return E[i];
    }

    //******************* Build graph operators *********************//

    // Insert a type into the graph.
    GraphNode<type>* insert(type* v)
    {
        GraphNode<type>* new_node = new GraphNode<type>(v, V.size());
        V.push_back(new_node);
        return new_node;
    }
    // Insert a node into the graph.
    GraphNode<type>* insert(GraphNode<type>* v)
    {
        v->id = V.size();
        V.push_back(v);
        return v;
    }
    // Insert an edge into the graph.
    GraphEdge<type>* insert(GraphEdge<type>* e)
    {
        if (e->first->uninserted())
            insert(e->first);
        if (e->second->uninserted())
            insert(e->second);
        E.push_back(e);
        VE_table.insert(e->first->id, e->second->id, e->id);
        VE_table.insert(e->second->id, e->first->id, e->id);
        return e;
    }
    // Find the equivalent class of type t in the union.
    GraphNode<type>* find(GraphNode<type>* v)
    {
        if (v == NULL)
            return NULL;
        else
            return v->find();
    }
    // Make the edge constructed by u and v.
    GraphEdge<type>* link(GraphNode<type>* u, GraphNode<type>* v)
    {
        if (VE_table.find(u->id, v->id))
            return E[VE_table.coeff(u->id, v->id)];
        GraphEdge<type>* new_edge = new GraphEdge<type>(u, v, E.size());
        return insert(new_edge);
    }
    // Make the edge constructed by V[i] and V[j].
    GraphEdge<type>* link(int i, int j)
    {
        if (i < 0 || j < 0)
            return NULL;
        return link(V[i], V[j]);
    }

    //******************** Find path algorithm **********************//

    // If node u and node v are in the same union.
    bool connected(GraphNode<type>* u, GraphNode<type>* v)
    {
        return u->find() == v->find() && u->find() != NULL;
    }
    // DFS visit from u to v.
    void DFS(GraphNode<type>* u, GraphNode<type>* v)
    {
        if (u == v)
            found = true;
        if (found)
            return;
        for (int i = 0; i < u->successor.size(); ++i)
            if (mark.find(u->succ(i)->id) == mark.end())
            {
                mark.insert(u->succ(i)->id);
                pi.insert(make_pair(u->succ(i)->id, u->id));
                DFS(u->succ(i), v);
            }
    }
    // BFS visit from u to v.
    void BFS(GraphNode<type>* u, GraphNode<type>* v)
    {
        queue<GraphNode<type>*> NEXT;
        NEXT.push(u);
        while (!NEXT.empty())
        {
            GraphNode<type>* next = NEXT.front();
            NEXT.pop();
            for (int i = 0; i < next->successor.size(); ++i)
                if (!mark.find(next->succ(i)->id))
                {
                    mark.insert(next->succ(i)->id);
                    pi.insert(make_pair(next->succ(i)->id, next->id));
                    NEXT.push(next->succ(i));
                }
        }
    }
    // Find path from node u to node v, return by sequence of nodes.
    list<GraphNode<type>*> path(GraphNode<type>* u, GraphNode<type>* v)
    {
        // Resulted new path.
        list<GraphNode<type>*> new_path;

        // If not connected, return empty path.
        if (!connected(u, v))
            return new_path;

        // Initialization for DFS/BFS.
        mark.clear();
        pi.clear();
        found = false;

        // DFS/BFS
        DFS(u, v);

        // Build the path from the pi-table.
        for (GraphNode<type>* w = v; w != u; w = V[pi[w->id]])
            new_path.push_front(w);
        new_path.push_front(u);
        return new_path;
    }
    // Show edges.
    void show_edge()
    {
        for (int i = 0; i < E.size(); ++i)
            cout << edge(i)->first->id << " " << edge(i)->second->id << endl;
    }
    // Show the pi array.
    void show_pi()
    {
        for (map<int, int>::iterator it = pi.begin(); it != pi.end(); ++it)
            cout << it->first << " " << it->second << endl;
    }
};

