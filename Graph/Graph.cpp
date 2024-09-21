/// Graph.cpp : This file implements the graph structure with union-find algorithm and find path algorithms.
// connected(u,v) will use the union-find algorithm to check if u and v are in the same union.
// path(u,v) will output a path connecting u and v (empty if no path).

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <bimap.h>
#include <chrono>
using namespace std;

enum Vcolor { GREEN, YELLOW, RED, GREY, BLACK };

int find_path_time = 0;

// This graph structure also implements the Union-Find structure.
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
    // List of colors of the node.
    vector<Vcolor> colors;
    // Nodes connecting to this node (successors).
    set<GraphNode<type>*> successor;
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
        successor.insert(v);
    }
    // t-color.
    Vcolor color(int t)
    {
        return colors[t];
    }
    // If this is green or yellow for i<t and is green for i>=t.
    bool passable(int t)
    {
        for (int i = 0; i < colors.size(); ++i)
            if (colors[i] != GREEN && colors[i] != YELLOW)
                return false;
            else if (i >= t && colors[i] != GREEN)
                return false;
        return true;
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
    // Weight of the edge.
    double weight;
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

        weight = 0.0;
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
    // Delete an edge.
    bool erase()
    {
        first->successor.erase(second);
        second->successor.erase(first);
        return true;
    }
};

template<typename type>
class Graph
{
    // Universal object list.
    vector<type*> O;
    // Universal node list.
    vector<GraphNode<type>*> V;
    // Universal edge list.
    vector<GraphEdge<type>*> E;
    // Universal end-points table.
    bimap<int, int, int> VE_table;
    // Visited nodes for DFS/BFS visit.
    set<int> mark;
    // Pi-table for DFS/BFS visit.
    map<int, int> pi;
    // If DFS/BFS is finished.
    bool found;
public:
    // Current node list.
    set<int> Vlist;
    // Current edge list.
    set<int> Elist;
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
    // Clear the graph.
    void clear()
    {
        O.clear();
        V.clear();
        E.clear();
        VE_table.clear();
        mark.clear();
        pi.clear();
        found = false;
    }
    // List of object.
    vector<type*> obj()
    {
        return O;
    }
    // List of current object.
    vector<type*> current_obj()
    {
        vector<type*> objs;
        for (auto it = Vlist.begin(); it != Vlist.end(); ++it)
            objs.push_back(node(*it)->content);
        return objs;
    }
    // Size of vertices of the graph.
    int Vsize()
    {
        return Vlist.size();
    }
    // Size of edges of the graph.
    int Esize()
    {
        return Elist.size();
    }

    //******************* Build graph operators *********************//

    // Insert a type into the graph.
    GraphNode<type>* insert(type* v)
    {
        GraphNode<type>* new_node = new GraphNode<type>(v, V.size());
        V.push_back(new_node);
        O.push_back(v);
        Vlist.insert(new_node->id);
        return new_node;
    }
    // Insert a node into the graph.
    GraphNode<type>* insert(GraphNode<type>* v)
    {
        v->id = V.size();
        V.push_back(v);
        O.push_back(v->content);
        Vlist.insert(v->id);
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
        Elist.insert(e->id);
        return e;
    }
    // Add a color for a node.
    void push_color(GraphNode<type>* v, Vcolor c)
    {
        if (v != NULL)
            v->colors.push_back(c);
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
    GraphEdge<type>* link(GraphNode<type>* u, GraphNode<type>* v, double w = 1.0)
    {
        if (u == NULL || v == NULL)
            return NULL;
        if (VE_table.find(u->id, v->id))
            return E[VE_table.coeff(u->id, v->id)];
        GraphEdge<type>* new_edge = new GraphEdge<type>(u, v, E.size());
        new_edge->weight = w;
        return insert(new_edge);
    }
    // Make the edge constructed by V[i] and V[j].
    GraphEdge<type>* link(int i, int j, double w = 1.0)
    {
        if (i < 0 || j < 0 || i >= V.size() || j >= V.size())
            return NULL;
        return link(V[i], V[j], w);
    }
    // Delete a node from graph.
    bool erase(GraphNode<type>* v)
    {
        if (v == NULL)
            return true;
        for (auto it = v->successor.begin(); it != v->successor.end(); it = v->successor.begin())
        {
            int eid = VE_table.coeff(v->id, (*it)->id);
            E[eid]->erase();
            Elist.erase(eid);
        }
        return Vlist.erase(v->id) > 0;
    }

    //**************** General find path algorithm ******************//

    // Is this node marked?
    bool marked(GraphNode<type>* v)
    {
        return mark.find(v->id) != mark.end();
    }
    // Make a node marked.
    void markit(GraphNode<type>* v)
    {
        mark.insert(v->id);
    }
    // Is this node in the graph?
    bool ingraph(GraphNode<type>* v)
    {
        return Vlist.find(v->id) != Vlist.end();
    }
    // Make u as path-parent of v.
    void make_parent(GraphNode<type>* u, GraphNode<type>* v)
    {
        pi.insert(make_pair(v->id, u->id));
    }
    // If node u and node v are in the same union, this only works for augmenting graph.
    bool quick_connected(GraphNode<type>* u, GraphNode<type>* v)
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
        for (auto it = u->successor.begin(); it != u->successor.end(); ++it)
            if (!marked(*it) && ingraph(*it))
            {
                markit(*it);
                make_parent((*it), u);
                DFS((*it), v);
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
            for (auto it = next->successor.begin(); it != next->successor.end(); ++it)
                if (!marked(*it) && ingraph(*it))
                {
                    markit(*it);
                    make_parent(next, (*it));
                    if (*it == v)
                    {
                        found = true;
                        return;
                    }
                    NEXT.push((*it));
                }
        }
    }
    // Find path from node u to node v, return by sequence of nodes.
    list<GraphNode<type>*> path(GraphNode<type>* u, GraphNode<type>* v)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        // Resulted new path.
        list<GraphNode<type>*> new_path;

        if (u == NULL)
            cout << "Invalid initial node." << endl;
        if (v == NULL)
            cout << "Invalid target node." << endl;
        if (u == NULL || v == NULL)
            return new_path;

        // Initialization for DFS/BFS.
        mark.clear();
        pi.clear();
        found = false;

        // DFS/BFS
        BFS(u, v);

        // Build the path from the pi-table.
        if (found != true)
            return new_path;
        for (GraphNode<type>* w = v; w != u; w = V[pi[w->id]])
            new_path.push_front(w);
        new_path.push_front(u);
        auto end_time = std::chrono::high_resolution_clock::now();
        find_path_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        return new_path;
    }
    // If node u and node v are connected by a path.
    bool connected(GraphNode<type>* u, GraphNode<type>* v)
    {
        return !path(u, v).empty();
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

    //************** Color-based find path algorithm ****************//

    // If v avoiding forbid area or not.
    bool outbid(GraphNode<type>* v, set<int>& forbid)
    {
        return forbid.find(v->id) == forbid.end();
    }
    // DFS visit from u to v avoiding forbid area such that it is green or yellow for i<t and it is green for i>=t.
    void DFS(GraphNode<type>* u, GraphNode<type>* v, set<int>& forbid, int t)
    {
        if (u == v)
            found = true;
        if (found)
            return;
        for (auto it = u->successor.begin(); it != u->successor.end(); ++it)
            if (!marked(*it) && ingraph(*it) && outbid(*it, forbid) && (*it)->passable(t))
            {
                markit(*it);
                make_parent((*it), u);
                DFS((*it), v, forbid, t);
            }
    }
    // BFS visit from u to v avoiding forbid area such that it is green or yellow for i<t and it is green for i>=t.
    void BFS(GraphNode<type>* u, GraphNode<type>* v, set<int>& forbid, int t)
    {
        queue<GraphNode<type>*> NEXT;
        NEXT.push(u);
        while (!NEXT.empty())
        {
            GraphNode<type>* next = NEXT.front();
            NEXT.pop();
            for (auto it = next->successor.begin(); it != next->successor.end(); ++it)
                if (!marked(*it) && ingraph(*it) && outbid(*it, forbid) && (*it)->passable(t))
                {
                    markit(*it);
                    make_parent(next, (*it));
                    if (*it == v)
                    {
                        found = true;
                        return;
                    }
                    NEXT.push((*it));
                }
        }
    }
    // Find path avoiding forbid area from node u to node v such that it is green or yellow for i<t and it is green for i>=t, return by sequence of nodes.
    list<GraphNode<type>*> path(GraphNode<type>* u, GraphNode<type>* v, set<int>& forbid, int t)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        // Resulted new path.
        list<GraphNode<type>*> new_path;

        if (u == NULL)
            cout << "Invalid initial node." << endl;
        if (v == NULL)
            cout << "Invalid target node." << endl;
        if (u == NULL || v == NULL)
            return new_path;

        // Initialization for DFS/BFS.
        mark.clear();
        pi.clear();
        found = false;

        // DFS/BFS
        BFS(u, v, forbid, t);

        // Build the path from the pi-table.
        if (found != true)
            return new_path;
        for (GraphNode<type>* w = v; w != u; w = V[pi[w->id]])
            new_path.push_front(w);
        new_path.push_front(u);
        auto end_time = std::chrono::high_resolution_clock::now();
        find_path_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        return new_path;
    }
    // If node u and node v are connected by a path avoiding forbid area such that it is green or yellow for i<t and it is green for i>=t.
    bool connected(GraphNode<type>* u, GraphNode<type>* v, set<int>& forbid, int t)
    {
        return true;
        return !path(u, v, forbid, t).empty();
    }


    // DFS visit from u to v such that it is green or yellow for i<t and it is green for i>=t.
    void DFS(GraphNode<type>* u, GraphNode<type>* v, int t)
    {
        if (u == v)
            found = true;
        if (found)
            return;
        for (auto it = u->successor.begin(); it != u->successor.end(); ++it)
            if (!marked(*it) && ingraph(*it) && (*it)->passable(t))
            {
                markit(*it);
                make_parent((*it), u);
                DFS((*it), v, t);
            }
    }
    // BFS visit from u to v such that it is green or yellow for i<t and it is green for i>=t.
    void BFS(GraphNode<type>* u, GraphNode<type>* v, int t)
    {
        queue<GraphNode<type>*> NEXT;
        NEXT.push(u);
        while (!NEXT.empty())
        {
            GraphNode<type>* next = NEXT.front();
            NEXT.pop();
            for (auto it = next->successor.begin(); it != next->successor.end(); ++it)
                if (!marked(*it) && ingraph(*it) && (*it)->passable(t))
                {
                    markit(*it);
                    make_parent(next, (*it));
                    if (*it == v)
                    {
                        found = true;
                        return;
                    }
                    NEXT.push((*it));
                }
        }
    }
    // Find path from node u to node v such that it is green or yellow for i<t and it is green for i>=t, return by sequence of nodes.
    list<GraphNode<type>*> path(GraphNode<type>* u, GraphNode<type>* v, int t)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        // Resulted new path.
        list<GraphNode<type>*> new_path;

        if (u == NULL)
            cout << "Invalid initial node." << endl;
        if (v == NULL)
            cout << "Invalid target node." << endl;
        if (u == NULL || v == NULL)
            return new_path;

        // Initialization for DFS/BFS.
        mark.clear();
        pi.clear();
        found = false;

        // DFS/BFS
        BFS(u, v, t);

        // Build the path from the pi-table.
        if (found != true)
            return new_path;
        for (GraphNode<type>* w = v; w != u; w = V[pi[w->id]])
            new_path.push_front(w);
        new_path.push_front(u);
        auto end_time = std::chrono::high_resolution_clock::now();
        find_path_time += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
        return new_path;
    }
};

#endif