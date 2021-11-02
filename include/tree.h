/*
 *  lib-io-tree-utils.h
 *  lib-io-tree
 *
 *  Created by defbond on 8/10/10.
 *  Copyright 2010 LIP/ENS-Lyon. All rights reserved.
 *
 */
#ifndef LIB_IO_TREE_UTILS_H
#define LIB_IO_TREE_UTILS_H

#ifdef __cplusplus

#include <ostream>
#include <iostream>
#include <stdio.h>
#include <list>

#include <vector>
#include <forward_list>
#include <assert.h>
#include <map>
#include "cluster.h"

#ifndef DEBUG_MEMUSAGE
#define DEBUG_MEMUSAGE 0
#endif
#ifndef VERBOSE
#define VERBOSE 0
#endif
#ifndef STRONG_ASSERT
#define STRONG_ASSERT 0
#endif

using namespace std;

#ifndef MAX_COMBI_SIZE
#define MAX_COMBI_SIZE 5
#endif

typedef enum {
    FURTHEST_NODE = 1,
    BEST_K_COMBI,
    BEST_FIT_ABS,
    FIRST_FIT_ABS,
    BEST_FIT,
    FIRST_FIT,
    BEST_INC_COMBI,
    BEST_COMBI,
    LARGEST_FIT,
    IMMEDIATELY
} io_method_t;

double u_wseconds(void);

class Task {
protected:
    bool cost_computed;
    double cost;
    double edge_weight;
    double node_weight;
    double MS_weight = 0;//assume execution time for any node is larger than 0
    double makespan_nocommu;
    bool makespan_computed = false;
    vector<Task *> *children;
    Task *parent;
    unsigned int parent_id;
    unsigned int id;
    bool broken = false;
    int label;
    double MS_sequentialPart, MS_parallelPart;
    double makespan_difference;
    unsigned int Qtree_id;

public :
    Task() {
        id = 0;
        Mpeak = 0;
        parent_id = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        children = new vector<Task *>();
    }

    Task(double nw, double ew, double mw) {
        id = 0;
        Mpeak = 0;
        parent_id = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        children = new vector<Task *>();

        edge_weight = ew;
        node_weight = nw;
        MS_weight = mw;
        makespan_nocommu = mw;
    }

    Task(unsigned int pparent_id, double nw, double ew, double mw) {
        id = 0;
        Mpeak = 0;
        Mavail = 0;
        parent = 0;
        cost_computed = false;
        makespan_computed = false;
        children = new vector<Task *>();

        edge_weight = ew;
        node_weight = nw;
        MS_weight = mw;
        makespan_nocommu = mw;
        parent_id = pparent_id;
    }

    ~Task() {
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
            delete *iter;
        }
        delete children;
    }

    void set_ms_diff(double slack) {
        makespan_difference = slack;
    }

    double get_ms_diff() {
        return makespan_difference;
    }

    void set_parent(Task *pparent) {
        this->parent = pparent;
    }

    void add_child(Task *pchild) {
        this->children->push_back(pchild);
        cost_computed = false;
    }

    vector<Task *> *get_children() {
        return children;
    }

    Task *get_child(unsigned int node_id) {
        return children->at(node_id);
    }

    Task *get_parent() {
        return parent;
    }

    bool is_leaf() const {
        return children->size() == 0;
    }

    bool is_root() const {
        return parent_id == 0;
    }
    
    double get_cost() {
        if (!cost_computed) {
            cost = edge_weight + node_weight;
            for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
                cost += (*iter)->get_edge_weight();
            }
            cost_computed = true;
        }
        return cost;
    }

    void set_parent_id(unsigned int pparent_id) {
        parent_id = pparent_id;
    }

    void set_id(unsigned int pid) {
        id = pid;
    }

    void set_label(int pid) {
        label = pid;
    }

    void set_edge_weight(double ew) {
        edge_weight = ew;
    }

    void set_node_weight(double nw) {
        node_weight = nw;
    }

    void set_makespan_weight(double mw) {
        MS_weight = mw;
    }

    unsigned int get_parent_id() const {
        return parent_id;
    }

    double get_edge_weight() const {
        return edge_weight;
    }

    double get_node_weight() const {
        return node_weight;
    }

    double get_makespan_weight() const {
        return MS_weight;
    }

    unsigned int get_id() const {
        return id;
    }

    int get_label() const {
        return label;
    }

    void print(ostream &out) const {
        out << max((unsigned int) 0, get_parent_id()) << " " << get_node_weight() << " " << get_edge_weight() << endl;
        for (vector<Task *>::iterator iter = children->begin(); iter != children->end(); iter++) {
            (*iter)->print(out);
        }

    }

    void break_edge() {
        broken = true;//break this edge
    }

    void restore_edge() {
        broken = false;//resotre this edge
    }

    bool is_broken() {
        if (broken == true) {
            return true;
        } else {
            return false;
        }
    }

    void update_makespan_cost() {
        makespan_nocommu = MS_sequentialPart + MS_parallelPart;
    }

    double get_sequential_part() {
        return MS_sequentialPart;
    }

    double get_parallel_part() {
        return MS_parallelPart;
    }

    double get_makespan_sequential(bool updateEnforce, double &MS_parallel) {
        if ((makespan_computed == true) & (updateEnforce == false)) {
            return MS_sequentialPart;
        }

        MS_sequentialPart = MS_weight;
        MS_parallelPart = 0;
        double temp;
        for (vector<Task *>::iterator iter = this->get_children()->begin(); iter != this->get_children()->end(); ++iter) {
            if ((*iter)->is_broken()) {
                //cout<<"edge "<<(*iter)->get_id()<<" broken"<<endl;
                temp = (*iter)->get_makespan_cost(true, updateEnforce);
                if (temp > MS_parallelPart) {
                    MS_parallelPart = temp;
                }
            } else {
                MS_sequentialPart += (*iter)->get_makespan_sequential(updateEnforce, MS_parallelPart);
                if (updateEnforce == true) {
                    (*iter)->update_makespan_cost();
                }
            }
        }

        if (MS_parallelPart > MS_parallel) {
            MS_parallel = MS_parallelPart;
        }

        return MS_sequentialPart;
    }

    double get_makespan_minus_comu() {
        if (Cluster::getFixedCluster()->isHomogeneous()) {
            return (makespan_nocommu - edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth());
        } else throw "Cluster not homogeneous";
    }

    double get_makespan_minus_w() {
        if (Cluster::getFixedCluster()->isHomogeneous()) {
            return (makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth() - MS_weight);
        } else throw "Cluster not homogeneous";
    }

    void set_makespan_uncomputed() {
        makespan_computed = false;
    }

    double get_makespan_cost(bool commulication = false, bool updateEnforce = false) {
        if (!Cluster::getFixedCluster()->isHomogeneous()) throw "Cluster not homogeneous";

        if ((makespan_computed == true) & (updateEnforce == false)) {
            if (commulication == true) {
                return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
            } else {
                return makespan_nocommu;
            }
        }

        MS_parallelPart = 0;
        MS_sequentialPart = this->get_makespan_sequential(updateEnforce, MS_parallelPart);//MS_parallelPart will be update here.
        makespan_nocommu = MS_sequentialPart + MS_parallelPart;

        makespan_computed = true;
        if (commulication == true) {
            //cout<<id<<"-"<<makespan_nocommu<<endl;//test
            return makespan_nocommu + edge_weight / Cluster::getFixedCluster()->getHomogeneousBandwidth();
        }

        //cout<<id<<"-"<<makespan_nocommu<<endl; //test
        return makespan_nocommu;
    }

    void set_other_side_id(unsigned int qtreeID) {
        Qtree_id = qtreeID;
    }

    unsigned int get_other_side_id() {
        return Qtree_id;
    }

    void remove_child(unsigned int childId) {
        for (vector<Task *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter) {
            if ((*iter)->get_id() == childId) {
                this->children->erase(iter);
                break;
            }
        }
    }

    void merge_to_parent() {
        this->get_parent()->set_makespan_weight(this->get_makespan_weight() + this->get_parent()->get_makespan_weight());
        this->get_parent()->remove_child(this->id);
        this->get_parent()->get_children()->insert(this->get_parent()->get_children()->end(), this->children->begin(),
                                                 this->children->end());
        //cout<<", children: ";
        for (vector<Task *>::iterator iter = this->children->begin(); iter != this->children->end(); ++iter) {
            //cout<<(*iter)->get_id()<<" ";
            (*iter)->set_parent(this->get_parent());
            (*iter)->set_parent_id(this->get_parent()->get_id());
        }
        //cout<<endl;
        this->children->clear();
        this->~Task();
    }

    unsigned int Ci;
    double Mpeak;
    double Mavail;

    double Sequence();

    double splt_subtrees(bool twolevel, list<Task *> &parallelRoots, unsigned long &sequentialLength);

    list<Task *> fillParallelRootsUntilBestMakespan(vector<double> &makespansOfSplittings,
                                                    unsigned long stepsUntilMinimalMakespan) const;


};

class Tree {
protected:
    vector<Task *> *nodes;
    unsigned int root_index;
    unsigned int root_count;
    unsigned int offset_id;
    unsigned int tree_id;

    Task * root;

    static Tree *originalTree;

public:

    Tree() {
        root_index = 0;
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        nodes = new vector<Task *>();
    }

    Tree(int N, int *prnts, double *nwghts, double *ewghts, double *mswghts) {
        root_index = 1;
        root_count = 0;
        offset_id = 0;
        tree_id = 1;
        nodes = new vector<Task *>();

        this->allocate_nodes(N);

        for (int i = 1; i < N + 1; i++) {
            //cout << "node id: " << i<<endl;
            Task *cur_node = this->get_node(i);
            cur_node->get_children()->clear();
            cur_node->set_edge_weight(ewghts[i]);
            cur_node->set_node_weight(nwghts[i]);
            cur_node->set_makespan_weight(mswghts[i]);
            cur_node->set_id(i);
            cur_node->set_label(i);
        }

        for (int i = 1; i < N + 1; i++) {
            Task *cur_node = this->get_node(i);

            if (prnts[i] > 0) {
                cur_node->set_parent_id(prnts[i]);
                cur_node->set_parent(this->get_node(prnts[i]));
                this->get_node(prnts[i])->add_child(cur_node);
            } else {
                cur_node->set_parent_id(0);
                this->set_root_id(i);
                this->set_tree_id(i);
            }
        }
    }

    Tree(vector<Task *> nodes, Tree *originalTree) {
        root_index = 1;
        root_count = 0;
        offset_id = 0;
        tree_id = 1;

        *(this->nodes) = nodes;
        this->originalTree = originalTree;
    }


    ~Tree() {
        if (root_index != 0 && nodes->size() > 0) {
            delete get_root();
        }

        delete nodes;
    }


    void print(ostream &out) const {
        out << nodes->size() << endl;

        for (vector<Task *>::iterator iter = nodes->begin(); iter != nodes->end(); iter++) {
            out << max((unsigned int) 0, (*iter)->get_parent_id()/*+1-offset_id*/) << " " << (*iter)->get_node_weight() << " "
                << (*iter)->get_edge_weight() << endl;
        }

    }


    void allocate_nodes(int new_node_count) {
        if (root_count > 0 && nodes->size() > 0) {
            delete get_root();
        }

        nodes->resize(new_node_count);

        unsigned int i = 0;
        for (vector<Task *>::iterator iter = nodes->begin(); iter != nodes->end(); iter++) {
            *iter = new Task();
            (*iter)->set_id(i++);
        }

        offset_id = nodes->front()->get_id();
    }
    void reverse_vector(){
        reverse(nodes->begin(),nodes->end());
    }

    void add_node(Task *newNode) {
        nodes->push_back(newNode);
    }

    void add_root(Task *newNode) {
        root_count++;
        assert(root_count == 1);
        nodes->push_back(newNode);
        root_index = nodes->size() - 1;
        this->root = newNode;
    }

    Task *get_root() const {
        assert (root_count ==1);
        return this->root;
    }

    unsigned int get_root_id() const {
        return this->get_root()->get_id();
    }

    void set_root_id(unsigned int root_id) {
        root_index = root_id;
    }

    void set_tree_id(unsigned int _id) {
        tree_id = _id;
    }

    Task *get_node(unsigned int node_id) const {
        Task * task =nodes->at(node_id - 1);
        if (task->get_id()!=node_id){
            throw runtime_error("Task not found for id " + to_string(node_id));
        }
        return task;
    }

    Task *get_node_by_position(unsigned int node_idx) const {
        assert(node_idx<this->get_nodes()->size());
        return nodes->at(node_idx);
    }

    const vector<Task *> *get_nodes() const {
        return nodes;
    }

   static void set_original_tree(Tree *origTree) {
        Tree::originalTree = origTree;
    }

    static Tree *get_original_tree() {
        return Tree::originalTree;
    }

    void print_broken_edges() {
        cout << "print broken edges" << endl;
        unsigned long treeSize = this->get_nodes()->size();
        for (unsigned int i = treeSize; i >= 1; --i) {
            Task *currentnode = this->get_node(i);
            if (currentnode->is_broken()) {
                cout << i << " ";
            }
        }
        cout << "End" << endl;


    }

    Tree *BuildQtree();

    unsigned int how_many_subtrees(bool quiet);

    bool memory_enough(Task *Qrootone, Task *Qroottwo, bool leaf, double memory_size, int *chstart, int *children);

    double improved_split(unsigned int number_processor, int *chstart, int *childrenID);

    double merge(unsigned int num_subtrees, int *chstart, int *childrenID, bool CheckMemory);

    double merge_v2(unsigned int num_subtrees, unsigned int processor_number, double const memory_size, int *chstart,
                   int *childrenID, bool CheckMemory);

    double split_again();

    double split_again_v2(unsigned int processor_number, unsigned int num_subtrees, std::map<int, int> &taskToPrc,
                        std::map<int, bool> &isProcBusy);

    vector<Task *> build_critical_path();
};


typedef list<int> schedule_t;

typedef map<unsigned int, double> io_map;
typedef pair<unsigned int, unsigned int> node_sche;
typedef pair<unsigned int, double> node_ew;


void parse_tree(const char *filename, int *N, int **prnts, double **nwghts, double **ewghts, double **mswghts);
Tree * read_tree(const char *filename);

extern "C"
{
void po_construct(const int N, const int *prnts, int **chstart, int **chend, int **children, int *root);
void poaux(const int *chstart, const int *children, int N, int r, int *por, int *label);
} /* closing brace for extern "C" */

double max_out_degree(Tree *tree, int quiet);

double max_out_degree(int N, double *nwghts, double *ewghts, int *chstart, int *children);

double io_counter(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                 double available_memory, bool divisible, int quiet, unsigned int &com_freq,
                 vector<unsigned int> *brokenEdges, io_method_t method);

double
io_counter_with_variable_mem(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                         vector<double> availableMemorySizesA2, int &currentProcessor,
                         std::map<int, int> &taskToPrc, std::map<int, bool> &isProcBusy, bool divisible, int quiet,
                         unsigned int &com_freq, vector<unsigned int> *brokenEdges, io_method_t method);

Tree *build_subtree(Tree *tree, Task *SubtreeRoot, unsigned int new_tree_size, int **prnts, double **ewghts,
                   double **timewghts, double **spacewghts, int *chstart, int *children);

void pop_smallest_roots_to_fit_to_cluster(list<Task *> &parallelRoots, unsigned long amountSubtrees);

void break_prepared_edges(Task *root, list<Task *> &parallelRoots);

double get_weight_pq(list<Task *> &parallelRoots, Task *currentNode);

double get_weight_surplus_from_smallest_node(list<Task *> &parallelRoots);

#endif
#endif
