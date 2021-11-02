#include "tree.h"
#include "cluster.h"

double io_counter(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                 double available_memory, bool divisible, int quiet, unsigned int &com_freq,
                 vector<unsigned int> *brokenEdges, io_method_t method);

double
io_counter_with_variable_mem(Tree *tree, int N, double *nwghts, double *ewghts, int *chstart, int *children, int *schedule,
                         Cluster *cluster, bool divisible, int quiet, unsigned int &com_freq,
                         vector<unsigned int> *brokenEdges, io_method_t method);
