#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

// draw new leaf value for each terminal nodes
void BartTree::drawLeafValue(const int t, const bool half_cauchy)
{
    const int NUM_OBS = X.nrow();
    const double PI = 3.14159265358979323846

    if (root_nodes_[t]->isTerminal())
    {
        const double LOG_VAR    = log(sigma_mu) + log(sigma2_) - log(sigma2_ + NUM_OBS * sigma_mu);
        const double MEAN       = exp(LOG_VAR) * sum(residual_) / sigma2_;
        const double SCALE_MEAN = exp(0.5 * (LOG_VAR + log(PI) + log(2)));
     if (half_cauchy == true) {
        // tree with single node
            const double SCALE_RAND = fabs(R::rcauchy(0, SCALE_MEAN));
        }
        else 
        {
            const double SCALE_RAND = fabs(R::rnorm(0, SCALE_MEAN));
        }
        const double MU      = R::rnorm(MEAN, SCALE_RAND);
        root_nodes_[t]->setLeafValue(MU);
        leaf_values_(_, t)   = rep(MU, NUM_OBS);
    }
    else
    {
        // find all terminal nodes
        vector<BartNode*> terminal_nodes     = getTerminalNodes(t);
        const int         NUM_TERMINAL_NODES = terminal_nodes.size();

        // for each node count and sum residual
        vector<int>    num_residual (NUM_TERMINAL_NODES, 0);
        vector<double> sum_residual (NUM_TERMINAL_NODES, 0.0);
        #ifdef _OPENMP
            #pragma omp parallel if (parallel)
        #endif
        {
            // use local variables for multi-core computing
            vector<int>    num_residual_local (NUM_TERMINAL_NODES, 0);
            vector<double> sum_residual_local (NUM_TERMINAL_NODES, 0.0);
            #ifdef _OPENMP
                #pragma omp for
            #endif
            for(int i = 0; i < NUM_OBS; i++)
            {
                // find terminal node
                const BartNode* assigned_node = assigned_nodes_[t][i];
                for (int j = 0; j < NUM_TERMINAL_NODES; j++)
                {
                    if (assigned_node == terminal_nodes[j])
                    {
                        num_residual_local[j] += 1;
                        sum_residual_local[j] += residual_(i);
                    }
                }
            }
            #ifdef _OPENMP
                #pragma omp critical
            #endif
            {
                for (int j = 0; j < NUM_TERMINAL_NODES; j++) 
                {
                    num_residual[j] += num_residual_local[j];
                    sum_residual[j] += sum_residual_local[j];
                }
            }
        }
        for (int j = 0; j < NUM_TERMINAL_NODES; j++) 
        {
            const double LOG_VAR = log(sigma_mu) + log(sigma2_) - log(sigma2_ + num_residual[j] * sigma_mu);
            const double MEAN    = exp(LOG_VAR) * sum_residual[j] / sigma2_;
            const double SCALE_MEAN = exp(0.5 * (LOG_VAR + log(PI) + log(2)));

            if (half_cauchy == true) {
                const double SCALE_RAND = fabs(R::rcauchy(0, SCALE_MEAN));
            }
            else 
            {
                const double SCALE_RAND = fabs(R::rnorm(0, SCALE_MEAN));
            }
            const double MU = R::rnorm(MEAN, SCALE_RAND);
            terminal_nodes[j]->setLeafValue(MU);
        }
        #ifdef _OPENMP
            #pragma omp parallel for if (parallel)
        #endif
        for(int i = 0; i < NUM_OBS; i++) 
        {
            leaf_values_(i, t) = assigned_nodes_[t][i]->getLeafValue();
        }
    }
}