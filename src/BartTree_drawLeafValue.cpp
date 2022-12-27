#include <Rcpp.h>
#include "BartTree.h"

using namespace Rcpp;
using namespace std;

// draw new leaf value for each terminal nodes
void BartTree::drawLeafValue(const int t, const bool half_cauchy)
{
    const int NUM_OBS = X.nrow();
    const double PI = 3.14159265358979323846;
    const double QN = 0.67448975019608170545; 
    // 표준편차 수정필요
    double ssd = 0.05;
    //double sigma_mu_temp1 = sigma_mu;

    //double SCALE_RAND;

    if (root_nodes_[t]->isTerminal())
    {
        const BartNode* current_node = root_nodes_[t];
        double MEAN_ETA = current_node->getLeafValue();
        //const double MEAN_ETA->getLeafValue(t);
        const double LOG_VAR    = log(sigma_mu) + log(sigma2_) - log(sigma2_ + NUM_OBS * sigma_mu);
        const double MEAN       = exp(LOG_VAR) * sum(residual_) / sigma2_;
        const double GIVEN_LIK  = -0.5 * LOG_VAR - 0.5 / exp(LOG_VAR) * (MEAN_ETA - MEAN);
        //const double GIVEN_LIK  = -0.5 * LOG_VAR - 0.5 / exp(LOG_VAR) * (current_node - MEAN);
        //const double SCALE_MEAN = exp(0.5 * (LOG_VAR + log(PI) + log(2)));
        
        if (half_cauchy == true) {
        // tree with single node
           const double GIVEN_HC = -log(1 + 0.25 * sigma_mu / pow(ssd,2));
           const double GIVEN    = GIVEN_LIK + GIVEN_HC;
           
           const double DRAW_SIGMA_MU = pow(fabs(R::rcauchy(0, 4*pow(ssd, 2))), 2);
           const double prop_LOG_VAR  = log(DRAW_SIGMA_MU) + log(sigma2_) - log(sigma2_ + NUM_OBS * DRAW_SIGMA_MU);
           const double prop_MEAN     = exp(prop_LOG_VAR) * sum(residual_) / sigma2_;
           const double prop_lik      = -0.5 * prop_LOG_VAR - 0.5 / exp(prop_LOG_VAR) * (MEAN_ETA - prop_MEAN);
           //const double prop_lik      = -0.5 * prop_LOG_VAR - 0.5 / exp(prop_LOG_VAR) * (current_node - MEAN);
           const double prop_HC       = -log(1 + 0.25*DRAW_SIGMA_MU / pow(ssd,2));
           const double prop          = prop_lik + prop_HC;

           double ratio = prop_lik - GIVEN_LIK;
           if (ratio > log(R::runif(0,1)))
           {
            sigma_mu = DRAW_SIGMA_MU;
           }       

        }
        else 
        {
            const double GIVEN_SIGMA = pow(ssd/QN, 2);
            const double GIVEN_HN    = -sigma_mu/(2*GIVEN_SIGMA);
            const double GIVEN       = GIVEN_LIK + GIVEN_HN;

            const double DRAW_SIGMA_MU = pow(fabs(R::rnorm(0, sqrt(GIVEN_SIGMA))), 2);
            const double prop_LOG_VAR  = log(DRAW_SIGMA_MU) + log(sigma2_) - log(sigma2_ + NUM_OBS * DRAW_SIGMA_MU);
            const double prop_MEAN     = exp(prop_LOG_VAR) * sum(residual_) / sigma2_;
            const double prop_lik      = -0.5 * prop_LOG_VAR - 0.5 / exp(prop_LOG_VAR) * (MEAN_ETA - prop_MEAN);
            const double prop_HC       = -DRAW_SIGMA_MU;
            const double prop          = prop_lik + prop_HC;
            
            double ratio = prop_lik - GIVEN_LIK;
            if (ratio > log(R::runif(0,1)))
           {
            sigma_mu = DRAW_SIGMA_MU;
           } 
        }
        const double MU      = R::rnorm(MEAN, sqrt(sigma_mu));
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

            const BartNode* current_node = root_nodes_[t];
            double MEAN_ETA = current_node->getLeafValue();
            //const double MEAN_ETA->getLeafValue(t);
            const double LOG_VAR    = log(sigma_mu) + log(sigma2_) - log(sigma2_ + num_residual[j] * sigma_mu);
            const double MEAN       = exp(LOG_VAR) * sum_residual[j] / sigma2_;
            const double GIVEN_LIK  = -0.5 * LOG_VAR - 0.5 / exp(LOG_VAR) * (MEAN_ETA - MEAN);
            //const double GIVEN_LIK  = -0.5 * LOG_VAR - 0.5 / exp(LOG_VAR) * (current_node - MEAN);
            //const double SCALE_MEAN = exp(0.5 * (LOG_VAR + log(PI) + log(2)));
            
            //const double LOG_VAR = log(sigma_mu) + log(sigma2_) - log(sigma2_ + num_residual[j] * sigma_mu);
            //const double MEAN    = exp(LOG_VAR) * sum_residual[j] / sigma2_;
            //const double SCALE_MEAN = exp(0.5 * (LOG_VAR + log(PI) + log(2)));

            if (half_cauchy == true) {
                const double GIVEN_HC = -log(1 + 0.25 * sigma_mu / pow(ssd,2));
                const double GIVEN    = GIVEN_LIK + GIVEN_HC;
           
                const double DRAW_SIGMA_MU = pow(fabs(R::rcauchy(0, 4*pow(ssd, 2))), 2);
                const double prop_LOG_VAR  = log(DRAW_SIGMA_MU) + log(sigma2_) - log(sigma2_ + NUM_OBS * DRAW_SIGMA_MU);
                const double prop_MEAN     = exp(prop_LOG_VAR) * sum_residual[j] / sigma2_;
                const double prop_lik      = -0.5 * prop_LOG_VAR - 0.5 / exp(prop_LOG_VAR) * (MEAN_ETA - prop_MEAN);
                //const double prop_lik      = -0.5 * prop_LOG_VAR - 0.5 / exp(prop_LOG_VAR) * (current_node - MEAN);
                const double prop_HC       = -log(1 + 0.25*DRAW_SIGMA_MU / pow(ssd,2));
                const double prop          = prop_lik + prop_HC;

                double ratio = prop_lik - GIVEN_LIK;
                if (ratio > log(R::runif(0,1)))
                {
                    sigma_mu = DRAW_SIGMA_MU;
                }       
            }
            else 
            {
                const double GIVEN_SIGMA = pow(ssd/QN, 2);
                const double GIVEN_HN    = -sigma_mu/(2*GIVEN_SIGMA);
                const double GIVEN       = GIVEN_LIK + GIVEN_HN;

                const double DRAW_SIGMA_MU = pow(fabs(R::rnorm(0, sqrt(GIVEN_SIGMA))), 2);
                const double prop_LOG_VAR  = log(DRAW_SIGMA_MU) + log(sigma2_) - log(sigma2_ + NUM_OBS * DRAW_SIGMA_MU);
                const double prop_MEAN     = exp(prop_LOG_VAR) * sum_residual[j] / sigma2_;
                const double prop_lik      = -0.5 * prop_LOG_VAR - 0.5 / exp(prop_LOG_VAR) * (MEAN_ETA - prop_MEAN);
                const double prop_HC       = -DRAW_SIGMA_MU;
                const double prop          = prop_lik + prop_HC;
            
                double ratio = prop_lik - GIVEN_LIK;
                if (ratio > log(R::runif(0,1)))
                {
                sigma_mu = DRAW_SIGMA_MU;
                } 
            }
            const double MU = R::rnorm(MEAN, sqrt(sigma_mu));
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