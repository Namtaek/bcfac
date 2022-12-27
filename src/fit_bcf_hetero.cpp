#include <Rcpp.h>
#include "BartTree.h"
#include "ProgressBar.h"

using namespace Rcpp;
using namespace R;
using namespace std;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
void fit_bcf_hetero(
        NumericVector&       Y1,
        NumericVector&       Y0,
        NumericMatrix&       var_count,
        NumericVector&       var_prob,
        double&              sigma2_exp,
        NumericVector&       sigma2_out_hist,
        NumericVector&       dir_alpha_hist,
        const NumericVector& Y,
        const NumericVector& trt,
        const NumericMatrix& X,
        const double         trt_treated,
        const double         trt_control,
        const int            chain_idx,
        const int            num_chain,
        const int            total_iter,
        const int            num_burn_in,
        const int            num_thin,
        const int            num_post_sample,
        const int            num_tree,
        const NumericVector& step_prob,
        const double         alpha,
        const double         beta,
        const double         nu,
        const double         alpha2,
        const double         beta2,
        const double         nu2,
        const int            num_tree_mod,
        const double         lambda_exp,
        const double         lambda_out,
        const double         lambda_mod, //
        const int            boot_size,
        const bool           is_binary_trt,
        const bool           parallel,
        const bool           verbose
) {
    
    const int NUM_OBS = X.nrow();
    const int NUM_VAR = X.ncol();
    const int TRT_IDX = X.ncol();

    //NumericVector PS(NUM_OBS);
    //NumericVector PS_vec(NUM_OBS);
    //NumericMatrix X1(NUM_OBS, NUM_VAR + 1);

    
    // unique value of potential confounders
    vector<NumericVector> Xcut (NUM_VAR);
    for (int j = 0; j < NUM_VAR; j++)
    {
        NumericVector temp;
        temp = unique(X(_, j));
        temp.sort();
        Xcut[j] = clone(temp);
    }
    /*
    for (int j = 0; j < NUM_VAR + 1; j++)
    {
        NumericVector temp;
        if (j == TRT_IDX)
        {
            // keep unique values of trt in Xcut
            temp = unique(trt);
            temp.sort();
            Xcut[j] = clone(temp);
        } 
        else 
        {
            temp = unique(X(_, j));
            temp.sort();
            Xcut[j] = clone(temp);
        }
    }
    */

    /*
    NumericMatrix temp_seq(NUM_OBS, 1);
    temp_seq(_, 0) = seq_len(NUM_OBS);
    for (int i = 0; i < NUM_OBS; i++) {
        float j = i;
        temp_seq(i, 0) = j / (NUM_OBS + 1);
    }
    */

    //X1 = cbind(X, temp_seq);
    //vector<NumericVector> Xcut1 (NUM_VAR + 1);
    /*
    for (int j = 0; j < NUM_VAR + 1; j++)
    {
        NumericVector temp;
        temp = unique(X1(_, j));
        temp.sort();            
        Xcut1[j] = clone(temp);
    }
    */


    
    NumericVector latent_variable;
    //NumericVector latent_variable_before;
    {
        const double MEAN        = R::qnorm(mean(trt), 0, 1, true, false);
        const double SD          = 1.0;
        latent_variable          = Rcpp::rnorm(NUM_OBS, MEAN, SD);
        //latent_variable_before   = Rcpp::rnorm(NUM_OBS, MEAN, SD);
    }
    
    double sigma2_out = sigma2_out_hist(0);
    double sigma2_mod = sigma2_out_hist(0);
    double dir_alpha  = dir_alpha_hist(0);
    
    // sigma_mu based on min/max of Y, Y (Tr=1) and Y (Tr=0)    
    double sigma_mu_exp = max(
        pow(min(latent_variable) / (-2*sqrt(num_tree)), 2),
        pow(max(latent_variable) / ( 2*sqrt(num_tree)), 2)
    );
    double sigma_mu_out = max(
        pow(min(Y) / (-2*sqrt(num_tree)), 2),
        pow(max(Y) / ( 2*sqrt(num_tree)), 2)
    );
    double sigma_mu_mod = max(
        pow(min(Y) / (-2*sqrt(num_tree)), 2),
        pow(max(Y) / ( 2*sqrt(num_tree)), 2)
    );


    // Initial values of R
    NumericVector residual_exp = clone(latent_variable);
    NumericVector residual_out = clone(Y);
    NumericVector residual_mod = clone(residual_out);
    //NumericVector residual_exp_before = clone(latent_variable);

    // Initial values for the selection probabilities
    //NumericVector post_dir_alpha = rep(1.0, NUM_VAR + 1); 
    NumericVector post_dir_alpha = rep(1.0, NUM_VAR);
    NumericVector mod_dir_alpha = rep(1.0, NUM_VAR);
    
    // Obtaining namespace of MCMCpack package
    // Then pick up rinvgamma() and rdirichlet() function from MCMCpack package
    Environment MCMCpack   = Environment::namespace_env("MCMCpack");
    Function    rinvgamma  = MCMCpack["rinvgamma"];
    Function    rdirichlet = MCMCpack["rdirichlet"];
    
    // initialize model
    BartTree exposure = BartTree(
        residual_exp, var_prob, sigma2_exp,   // mutable variables
        1, trt, X, Xcut, step_prob, num_tree, // const variables 1
        alpha, beta, sigma_mu_exp, parallel   // const variables 2
    );

    /*
    BartTree exposure_before = BartTree(   
        residual_exp_before, var_prob, sigma2_exp,   // mutable variables
        1, trt, X, Xcut, step_prob, num_tree, // const variables 1
        alpha, beta, sigma_mu_exp, parallel   // const variables 2
    );
    exposure_before.updateLatentVariable(latent_variable_before, is_binary_trt);
    exposure_before.step(latent_variable_before, is_binary_trt, false);

    PS_vec = exposure_before.getFittedValues();
    for (int l = 0; l < NUM_OBS; l++) {
        PS[l] = R::pnorm(PS_vec[l], 0 ,1, true, false);
    }
    */
    //X1 = cbind(X, PS);

    
    BartTree outcome  = BartTree(
        residual_out, var_prob, sigma2_out,   // mutable variables
        0, trt, X, Xcut, step_prob, num_tree, // const variables 1
        alpha, beta, sigma_mu_out, parallel   // const variables 2
    );

    
    BartTree modifier = BartTree(
        residual_mod, var_prob, sigma2_mod,   // mutable variables
        0, trt, X, Xcut, step_prob, num_tree, // const variables 1
        alpha2, beta2, sigma_mu_mod/10, parallel   // const variables 2
    );
    
    
    
    IntegerVector boot_idx = sample(NUM_OBS, boot_size, true) - 1;
    
    int thin_count = 0, res_idx = 0;

    // init ProgressBar
    ProgressBar progress_bar = ProgressBar(chain_idx, num_chain, total_iter, 40); 
    
    // run MCMC
    for (int iter = 1; iter < total_iter + 1; iter++)
    {
        // check for break from R
        Rcpp::checkUserInterrupt();
        if (verbose)
            progress_bar.print(iter);
        
        // update latent_variable
        exposure.updateLatentVariable(latent_variable, is_binary_trt);
        //std::cout << "\n" << "Latent Variable updated " << "\n";
        
        // update tree
        exposure.step(latent_variable, is_binary_trt, false, false);
        //std::cout << "exposure Variable updated " << "\n";

        outcome.step(Y, is_binary_trt, true, false);
        //std::cout << "outcome updated " << "\n";

        // update sigma
        if (!is_binary_trt)
            exposure.updateSigma2(rinvgamma, Y, nu, lambda_exp);
        outcome.updateSigma2(rinvgamma,  Y, nu, lambda_out);
        modifier.updateSigma2(rinvgamma,  residual_out, nu2, lambda_mod);
        sigma2_out_hist(iter) = outcome.getSigma2();
        //std::cout << "sigma2 updated " << "\n";

        // count included
        NumericVector var_count_exp = exposure.countSelectedVariables();
        NumericVector var_count_out =  outcome.countSelectedVariables();
        
        // MH to update dir_alpha
        exposure.updateDirAlpha(dir_alpha);
        //std::cout << "dirichlet alpha updated " << "\n";
        dir_alpha_hist(iter) = dir_alpha;
        
        // then update post_dir_alpha
        //post_dir_alpha = rep(dir_alpha / (NUM_VAR + 1), NUM_VAR + 1);
        //post_dir_alpha = rep(dir_alpha, NUM_VAR);
        post_dir_alpha = rep(dir_alpha / NUM_VAR, NUM_VAR) + var_count_exp + var_count_out;

        // MH algorithm to update inclusion probabilities
        /*
        exposure.updateVarProb(
            rdirichlet, post_dir_alpha, var_count_exp, var_count_out
        );
        */
        var_prob = rdirichlet(1, post_dir_alpha);

        modifier.step(residual_out, is_binary_trt, false, true);
        
        // sample E[Y(1) - Y(0)]
        if (iter > num_burn_in) 
        {
            if (thin_count == num_thin)
            {
                // record selected variables
                var_count(res_idx, _) = var_count_out;
                
                // predict effect and potential outcomes
                Y1(res_idx) = outcome.predict(boot_idx, trt_treated) + modifier.predict(boot_idx, trt_treated);
                Y0(res_idx) = outcome.predict(boot_idx, trt_control);
                if (res_idx == num_post_sample)
                    break;
                res_idx++;
                thin_count = 0;
            } 
            else 
            {
                thin_count++;
            }
        }
    // end of an MCMC iteration
    }
}
