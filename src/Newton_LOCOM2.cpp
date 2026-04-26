#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::cube CalculateXX(const arma::mat & X){
    
    const int K = X.n_cols;
    const int n_sam = X.n_rows;
    arma::cube XX(K, K, n_sam);
    
    arma::rowvec xs;
    for (int s = 0; s < n_sam; ++s){
        xs = X.row(s);
        XX.slice(s) = xs.t() * xs;
    }
    return XX;
}


arma::cube UpdateXX(const arma::cube & XX, 
                    const arma::mat & X, 
                    const arma::mat & Yr){
    
    const int K = XX.n_rows;
    const int n_sam = XX.n_slices;
    const int n_trait = Yr.n_cols;
    arma::cube XX_perm = XX;
    
    const arma::mat Yt = Yr.t();
    const arma::mat Xt_part = X.cols(n_trait, K-1).t();
    
    for (int s = 0; s < n_sam; ++s){
        const arma::subview_col<double> y = Yt.col(s);
        const arma::subview_col<double> x = Xt_part.col(s);
        
        arma::mat & slice = XX_perm.slice(s);
        slice.submat(0, 0, n_trait-1, n_trait-1) = y * y.t();
        slice.submat(0, n_trait, n_trait-1, K-1) = y * x.t();
        slice.submat(n_trait, 0, K-1, n_trait-1) = x * y.t();
    }
    return XX_perm;
}


arma::mat compute_weighted_XtWX(const arma::mat& X, const arma::vec& w) {
    return X.t() * (X.each_col() % w);
}


// [[Rcpp::export]]
List Newton_LOCOM2(arma::mat & rel_table,  // P_ij: relative abundance table
                   arma::mat & rel_table_sum,  // M_ij: rel_table_sum = rel_table + rel_table[,ref.otu]
                   arma::mat & X, 
                   arma::cube & XX,
                   arma::mat & beta_init, 
                   arma::mat & weight,
                   double tol, 
                   int iter_max, 
                   double Firth_thresh, 
                   bool robust_var,
                   arma::vec & prop_presence, 
                   bool get_var, 
                   double step_rate = 1.0) {
    
    const int n_otu = rel_table.n_cols;
    const int n_sam = rel_table.n_rows;
    int K = X.n_cols;
    
    arma::mat beta = beta_init;
    
    arma::vec z(n_sam), u(n_sam), DV(n_sam), DV_no_w(n_sam);
    arma::vec J_temp(n_sam), J_temp_1(n_sam), d_temp_1(n_sam);
    arma::vec H_temp(n_sam), d_temp(n_sam);
    arma::vec diag_quad_Sigma(n_sam), diag_quad_Jinv(n_sam);
    
    arma::mat J(K, K), J_inv(K, K), Sigma(K, K), e(K, K);
    arma::vec step(K), S(K);
    arma::vec n_steps(n_otu, fill::zeros);
    arma::mat S_all(K, n_otu);
    
    arma::cube variance;
    if (get_var) variance = cube(K, K, n_otu, fill::zeros);
    
    for (int j = 0; j < n_otu; j++){
        
        beta.col(j) = beta_init.col(j);
        const arma::vec w_j = weight.col(j);
        const arma::vec M_j = rel_table_sum.col(j);
        const arma::vec P_j = rel_table.col(j);
        
        int i = 0;
        for (; i < iter_max; i++){
            
            z = exp(X * beta.col(j));
            u = z / (1.0 + z);
            
            double weight_sum = dot(w_j, M_j);
            double C = 1.0 / weight_sum;
            
            for(int s=0; s<n_sam; ++s) {
                double M_ij = M_j[s];
                if(M_ij > 0) {
                    DV_no_w[s] = (P_j[s] / M_ij) - u[s];
                } else {
                    DV_no_w[s] = 0.0;
                }
            }
            
            DV = (P_j - M_j % u) % w_j * C;  
            S = n_sam * X.t() * DV;
            J_temp = - (u % (1.0 - u)) % w_j % M_j * C;
            J = compute_weighted_XtWX(X, J_temp);
            J_inv = inv(J);
            
            if (prop_presence(j) < Firth_thresh){
                
                if (robust_var) {
                    arma::vec sigma_weights = square(DV_no_w) % w_j % M_j * C;
                    
                    e = compute_weighted_XtWX(X, sigma_weights); 
                    Sigma = J_inv * e * J_inv;
                    
                } else {
                    Sigma = -J_inv;
                }
                
                J_temp_1 = J_temp % (1.0 - 2.0 * u);
                d_temp_1 = J_temp % DV_no_w;
                
                // Efficient way to compute diag(X * M * X'): sum( (X*M) % X, dim=1 )
                diag_quad_Sigma = sum((X * Sigma) % X, 1);
                diag_quad_Jinv  = sum((X * J_inv) % X, 1);
                
                for (int k = 0; k < K; k++) {
                    
                    arma::vec H_weights_k = J_temp_1 % X.col(k);
                    double trace_H_Sigma = dot(H_weights_k, diag_quad_Sigma);
                    
                    arma::vec d_weights_k = d_temp_1 % X.col(k);
                    double trace_Jinv_d = dot(d_weights_k, diag_quad_Jinv);
                    
                    double trace_H_Jinv = dot(H_weights_k, diag_quad_Jinv);
                    
                    S(k) = S(k) - 0.5 * trace_H_Sigma + trace_Jinv_d + (0.5 / n_sam) * trace_H_Jinv;;
                    
                }
            }
            
            step = (step_rate / n_sam) * J_inv * S;
            
            if (sum(abs(step)) > 10*K){
                beta.col(j).fill(0);
            } else {
                beta.col(j) -= step;
            }
            
            if (sum(abs(step) < tol) == K){
                break;
            }
        } // end iter
        
        n_steps(j) = i;
        S_all.col(j) = S;
        
        if (get_var) {
            if (!robust_var) {
                arma::vec final_weights = square(DV); 
                e = compute_weighted_XtWX(X, final_weights);
                Sigma = J_inv * e * J_inv;
            }
            variance.slice(j) = Sigma;
        }
        
    } // otu
    
    List L = List::create(beta, variance, n_steps, Sigma, e, S_all);
    return L;
    
} // Newton_LOCOM2()


// [[Rcpp::export]]
arma::cube perm_Newton_LOCOM2(arma::mat rel_table, arma::mat rel_table_sum, 
                              arma::mat Yr, arma::mat X, arma::cube XX, //(LR) replaced freq_table by rel_table
                              arma::mat beta_init, arma::mat weight,
                              arma::umat perm,
                              double tol, int iter_max, double Firth_thresh, bool robust_var,
                              arma::vec prop_presence, bool get_var, double step_rate = 1) {
    
    int n_trait = Yr.n_cols;
    int n_otu = rel_table.n_cols;
    int n_perm = perm.n_cols; 
    
    arma::cube beta_est_temp(n_trait, n_otu, n_perm);
    
    arma::mat Yr_perm;
    arma::mat X_perm = X;
    uvec o;
    
    
    for (int i_perm = 0; i_perm < n_perm; i_perm++){
        
        o = perm.col(i_perm);
        
        Yr_perm = Yr.rows(o);             
        X_perm.head_cols(n_trait) = Yr_perm;
        XX = UpdateXX(XX, X, Yr_perm);
        
        List res = Newton_LOCOM2(rel_table, rel_table_sum, X_perm, XX, beta_init, weight, tol, iter_max, Firth_thresh, robust_var, prop_presence, get_var, step_rate);
        arma::mat res_beta = res[0];
        beta_est_temp.slice(i_perm) = res_beta.head_rows(n_trait);
    }
    
    return beta_est_temp;
    
} // perm_Newton_LR()






