#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
arma::mat scaleX_cpp(const arma::mat& X,
                     const double& tol = 1e-10){
  
  arma::mat Xc = X;
  
  // Column-center:
  Xc.each_row() -= mean(Xc, 0);
  
  // Column-scale:
  arma::rowvec s = sqrt(mean(Xc % Xc, 0));
  s.transform( [tol](double val) { return (val > tol) ? val : std::numeric_limits<double>::infinity(); } );

  Xc.each_row() %= (1 / s);
  
  // Return:
  return(Xc);
  
}

//' @export
// [[Rcpp::export]]
Rcpp::List scaleXZ_cpp(const arma::mat& X,
                       const arma::mat& Z,
                       const double& tol = 1e-10){
  
  arma::mat Xc = X;
  arma::mat Zc = Z;
  
  // Column-center:
  arma::rowvec m = mean(Xc, 0);
  Xc.each_row() -= m;
  Zc.each_row() -= m;
  
  // Column-scale:
  arma::rowvec s = sqrt(mean(Xc % Xc, 0));
  s.transform( [tol](double val) { return (val > tol) ? val : std::numeric_limits<double>::infinity(); } );
  s = 1 / s;
  
  Xc.each_row() %= s;
  Zc.each_row() %= s;
  
  // Return:
  return Rcpp::List::create(Rcpp::Named("Xc") = Xc,
                            Rcpp::Named("Zc") = Zc);
  
}

//' @export
// [[Rcpp::export]]
arma::uvec modulo(const arma::uvec& x, const int& m){
  
  arma::uvec y = x - arma::floor(x / m) * m;
  return(y);
  
}

//' @export
// [[Rcpp::export]]
arma::colvec Thomas_algorithm(const arma::colvec& v,
                              const int& method = 2){
  
  int m = v.n_elem;
  arma::colvec x(m); 
  if(method == 1){
    
    arma::colvec w = arma::cumsum(v);
    x(m - 1) = w(m - 1);
    for(int i = m - 2; i >= 0; i--){
      
      x(i) = x(i + 1) + w(i);
      
    }
    return(x);
    
  }else if(method == 2){
    
    arma::colvec w(m);
    arma::colvec a = arma::regspace(1, m) / arma::regspace(2, m + 1);
    w(0) = v(0) / 2;
    for(int i = 1; i < m; i++){
      
      w(i) = (v(i) + w(i - 1)) * a(i);
      
    }
    
    x(m - 1) = w(m - 1);
    for(int i = m - 2; i >= 0; i--){
      
      x(i) = w(i) + x(i + 1) * a(i);
      
    }
    return(x);
    
  }else if(method == 3){
    
    // arma::colvec w = arma::cumsum(reverse(v));
    arma::colvec w(m);
    w(0) = v(m - 1);
    for(int i = 1; i < m; i++){
      
      w(i) = w(i - 1) + v(m - i - 1);
      
    }
    x(m - 1) = w(m - 1);
    for(int i = m - 2; i >= 0; i--){
      
      x(i) = x(i + 1) + w(i);
      
    }
    return(reverse(x));
    
  }else{
    
    return(v);
    
  }
  
}

//' @export
// [[Rcpp::export]]
arma::colvec hW_solver_Thomas(const arma::uvec& W,
                              const arma::colvec& c){
  
  arma::uword maxW = c.n_elem - 1;
  int LW = W.n_elem;
  arma::colvec cW = -1 * c(W);
  arma::colvec hW(LW);
  arma::uvec ij;
  
  arma::mat X;
  arma::uvec diffW = arma::diff(W);
  
  int i = 0;
  int j = 0;
  
  while(i < LW){
    
    if(j == (LW - 1)){
      
      if(i == j){
        
        if(W(i) == maxW){
          
          hW(i) = cW(i);
          
        }else{
          
          hW(i) = cW(i) / 2;
          
        }
        
      }else{
        
        ij = arma::regspace<arma::uvec>(i, j);
        if(W(i) == 0){
          
          hW(ij) = Thomas_algorithm(cW(ij), 1);
          
        }else if(W(j) == maxW){
          
          hW(ij) = Thomas_algorithm(cW(ij), 3);
          
        }else{
          
          hW(ij) = Thomas_algorithm(cW(ij), 2);
          
        }
        
      }
      
      break;
      
    }else{
      
      if(diffW(j) == 1){
        
        j++;
        
      }else{
        
        if(i == j){
          
          if(W(i) == 0){
            
            hW(i) = cW(i);
            
          }else{
            
            hW(i) = cW(i) / 2;
            
          }
          
        }else{
          
          ij = arma::regspace<arma::uvec>(i, j);
          if(W(i) == 0){
            
            hW(ij) = Thomas_algorithm(cW(ij), 1);
            
          }else{
            
            hW(ij) = Thomas_algorithm(cW(ij), 2);
            
          }
          
        }
        
        j++;
        i = j;
        
      }
      
    }
    
  }
  
  return(hW);
  
}

//' @export
// [[Rcpp::export]]
arma::mat Custom_Active_Set(const arma::mat& Yhat,
                            const arma::mat& L,
                            const arma::mat& U,
                            const double& eps = 1e-8){
  
  // Create constraint matrix:
  arma::mat C = (arma::join_rows(Yhat, U) - arma::join_rows(L, Yhat)).t();
  // arma::mat Ct = C.t();
  double neps = -1 * eps;
  
  // Extract dimensions:
  arma::uword m = Yhat.n_cols;
  int n = Yhat.n_rows;
  
  if(C.min() > neps){
    
    return(arma::mat(n, m + 1, arma::fill::zeros));
    
  }
  
  // Initialize objects for algorithm:
  arma::mat H(n, m + 1);
  arma::colvec h(m + 1);
  arma::colvec h_old(m + 1);
  arma::colvec p(m + 1);
  arma::colvec hp_ratio;
  arma::colvec hW;
  arma::colvec hW2;
  arma::colvec mu;
  arma::colvec a(m + 2);
  arma::colvec b(m);
  arma::colvec c(m + 1);
  arma::uvec B;
  arma::uvec W;
  arma::uvec Wc = arma::regspace<arma::uvec>(0, m);
  arma::uvec ind1;
  arma::uvec ind2;
  arma::uvec j;
  arma::uvec Wcj;
  arma::uvec insertj;
  int LW = 0;
  bool keep_going;
  
  // full_indices = 1:m1
  for(int i = 0; i < n; i++){
    
    // Reset or zero-out objects for next loop iteration:
    h.zeros();
    h_old.zeros();
    hW.reset();
    mu.reset();
    W.reset();
    LW = 0;
    Wc = arma::regspace<arma::uvec>(0, m);
    c = C.col(i);
    a.zeros();
    b.zeros();
    p.zeros();
    hp_ratio.reset();
    insertj.reset();
    
    // Set variable to continue iterating:
    keep_going = true;
    
    while(keep_going){
      
      // Count how many active Lagrange indices there are:
      LW = W.n_elem;
      
      // Save Lagrange multiplier value:
      h_old = h;
      
      if(LW == 0){}else{
        
        if(LW == 1){
          
          if((W(0) > 0) & (W(0) < m)){
            
            hW = -c(W) * 0.5;
            
          }else{
            
            hW = -c(W);
            
          }
          
        }else{
          
          hW = hW_solver_Thomas(W, c);
          
        }
        
      }
      
      if(h.min() > neps){
        
        h(W) = hW;
        
        if(LW == 0){
          
          mu = c;
          
        }else{
          
          b.zeros();
          ind1 = find(W < m);
          ind2 = find(W > 0);
          
          b(W.elem(ind1)) = hW(ind1);
          b(W.elem(ind2) - 1) -= hW(ind2);
          
          a.zeros();
          a.rows(1, m) = b;
          a(modulo(Wc + 1, m + 1)) -= a(Wc);
          // a(modulo(Wc + 1, m + 2)) -= a(Wc);   // This is the correct one!
          
          mu = a(modulo(Wc + 1, m + 2)) + c(Wc);  
          // mu = a(Wc + 1) + c(Wc);  // This is the correct one!
          
        }
        
        if(mu.min() > neps){
          
          keep_going = false;
          
        }else{
          
          j = mu.index_min();
          Wcj = Wc(j);
          
          if(LW == 0){
            
            W = Wcj;
            Wc.shed_rows(j);
            
          }else{
            
            insertj.reset();
            insertj = find(W > Wcj(0));
            if(insertj.n_elem == 0){
              
              insertj = W.n_elem;
              
            }
            
            W.insert_rows(insertj(0), Wcj);
            Wc.shed_rows(j);
            
          }
          
        }
        
      }else{
        
        p = h - h_old;
        B = W(arma::find(hW <= neps));
        hp_ratio = -1 * h_old(B) / p(B);
        j = B(hp_ratio.index_min());
        insertj = find(Wc > j);
        if(insertj.n_elem == 0){
          
          insertj = Wc.n_elem;
          
        }
        
        Wc.insert_rows(insertj(0), j);
        W.shed_rows(arma::find(W == j));
        h = h_old - (h_old(j) / p(j)) * p;
        
      }
      
    }
    
    h(Wc).zeros();
    H.row(i) = h.t();
    
  }
  
  return(H);
  
}

//' @export
// [[Rcpp::export]]
arma::mat XDXt(const arma::mat& X,
               const arma::colvec& d){
  
  arma::mat Xd = X;
  Xd.each_row() %= d.t();
  return(X * Xd.t());
  
}

//' @export
// [[Rcpp::export]]
arma::colvec projA0(const arma::colvec& v,
                    const arma::colvec& h){
  
  // The vector h corresponds to the associated i^th Lagrange multiplier,
  // which has positive entries or zero entries. The entries which are positive
  // correspond to constraints which are exactly zero, hence we preserve the
  // associated columns of A as A0 = A( , bool(h > 0)). This function is for
  // projecting onto ColSp(A0).
  // Let n = sum(h > 0). Conducting this operation with the usual projection
  // equation [ A0 (A0' A0)^{-1} A0' v ], perhaps through SVD(A0), would be at
  // least O(n * m^2) complexity.  However, due to the special structure of A,
  // and therefore A0, the most involved operations scale with complexity
  // ~O(n log n), with ~O(m) memory access operations.
  
  // Get indices where h is positive:
  arma::uvec w = find(h);
  arma::uvec dw = arma::diff(w);
  
  // Get dimension of input/output vector:
  double m = v.n_elem;
  arma::colvec x(m, arma::fill::zeros);
  int i = 0;
  int j = 0;
  int LW = w.n_elem;
  arma::uvec ij;
  
  while(i < LW){
    
    if(j == (LW - 1)){
      
      if(w(i) == 0){
        
        ij = arma::regspace<arma::uvec>(i, j);
        x(w(ij)) = v(w(ij));
        
      }else if(w(j) == m){
        
        ij = arma::regspace<arma::uvec>(i, j);
        x(w(ij) - 1) = v(w(ij) - 1);
        
      }else{
        
        ij = arma::regspace<arma::uvec>(w(i) - 1, w(j));
        x(ij) = v(ij);
        x(ij) -= mean(x(ij));
        
      }
      
      j++;
      i = j;
      
    }else{
      
      if(dw(j) == 1){
        
        j++;
        
      }else{
        
        if(w(i) == 0){
          
          ij = arma::regspace<arma::uvec>(i, j);
          x(w(ij)) = v(w(ij));
          
        }else{
          
          ij = arma::regspace<arma::uvec>(w(i) - 1, w(j));
          x(ij) = v(ij);
          x(ij) -= mean(x(ij));
          
        }
        
        j++;
        i = j;
        
      }
      
    }
    
  }
  
  return(x);
  
}

//' @export
// [[Rcpp::export]]
Rcpp::List FRiSO_GSD(const arma::mat& X,
                     const arma::mat& Y,
                     const arma::colvec& gamma_init,
                     const arma::colvec& tauseq,
                     const arma::colvec& nudge,
                     const double& lower = 0,
                     const double& upper = 1,
                     const double& alpha = 1,
                     const double& eps = 0.0001,
                     const double& max_theta = 0.78539816,
                     const double& impulse = 1,
                     const int& maxIter = 1000){
  
  // Grab dimensions:
  double n = Y.n_rows;
  double m = Y.n_cols;
  double p = X.n_cols;
  
  // Initialize common double values:
  double tau;       // (current) tau value
  double tau_root;  // square root of tau value
  double tau_sq;    // square of tau value
  double der1;      // d / d\theta
  double der2;      // d2 / d\theta2
  double normv;     // norm of negative tangent gradient
  double theta;     // angle to rotate
  double max_theta_j;
  arma::colvec t1(p, arma::fill::zeros);
  arma::colvec t2(p, arma::fill::zeros);
  arma::colvec tc(p, arma::fill::zeros);
  
  // Initialize common integers:
  int n_tau = tauseq.n_elem; // number of tau values to evaluate at
  
  // Initialize common vectors for algorithm:
  arma::uvec w;                        // vector for defining A0
  arma::colvec gamma_new = gamma_init; // gamma vector for looping
  arma::colvec gamma_old = gamma_init; // other gamma vector for looping
  arma::colvec u(p);                   // normal vector for rotations
  arma::colvec v(p);                   // tangent vector for rotations
  arma::colvec sum_C(n);               // vector samples w/active constraints
  arma::colvec dn(p);                  // diagonal of N matrix
  arma::colvec grad(p);                // gradient vector
  arma::colvec grad_tan(p);            // tangent component of gradient vector
  arma::colvec L(n, arma::fill::value(lower)); // lower box constraint vector
  arma::colvec U(n, arma::fill::value(upper)); // upper box constraint vector
  
  // Set up matrices and vectors for function:
  arma::mat Xn = X / sqrt(n); // scaled X matrix (i.e. tilX)
  arma::mat Yt = Y.t();       // transposed Y (calculate once)
  arma::mat Xnt = Xn.t();     // transposed scaled X matrix (calculate once)
  arma::mat H(n, m + 1);      // Lagrange multiplier (capital Eta)
  arma::mat Yhat(n, m);       // unconstrained quantile (embedded) solution
  arma::mat Q(n, m);          // constrained quantile (embedded) solution
  arma::mat E(n, m);          // residuals matrix
  arma::mat tilE(m, n);       // adjusted residuals matrix
  arma::mat C(n, m - 1);      // active constraint matrix
  arma::mat V(m, n);          // matrix for second derivative
  arma::colvec Vi(m + 2);     // vector for holding (0, V[ , i]', 0)'
  arma::colvec AAAVi;         // vector for holding (A[ , w]' A[ , w])^{-1}A[ , w]' V[ , i]
  arma::colvec Vnew(m);       // vector for holding Proj{A[ , w]} V[ , i]
  arma::uvec mw;              // index vector for modulo(w) operations
  arma::mat W(m, n);          // other matrix for second derivative
  arma::uvec ind1;
  arma::uvec ind2;
  
  // Set up objects to store outputs:
  arma::mat GAMMA(p, tauseq.n_elem); // matrix to store all FRiSO solutions (as
  // gammas)
  
  // Two routes depending on how p and n compare:
  if(p > (1.2 * n)){
    
    // Initialize for (p > 1.2 * n) case:
    arma::mat YbarY = arma::mat(n, 1, arma::fill::ones) * (arma::mat(1, n, arma::fill::value(1 / n)) * Y) + Y;
    arma::mat Ginv(n, n);       // Xn D Xn' + I
    arma::mat GX(n, p);         // solve(Xn D Xn' + I, Xn)
    arma::mat GY(n, m);         // solve(Xn D Xn' + I, Y)
    arma::mat XtG(p, n);        // solve(Xn D Xn' + I, Xn)'
    arma::mat G_XnY_(n, p + m); // solve(Xn D Xn' + I, [Xn, Y])
    arma::mat DuvXtG(p, n);     // solve(Xn D Xn' + I, Xn Du Dv)'
    arma::mat YtGX(m, p);       // Y' solve(Xn D Xn' + I, Xn)
    arma::mat XtGY(p, m);       // Xn' solve(Xn D Xn' + I, Y)
    
    // Loop through tau:
    for(int t = 0; t < n_tau; t++){
      
      /////////////////////////////////////////////////////////////////////////
      // Setup with new tau:
      //
      // Set current tau value:
      tau = tauseq(t);
      tau_root = sqrt(tau);
      tau_sq = tau * tau;
      
      // Zero out tangent vectors:
      t1.zeros();
      t2.zeros();
      tc.zeros();
      
      // Reset maximum angle of rotation:
      max_theta_j = max_theta;
      
      // Update gamma value, 'nudging' it off the boundary:
      gamma_new += nudge(t);
      gamma_new *= tau_root / arma::norm(gamma_new);
      gamma_old = gamma_new;
      u = gamma_old / tau_root;
      
      // While loop:
      int j = 0;
      bool continueLoop = true;
      
      while(continueLoop & (j < maxIter)){
        
        ///////////////////////////////////////////////////////////////////////////
        // Set old values:
        gamma_old = gamma_new;
        u = gamma_old / tau_root;
        
        ///////////////////////////////////////////////////////////////////////////
        // Global solution:
        Ginv = XDXt(Xn, gamma_old % gamma_old);
        Ginv.diag() += 1;
        G_XnY_ = arma::solve(Ginv, arma::join_rows(Xn, Y), arma::solve_opts::likely_sympd);
        XtG = G_XnY_.head_cols(p).t();
        GY = G_XnY_.tail_cols(m);
        Yhat = YbarY - GY;
        XtGY = Xnt * GY;
        
        // Lagrange multiplier and constrained solution:
        H = Custom_Active_Set(Yhat, L, U);
        Q = Yhat + (H.head_cols(m) - H.tail_cols(m));
        sum_C = sum(H, 1);
        
        // Residuals matrix:
        E = Q - Y;
        
        ////////////////////////////////////////////////////////////////////////
        // Step 1: Calculate YtGX:
        // XtGY = Xnt * GY;
        // YtGX = Y.t() * GX;
        
        // Step 3: Loop through i = {1 ... n} to calculate tilE:
        tilE = E.t();
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              tilE.col(i) -= projA0(tilE.col(i), H.row(i).t());
              
            }
            
          }
        }
        
        // Step 4: Calculate dn:
        dn = sum(XtG % (XtGY * tilE), 1);
        
        // Step 4.25: Calculate gradient:
        grad = 4 * gamma_old % dn;
        
        // Step 4.5: Calculate tangential gradient:
        grad_tan = grad - u * arma::dot(u, grad);
        
        // Step 4.75: Calculate v:
        normv = arma::norm(grad_tan);
        v = -grad_tan / normv;
        
        // Step 5: Calculate GXDuDv.t():
        DuvXtG = XtG;
        DuvXtG.each_col() %= (u % v);
        
        // Step 6: Calculate V:
        V = XtGY.t() * DuvXtG;
        
        // Step 7: Calculate W:
        W = V * (Xn * DuvXtG);
        
        // Step 8: Calculate sum(W % tilE):
        der2 = -16 * tau_sq * arma::accu(W % tilE);
        
        // Step 9: Loop through i = {1 ... n} to calculate tilV and sum(V % tilV):
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              V.col(i) -= projA0(V.col(i), H.row(i).t());
              
            }
            
          }
        }
        
        der2 += 8 * tau_sq * arma::accu(V % V);
        
        //////////////////////////////////////////////////////////////////////
        // Finish first derivative:
        der1 = -tau_root * normv;
        
        // Finish second derivative:
        der2 += 4 * tau * sum(dn % (v % v - u % u));
        
        // Calculate angle theta:
        if(der2 == 0){
          
          theta = max_theta_j;
          
        } else {
          
          theta = alpha * abs(der1 / der2);
          if(theta > max_theta_j){
            
            theta = max_theta_j;
            
          }
          
        }
        
        if(false){ // without momentum:
          
          // Rotate gamma_old:
          // Rcout << theta << arma::endl;
          gamma_new = cos(theta) * gamma_old + sin(theta) * v * tau_root;
          gamma_new = gamma_new * tau_root / arma::norm(gamma_new);
          
        }else{ // with possible momentum:
          
          // t1 is the new line segment determined by the local gradient:
          t1 = v * theta * tau_root;
          
          // t2 is the old line segment determined by the old gradient.
          // tc is the linear combination of these.
          //  - If impulse = 1, then we are entirely using the local gradient.
          //  - If impulse = 0, then we are an "immovable object" (though the
          //    initial momentum from t2 is zero, so we are in fact an
          //    "unmoving object").
          tc = impulse * t1 + (1 - impulse) * t2;
          
          // The corresponding angle of rotation is given by theta = L / r:
          theta = arma::norm(tc) / tau_root;
          
          // We now rotate gamma_old in the direction of tc by theta:
          // Rcout << theta << arma::endl;
          gamma_new = cos(theta) * gamma_old + tc * (sin(theta) * tau_root / arma::norm(tc));
          gamma_new = gamma_new * tau_root / arma::norm(gamma_new);
          
          // We have rotated to gamma_new value; we need to calculate the
          // tangent line segment as if we had rotated it as well.
          t2 = theta * (cos(theta) * tc * tau_root - sin(theta) * gamma_old);
          
        }
        
        
        // Check if gamma_new is close to gamma_old:
        if(abs(gamma_new - gamma_old).max() < eps){
          
          continueLoop = false;
          
        }
        
        // Increment counter, max_theta_j:
        j++;
        max_theta_j = max_theta * 2 / (2 + j);
        
      }
      
      // Set the solution for the current tau:
      GAMMA.col(t) = gamma_new;
      // Rcout << j << arma::endl;
      
    }
    
  } else {
    
    // Initialize for (p <= 1.2 * n) case:
    arma::mat Ybar = arma::mat(n, 1, arma::fill::ones) * (arma::mat(1, n, arma::fill::value(1 / n)) * Y);
    arma::mat XtY = Xnt * Y;                      // Xn' Y
    arma::mat XtXty = arma::join_rows(Xnt, XtY); // [Xn', Xn'Y]
    arma::mat Sigma = Xnt * Xn;                   // Xn' Xn (estimated cov(X))
    arma::mat Ftinv(p, p);                        // Sigma D + I
    arma::mat FtXt(p, n);                         // solve(Sigma D + I, Xn')
    arma::mat FtXtY(p, m);                        // solve(Sigma D + I, Xn' Y)
    arma::mat Ft_XtXty_(p, n + m);                 // solve(Sigma D + I, [Xn', Xn' Y])
    arma::mat DFtXtY(p, m);                       // D solve(Sigma D + I, Xn' Y)
    arma::mat DuvFtXt(p, n);                      // Du Dv solve(Sigma D + I, Xn')
    
    // Loop through tau:
    for(int t = 0; t < n_tau; t++){
      
      /////////////////////////////////////////////////////////////////////////
      // Setup with new tau:
      //
      // Set current tau value:
      tau = tauseq(t);
      tau_root = sqrt(tau);
      tau_sq = tau * tau;
      
      // Reset maximum angle of rotation:
      max_theta_j = max_theta;
      
      // Update gamma value, 'nudging' it off the boundary:
      gamma_new += nudge(t);
      gamma_new *= tau_root / arma::norm(gamma_new);
      gamma_old = gamma_new;
      u = gamma_old / tau_root;
      
      // While loop:
      int j = 0;
      bool continueLoop = true;
      
      while(continueLoop & (j < maxIter)){
        
        ///////////////////////////////////////////////////////////////////////////
        // Set old values:
        gamma_old = gamma_new;
        u = gamma_old / tau_root;
        
        ////////////////////////////////////////////////////////////////////////
        // Step 1: Calculate FtXt:
        Ftinv = Sigma;
        Ftinv.each_row() %= (gamma_old % gamma_old).t();
        Ftinv.diag() += 1;
        
        // below: [p, p] [p, n + m]
        Ft_XtXty_ = arma::solve(Ftinv, XtXty, arma::solve_opts::fast);
        FtXt = Ft_XtXty_.head_cols(n);
        FtXtY = Ft_XtXty_.tail_cols(m);
        // FtXt = arma::solve(Ftinv, Xnt, arma::solve_opts::fast);
        
        // XF = FtXt.t(), [n, p] [p, p]
        
        // Step 2: Calculate FtXtY:
        // FtXtY = arma::solve(Ftinv, XtY, arma::solve_opts::fast);
        
        // Global solution:
        DFtXtY = FtXtY;
        DFtXtY.each_col() %= (gamma_old % gamma_old);
        Yhat = Ybar + Xn * DFtXtY;
        
        // Lagrange multiplier and constrained solution:
        H = Custom_Active_Set(Yhat, L, U);
        Q = Yhat + (H.head_cols(m) - H.tail_cols(m));
        sum_C = sum(H, 1);
        
        // Residuals matrix:
        E = Q - Y;
        
        // Step 3: Loop through i = {1 ... n} to calculate tilE:
        tilE = E.t();
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              tilE.col(i) -= projA0(tilE.col(i), H.row(i).t());
              
            }
            
          }
        }
        
        // Step 4: Calculate dn:
        // [p, n] % ( [p, m] * [m, n] )
        dn = sum(FtXt % (FtXtY * tilE), 1);
        
        // Step 4.25: Calculate gradient:
        grad = 4 * gamma_old % dn;
        
        // Step 4.5: Calculate tangential gradient:
        grad_tan = grad - u * arma::dot(u, grad);
        
        // Step 4.75: Calculate v:
        normv = arma::norm(grad_tan);
        v = -grad_tan / normv;
        
        // Step 5: Calculate DuvFtXt:
        DuvFtXt = FtXt;
        DuvFtXt.each_col() %= (u % v);
        
        // Step 6: Calculate V:
        V = FtXtY.t() * DuvFtXt;
        
        // Step 7: Calculate W:
        W = (V * Xn) * DuvFtXt;
        
        // Step 8: Calculate sum(W % tilE):
        der2 = -16 * tau_sq * arma::accu(W % tilE);
        
        // Step 9: Loop through i = {1 ... n} to calculate tilV and sum(V % tilV):
        if(sum(sum_C) == 0){}else{
          for(int i = 0; i < n; i++){
            
            // If sum_C(i) is non-zero, we need to adjust this row of Eadj:
            if(sum_C(i) > 0){
              
              V.col(i) -= projA0(V.col(i), H.row(i).t());
              
            }
            
          }
        }
        
        der2 += 8 * tau_sq * arma::accu(V % V);
        
        //////////////////////////////////////////////////////////////////////
        // Finish first derivative:
        der1 = -tau_root * normv;
        
        // Finish second derivative:
        der2 += 4 * tau * sum(dn % (v % v - u % u));
        
        // Calculate angle theta:
        if(der2 == 0){
          
          theta = max_theta_j;
          
        } else {
          
          theta = alpha * abs(der1 / der2);
          if(theta > max_theta_j){
            
            theta = max_theta_j;
            
          }
          
        }
        
        if(false){ // without momentum:
          
          // Rotate gamma_old:
          // Rcout << theta << arma::endl;
          gamma_new = cos(theta) * gamma_old + sin(theta) * v * tau_root;
          gamma_new = gamma_new * tau_root / arma::norm(gamma_new);
          
        }else{ // with possible momentum:
          
          // t1 is the new line segment determined by the local gradient:
          t1 = v * theta * tau_root;
          
          // t2 is the old line segment determined by the old gradient.
          // tc is the linear combination of these.
          //  - If impulse = 1, then we are entirely using the local gradient.
          //  - If impulse = 0, then we are an "immovable object" (though the
          //    initial momentum is zero).
          tc = impulse * t1 + (1 - impulse) * t2;
          
          // The corresponding angle of rotation is given by theta = L / r:
          theta = arma::norm(tc) / tau_root;
          
          // We now rotate gamma_old in the direction of tc by theta:
          // Rcout << theta << arma::endl;
          gamma_new = cos(theta) * gamma_old + tc * (sin(theta) * tau_root / arma::norm(tc));
          gamma_new = gamma_new * tau_root / arma::norm(gamma_new);
          
          // We have rotated to gamma_new value; we need to calculate the
          // tangent line segment as if we had rotated it as well.
          t2 = theta * (cos(theta) * tc * tau_root - sin(theta) * gamma_old);
          
        }
        
        // Check if gamma_new is close to gamma_old:
        if(abs(gamma_new - gamma_old).max() < eps){
          
          continueLoop = false;
          
        }
        
        // Increment counter, max_theta_j:
        j++;
        max_theta_j = max_theta * 2 / (2 + j);
        
      }
      
      // Set the solution for the current tau:
      GAMMA.col(t) = gamma_new;
      // Rcout << j << arma::endl;
      
    }
    
  }
  
  // Return output list:
  return Rcpp::List::create(Rcpp::Named("LAMBDA") = GAMMA % GAMMA);
  
}