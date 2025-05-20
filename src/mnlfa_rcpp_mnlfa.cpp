//// File Name: mnlfa_rcpp_mnlfa.cpp
//// File Version: 0.409



// [[Rcpp::depends(RcppArmadillo)]]

// #include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//*** constants
const double pi1 = 3.14159265358979;
const double sqpi2 = 0.398942280401433;


///********************************************************************
///** mnlfa_rcpp_plogis
// [[Rcpp::export]]
double mnlfa_rcpp_plogis( double x)
{
    double z=0;
    if (x<0){
        z = std::exp(x);
        z = z / ( 1 + z );
    } else {
        z = std::exp(-x);
        z = 1 / ( 1 + z);
    }
    //--- OUTPUT
    return z;
}
///********************************************************************


///********************************************************************
///** mnlfa_rcpp_calc_probs_2pl
// [[Rcpp::export]]
Rcpp::NumericMatrix mnlfa_rcpp_calc_probs_2pl(
    Rcpp::NumericMatrix a, Rcpp::NumericMatrix b, Rcpp::NumericMatrix theta,
    Rcpp::IntegerVector y, Rcpp::LogicalVector y_resp )
{
    int N = a.nrow();
    int TP = theta.nrow();
    Rcpp::NumericMatrix like(N,TP);
    like.fill(1);
    double temp=0;
    for (int nn=0;nn<N; nn++){
        if (y_resp[nn] ){
            for (int tt=0; tt<TP; tt++){
                temp = mnlfa_rcpp_plogis( a(nn,0)*theta(tt,0) - b(nn,0) );
                if (y[nn] == 1){
                    like(nn,tt) = temp;
                } else {
                    like(nn,tt) = 1 - temp;
                }
            }
        }
    }
    //--- OUTPUT
    return like;
}
///********************************************************************



///********************************************************************
///** mnlfa_rcpp_calc_dnorm
// [[Rcpp::export]]
Rcpp::NumericMatrix mnlfa_rcpp_calc_dnorm(
    Rcpp::NumericMatrix a, Rcpp::NumericMatrix b, Rcpp::NumericMatrix psi,
    Rcpp::NumericMatrix theta, Rcpp::NumericVector y, Rcpp::LogicalVector y_resp )
{
    int N = a.nrow();
    int TP = theta.nrow();
    Rcpp::NumericMatrix like(N,TP);
    like.fill(1);
    double temp=0;
    double t2=0;
    double ps=0;
    for (int nn=0;nn<N; nn++){
        if (y_resp[nn] ){
            ps = psi(nn,0) + 1e-30;
            t2 = sqpi2 / ps;
            for (int tt=0; tt<TP; tt++){
                temp = ( y[nn]-a(nn,0)*theta(tt,0)-b(nn,0) ) / ps;
                like(nn,tt) = t2 * std::exp( -0.5*temp*temp );
            }
        }
    }
    //--- OUTPUT
    return like;
}
///********************************************************************


///********************************************************************
///** mnlfa_rcpp_calc_probs_gpcm
// [[Rcpp::export]]
Rcpp::NumericMatrix mnlfa_rcpp_calc_probs_gpcm(
    Rcpp::NumericMatrix a, Rcpp::NumericMatrix b, Rcpp::NumericMatrix theta,
    Rcpp::NumericVector y, Rcpp::LogicalVector y_resp)
{
    int N = a.nrow();
    int TP = theta.nrow();
    int K = b.ncol();
    Rcpp::NumericMatrix like(N,TP);
    like.fill(1);
    Rcpp::NumericVector unnorm(K+1);
    unnorm[0]=1.0;
    double temp=0.0;
    double sump=1.0;

    for (int nn=0;nn<N; nn++){
        if (y_resp[nn] ){
            for (int tt=0; tt<TP; tt++){
                // KM <- sirt::sirt_matrix2(1:K, nrow=N)
                // M1 <- a[,1]*theta[tt,1]*KM - b
                // M1 <- cbind(1, exp( M1 ) )
                // M1 <- M1 / rowSums( M1 )
                // ind1 <- cbind( 1:N, y+1 )
                // like_ii[,tt] <- M1[ind1]
                sump=1.0;
                for (int kk=0; kk<K; kk++){
                    temp=a(nn,0)*theta(tt,0)*(kk+1)-b(nn,kk);
                    unnorm[kk+1] = std::exp(temp);
                    sump += unnorm[kk+1];
                }
                like(nn,tt) = unnorm[ y[nn] ] / sump;
            } // end tt
        } // end y_resp
    }  // end nn

    //--- OUTPUT
    return like;
}
///********************************************************************



///********************************************************************
///** mnlfa_rcpp_calc_probs_gpcm_deriv
// [[Rcpp::export]]
Rcpp::List mnlfa_rcpp_calc_probs_gpcm_deriv(
    Rcpp::NumericVector a, Rcpp::NumericMatrix b, Rcpp::NumericMatrix theta,
    Rcpp::NumericVector y, Rcpp::LogicalVector y_resp,
    Rcpp::NumericMatrix Xdes_int, Rcpp::NumericMatrix Xdes_slo,
    int ind, bool is_b, int uu )
{
    int N = Xdes_int.nrow();
    int TP = theta.nrow();
    int K = b.ncol();
    Rcpp::NumericMatrix H1(N,TP);
    Rcpp::NumericMatrix H2(N,TP);

    Rcpp::NumericVector unnorm(K+1);
    unnorm[0]=1.0;
    Rcpp::NumericVector fj(K+1);
    double temp=0.0;
    double sump=1.0;
    double Efj=0.0;
    double Vfj=0.0;
    double g1=0.0;
    double g2=0.0;
    double z=0;

    for (int nn=0;nn<N; nn++){
        if (y_resp[nn] ){
            for (int tt=0; tt<TP; tt++){
                sump=1.0;
                for (int kk=0; kk<K; kk++){
                    temp=a[nn]*theta(tt,0)*(kk+1)-b(nn,kk);
                    unnorm[kk+1] = std::exp(temp);
                    sump += unnorm[kk+1];
                }
                for (int kk=0; kk<K+1; kk++){
                    fj[kk] = unnorm[kk] / sump;
                }
                if ( (ind>K) | (!is_b) ){
                    Efj=0;
                    Vfj=0;
                    for (int kk=1; kk<K+1; kk++){
                        Efj += fj[kk]*kk;
                        Vfj += fj[kk]*kk*kk;
                    }
                    Vfj = Vfj - Efj*Efj;
                }

                if ( (ind<=K) & is_b){
                    //        g1 <- fj[,uu+1]
                    //        g2 <- ifelse(y==uu, g1-1, g1)
                    g1 = fj[uu];
                    if (y[nn]==uu){
                        g2=g1-1.0;
                    } else {
                        g2=g1;
                    }
                    H1(nn,tt)=g2;
                    H2(nn,tt)=-g1*(1-g1);
                }

                if ( (ind>K) & is_b){
                    // z <- Xdes_int[,ind-K+1]
                    z=Xdes_int(nn,ind-K);
                    // g2 <- -z*(y-Efj)
                    H1(nn,tt)=-z*(y[nn]-Efj);
                    // h2 <- -z^2*Vfj
                    H2(nn,tt)=-z*z*Vfj;
                }

                if (!is_b){
                    // z <- Xdes_slo[,ind]
                    z=Xdes_slo(nn,ind-1);
                    temp=theta(tt,0)*z;
                    // g2 <- theta[tt,1]*z*(y-Efj)
                    H1(nn,tt)=temp*(y[nn]-Efj);
                    // h2 <- -theta[tt,1]^2*z^2*Vfj
                    H2(nn,tt)=-temp*temp*Vfj;
                }

            } // end tt
        } // end y_resp
    }  // end nn

    //--- OUTPUT
    return Rcpp::List::create(
            Rcpp::Named("H1") = H1,
            Rcpp::Named("H2") = H2
        );
}
///********************************************************************


///********************************************************************
///** mnlfa_rcpp_mstep_trait_unidim_fun
// [[Rcpp::export]]
double mnlfa_rcpp_mstep_trait_unidim_fun( Rcpp::NumericMatrix theta,
    Rcpp::NumericMatrix mu_p, Rcpp::NumericMatrix sigma_p, Rcpp::NumericMatrix post,
    Rcpp::NumericVector weights)
{
    double theta_tt=0;
    double tmp=0;
    double tmp1=0;
    int N = mu_p.nrow();
    int TP = theta.nrow();
    double val=0;
    double const1 = -std::log(2*pi1)/2;
    double eps=1e-10;;
    for (int tt=0; tt<TP; tt++){
        theta_tt = theta(tt,0);
        for (int nn=0; nn<N; nn++){
            tmp = ( theta_tt - mu_p(nn,0) ) / ( sigma_p(nn,0)+eps);
            tmp1 = const1 - std::log(sigma_p(nn,0)+eps) - tmp*tmp/2;
            val += weights[nn]*post(nn,tt)*tmp1;
        }
    }
    //--- OUTPUT
    return val;
}
///********************************************************************


///********************************************************************
///** mnlfa_rcpp_mstep_trait_unidim_grad
// [[Rcpp::export]]
Rcpp::NumericVector mnlfa_rcpp_mstep_trait_unidim_grad( Rcpp::NumericMatrix theta,
        Rcpp::NumericVector mu_p, Rcpp::NumericVector sigma_p, Rcpp::NumericMatrix post,
        int IM, int IS, Rcpp::NumericMatrix Xdes_mean, Rcpp::NumericMatrix Xdes_sd,
        Rcpp::NumericVector weights)
{
    int NP = IM+IS;
    int N = mu_p.size();
    int TP = theta.nrow();
    Rcpp::NumericVector grad(NP);
    Rcpp::NumericVector H1b(N);
    Rcpp::NumericVector H2b(N);
    double sp_temp=0;

    // H1 <- post * (thetaM - mu_p) / sigma_p^2
    // for (hh in 1:IM){
    //    grad[hh] <- - sum( H1 * Xdes_mean[,hh] )
    // }

    for (int nn=0; nn<N; nn++){
        sp_temp = sigma_p[nn]*sigma_p[nn];
        for (int tt=0; tt<TP; tt++){
            H1b[nn] += post(nn,tt)*( theta(tt,0) - mu_p[nn]) / sp_temp;
        }
    }

    for (int hh=0;hh<IM; hh++){
        for (int nn=0; nn<N; nn++){
            grad[hh] -= weights[nn]*H1b[nn]*Xdes_mean(nn,hh);
        }
    }

    // ## derivative with respect to sigma is -1/sigma + (x-mu)^2/sigma^3
    // ## multiply it with inner derivative sigma
    // H2 <- post * ( - 1 + (thetaM - mu_p)^2 / sigma_p^2    )
    // for (hh in 1:IS){
    //    grad[IM+hh] <- - sum( H2  * Xdes_sd[,hh] )
    // }

    for (int nn=0; nn<N; nn++){
        for (int tt=0; tt<TP; tt++){
            sp_temp = ( theta(tt,0) - mu_p[nn]) / sigma_p[nn];
            H2b[nn] += post(nn,tt)*(sp_temp*sp_temp-1);
        }
    }

    for (int hh=0;hh<IS; hh++){
        for (int nn=0; nn<N; nn++){
            grad[IM+hh] -= weights[nn]*H2b[nn]*Xdes_sd(nn,hh);
        }
    }

    //--- OUTPUT
    return grad;
}
///********************************************************************
