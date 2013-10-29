
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double twoTumors_CaseB_TA(int n1, int n2, double delta, 
                          double gamma1, double gamma2, 
                          double tau,
                          int max_samples = 10) 
{
    int samples = 0;
    
    // Get x (see Wiuf's note)
    double x = 0.0;
    if (gamma1 != 1.0) {
        x = -n1 * gamma1 + sqrt((n1*gamma1)*(n1*gamma1) + 4*(1-gamma1))/(2*(1-gamma1));
    } else {
        x = 1.0 / float(n1);
    }
    
    // M -- see Wiuf's note.
    double M = (x*pow(1-x, n1-1))/pow(1-(1-gamma1)*x, n1+1);
    
    // Rejection sampling following Wiuf's note.
    for (;;) {
        double t = 0.0;
        double s = 0.0;
        for (;;) {
            double U = runif(1, 0.0, 1.0)[0];
            double X = 1/pow(U, 1.0/float(n2)) - 1;
            t = log(1 + gamma2/X) / delta;
            s = delta*(t - tau);
            if (s >= 0) break;
            
            if (samples++ == max_samples)
              return NA_REAL;
        }
        double accept = runif(1, 0.0, 1.0)[0];
        double accept_prob_nom = (1.0/M) * exp(-s)*pow(1-exp(-s), n1-1);
        double accept_prob_denom = pow(1-(1-gamma1)*exp(-s), n1+1);
        double accept_prob = accept_prob_nom / accept_prob_denom;
        if (accept <= accept_prob)
            return t;
        if (samples++ == max_samples)
            return NA_REAL;
    }
    
}

// [[Rcpp::export]]
double twoTumors_CaseB_TA_alt(int n1, int n2, double delta, 
                              double gamma1, double gamma2, 
                              double tau,
                              int max_samples = 10)
{
    int samples = 0;
    
    // Get x (see Wiuf's note)
    double x = 0.0;
    if (gamma2 != 1.0) {
        x = -n2 * gamma2 + sqrt((n2*gamma2)*(n2*gamma2) + 4*(1-gamma2)) / (2*(1-gamma2));
    } else {
        x = 1.0 / float(n2);
    }
    
    // M -- see Wiuf's note.
    double M = (x*pow(1-x, n2-1))/pow(1-(1-gamma2)*x, n2+1);
    
    // Rejection sampling following Wiuf's note.
    for (;;) {
        double U = runif(1, 0.0, 1.0)[0];
        double XA = 1.0/(pow(U,1.0/n1)) - 1.0;
        double t = log(1+gamma1/XA) / delta;
        double s = delta * (t + tau);
        
        double accept = runif(1, 0.0, 1.0)[0];
        double accept_prob_nom = (1.0/M) * exp(-s)*pow(1-exp(-s), n2-1);
        double accept_prob_denom = pow(1-(1-gamma2)*exp(-s), n2+1);
        double accept_prob = accept_prob_nom / accept_prob_denom;
        if (accept <= accept_prob)
            return t + tau;
            
        if (samples++ == max_samples)
            return NA_REAL;
    }
    
}


// [[Rcpp::export]]
double twoTumors_CaseB_TA_hybrid(int n1, int n2, double delta, 
                          double gamma1, double gamma2, 
                          double tau,
                          int max_samples = 10) 
{
    double t = NA_REAL;
    int samples = 0;
    for (int i = 0; i < max_samples; ++i) {
        if (i % 2 == 0) {
          t = twoTumors_CaseB_TA(n1, n2, delta, gamma1, gamma2, tau, 1);
        } else {
          t = twoTumors_CaseB_TA_alt(n1, n2, delta, gamma1, gamma2, tau, 1);
        }
        if (! R_IsNA(t))
          return t;
    }
    return NA_REAL;
}