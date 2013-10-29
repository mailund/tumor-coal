
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double twoTumors_CaseB_TA(int n1, int n2, double delta, 
                          double gamma1, double gamma2, double tau) {
    
    // Get x (see Wiuf's note)
    double x = 0.0;
    if (gamma1 != 1.0) {
        x = -n1 * gamma1 + sqrt((n1*gamma1)*(n1*gamma1) + 4*(1-gamma1))/(2*(1-gamma1));
    } else {
        x = 1.0 / float(n1);
    }
    
    // M -- see Wiuf's note.
    double M = (x*pow(1-x, n1-1))/pow(1-(1-gamma2)*x, n1+1);
    
    // Rejection sampling following Wiuf's note.
    for (;;) {
        double t = 0.0;
        for (;;) {
            double U = runif(1, 0.0, 1.0)[0];
            double X = 1/pow(U, 1.0/float(n2)) - 1;
            t = 1/delta * log(1 + gamma2/X);
            double s = delta*(t - tau);
            if (s >= 0) break;
        }
        double accept = runif(1, 0.0, 1.0)[0];
        double accept_prob_nom = (1.0/M) * exp(-t)*pow(1-exp(-t), n1-1);
        double accept_prob_denom = pow(1-(1-gamma1)*exp(-t), n1+1);
        double accept_prob = accept_prob_nom / accept_prob_denom;
        if (accept <= accept_prob)
            return t;
    }
    
}
