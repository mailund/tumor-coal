
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double twoTumors_CaseB_TA(int n, double delta, double gamma, double tau) {
    // Get x (see Wiuf's note)
    double x = 0.0;
    if (gamma != 1.0) {
        x = -n * gamma + sqrt((n*gamma)*(n*gamma) + 4*(1-gamma))/(2*(1-gamma));
    } else {
        x = 1.0 / float(n);
    }
    
    // M -- see Wiuf's note.
    double M = (x*pow(1-x, n-1))/pow(1-(1-gamma)*x, n+1);
    
    // Rejection sampling following Wiuf's note.
    for (;;) {
        double t = 0.0;
        for (;;) {
            double U = runif(1, 0.0, 1.0)[0];
            double X = 1/pow(U, 1.0/float(n)) - 1;
            t = 1/delta * log(1 + gamma/X);
            double s = delta*(t - tau);
            if (s >= 0) break;
        }
        double accept = runif(1, 0.0, 1.0)[0];
        double accept_prob_nom = (1.0/M) * exp(-t)*pow(1-exp(-t), n-1);
        double accept_prob_denom = pow(1-(1-gamma)*exp(-t), n+1);
        double accept_prob = accept_prob_nom / accept_prob_denom;
        if (accept <= accept_prob)
            return t;
    }
    
}
