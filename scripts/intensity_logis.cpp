// [[Rcpp::plugins(cpp11)]]
#include <cmath>  // std::pow
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <roptim.h>
// [[Rcpp::depends(roptim)]]
using namespace roptim;

arma::vec logistic(arma::vec x, double b){
  return 1/(1+exp(-b*x));
}

arma::vec dlogistic(arma::vec x, double b){
  arma::vec t1 = exp(-b*x);
  return t1/(1+t1)%(1+t1);
}

arma::vec ddlogis(arma::vec x){
  arma::vec t1 = exp(x);
  arma::vec t2 = exp(-x);
  arma::vec den = (t1+2.0+t2);
  return -(t1-t2)/(den%den);
}

// [[Rcpp::export]]
arma::vec d3logis(arma::vec x, arma::vec y, arma::vec z,
                  double bx, double by, double bz){
  return dlogistic(x,bx)%dlogistic(y,by)%dlogistic(y,bz);
}

arma::vec p3logis(arma::vec x, arma::vec y, arma::vec z,
                  double bx, double by, double bz){
  return logistic(x,bx)%logistic(y,by)%logistic(y,bz);
}


class LLL : public Functor {
private:
  arma::vec X_;
  arma::vec Y_;
  arma::vec Z_;
  arma::vec Tau_;
  int N = X_.n_rows;
  
public:
  LLL(const arma::vec &X, const arma::vec &Y, const arma::vec &Z, const arma::vec &Tau) : X_(X), Y_(Y), Z_(Z), Tau_(Tau) {
  }
  
  double operator()(const arma::vec &par) override {
    double a = exp(arma::as_scalar(par.row(0)));
    double bx = exp(arma::as_scalar(par.row(1)));
    double by = exp(arma::as_scalar(par.row(2)));
    double bz = exp(arma::as_scalar(par.row(3)));
    arma::vec loglambda(N);
    for(int i=0; i<N; i++){
      arma::vec Xmi = X_;
      arma::vec Ymi = Y_;
      arma::vec Zmi = Z_;
      Xmi.shed_row(i);
      Ymi.shed_row(i);
      Zmi.shed_row(i);
      loglambda.row(i) = log(sum(a*bx*by*bz*d3logis((arma::as_scalar(X_.row(i))-Xmi),(arma::as_scalar(Y_.row(i))-Ymi),(arma::as_scalar(Z_.row(i))-Zmi),bx,by,bz)));
    }
    return -sum(loglambda) + 
      sum(a*p3logis(arma::as_scalar(Tau_.row(0))-X_,arma::as_scalar(Tau_.row(1))-Y_,arma::as_scalar(Tau_.row(2))-Z_,bx,by,bz));
  }
  
  void Gradient(const arma::vec &par, arma::vec &gr) override {
    gr = arma::zeros<arma::vec>(4);
    arma::vec lambda(N);
    double a = exp(arma::as_scalar(par.row(0)));
    double bx = exp(arma::as_scalar(par.row(1)));
    double by = exp(arma::as_scalar(par.row(2)));
    double bz = exp(arma::as_scalar(par.row(3)));
    arma::vec dG1(N);
    arma::vec dG2(N);
    arma::vec dG3(N);
    for(int i=0; i<N; i++){
      arma::vec Xmi = X_;
      arma::vec Ymi = Y_;
      arma::vec Zmi = Z_;
      Xmi.shed_row(i);
      Ymi.shed_row(i);
      Zmi.shed_row(i);
      lambda.row(i) = sum(a*bx*by*bz*d3logis((arma::as_scalar(X_.row(i))-Xmi),(arma::as_scalar(Y_.row(i))-Ymi),(arma::as_scalar(Z_.row(i))-Zmi),bx,by,bz));
      dG1.row(i) = arma::as_scalar(lambda.row(i))+ 
        sum(a*bx*by*bz*(arma::as_scalar(X_.row(i))-Xmi)%ddlogis(bx*(arma::as_scalar(X_.row(i))-Xmi))%dlogistic((arma::as_scalar(Y_.row(i))-Ymi),by)%dlogistic((arma::as_scalar(Z_.row(i))-Zmi),bz));
      dG2.row(i) = arma::as_scalar(lambda.row(i))+ 
        sum(a*bx*by*bz*(arma::as_scalar(X_.row(i))-Xmi)%dlogistic((arma::as_scalar(X_.row(i))-Xmi),bx)%ddlogis((arma::as_scalar(Y_.row(i))-Ymi)*by)%dlogistic((arma::as_scalar(Z_.row(i))-Zmi),bz));
      dG3.row(i) = arma::as_scalar(lambda.row(i))+ 
        sum(a*bx*by*bz*(arma::as_scalar(X_.row(i))-Xmi)%dlogistic((arma::as_scalar(X_.row(i))-Xmi),bx)%dlogistic((arma::as_scalar(Y_.row(i))-Ymi),by)%ddlogis(bz*(arma::as_scalar(Z_.row(i))-Zmi)));
    }
    gr.row(0) = N/a - sum(a*p3logis(arma::as_scalar(Tau_.row(0))-X_,arma::as_scalar(Tau_.row(1))-Y_,arma::as_scalar(Tau_.row(2))-Z_,
                          bx,by,bz));
    gr.row(1) = -sum(dG1/lambda) + sum(a*bx*(arma::as_scalar(Tau_.row(0))-X_)%dlogistic(arma::as_scalar(Tau_.row(0))-X_,bx)%logistic(arma::as_scalar(Tau_.row(1))-Y_,by)%logistic(arma::as_scalar(Tau_.row(2))-Z_,bz));
    gr.row(2) = -sum(dG2/lambda) + sum(a*by*(arma::as_scalar(Tau_.row(1))-Y_)%dlogistic(arma::as_scalar(Tau_.row(1))-Y_,by)%logistic(arma::as_scalar(Tau_.row(0))-X_,bx)%logistic(arma::as_scalar(Tau_.row(2))-Z_,bz));
    gr.row(3) = -sum(dG3/lambda) + sum(a*bz*(arma::as_scalar(Tau_.row(2))-Z_)%dlogistic(arma::as_scalar(Tau_.row(2))-Z_,bz)%logistic(arma::as_scalar(Tau_.row(0))-X_,bx)%logistic(arma::as_scalar(Tau_.row(1))-Y_,by));
  }
};


// [[Rcpp::export]]
Rcpp::List intensity_est_logis_bfgs(arma::vec X, arma::vec Y, arma::vec Z, arma::vec Tau, arma::vec par) {
  LLL lll(X,Y,Z,Tau);
  Roptim<LLL> opt("BFGS");
  
  opt.minimize(lll, par);
  
  return Rcpp::List::create(
    Rcpp::Named("par") = opt.par(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("convergence") = opt.convergence()
  );
}
