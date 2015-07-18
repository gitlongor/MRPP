#include <RcppArmadilloExtensions/sample.h>
// #include <omp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

NumericVector csample_num(NumericVector x,int size,bool replace,
                          NumericVector prob = NumericVector::create()){
  NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

// function that calculate the Euclidean norm of two vectors
// y1 and y2 must have same length
// [[Rcpp::export]]
double DistEUC(NumericVector y1, NumericVector y2){
  int p = y1.size();
  if (p != y2.size())
  {
    cout << "y1 and y2 have different size!" << endl;
    return -1;
  }
  double temp = 0.0;
  for (int i = 0; i < p; i++){
    temp += pow(y1[i] - y2[i], 2.0);
  }
  return sqrt(temp);
}

double DistPLLR(NumericVector y1, NumericVector y2, int ls1, int ls2){
  int p = y1.size();
  if (p != y2.size())
  {
    cout << "y1 and y2 have different size!" << endl;
    return -1;
  }
  double temp1 = 0.0;
  for (int i = 0; i < p; i++){
    double temp = 0.0;
    if (y1[i] == 0 && y2[i] == 0) temp= 0.0;
    else if (y1[i] == 0 && y2[i] != 0) temp = 2.0*log(double(ls1+ls2)/double(ls2))*double(y2[i]);
    else if (y1[i] != 0 && y2[i] == 0) temp = 2.0*log(double(ls1+ls2)/double(ls1))*double(y1[i]);
    else temp = -2.0*(double(y1[i]+y2[i])*log(double(y1[i]+y2[i]))-double(y1[i])*log(double(y1[i]))-double(y2[i])*log(double(y2[i]))+double(y1[i])*log(double(ls1))+double(y2[i])*log(double(ls2))-double(y1[i]+y2[i])*log(double(ls1+ls2)));
    temp1 += temp;    
  }
  return temp1;
}

// @param method the distance measure, "EUC" means Euclidean norm, "PLLR" means
//               distance measure based on Poisson log-likelihood ratio statistic 
double Dist(NumericVector y1, NumericVector y2, int ls1, int ls2,
              String method) {
  if (method == "EUC") return DistEUC(y1, y2);
  else if (method == "PLLR") return DistPLLR(y1, y2, ls1, ls2);
  else {
    cout << "Wrong distance measure!" << endl;
    return -1;
  }
}

NumericMatrix transpose(NumericMatrix x) {
  arma::mat y = as<arma::mat>(x) ;
  return wrap(trans(y));
}

// @param y1 and @param y2 transformed data
// @param RNA_Seq TRUE or FALSE, indicate the data is from RNA_Seq or not
// [[Rcpp::export]]
NumericMatrix DistMatrixEUC(NumericMatrix y1, NumericMatrix y2, bool RNA_Seq) {
  if (RNA_Seq && y1.ncol() > 1 && y1.nrow() > 1) {
    y1 = transpose(y1);
    y2 = transpose(y2);
  }
  int n1 = y1.nrow();
  int n2 = y2.nrow();
  int n = n1 + n2;
  NumericMatrix result(n,n);
  for (int i = 0; i < n1; i++) {
    result(i,i) = DistEUC(y1(i,_), y1(i,_));
  }
  for (int i = n1; i < n; i++) {
    result(i,i) = DistEUC(y2(i - n1,_), y2(i - n1,_));
  }
  for (int i = 0; i < n1; i++) {
    for (int j = (i + 1); j < n1; j++) {
      result(i,j) = DistEUC(y1(i,_), y1(j,_));
    }
  }
  for (int i = n1; i < n; i++) {
    for (int j = (i + 1); j < n; j++) {
      result(i,j) = DistEUC(y2(i - n1,_), y2(j - n1,_));
    }
  }
  for (int i = 0; i < n1; i++) {
    for (int j = n1; j < n; j++) {
      result(i,j) = DistEUC(y1(i,_), y2(j - n1,_));
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = (i + 1); j < n; j++) {
      result(j,i) = result(i,j);
    }
  }
  return result;
}

// [[Rcpp::export]]
double deltaZ(NumericMatrix disMatrix, NumericVector perm, int n1) {
  int n = disMatrix.nrow();
  double z, temp1, temp2, c1, c2;
  temp1 = 0.0;
  temp2 = 0.0;
  c1 = (double(n1) - 1.0)/(double(n) - 2.0);
  c2 = 1.0 - c1;
  int n2 = n - n1;
  for (int i = 0; i < n1; i++){
    for (int j = (i+1); j < n1; j++){
      temp1 += disMatrix(perm[i]-1,perm[j]-1);
    }
  }
  for (int i = n1; i < n; i++){
    for (int j = (i+1); j < n; j++){
      temp2 += disMatrix(perm[i]-1,perm[j]-1);
    }
  }
  z = 2.0*c1*temp1/(double(n1)*(double(n1)-1.0)) + 2.0*c2*temp2/(double(n2)*(double(n2)-1.0));
  return z;
}

// [[Rcpp::export]]
double MRPP_pvalue_EUC(NumericMatrix y1, NumericMatrix y2, bool RNA_Seq, int B) {
  if (RNA_Seq && y1.ncol() > 1 && y1.nrow() > 1) {
    y1 = transpose(y1);
    y2 = transpose(y2);
  }
  int n1 = y1.nrow();
  int n2 = y2.nrow();
  int n = n1 + n2;
  NumericMatrix disMatrix = DistMatrixEUC(y1, y2, RNA_Seq);
  double z0;
  int b;
  NumericVector d(n);
  for (int i = 0; i < n; i++){
    d[i] = i + 1;
  }
  z0 = deltaZ(disMatrix, d, n1);

  bool replace1=FALSE;

  NumericVector output(B);
  
//  int max=omp_get_max_threads();
//  omp_set_num_threads(max-4);

//  #pragma omp parallel for
  for (b = 0; b < B; b++) {
    NumericVector perm = csample_num(d,n,replace1);
    output[b] = deltaZ(disMatrix, perm, n1);
  }

  double result = 0.0;
  for(int b=0; b<B; b++){
    if(output[b] < z0){
      result += 1.0;
    }
  }
  return (result/double(B));
}

NumericVector rdeltaZ_EUC(NumericMatrix y1, NumericMatrix y2, NumericMatrix disMatrix,
                      NumericVector perm, bool RNA_Seq) {
  if (RNA_Seq && y1.ncol() > 1 && y1.nrow() > 1) {
    y1 = transpose(y1);
    y2 = transpose(y2);
  }  
  int n1 = y1.nrow();
  int n2 = y2.nrow();
  int p = y1.ncol();
  int n = n1 + n2;
  
  NumericMatrix c(n,p);
  for (int i = 0; i < n1; i ++) {
    c(i,_) = y1(i,_);
  }
  for (int i = n1; i < n; i++) {
    c(i,_) = y2(i - n1,_);
  }
  NumericVector rz(p);
  double temp1, temp2, c1, c2;
  c1 = (double(n1) - 1.0)/(double(n) - 2.0);
  c2 = 1.0 - c1;
  for (int r = 0; r < p; r++) {
    temp1 = 0.0;
    temp2 = 0.0;
    for (int i = 0; i < (n1 - 1); i++){
      for (int j = (i + 1); j < n1; j++){
        temp1 += pow(c(perm[i] - 1,r) - c(perm[j] - 1,r), 2.0) / (2 * disMatrix(perm[i] - 1, perm[j] - 1));
      }
    }
    for (int i = n1; i < (n - 1); i++){
      for (int j = (i + 1); j < n; j++){
        temp2 += pow(c(perm[i] - 1,r) - c(perm[j] - 1,r), 2.0) / (2 * disMatrix(perm[i] - 1, perm[j] - 1));
      }
    }
    rz[r] = 2.0 * c1 * temp1 / (double(n1) * (double(n1) - 1.0)) + 2.0 * c2 * temp2 / (double(n2) * (double(n2) - 1.0));
  }
  return rz;
}

// [[Rcpp::export]]
NumericVector Importance(NumericMatrix y1, NumericMatrix y2, bool RNA_Seq, int B) {
  if (RNA_Seq && y1.ncol() > 1 && y1.nrow() > 1) {
    y1 = transpose(y1);
    y2 = transpose(y2);
  }
  int n1 = y1.nrow();
  int n2 = y2.nrow();
  int p = y1.ncol();
  int n = n1 + n2;
  NumericMatrix disMatrix = DistMatrixEUC(y1, y2, RNA_Seq);
  NumericVector rz0;
  int b;
  NumericVector d(n);
  for (int i = 0; i < n; i++){
    d[i] = i + 1;
  }
  rz0 = rdeltaZ_EUC(y1, y2, disMatrix, d, RNA_Seq);

  bool replace1=FALSE;
  NumericVector output(p);
  for(int i = 0; i < p; i++) {
    output[i] = 0.0;
  }
  NumericVector temp(p);
  for (b = 0; b < B; b++) {
    NumericVector perm = csample_num(d, n, replace1);
    temp = rdeltaZ_EUC(y1, y2, disMatrix, perm, RNA_Seq);
    for (int r = 0; r < p; r++) {
      output[r] += temp[r];
    }
  }

  NumericVector result(p);
  for (int r = 0; r < p; r++){
    result[r] = rz0[r] - output[r]/double(B);
  }
  return (result);
}

// [[Rcpp::export]]
NumericVector Importance1(NumericMatrix y1, NumericMatrix y2, bool RNA_Seq) {
  if (RNA_Seq && y1.ncol() > 1 && y1.nrow() > 1) {
    y1 = transpose(y1);
    y2 = transpose(y2);
  }  
  int n1 = y1.nrow();
  int n2 = y2.nrow();
  int p = y1.ncol();
  int n = n1 + n2;
  
  NumericMatrix c(n,p);
  for (int i = 0; i < n1; i ++) {
    c(i,_) = y1(i,_);
  }
  for (int i = n1; i < n; i++) {
    c(i,_) = y2(i - n1,_);
  }
  
  NumericMatrix disMatrix = DistMatrixEUC(y1, y2, RNA_Seq);
  
  NumericVector rz(p);
  double temp1, temp2, temp3, c1, c2;
  c1 = (double(n1) - 1.0)/(double(n) - 2.0);
  c2 = 1.0 - c1;
  for (int r = 0; r < p; r++) {
    NumericMatrix temp(n,n);
    for (int i = 0; i < n; i++) {
      temp(i,i) = 0.0;
      for (int j = (i + 1); j < n; j++) {
        temp(i,j) = disMatrix(i,j) - sqrt(pow(disMatrix(i,j), 2.0) - pow(c(i,r) - c(j,r), 2.0));
      }
    }
    temp1 = 0.0;
    temp2 = 0.0;
    temp3 = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = (i + 1); j < n; j++) {
        temp1 += temp(i,j);
      }
    }
    for (int i = 0; i < (n1 - 1); i++){
      for (int j = (i + 1); j < n1; j++){
        temp2 += temp(i,j);
      }
    }
    for (int i = n1; i < (n - 1); i++){
      for (int j = (i + 1); j < n; j++){
        temp3 += temp(i,j);
      }
    }
    rz[r] = 2.0*c1*temp2/(double(n1)*(n1 - 1.0)) + 2.0*c2*temp3/(double(n2)*(n2 - 1.0)) - 2.0*temp1/(double(n)*(n - 1.0));
  }
  return rz;
}

// [[Rcpp::export]]
NumericVector JackknifeImp(NumericMatrix y1, NumericMatrix y2, bool RNA_Seq) {
  if (RNA_Seq && y1.ncol() > 1 && y1.nrow() > 1) {
    y1 = transpose(y1);
    y2 = transpose(y2);
  }  
  int n1 = y1.nrow();
  int n2 = y2.nrow();
  int p = y1.ncol();

  NumericMatrix jack(n1*n2, p);
  
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      NumericMatrix y1Temp(n1-1,p);
      NumericMatrix y2Temp(n2-1,p);
      for (int k = 0; k < i; k++) {
        y1Temp(k,_) = y1(k,_);
      }
      for (int k = (i + 1); k < n1; k++) {
        y1Temp(k - 1,_) = y1(k,_);
      }
      for (int k = 0; k < j; k++) {
        y2Temp(k,_) = y2(k,_);
      }
      for (int k = (j + 1); k < n2; k++) {
        y2Temp(k - 1,_) = y2(k,_);
      }
      jack(i*n2 + j,_) = Importance1(y1Temp, y2Temp, RNA_Seq);
    }
  }
  
  NumericVector rz(p);
  for (int i = 0; i < p; i++) {
    rz[i] = 0.0;
    for (int j = 0; j < (n1*n2); j++) {
      rz[i] += jack(j,i);
    }
    rz[i] = rz[i]/double(n1*n2);
  }
  return rz;
}
