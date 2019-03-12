#include <iostream>
#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>

using namespace std;
using namespace Rcpp;
using namespace arma;

#define ARMA_DONT_PRINT_ERRORS

  //*********************//
 //  Generate K matrix  //
//*********************//
double GenKElem(vec v1, vec v2) {

	vec m = v1 % v2;
	double sum1 = sum(m);
	double sum2 = 0;
	for (int i = 0; i < (m.n_elem - 1); i++)
		for (int j = i + 1; j < m.n_elem; j++) 
			sum2 = sum2 + m[i] * m[j];
	double sum = 1 + sum1 + sum2;
	return sum;
}

mat GenK(mat var) {

	int n = var.n_rows;
	mat K = zeros<mat>(n, n);

	for (int i = 0; i < n; i++) 
		for (int j = i; j < n; j++) 
			K(i, j) = K(j, i) = GenKElem(var.row(i).t(), var.row(j).t());
	return K;
}

  //*************************//
 //  Generate K del matrix  //
//*************************//
mat GenKDel(mat var, vec testVar) {

	int n = var.n_rows;
	//KDel1
	mat m1 = zeros<mat>(n, n);
	for (int i = 0; i < n; i++) 
		for (int j = i; j < n; j++) 
			m1(i, j) = m1(j, i) = testVar[i] * testVar[j];
	//KDel2
	mat m2 = zeros<mat>(n, n);
	for (int i = 0; i < n; i++) 
		for (int j = i; j < n; j++)
			m2(i, j) = m2(j, i) = 1 + sum(var.row(i) % var.row(j));
	mat delMat = m1 % m2;
	return delMat;
}

  //***********************//
 //      Optimize tau     //
//***********************//
double tauOptim(double tau, mat K, vec DVec, vec residTilde, mat covar, int iter) {

	mat DInv = diagmat(1 / DVec);
	mat V = DInv + tau * K;
	mat VInv = inv(V);
	mat P = VInv - VInv * covar * inv(covar.t() * VInv * covar) * covar.t() * VInv;
	double delta = -as_scalar(residTilde.t() * VInv * K * VInv * residTilde - trace(P*K)) / 2;
	double info = trace(P * K * P * K) / 2;
	//pow(0.5, iter) *
	double tauOp = tau - pow(0.5, iter) * delta / info;

	return tauOp;
}

  //*************************//
 // Optimize alpha and beta //
//*************************//
vec coefOptim(double tau, mat K, mat D, vec yTilde, mat covar, int iter) {

	int numSample = covar.n_rows;
	int numCovar = covar.n_cols;
	// mat H = zeros<mat>(numSample + numCovar, numSample + numCovar);
	mat Ha = covar.t() * D * covar;
	mat Hb = covar.t() * D * K;
	mat Hc = D * covar;
	mat Hd = (1 / tau) * eye(numSample, numSample) + D * K;
	mat H = join_cols(join_rows(Ha, Hb), join_rows(Hc, Hd));

	vec deltaa = covar.t() * D * yTilde;
	vec deltab = D * yTilde;
	vec delta = join_cols(deltaa, deltab);

	vec coefOp = solve(HInv, delta);

	return coefOp;
}

  //******************************//
 // Estimate tau, alpha and beta //
//******************************//
int est(vec y, mat covar, mat K, double tauIn, vec betaIn, vec alphaIn,
	double &tauHat, vec &betaHat, vec &alphaHat) {

	int numSample = covar.n_rows;
	int numCovar = covar.n_cols;

	tauHat = tauIn;
	betaHat = betaIn;
	alphaHat = alphaIn;

	int niter = 20;
	double convCiter = 1e-3;
	double tauOld = tauHat;
	vec betaOld = betaHat;
	vec alphaOld = alphaHat;

	// Newton method
	for (int iter = 0; iter < niter; iter++) {

		vec eff = covar * betaHat + K * alphaHat;
		vec gamma = 1 / (1 + exp(-eff));
		vec DVec = gamma % (1 - gamma);
		mat D = diagmat(DVec);
		mat DInv = diagmat(1 / DVec);
		vec yTilde = eff + DInv * (y - gamma);
		vec residTilde = yTilde - covar * betaHat;
		// cout << gamma << endl;
		// cout << DVec << endl;
		// Optimize beta and alpha
		// cout << "iter " << iter << ": optimizing beta and alpha.";
		vec coef = coefOptim(tauHat, K, D, yTilde, covar, iter);
		betaHat = coef.subvec(0, numCovar - 1);
		alphaHat = coef.subvec(numCovar, numSample + numCovar-1);
		cout << "betaHat: " << betaHat << endl;
		// Optimize tau
		// cout << "iter " << iter << ": optimizing tau.";
		tauHat = tauOptim(tauHat, K, DVec, residTilde, covar, iter);
		cout << "tauHat: " << tauHat << endl;
		// if (tauHat < 0) { tauHat = tauIn; iter = 0; }
		// 
		double alphaDiff = dot(alphaOld - alphaHat, alphaOld - alphaHat);
		double betaDiff = dot(betaOld - betaHat, betaOld - betaHat);
		double tauDiff = (tauOld - tauHat)*(tauOld - tauHat);
		double citer = sqrt(alphaDiff) + sqrt(betaDiff) + sqrt(tauDiff);
		cout << "citer: " << citer << endl;
		if (citer < convCiter) {

			cout << "iter " << iter << ": algorithm converage.";
			break;
		}
		else {

			alphaOld = alphaHat;
			betaOld = betaHat;
			tauOld = tauHat;
		}

		// if (iter == niter) cout << "Not converage!" << endl;
	}
	return 0;
}

  //********************//
 //     Test delta     //
//********************//
// [[Rcpp::export]]
SEXP PEA(SEXP yIn, SEXP covarIn, SEXP varIn, SEXP testVarIn, SEXP permIn) {

	try {
		vec y = as<vec>(yIn);  // *dim = numSample x 1
		mat covar = as<mat>(covarIn);  // *dim = numSample x numCovar
		mat var = as<mat>(varIn);  // *dim = numSample x numVar
		vec testVar = as<vec>(testVarIn); // *dim = numSample x 1
		int perm = as<int>(permIn);

		  ////////////////////////////////
		 //   Generate kernel matrix   //
		////////////////////////////////
		cout << "Kernel matrix ... " << endl;
		mat K = GenK(var);
		mat KDel = GenKDel(var, testVar);

		  ///////////////////////////
		 // Parameters Estimation //
		///////////////////////////
		cout << "Estimating ..." << endl;
		int numSample = covar.n_rows;
		int numCovar = covar.n_cols;
		double tauInitial = 0.01;
		vec betaInitial = zeros<vec>(numCovar);
		vec alphaInitial = zeros<vec>(numSample);
		double tauHat;
		vec betaHat = zeros<vec>(numCovar);
		vec alphaHat = zeros<vec>(numSample);
		double tstart1 = clock();
		est(y, covar, K, tauInitial, betaInitial, alphaInitial, tauHat, betaHat, alphaHat);
		double time_est = (clock() - tstart1) / (double(CLOCKS_PER_SEC));
		// cout << "time_est: " << time_est << endl;

		////////////////////////////
		// Calculating Statistics //
		////////////////////////////
		cout << "Testing..." << endl;
		vec eff = covar * betaHat + K * alphaHat;
		vec gamma = 1 / (1 + exp(-eff));
		vec DVec = gamma % (1 - gamma);
		mat D = diagmat(DVec);
		mat DInv = diagmat(1 / DVec);
		vec yTilde = eff + DInv * (y - gamma);
		vec residTilde = yTilde - covar * betaHat;
		// cout << "residTilde: " << residTilde << endl;
		mat V = DInv + tauHat * K;
		mat VInv = inv(V);
		// cout << "VInv: " << VInv << endl;
		vec effDiff = exp(eff) - exp(-eff);
		vec KDel_alpha = KDel * alphaHat;
		mat DInvDel = diagmat((exp(eff) - exp(-eff)) % (KDel * alphaHat));
		mat VDel = DInvDel + tauHat * KDel;
		mat VTau = K;
		// cout << "KDel: " << KDel << endl;
		// Score statistics
		double U = as_scalar(residTilde.t() * VInv * VDel * VInv * residTilde) / 2;

		//permutation
		vec U_perm = zeros<vec>(perm);
		for (int i = 0; i < perm; i++) {

			vec testVar_perm = shuffle(testVar);
			mat KDel_perm = GenKDel(var, testVar_perm);
			mat DInvDel_perm = diagmat((exp(eff) - exp(-eff)) % (KDel_perm * alphaHat));
			mat VDel_perm = DInvDel_perm + tauHat * KDel_perm;
			U_perm[i] = as_scalar(residTilde.t() * VInv * VDel_perm * VInv * residTilde) / 2;
		}
		double P = 2 * min(sum(U_perm > U), perm - sum(U_perm > U)) / (double) perm;
		return List::create(Named("U") = U, Named("P") = P, Named("betaHat") = betaHat);
	}
	catch (std::exception &ex) {
		forward_exception_to_r(ex);
	}
	catch (...) {
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
}

  ///////////////////////////////////////
 ///  		CODE END HERE        ///
///////////////////////////////////////

