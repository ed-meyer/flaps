//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef QR_H
#define QR_H

// minimum 2-norm solutions of underdetermined linear systems
// as described in Algorithm 5.6.2 in
//   \bibitem[Golub and Van~Loan(2013)]{golub2013matrix}
//   Golub, G.H.; Van~Loan, C.F.
//   \newblock {\em Matrix Computations}, 4th ed.; The Johns Hopkins University
//   Press: Baltimore,  MD, USA, 2013.

#include <string>
#include <vector>

#include "message.h"
#include "svd.h"

// QR factorizes the transpose of an (nr,nc) matrix where
// nc >= nr. The factorization is in the (nc,nr) matrix af
class QR {
public:
	size_t nr, nc;
	std::vector<double> af;			// nc by nr matrix of factors
	std::vector<int> pivots;
	std::vector<double> tau;		// tau array for QR factorization
	double rcond;
	int rankdef;
	SVD* svd;

	// constructor(s)
	QR() : nr(0), nc(0), rcond(0.0), rankdef(0), svd(nullptr) {}
	QR (size_t nrow, size_t ncol, double const* A); // factorize A^t
	~QR();

	void solve (std::vector<double> const& b, std::vector<double>& x);
	int performanceIndex(const double* A, const double* x, const double* b);
	// augJacDet: compute the determinant of the Jacobian augmented
	// with the tangent vector (see eqns 10,11 in ref. in bifurcation())
	// The "jac" arg is only if it is rank-deficient and an SVD was used
	double augJacDet(std::vector<double> const& jac,
			std::vector<double> const& tan, double& expon) const;
	std::vector<double> nullProj (std::vector<double> const& proj, double scale);

	// force svd
	static bool usesvd(int act=-1);
	// compute the determinant
	double determinant(double& expon) const;
};

// for checking augJacDet:
double
augeigdet(size_t nr, size_t nc, std::vector<double> const& J,
		std::vector<double> const& tan, double& expon);

// for determining rank deficiencies
double
rcondeps();


#endif	/* QR_H */
