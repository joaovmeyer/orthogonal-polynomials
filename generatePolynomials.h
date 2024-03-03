#ifndef GENERATE_POLYNOMIALS_H
#define GENERATE_POLYNOMIALS_H

#include <vector>


using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

std::ostream& operator << (std::ostream& os, const Vec& v) {
	for (size_t i = 0; i < v.size(); ++i) {
		os << v[i] << ", ";
	}
	os << "\n";

	return os;
}

std::ostream& operator << (std::ostream& os, const Mat& m) {
	for (size_t i = 0; i < m.size(); ++i) {
		os << m[i];
	}

	return os;
}


// uses Horner's method for evaluating polynomials
double evaluatePolynomial(const Vec& p, double x) {
	size_t n = p.size() - 1; // degree of polynomial

	double b = p[n];
	for (size_t i = n; i > 0; --i) {
		b = p[i - 1] + b * x;
	}

	return b;
}

// without pre-evaluation
double polynomialInnerProduct(const Vec& p1, const Vec& p2, const Vec& x) {
	double res = 0;

	for (size_t i = 0; i < x.size(); ++i) {
		res += evaluatePolynomial(p1, x[i]) * evaluatePolynomial(p2, x[i]);
	}

	return res;
}

// with pre-evaluated polynomials
double polynomialInnerProduct(const Vec& e1, const Vec& e2) {
	double res = 0;

	for (size_t i = 0; i < e1.size(); ++i) {
		res += e1[i] * e2[i];
	}

	return res;
}

// returns n polynomial equations such that p_i(x) * p_j(x) = 0 only if i != j, given that x is in given x list
Mat getOrthogonalPolynomials(int n, const Vec& x) {

	Mat polynomials(n, Vec(n));

	// precompute some values so we don't need to do it multiple times.
	// this is about 5x faster than if we wouldn't precompute evaluations
	// if we didn't even cache the inner products, than this version would
	// be about 8.5x faster!
	Vec inners(n);
	Mat evaluations(n, Vec(x.size()));

	for (int i = 0; i < n; ++i) {
		polynomials[i][i] = 1; // make a i degree polynomial

		// evaluate our new polynomial in all the x points
		for (size_t j = 0; j < x.size(); ++j) {
			evaluations[i][j] = evaluatePolynomial(polynomials[i], x[j]);
		}

		// accumulate the projections
		Vec projectionsSum = Vec(n);
		for (int j = 0; j < i; ++j) {
			double s = polynomialInnerProduct(evaluations[i], evaluations[j]) / inners[j];

			for (int k = 0; k < n; ++k) {
				projectionsSum[k] += polynomials[j][k] * s;
			}
		}

		// make i-th polynomial orthogonal to the rest
		for (int k = 0; k < n; ++k) {
			polynomials[i][k] -= projectionsSum[k];
		}

		// store everything we will need in the future
		for (size_t j = 0; j < x.size(); ++j) {
			evaluations[i][j] = evaluatePolynomial(polynomials[i], x[j]);
		}

		inners[i] = polynomialInnerProduct(evaluations[i], evaluations[i]);
	}

	return polynomials;
}


#endif
