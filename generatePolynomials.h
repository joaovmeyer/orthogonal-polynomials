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

double polynomialInnerProduct(const Vec& p1, const Vec& p2, const Vec& x) {
	double res = 0;

	for (size_t i = 0; i < x.size(); ++i) {
		res += evaluatePolynomial(p1, x[i]) * evaluatePolynomial(p2, x[i]);
	}

	return res;
}

// returns n polynomial equations such that p_i(x) * p_j(x) = 0 only if i != j, given that x is in given x list
Mat getOrthogonalPolynomials(int n, const Vec& x) {
	Mat polynomials(n);
	Vec inners(n);

	for (int i = 0; i < n; ++i) {
		Vec polynomial(n); polynomial[i] = 1; // make a i degree polynomial


		for (int j = 0; j < i; ++j) {
			double s = polynomialInnerProduct(polynomial, polynomials[j], x) / inners[j];

			for (int k = 0; k < n; ++k) {
				polynomial[k] -= polynomials[j][k] * s;
			}
		}

		inners[i] = polynomialInnerProduct(polynomial, polynomial, x);
		polynomials[i] = polynomial;
	}

	return polynomials;
}


#endif
