#include <iostream>
#include <iomanip>

#include "generatePolynomials.h"

using namespace std;


int main() {

	Vec x = { -1, -0.5, 0, 0.5, 1 };
	Mat polynomials = getOrthogonalPolynomials(3, x);

	cout << polynomials << "\n\n";

	cout << std::fixed << std::setprecision(2);

	// checking for orthogonality (only diagonal elements should be non-zero)
	for (size_t i = 0; i < polynomials.size(); ++i) {
		for (size_t j = 0; j < polynomials.size(); ++j) {
			cout << polynomialInnerProduct(polynomials[i], polynomials[j], x) << ", ";
		}
		cout << "\n";
	}

	return 0;
}
