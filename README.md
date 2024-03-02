# Orthogonal Polynomials
code that generates n polynomials orthogonal to each other 

Given a set of m points x_1, ..., x_m, and the inner product of two functions f and g given by the formula <f, g> = f(x_1) * g(x_1) + ... + f(x_m) * g(x_m), we want to find n polynomials p_1, ..., p_n (degrees 0, ..., n - 1) such that <p_i, p_j> = 0 <=> i != j. For this we use the Gram-Schmidt process. Generating orthogonal polynomials can be usefull in a lot of situations, but my interest comes from it's use in orthogonal polynomial regression, where it mitigates multicollinearity issues and is more numerically stable than polynomial regression with non-orthogonal polynomials.
