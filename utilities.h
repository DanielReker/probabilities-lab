#pragma once

#include <cmath>
#include <numbers>



//std::ofstream ostream("output.txt");
std::ostream& ostream = std::cout;


template<class Type>
Type factorial(Type n) {
	static std::vector<Type> factorials = { 1 };
	while (n >= factorials.size()) {
		factorials.push_back(factorials.back() * static_cast<Type>(factorials.size()));
	}
	return factorials[static_cast<size_t>(n)];
}


template<class Type>
Type binomialCoefficient(Type n, Type k) {
	return factorial(n) / (factorial(k) * factorial(n - k));
}



template<class Type>
inline Type phi(Type x) {
	return 0.5 * (1 + std::erf(x / std::numbers::sqrt2_v<Type>));
}


template<class Type>
inline Type normalCDF(Type x, Type mean, Type standardDeviation) {
	return phi<Type>((x - mean) / standardDeviation);
}

