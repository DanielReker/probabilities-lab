#include <iostream>
#include <ranges>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>
#include <random>
#include <functional>
#include <map>
#include <format>
#include <fstream>

#include <matplot/matplot.h>

#include "utilities.h"
#include "DiscreteDistribution.h"


using FloatType = long double;
using IntType = int64_t;


unsigned int seed;
FloatType precision, lambda;
IntType experiments;



FloatType poissonDistributionPMF(IntType value, FloatType lambda) {
	return std::pow(lambda, value) * std::exp(-lambda) / factorial<FloatType>(value);
}

FloatType hypergeometricDistributionPMF(IntType value, IntType N, IntType M, IntType n) {
	return binomialCoefficient(value, M) * binomialCoefficient(n - value, N - M) / binomialCoefficient(n, N);
}


template<class Type>
void addBarLabels(auto bar, const auto& x, const auto& y, int fontSize, bool isFloating) {
	double offsetY = *std::max_element(y.begin(), y.end()) * 0.04;

	std::vector<double> label_x;
	std::vector<double> label_y;
	std::vector<std::string> labels;
	for (size_t j = 0; j < x.size(); ++j) {
		label_x.emplace_back(bar->x_end_point(0, j) - 0.5);
		label_y.emplace_back(y[j] + offsetY);
		if(!isFloating) labels.emplace_back(matplot::num2str(static_cast<Type>(y[j])));
		else labels.emplace_back(std::format("{:03}", static_cast<Type>(y[j])));
	}
	matplot::hold(true);
	matplot::text(label_x, label_y, labels)->font_size(static_cast<float>(fontSize));
}


void printSampleInfo (const std::map<IntType, IntType>& sample) {
	std::vector<double> values, counts;

	static const int width = 25;
	ostream << "Generated values:\n";
	ostream << std::setw(width) << "Value" << std::setw(width) << "Absolute frequency" << std::setw(width) << "Relative frequency" << '\n';
	for (const auto& [value, count] : sample) {
		values.push_back(static_cast<double>(value));
		counts.push_back(static_cast<double>(count));

		ostream << std::setprecision(17) <<
			std::fixed << std::setw(width) << value <<
			std::scientific << std::setw(width) << count <<
			std::fixed << std::setw(width) << static_cast<double>(count) / experiments << '\n';
	}


	matplot::grid(true);

	matplot::ylim({ 0, 1.3 * *std::max_element(counts.begin(), counts.end()) });
	auto absoluteBar = matplot::bar(values, counts);
	matplot::title(std::format("{} experiments, lambda = {}", experiments, lambda, seed));
	addBarLabels<IntType>(absoluteBar, values, counts, 7, false);

	matplot::show();
	matplot::save(std::format("img/poisson_exps{}_lambda{}_seed{}/absolute.pdf", experiments, lambda, seed));


	for (auto& count : counts) {
		count /= experiments;
	}

	matplot::hold(false);

	matplot::ylim({ 0, 1.3 * *std::max_element(counts.begin(), counts.end()) });
	matplot::bar(values, counts);
	matplot::title(std::format("{} experiments, lambda = {}", experiments, lambda, seed));
	addBarLabels<double>(absoluteBar, values, counts, 7, true);


	matplot::show();
	matplot::save(std::format("img/poisson_exps{}_lambda{}_seed{}/relative.pdf", experiments, lambda, seed));
}


int main()
{
	std::ifstream config("config.txt");

	ostream << "Poisson distribution\n";

	
	config >> precision >> lambda >> seed;
	if (seed == 0) seed = std::random_device{}();
	ostream << std::format("precision = {}\nlambda = {}\nseed = {}\n", precision, lambda, seed);


	ostream << "number of experiments = ";
	std::cin >> experiments;


	auto poissonDistributionRange = std::views::iota(0) |
		std::views::transform([&](IntType val) { return std::pair(val, poissonDistributionPMF(val, lambda)); });

	DiscreteDistributuion<IntType, FloatType> poissonDist(poissonDistributionRange, seed, precision);

	ostream << '\n';
	poissonDist.printInfo();
	//poissonDist.printPrecisionDistribution(experiments, 5);


	std::map<IntType, IntType> histogram;


	for (IntType i = 0; i < experiments; i++) {
		histogram[poissonDist.generateValue()]++;
	}

	

	return 0;
}
