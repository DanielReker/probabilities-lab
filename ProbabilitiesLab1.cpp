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
FloatType precision;
IntType experiments, N, M, n;


FloatType hypergeometricDistributionPMF(FloatType value, FloatType N, FloatType M, FloatType n) {
	return binomialCoefficient(M, value) * binomialCoefficient(N - M, n - value) / binomialCoefficient(N, n);
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
	matplot::title(std::format("{} experiments, N = {}, M = {}, n = {}", experiments, N, M, n));
	addBarLabels<IntType>(absoluteBar, values, counts, 6, false);

	matplot::show();
	matplot::save(std::format("img/hypergeometric_exps{}_N{}_M_{}_n{}_seed{}/absolute.pdf", experiments, N, M, n, seed));


	for (auto& count : counts) {
		count /= experiments;
	}

	matplot::hold(false);

	matplot::ylim({ 0, 1.3 * *std::max_element(counts.begin(), counts.end()) });
	matplot::bar(values, counts);
	matplot::title(std::format("{} experiments, N = {}, M = {}, n = {}", experiments, N, M, n));
	addBarLabels<double>(absoluteBar, values, counts, 7, true);


	matplot::show();
	matplot::save(std::format("img/hypergeometric_exps{}_N{}_M_{}_n{}_seed{}/relative.pdf", experiments, N, M, n, seed));
}


int main()
{
	std::ifstream config("config.txt");

	ostream << "Hypergeometric distribution\n";

	
	config >> N >> M >> n >> seed;
	if (seed == 0) seed = std::random_device{}();
	ostream << std::format("N = {}\nM = {}\nn = {}\nseed = {}\n", precision, N, M, n, seed);


	ostream << "number of experiments = ";
	std::cin >> experiments;


	IntType minm = std::max(static_cast<IntType>(0), M - N + n), maxm = std::min(M, n);
	auto distributionValuesProbabilities = std::views::iota(minm, maxm + 1) |
		std::views::transform([&](IntType val) { return std::pair(val, hypergeometricDistributionPMF(val, N, M, n)); });

	DiscreteDistributuion<IntType, FloatType> distribution(distributionValuesProbabilities, seed);



	ostream << '\n';
	distribution.printInfo();

	std::map<IntType, IntType> sample;
	for (IntType i = 0; i < experiments; i++) {
		sample[distribution.generateValue()]++;
	}

	printSampleInfo(sample);

	return 0;
}
