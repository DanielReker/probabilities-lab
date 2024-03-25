#pragma once



template <class ValueType, class FloatType>
class DiscreteDistributuion {
private:
	std::map<ValueType, FloatType> m_probabilityFunction;
	std::vector<std::pair<FloatType, ValueType>> m_cumulativeProbabilityFunction;

	std::mt19937_64 m_randomGenerator;
	std::uniform_real_distribution<FloatType> m_uniformDistribution;



	static inline FloatType probabiliyConstant(FloatType probability) {
		return 0.5 * std::pow(10, std::ceil(std::log10(probability))) / std::sqrt(2 * probability * (1 - probability));
	}

	static inline FloatType atLeastSignsProbability(FloatType probability, int64_t experiments, int64_t signs) {
		return std::erf(probabiliyConstant(probability) * std::sqrt(experiments) / std::pow(static_cast<FloatType>(10), signs));
	}

	static inline FloatType signsProbability(FloatType probability, int64_t experiments, int64_t signs) {
		return atLeastSignsProbability(probability, experiments, signs) - atLeastSignsProbability(probability, experiments, signs + 1);
	}

	static inline FloatType signsExpectation(FloatType probability, int64_t experiments) {
		FloatType sum = 0;
		for (int64_t signs = 1; signs <= 100; signs++) {
			sum += atLeastSignsProbability(probability, experiments, signs);
		}
		return sum;
	}


public:
	template<class ProbabilityDensity>
	DiscreteDistributuion(const ProbabilityDensity& probabilityDensity, unsigned int seed,
		FloatType precision = std::numeric_limits<FloatType>::quiet_NaN()) :
		m_uniformDistribution{ 0, 1 }, m_randomGenerator{ seed } {

		FloatType currentSum = 0;
		for (const auto& [value, probability] : probabilityDensity) {
			currentSum += probability;
			//m_cumulativeProbabilityFunction[value] = currentSum;
			m_cumulativeProbabilityFunction.push_back({ currentSum, value });

			m_probabilityFunction[value] = probability;

			if (!std::isnan(precision) && std::abs(currentSum - 1) < precision) break;
		}
	}

	ValueType generateValue() {
		auto foundPairIt = std::upper_bound(m_cumulativeProbabilityFunction.begin(), m_cumulativeProbabilityFunction.end(),
			std::pair<FloatType, ValueType>({ m_uniformDistribution(m_randomGenerator), m_cumulativeProbabilityFunction.front().second }));

		return (*foundPairIt).second;
	}

	const auto& getPDF() { return m_probabilityFunction; }

	void printInfo() {
		static const int64_t width = 27;
		static const int64_t outputPrecision = 17;


		ostream << "Distribution information:\n";
		ostream << std::setw(width) << "Value" << std::setw(width) << "Probability" << std::setw(width) << "CDF" << '\n';
		for (const auto& [probability, value] : m_cumulativeProbabilityFunction) {
			ostream << std::setprecision(outputPrecision) <<
				std::fixed << std::setw(width) << value <<
				std::fixed << std::setw(width) << m_probabilityFunction[value] <<
				std::fixed << std::setw(width) << probability << '\n';
		}
	}





	void printPrecisionDistribution(int64_t experiments, int64_t maxPrecision) {
		static const int64_t width = 7;
		static const int64_t outputPrecision = 2;

		ostream << std::setw(width) << "Val" << std::setw(width + 10) << "Prob";
		for (int64_t precision = 0; precision <= maxPrecision; precision++) {
			ostream << std::setw(width) << std::format("{} dec", precision);
		}
		ostream << '\n';

		for (const auto& [value, exactProbability] : m_probabilityFunction) {
			//FloatType mean = experiments * exactProbability;
			//FloatType standardDeviation = std::sqrt(experiments * exactProbability * (1 - exactProbability));

			ostream << std::setw(width) << value << std::setw(width + 10) << std::setprecision(outputPrecision + 10) << exactProbability;

			for (int64_t precision = 0; precision <= maxPrecision; precision++) {
				//FloatType maxDelta = 0.5 * std::pow(10, std::ceil(std::log10(exactProbability)) - precision);

				FloatType precisionProbability = signsProbability(exactProbability, experiments, precision);
					//normalCDF<FloatType>(std::ceil((exactProbability + maxDelta) * experiments) - 0.5, mean, standardDeviation) -
					//normalCDF<FloatType>(std::floor((exactProbability - maxDelta) * experiments) + 0.5, mean, standardDeviation);

				ostream << std::setw(width) << std::setprecision(outputPrecision) << precisionProbability * 100;
			}
			ostream << std::setw(width) << std::setprecision(outputPrecision) << signsExpectation(exactProbability, experiments);
			ostream << '\n';
		}
	}
};