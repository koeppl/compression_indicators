#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>
#include "/scripts/code/dcheck.hpp"

using namespace std;

template<class T>
double double_div(const T& a, const T& b) {
	return static_cast<double>(a)/static_cast<double>(b);
}

double entropy_div(const std::vector<size_t>& counts) {
	const size_t length = std::accumulate(counts.begin(), counts.end(), 0);
	DCHECK_EQ(counts.size(), std::numeric_limits<uint8_t>::max()+1);
	double ret = 0;
	for(size_t i = 0; i < std::numeric_limits<uint8_t>::max()+1; ++i) {
		if(counts[i] == 0) continue;
		ret += counts[i] * std::log2(double_div(length, counts[i]));
	}

	return ret;
}

size_t effective_alphabet(std::istream& is) {
	constexpr size_t counts_length = std::numeric_limits<uint8_t>::max()+1;
	uint8_t counts[counts_length];
	std::fill(counts, counts+counts_length, 0);
	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		DCHECK_LT(read_char, counts_length);
		counts[ read_char ] = 1;
	}
	return std::accumulate(counts, counts+counts_length, 0ULL);
}


int main(const int argc, const char *const argv[]) {
	if(argc != 2) {
		std::cerr << "Usage: " << argv[0] << " filename\ncomputes the effective alphabet size of `filename`." << std::endl;
		return 1;
	}

	ifstream file(argv[1]);
	if(!file.good()) {
		std::cerr << "file " << argv[1] << " is not readable." << std::endl;
		return 1;
	}
	std::cout << effective_alphabet(file) << std::endl;

	return 0;
}
