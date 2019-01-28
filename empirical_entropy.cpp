#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>
#include "dcheck.hpp"

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

double zeroth_entropy(std::istream& is) {
	constexpr size_t counts_length = std::numeric_limits<uint8_t>::max()+1;
	size_t counts[counts_length];
	std::fill(counts, counts+counts_length, 0);
	size_t length = 0;
	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		DCHECK_LT(read_char, counts_length);
		++counts[ read_char ];
		++length;
	}
	DCHECK_EQ(std::accumulate(counts,counts+counts_length,0ULL), length);
	double ret = 0;
	for(size_t i = 0; i < counts_length; ++i) {
		if(counts[i] == 0 || length == counts[i]) continue;
		ret += counts[i] * std::log2(double_div(length, counts[i]));
	}

	return ret/length;
}

double kth_entropy(std::istream& is, const size_t num_k) {
	if(num_k == 0) return zeroth_entropy(is);
	size_t length = 0;
	std::string ringbuffer;
	std::unordered_map<std::string, std::vector<size_t>> dict;
	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		ringbuffer.push_back(read_char);
		if(ringbuffer.size() < num_k+1) {
			continue;
		}
		DCHECK_EQ(ringbuffer.size(), num_k+1);
		const uint8_t after_char = ringbuffer[num_k];
		const std::string key = ringbuffer.substr(0, num_k);
		auto it = dict.find(key);
		if(it == dict.end()) {
			dict[key] = std::vector<size_t>(std::numeric_limits<uint8_t>::max()+1, 0);
			it = dict.find(key);
		}
		++it->second[after_char];
		ringbuffer.erase(0,1); //remove front
	}
	double ret = 0;
	for(const auto& element : dict) {
		ret += entropy_div(element.second);
	}
	return ret/length;
}


int main(const int argc, const char *const argv[]) {
	if(argc != 3) {
		std::cerr << "Usage: " << argv[0] << " filename k\ncomputes the `k`-th entropy of `filename`." << std::endl;
		return 1;
	}
	const size_t num_k = strtoul(argv[2], NULL, 10);

	{
		const double required = std::pow(static_cast<double>(std::numeric_limits<uint8_t>::max()+1), num_k);
		if(required > static_cast<double>(std::numeric_limits<size_t>::max())) {
			std::cerr << "Chosen k = " << num_k << " too large!" << std::endl;
			return 1;
		}
		void* mem = malloc(required);
		if(mem == nullptr) {
			std::cerr << "no memory available to compute " << num_k << "-th entropy!" << std::endl;
			return 1;
		}
		free(mem);
	}


	ifstream file(argv[1]);
	if(!file.good()) {
		std::cerr << "file " << argv[1] << " is not readable." << std::endl;
		return 1;
	}
	std::cout << kth_entropy(file, num_k) << std::endl;

	return 0;
}
