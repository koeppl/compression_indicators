#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>
#include <utility>
#include "dcheck.hpp"
#include "separate_chaining_table.hpp"
#include "chaining_bucket.hpp"
#include <tudocomp/ds/uint_t.hpp>

using namespace std;

template<class T>
double double_div(const T& a, const T& b) {
	return static_cast<double>(a)/static_cast<double>(b);
}

double entropy_div(const std::vector<size_t>& counts) {
	const size_t length = std::accumulate(counts.begin(), counts.end(), 0ULL);
	DCHECK_EQ(counts.size(), std::numeric_limits<uint8_t>::max()+1);
	double ret = 0;
	for(size_t i = 0; i < std::numeric_limits<uint8_t>::max()+1; ++i) {
		if(counts[i] == 0) continue;
		ret += counts[i] * std::log2(double_div(length, static_cast<size_t>(counts[i])));
	}

	return ret;
}

template<class T>
double entropy_div_map(const T& counts) {
	const size_t length = std::accumulate(counts.cbegin(), counts.cend(), 0ULL, [] (const auto& a, const auto& b) { return a + b.second; });
	DCHECK_LE(counts.size(), std::numeric_limits<uint8_t>::max()+1);
	double ret = 0;

	for(const auto& count : counts) {
		if(static_cast<size_t>(count.second) == 0) continue;
		ret += static_cast<size_t>(count.second) * std::log2(double_div(length, static_cast<size_t>(count.second)));
	}

	return ret;
}

double zeroth_entropy(std::istream& is, const size_t maxlength) {
	constexpr size_t counts_length = std::numeric_limits<uint8_t>::max()+1;
	size_t counts[counts_length];
	std::fill(counts, counts+counts_length, 0);
	size_t length = 0;
	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) { break; }
		DCHECK_LT(read_char, counts_length);
		++counts[ read_char ];
		++length;
		if(length == maxlength) { break; }
	}
	DCHECK_EQ(std::accumulate(counts,counts+counts_length,0ULL), length);
	double ret = 0;
	for(size_t i = 0; i < counts_length; ++i) {
		if(counts[i] == 0 || length == counts[i]) continue;
		ret += counts[i] * std::log2(double_div(length, counts[i]));
	}

	return ret/length;
}

double first_entropy(std::istream& is, const size_t maxlength) {
	size_t length = 0;
	std::unordered_map<uint8_t, std::vector<size_t>> dict;
	uint8_t oldchar;
	uint8_t newchar = is.get();
	++length;
	while(!is.eof()) {
		if(is.eof()) break;
		++length;
		oldchar = newchar;
		newchar = is.get();
		auto it = dict.find(oldchar);
		if(it == dict.end()) {
			dict[oldchar] = std::vector<size_t>(std::numeric_limits<uint8_t>::max()+1, 0);
			it = dict.find(oldchar);
		}
		DCHECK(it != dict.end());
		++it->second[newchar];
		if(length == maxlength) { break; }
	}
	double ret = 0;
	for(const auto& element : dict) {
		ret += entropy_div(element.second);
	}
	return ret/length;
}

double kth_entropy_naive(std::istream& is, const size_t num_k, const size_t maxlength) {
	if(num_k == 0) return zeroth_entropy(is, maxlength);
	if(num_k == 1) return first_entropy(is, maxlength);
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
		if(length == maxlength) { break; }
	}
	double ret = 0;
	for(const auto& element : dict) {
		ret += entropy_div(element.second);
	}
	return ret/length;
}

double kth_entropy(std::istream& is, const size_t num_k, const size_t maxlength) {
	using namespace separate_chaining;

	if(num_k == 0) return zeroth_entropy(is, maxlength);
	// if(num_k == 1) return first_entropy(is, maxlength);
	size_t length = 0;
	uint64_t buffer = 0;
	using bucket_type = chaining_bucket<avx2_key_bucket<uint8_t>, plain_key_bucket<tdc::uint_t<40>>, incremental_resize>;

	separate_chaining_map<varwidth_key_bucket, class_key_bucket<bucket_type>, xorshift_hash, incremental_resize> dict(63);
	//separate_chaining_map<avx2_key_bucket<size_t>, class_key_bucket<bucket_type>, hash_mapping_adapter<uint64_t, SplitMix>, incremental_resize> dict;

#ifndef NDEBUG
	std::unordered_map<size_t, std::vector<size_t>> dict_vector;
#endif

	while(!is.eof()) {
		const uint64_t key = buffer & (-1ULL >> (64- 8*num_k));
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		buffer <<= 8;
		buffer |= read_char;
		if(length < num_k+1) {
			continue;
		}
#ifndef NDEBUG
		if(dict_vector.find(key) == dict_vector.end()) {
			dict_vector[key] =  std::vector<size_t>(std::numeric_limits<uint8_t>::max()+1, 0);
		}
		++dict_vector[key][read_char];
#endif

		auto& bucket = dict[key];
		++bucket.find_or_insert(read_char, 0).value_ref();
		DCHECK_EQ(static_cast<size_t>(bucket[read_char]), static_cast<size_t>(dict[key][read_char]));
		DCHECK_EQ(static_cast<size_t>(bucket[read_char]), static_cast<size_t>(dict_vector[key][read_char]));
		if(length == maxlength) { break; }
	}

#ifndef NDEBUG
	double ret2 = 0;
	for(const auto& element : dict_vector) {
		ret2 += entropy_div(element.second);
	}
#endif
	double ret = 0;
	for(auto it = dict.begin_nav(); it != dict.end_nav(); ++it) {
#ifndef NDEBUG
	    const auto& vec = dict_vector[it.key()];
	    for(size_t i = 0; i < vec.size(); ++i) {
		if(it.value().find(i) != it.value().cend())
		    DCHECK_EQ(vec[i], it.value()[i]);
	    }
#endif
	    const double ret0 = entropy_div_map(it.value());
	    DCHECK_LT(ret0-0.002, entropy_div(vec)+0.002);
	    DCHECK_GT(ret0+0.002, entropy_div(vec)-0.002);

	    ret += ret0;
	}
	DCHECK_LT(ret-0.002, ret2+0.002);
	DCHECK_GT(ret+0.002, ret2-0.002);

	return ret/length;
}

double kth_entropy_compact(std::istream& is, const size_t num_k, const size_t maxlength) {
	using namespace separate_chaining;

	if(num_k == 0) return zeroth_entropy(is, maxlength);
	// if(num_k == 1) return first_entropy(is, maxlength);
	size_t length = 0;
	size_t buffer = 0;

	separate_chaining_map<avx2_key_bucket<size_t>, plain_key_bucket<tdc::uint_t<40>>, hash_mapping_adapter<uint64_t, SplitMix>, incremental_resize> dict;


	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		buffer <<= 8;
		buffer |= read_char;
		if(length < num_k+1) {
			continue;
		}
		const size_t composedkey = buffer & (-1ULL >> (64- 8*(num_k+1)));


		++dict.find_or_insert(composedkey, 0).value_ref();
		if(length == maxlength) { break; }
	}

	double ret = 0;

	std::vector<size_t> counts(std::numeric_limits<uint8_t>::max()+1, 0);
	while(!dict.empty()) {
		const auto itbegin = dict.cbegin();
		DCHECK_EQ(itbegin.key(), itbegin.key() &  (-1ULL >> (64- 8*(num_k+1))));
		const size_t kmer_base = itbegin.key() & (-1ULL << 8);

		std::fill(counts.begin(), counts.end(), 0);
		for(size_t c = 0; c < std::numeric_limits<uint8_t>::max()+1; ++c) {
			const size_t kmer_key = kmer_base | c;
			auto it = dict.find(kmer_key);
			if(it == dict.cend()) { continue; }
			counts[c] = it.value();
			dict.erase(kmer_key);
		}
		ret += entropy_div(counts);
	}
		return ret/length;
}


int main(const int argc, const char *const argv[]) {
	if(argc < 3) {
		std::cerr << "Usage: " << argv[0] << " filename k [prefix]\ncomputes the `k`-th entropy of `filename`, optionally for the prefix with `prefix` characters." << std::endl;
		return 1;
	}
	const size_t length = argc >= 4 ? strtoul(argv[3], NULL, 10) : 0;
	const size_t num_k = strtoul(argv[2], NULL, 10);

	// {
	// 	const double required = std::pow(static_cast<double>(std::numeric_limits<uint8_t>::max()+1), num_k);
	// 	if(required > static_cast<double>(std::numeric_limits<size_t>::max())) {
	// 		std::cerr << "Chosen k = " << num_k << " too large!" << std::endl;
	// 		return 1;
	// 	}
	// 	void* mem = malloc(required);
	// 	if(mem == nullptr) {
	// 		std::cerr << "no memory available to compute " << num_k << "-th entropy!" << std::endl;
	// 		return 1;
	// 	}
	// 	free(mem);
	// }
    //

	{
		ifstream file(argv[1]);
		if(!file.good()) {
			std::cerr << "file " << argv[1] << " is not readable." << std::endl;
			return 1;
		}
		std::cout << kth_entropy(file, num_k, length) << std::endl;
	}
	{
		ifstream file(argv[1]);
		if(!file.good()) {
			std::cerr << "file " << argv[1] << " is not readable." << std::endl;
			return 1;
		}
		std::cout << kth_entropy_compact(file, num_k, length) << std::endl;
	}
	{
		ifstream file(argv[1]);
		if(!file.good()) {
			std::cerr << "file " << argv[1] << " is not readable." << std::endl;
			return 1;
		}
		std::cout << kth_entropy_naive(file, num_k, length) << std::endl;
	}


	return 0;
}
