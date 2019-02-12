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


bool kVerbose = false;

inline void display_read_process(const size_t length, const size_t maxlength) {
    if(kVerbose && maxlength != 0 && (length-1)*100/maxlength < (length)*100/maxlength) { 
	std::cout << ((length)*100/maxlength) << "% of text read..." << std::endl; 
    }
}
template<class T>
inline void display_read_kmers(const size_t counter, const T& map) {
    if(kVerbose) {
	const size_t mapsize = map.size();
	if((counter-1)*100/mapsize < (counter)*100/mapsize) { 
	    std::cout << ((counter)*100/mapsize) << "% of k-mers processed..." << std::endl; 
	}
    }
}
template<class T>
inline void display_remain_kmers(const size_t counter, const T& map) {
    if(kVerbose) {
	const size_t mapsize = map.size()+counter;
	if((counter-1)*100/mapsize < (counter)*100/mapsize) { 
	    std::cout << ((counter)*100/mapsize) << "% of k-mers processed..." << std::endl; 
	}
    }
}

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
		display_read_process(length, maxlength);
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


double kth_entropy_naive(std::istream& is, const size_t num_k, const size_t maxlength) {
	if(num_k == 0) return zeroth_entropy(is, maxlength);
	size_t length = 0;
	std::string ringbuffer;
	std::unordered_map<std::string, std::vector<size_t>> dict;
	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		display_read_process(length, maxlength);
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
	size_t stat_counter = 0;
	for(const auto& element : dict) {
		ret += entropy_div(element.second);
		++stat_counter;
		display_read_kmers(stat_counter, dict);

	}
	return ret/length;
}

double kth_entropy(std::istream& is, const size_t num_k, const size_t maxlength) {
	using namespace separate_chaining;

	if(num_k == 0) return zeroth_entropy(is, maxlength);
	size_t length = 0;
	uint64_t buffer = 0;
	using bucket_type = chaining_bucket<avx2_key_bucket<uint8_t>, plain_key_bucket<tdc::uint_t<40>>, incremental_resize>;

	separate_chaining_map<varwidth_key_bucket, class_key_bucket<bucket_type>, xorshift_hash, incremental_resize> dict(num_k*8);
	//separate_chaining_map<avx2_key_bucket<size_t>, class_key_bucket<bucket_type>, hash_mapping_adapter<uint64_t, SplitMix>, incremental_resize> dict;

#ifndef NDEBUG
	std::unordered_map<size_t, std::vector<size_t>> dict_vector;
#endif

	while(!is.eof()) {
		const uint64_t key = buffer & (-1ULL >> (64- 8*num_k));
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		display_read_process(length, maxlength);
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
	size_t stat_counter = 0;
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
	    ++stat_counter;
	    display_read_kmers(stat_counter, dict);
	}
	DCHECK_LT(ret-0.002, ret2+0.002);
	DCHECK_GT(ret+0.002, ret2-0.002);

	return ret/length;
}

// template<class T>
// void kth_entropy_compact(std::istream& is, const size_t num_k, const size_t maxlength, T& dict, size_t& length, size_t& buffer) {
// 	while(!is.eof()) {
// 		const uint8_t read_char = is.get();
// 		if(is.eof()) break;
// 		++length;
// 		display_read_process(length, maxlength);
// 		buffer <<= 8;
// 		buffer |= read_char;
// 		const size_t composedkey = buffer & (-1ULL >> (64- 8*(num_k+1)));
//
// 		if( ++dict.find_or_insert(composedkey, 0).value_ref() == std::numeric_limits<typename T::value_type>::max() ) {
// 		    return;
// 		}
// 		
// 		if(length == maxlength) { break; }
// 	}
// }
// template<class A, class B>
// void copy_dict(const A& dictA, B& dictB) {
//
// }

double kth_entropy_compact(std::istream& is, const size_t num_k, const size_t maxlength) {
	using namespace separate_chaining;

	if(num_k == 0) return zeroth_entropy(is, maxlength);
	size_t length = 0;
	size_t buffer = 0;

	//separate_chaining_map<avx2_key_bucket<size_t>, plain_key_bucket<tdc::uint_t<40>>, hash_mapping_adapter<uint64_t, SplitMix>, incremental_resize> dict;
	// separate_chaining_map<varwidth_key_bucket, plain_key_bucket<uint8_t>, xorshift_hash, incremental_resize> dict8(64);
	// separate_chaining_map<varwidth_key_bucket, plain_key_bucket<uint16_t>, xorshift_hash, incremental_resize> dict16(64);
	// separate_chaining_map<varwidth_key_bucket, plain_key_bucket<tdc::uint_t<24>>, xorshift_hash, incremental_resize> dict24(64);
	// separate_chaining_map<varwidth_key_bucket, plain_key_bucket<uint32_t>, xorshift_hash, incremental_resize> dict32(64);
	separate_chaining_map<varwidth_key_bucket, plain_key_bucket<tdc::uint_t<40>>, xorshift_hash, incremental_resize> dict((num_k+1)*8);


	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) return 0;
		++length;
		buffer <<= 8;
		buffer |= read_char;
		if(length == num_k) {
		    break;
		}
	}

	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		display_read_process(length, maxlength);
		buffer <<= 8;
		buffer |= read_char;
		const size_t composedkey = buffer & (-1ULL >> (64- 8*(num_k+1)));

		++dict.find_or_insert(composedkey, 0).value_ref();
		if(length == maxlength) { break; }
	}

	double ret = 0;
	size_t stat_counter = 0;

	std::vector<size_t> counts(std::numeric_limits<uint8_t>::max()+1, 0);
	auto itbegin = dict.rbegin_nav();
	while(!dict.empty()) {
		DCHECK_EQ(itbegin.key(), itbegin.key() &  (-1ULL >> (64- 8*(num_k+1))));
		const size_t kmer_base = itbegin.key() & (-1ULL << 8);

		std::fill(counts.begin(), counts.end(), 0);
		for(size_t c = 0; c < std::numeric_limits<uint8_t>::max()+1; ++c) {
			const size_t kmer_key = kmer_base | c;
			auto it = dict.find(kmer_key);
			if(it == dict.cend()) { continue; }
			counts[c] = it.value();
			dict.erase(kmer_key);
			++stat_counter;
			display_remain_kmers(stat_counter, dict);
		}
		ret += entropy_div(counts);
		if(dict.bucket_size(itbegin.bucket()) == 0) {
		    dict.shrink_to_fit(itbegin.bucket());
		}
		--itbegin;
	}
		return ret/length;
}


#include <unistd.h>


#if __GNUC__ <= 7
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif//__GNUC__




void printUsage(char** argv) {
    std::cerr << "Usage: " << argv[0] << " -f filename -k k [-p prefix] [-v] [-m method]\ncomputes the `k`-th entropy of `filename`, optionally for the prefix with `prefix` characters with method m = {0,1,2}, optionally verbose (-v)." << std::endl;
}

int main(const int argc, char** argv) {

    size_t num_k = 0;
    size_t method = 0;
    size_t prefixlength = 0;
    std::string filename;
    int c;
    while((c = getopt (argc, argv, "k:vm:p:f:")) != -1) {
	switch(c) {
	    case 'k':
		num_k = strtoul(optarg, NULL, 10);
		break;
	    case 'v':
		kVerbose = 1;
		break;
	    case 'm':
		method = strtoul(optarg, NULL, 10);
		break;
	    case 'p':
		prefixlength = strtoul(optarg, NULL, 10);
		break;
	    case 'f':
		filename = optarg;
		break;
	    default:
		printUsage(argv);
		abort();
		break;
	}
    }
    if(kVerbose) {
	std::cout << "Parameters: \n - filname = " << filename << "\n - k = " << num_k << "\n - prefix = " << "\n - method = " << method << std::endl;
    }

    if(filename.empty()) {
	printUsage(argv);
	return 1;
    }
    if(prefixlength == 0) {
	try {
	    prefixlength = fs::file_size(filename); 
	} catch(fs::filesystem_error& e) {
	    std::cerr << "file " << filename << " is not readable." << std::endl;
	    std::cerr << e.what() << std::endl;
	    return 1;
	}
    }

	ifstream file(filename);
	if(!file.good()) {
	    std::cerr << "file " << filename << " is not readable." << std::endl;
	    return 1;
	}
	double entropy;
	switch(method) {
	    case 1:
		if(num_k > 8) {
		    std::cerr << "Method works only up to k = 8" << std::endl;
		    return 1;
		}
		if(kVerbose) { std::cout << "Running kth_entropy_compact..." << std::endl; }
		entropy = kth_entropy_compact(file, num_k, prefixlength);
		break;
	    case 2:
		if(kVerbose) { std::cout << "Running kth_entropy_naive..." << std::endl; }
		entropy = kth_entropy_naive(file, num_k, prefixlength);
		break;
	    default:
		if(num_k > 7) {
		    std::cerr << "Method works only up to k = 7" << std::endl;
		    return 1;
		}
		if(kVerbose) { std::cout << "Running kth_entropy..." << std::endl; }
		entropy = kth_entropy(file, num_k, prefixlength);
		break;

	}
	if(kVerbose) { std::cout << num_k << "-th empirical entropy of file " << filename << " is "; }
	std::cout << entropy << std::endl;

	return 0;
}
