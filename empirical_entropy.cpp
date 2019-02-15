#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>
#include <utility>
#include "dcheck.hpp"

#include <tudocomp/ds/uint_t.hpp>

#include <separate/separate_chaining_table.hpp>
#include <separate/bucket.hpp>
#include <separate/bucket_table.hpp>

#include "common.hpp"

using namespace std;


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



double kth_entropy_naive(std::istream& is, const size_t num_k, const size_t maxlength) {
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
		ret += entropy_in_container(element.second);
		++stat_counter;
		display_read_kmers(stat_counter, dict);

	}
	return ret/length;
}

double kth_entropy(std::istream& is, const size_t num_k, const size_t maxlength) {
	using namespace separate_chaining;

	size_t length = 0;
	uint64_t buffer = 0;
	using bucket_type = bucket_table<avx2_bucket<uint8_t>, plain_bucket<tdc::uint_t<40>>, incremental_resize>;

	separate_chaining_map<varwidth_bucket, class_bucket<bucket_type>, xorshift_hash<uint64_t>, incremental_resize> dict(num_k*8);
	//separate_chaining_map<avx2_bucket<size_t>, class_bucket<bucket_type>, hash_mapping_adapter<uint64_t, SplitMix>, incremental_resize> dict;

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
		ret2 += entropy_in_container(element.second);
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
	    DCHECK_LT(ret0-0.002, entropy_in_container(vec)+0.002);
	    DCHECK_GT(ret0+0.002, entropy_in_container(vec)-0.002);

	    ret += ret0;
	    ++stat_counter;
	    display_read_kmers(stat_counter, dict);
	}
	DCHECK_LT(ret-0.002, ret2+0.002);
	DCHECK_GT(ret+0.002, ret2-0.002);

	return ret/length;
}

template<class T>
void read_dict(std::istream& is, const size_t num_k, const size_t maxlength, T& dict, size_t& length, size_t& buffer) {
    while(!is.eof()) {
	const uint8_t read_char = is.get();
	if(is.eof()) break;
	++length;
	display_read_process(length, maxlength);
	buffer <<= 8;
	buffer |= read_char;
	const size_t composedkey = buffer & (-1ULL >> (64- 8*(num_k+1)));

	if(++dict.find_or_insert(composedkey, 0).value_ref() == std::numeric_limits<typename T::value_type>::max() ) {
	    return;
	}

	if(length == maxlength) { return; }
    }
}

template<class A, class B>
void move_dict(A& dictA, B& dictB) {
    DCHECK_LE(sizeof(typename A::value_type), sizeof(typename B::value_type));
    DCHECK_LE(dictA.key_bit_width(), dictB.key_bit_width());
    DCHECK_LE(sizeof(typename A::key_type), sizeof(typename B::key_type));
    dictB.reserve(dictA.bucket_count());
    for(auto it = dictA.rbegin_nav(); it != dictA.rend_nav(); --it) {
	dictB[it.key()] = static_cast<size_t>(it.value());
	dictA.erase(it);
	// if(it.position() == 0) {
	//     dictA.shrink_to_fit(it.bucket());
	// }
    }
    DCHECK(dictA.empty());
}

template<class dict_t>
double compute_entropy(dict_t& dict, [[maybe_unused]] const size_t num_k, const size_t length) {

	if(dict.empty()) return 0;
	double ret = 0;
	size_t stat_counter = 0;

	std::vector<size_t> counts(std::numeric_limits<uint8_t>::max()+1, 0);
	auto itbegin = dict.rbegin_nav();
	while(!dict.empty()) {
		DCHECK_EQ(itbegin.key(), itbegin.key() & (-1ULL >> (64- 8*(num_k+1))));
		const size_t kmer_base = itbegin.key() & (-1ULL << 8);

		std::fill(counts.begin(), counts.end(), 0);
		for(size_t c = 0; c < std::numeric_limits<uint8_t>::max()+1; ++c) {
			const size_t kmer_key = kmer_base | c;
			auto it = dict.find(kmer_key);
			if(it == dict.cend()) { continue; }
			counts[c] = it.value();
			dict.erase(it);
			//dict.erase(kmer_key);
			++stat_counter;
			display_remain_kmers(stat_counter, dict);
		}
		ret += entropy_in_container(counts);
		if(dict.bucket_size(itbegin.bucket()) == 0) {
		    dict.shrink_to_fit(itbegin.bucket());
		}
		--itbegin;
	}
	return ret/length;
}

namespace value_switch {
    using namespace separate_chaining;
    template<typename value_t> using map_type = separate_chaining_map<varwidth_bucket, plain_bucket<value_t>, xorshift_hash<uint64_t>, incremental_resize>;

    double kth_entropy_switch(std::istream& is, const size_t num_k, const size_t maxlength) {

	size_t length = 0;
	size_t buffer = 0;

	map_type<uint8_t> dict8((num_k+1)*8);

	while(!is.eof()) {
	    const uint8_t read_char = is.get();
	    if(is.eof()) return 0;
	    ++length;
	    buffer <<= 8;
	    buffer |= read_char;
	    if(length == num_k) {
		break;
	    }
	    if(length == maxlength) { return 0; }
	}

	read_dict(is, num_k, maxlength, dict8, length, buffer);
	if(length == maxlength || is.eof()) { return compute_entropy(dict8, num_k, length); }

	map_type<uint16_t> dict16((num_k+1)*8);
	if(kVerbose) {
	    std::cout << "Copy to dictionary with 16-bit values..."  << std::endl; 
	}
	move_dict(dict8, dict16);
	read_dict(is, num_k, maxlength, dict16, length, buffer);
	if(length == maxlength  || is.eof()) { return compute_entropy(dict16, num_k, length); }

	map_type<tdc::uint_t<24>> dict24((num_k+1)*8);
	if(kVerbose) {
	    std::cout << "Copy to dictionary with 24-bit values..."  << std::endl; 
	}
	move_dict(dict16, dict24);
	read_dict(is, num_k, maxlength, dict24, length, buffer);
	if(length == maxlength  || is.eof()) { return compute_entropy(dict24, num_k, length); }

	map_type<uint32_t> dict32((num_k+1)*8);
	if(kVerbose) {
	    std::cout << "Copy to dictionary with 32-bit values..."  << std::endl; 
	}
	move_dict(dict24, dict32);
	read_dict(is, num_k, maxlength, dict32, length, buffer);
	if(length == maxlength  || is.eof()) { return compute_entropy(dict32, num_k, length); }

	map_type<tdc::uint_t<40>> dict40((num_k+1)*8);
	if(kVerbose) {
	    std::cout << "Copy to dictionary with 40-bit values..."  << std::endl; 
	}
	move_dict(dict32, dict40);
	read_dict(is, num_k, maxlength, dict40, length, buffer);
	DCHECK(length == maxlength  || is.eof());
	return compute_entropy(dict40, num_k, length);
    }
}//ns value_switch

namespace keyvalue_switch {

    template<class storage_width> class nextkeywidth { };

    template<> struct nextkeywidth<uint64_t> { using storage_type = uint32_t; };
    template<> struct nextkeywidth<uint32_t> { using storage_type = uint16_t; };
    template<> struct nextkeywidth<uint16_t> { using storage_type = uint8_t; };
    template<> struct nextkeywidth<uint8_t> { using storage_type = uint8_t; }; // cannot recurse further

    template<class dict_t>
	bool next_storage_keywidth_available(const size_t num_k, const dict_t& dict) {
	    if(sizeof(typename dict_t::storage_type) == 8 && ((num_k+1)*8) <= dict.bucket_count_log2() + 32) { return true; }
	    if(sizeof(typename dict_t::storage_type) == 4 && ((num_k+1)*8) <= dict.bucket_count_log2() + 32 - 16) { return true; }
	    if(sizeof(typename dict_t::storage_type) == 2 && ((num_k+1)*8) <= dict.bucket_count_log2() + 32 - 16 -  8) { return true; }
	    return false;
	}

    template<class T>
	void read_dictfast(std::istream& is, const size_t num_k, const size_t maxlength, T& dict, size_t& length, size_t& buffer) {
	    while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		display_read_process(length, maxlength);
		buffer <<= 8;
		buffer |= read_char;
		const size_t composedkey = buffer & (-1ULL >> (64- 8*(num_k+1)));

		if(++dict.find_or_insert(composedkey, 0).value_ref() == std::numeric_limits<typename T::value_type>::max() ) { return; }
		if(next_storage_keywidth_available(num_k, dict)) { return; }

		if(length == maxlength) { return; }
	    }
	}

    using namespace separate_chaining;
    template<class storage_type, class value_type> 
	using map_type = separate_chaining::separate_chaining_map<avx2_bucket<storage_type>, plain_bucket<value_type>, xorshift_hash<uint64_t, storage_type>, incremental_resize>;



    template<class T>
	double kth_entropy_fastswitch(std::istream& is, const size_t num_k, const size_t maxlength, T& dict, size_t& length, size_t& buffer) {
	    using storage_type = typename T::storage_type;

	    read_dictfast(is, num_k, maxlength, dict, length, buffer);
	    if(length == maxlength || is.eof()) { return compute_entropy(dict, num_k, length); }
	    if(next_storage_keywidth_available(num_k, dict)) {
		if(kVerbose) {
		    std::cout << "Copy to dictionary with " << (sizeof(typename nextkeywidth<storage_type>::storage_type)*8) << " bit storage keys..."  << std::endl; 
		}
		map_type<typename nextkeywidth<storage_type>::storage_type, typename std::remove_reference<decltype(dict)>::type::value_type> newdict((num_k+1)*8);
		move_dict(dict, newdict);
		return kth_entropy_fastswitch(is, num_k, maxlength, newdict, length, buffer);
	    }


	    map_type<storage_type, uint16_t> dict16((num_k+1)*8);
	    if(sizeof(typename T::value_type) < 2) {
		if(kVerbose) {
		    std::cout << "Copy to dictionary with 16-bit values..."  << std::endl; 
		}
		move_dict(dict, dict16);
		read_dict(is, num_k, maxlength, dict16, length, buffer);
		if(length == maxlength  || is.eof()) { return compute_entropy(dict16, num_k, length); }
		if(next_storage_keywidth_available(num_k, dict16)) {
		    if(kVerbose) {
			std::cout << "Copy to dictionary with " << (sizeof(typename nextkeywidth<storage_type>::storage_type)*8) << " bit storage keys..."  << std::endl; 
		    }
		    map_type<typename nextkeywidth<storage_type>::storage_type, typename decltype(dict16)::value_type> newdict((num_k+1)*8);
		    move_dict(dict16, newdict);
		    return kth_entropy_fastswitch(is, num_k, maxlength, newdict, length, buffer);
		}
	    }


	    map_type<storage_type, tdc::uint_t<24>> dict24((num_k+1)*8);
	    if(sizeof(typename T::value_type) < 3) {
		if(kVerbose) {
		    std::cout << "Copy to dictionary with 24-bit values..."  << std::endl; 
		}
		if(sizeof(typename T::value_type) == 2) { 
		    DCHECK(dict16.empty());
		    move_dict(dict, dict24);
		} else {
		    DCHECK(dict.empty());
		    move_dict(dict16, dict24);
		}
		read_dict(is, num_k, maxlength, dict24, length, buffer);
		if(length == maxlength  || is.eof()) { return compute_entropy(dict24, num_k, length); }
		if(next_storage_keywidth_available(num_k, dict24)) {
		    if(kVerbose) {
			std::cout << "Copy to dictionary with " << (sizeof(typename nextkeywidth<storage_type>::storage_type)*8) << " bit storage keys..."  << std::endl; 
		    }
		    map_type<typename nextkeywidth<storage_type>::storage_type, typename decltype(dict24)::value_type> newdict((num_k+1)*8);
		    move_dict(dict24, newdict);
		    return kth_entropy_fastswitch(is, num_k, maxlength, newdict, length, buffer);
		}
	    }

	    map_type<storage_type, uint32_t> dict32((num_k+1)*8);
	    if(sizeof(typename T::value_type) < 4) {
		if(kVerbose) {
		    std::cout << "Copy to dictionary with 32-bit values..."  << std::endl; 
		}
		if(sizeof(typename T::value_type) == 3) { 
		    DCHECK(dict24.empty());
		    move_dict(dict, dict32);
		} else {
		    DCHECK(dict.empty());
		    move_dict(dict24, dict32);
		}
		read_dict(is, num_k, maxlength, dict32, length, buffer);
		if(length == maxlength  || is.eof()) { return compute_entropy(dict32, num_k, length); }
		if(next_storage_keywidth_available(num_k, dict32)) {
		    if(kVerbose) {
			std::cout << "Copy to dictionary with " << (sizeof(typename nextkeywidth<storage_type>::storage_type)*8) << " bit storage keys..."  << std::endl; 
		    }
		    map_type<typename nextkeywidth<storage_type>::storage_type, typename decltype(dict32)::value_type> newdict((num_k+1)*8);
		    move_dict(dict32, newdict);
		    return kth_entropy_fastswitch(is, num_k, maxlength, newdict, length, buffer);
		}
	    }

	    map_type<storage_type, tdc::uint_t<40>> dict40((num_k+1)*8);
	    if(kVerbose) {
		std::cout << "Copy to dictionary with 40-bit values..."  << std::endl; 
	    }
	    if(sizeof(typename T::value_type) == 4) { 
		DCHECK(dict32.empty());
		move_dict(dict, dict40);
	    } else {
		DCHECK(dict.empty());
		move_dict(dict32, dict40);
	    }
	    read_dict(is, num_k, maxlength, dict40, length, buffer);
	    if(next_storage_keywidth_available(num_k, dict40)) {
		if(kVerbose) {
		    std::cout << "Copy to dictionary with " << (sizeof(typename nextkeywidth<storage_type>::storage_type)*8) << " bit storage keys..."  << std::endl; 
		}
		map_type<typename nextkeywidth<storage_type>::storage_type, typename decltype(dict40)::value_type> newdict((num_k+1)*8);
		move_dict(dict40, newdict);
		return kth_entropy_fastswitch(is, num_k, maxlength, newdict, length, buffer);
	    }


	    DCHECK(length == maxlength  || is.eof());
	    return compute_entropy(dict40, num_k, length);
	}

    double kth_entropy_fastswitch(std::istream& is, const size_t num_k, const size_t maxlength) {
	using namespace separate_chaining;

	size_t length = 0;
	size_t buffer = 0;

	separate_chaining_map<avx2_bucket<uint64_t>, plain_bucket<uint8_t>, xorshift_hash<uint64_t>, incremental_resize> dict8((num_k+1)*8);

	while(!is.eof()) {
	    const uint8_t read_char = is.get();
	    if(is.eof()) return 0;
	    ++length;
	    buffer <<= 8;
	    buffer |= read_char;
	    if(length == num_k) {
		break;
	    }
	    if(length == maxlength) { return 0; }
	}
	return kth_entropy_fastswitch(is, num_k, maxlength, dict8, length, buffer);

    }
}//ns keyvalue

template<class bucket_type>
double kth_entropy_compact(std::istream& is, const size_t num_k, const size_t maxlength) {
	using namespace separate_chaining;

	size_t length = 0;
	size_t buffer = 0;

	separate_chaining_map<bucket_type, plain_bucket<tdc::uint_t<40>>, xorshift_hash<uint64_t>, incremental_resize> dict((num_k+1)*8);

	//separate_chaining_map<avx2_bucket<uint64_t>, plain_bucket<uint64_t>, xorshift_hash, incremental_resize> dict((num_k+1)*8);
	//separate_chaining_map<avx2_bucket<uint64_t>, plain_bucket<uint64_t>, hash_mapping_adapter<uint64_t, SplitMix>, incremental_resize> dict(64);

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

	read_dict(is, num_k, maxlength, dict, length, buffer);
	DCHECK(length == maxlength  || is.eof());

	// while(!is.eof()) {
	// 	const uint8_t read_char = is.get();
	// 	if(is.eof()) break;
	// 	++length;
	// 	display_read_process(length, maxlength);
	// 	buffer <<= 8;
	// 	buffer |= read_char;
	// 	const size_t composedkey = buffer & (-1ULL >> (64- 8*(num_k+1)));
        //
	// 	++dict.find_or_insert(composedkey, 0).value_ref();
	// 	if(length == maxlength) { break; }
	// }
	return compute_entropy(dict, num_k, length);
}



using function_type = double (*)(std::istream& is, const size_t num_k, const size_t maxlength);
constexpr function_type func[] { kth_entropy, kth_entropy_naive, kth_entropy_compact<separate_chaining::varwidth_bucket>, kth_entropy_compact<separate_chaining::avx2_bucket<uint64_t>>, value_switch::kth_entropy_switch, keyvalue_switch::kth_entropy_fastswitch };
const char*const func_label[] = { "default", "naive", "compact", "fast", "switch", "fastswitch" };
constexpr size_t func_upper_k[] = { 8, 64, 7, 7, 7, 7};
constexpr size_t func_length = sizeof(func)/sizeof(function_type);
static_assert(func_length == sizeof(func_label)/sizeof(char*const));
static_assert(func_length == sizeof(func_upper_k)/sizeof(size_t));

struct extended_parameters : public common_parameters {
    using super_class = common_parameters;
    size_t m_k = 2;
    size_t m_method = 0;


    protected:
    virtual void fallback_option(const int c, const int, char** argv) override {
	switch(c) {
	    case 'k':
		m_k = strtoul(optarg, NULL, 10);
		break;
	    case 'm':
		m_method = strtoul(optarg, NULL, 10);
		if(m_method >= func_length) {
		    std::cerr << "method must be an integer between 0 and " << (func_length-1) << std::endl;
		    abort();
		}
		break;
	    default:
		printUsage(argv);
		abort();
		break;
	}
    }


    public:
    extended_parameters(const char*const program_description) : super_class(program_description) {}

    size_t num_k() const { return m_k; }
    size_t method() const { return m_method; }

    void printUsage(char** argv) {
	super_class::printUsage(argv);
	std::cerr <<
	    "-k k\t compute the k-th entropy, k >= 1 [default: 2]"
	    "-m method \t use one of the following methods [default: 0]:"
	    ;
	for(size_t i = 0; i < func_length; ++i) {
	    std::cerr << "\t" <<  i  << "\t" << func_label[i] << std::endl;
	}
    }
};


int main(const int argc, char** argv) {
    extended_parameters parm("computes the k-th entropy");
    parm.getopts(argc, argv, "k:vm:p:f:");
    if(kVerbose) {
	std::cout << "k : " << parm.num_k() << std::endl;
	std::cout << "method: " << parm.method() << " (" << func_label[parm.method()] << ")" << std::endl;
    }
    if(parm.num_k() < 1) {
	std::cerr << "program works only for 1 <= k." << std::endl;
	return 1;
    }


    double entropy;

    if(parm.num_k() > func_upper_k[parm.method()]) {
	std::cerr << "Method works only for k <= "<< func_upper_k[parm.method()] << std::endl;
	return 1;
    }
    if(kVerbose) { std::cout << "Running " << func_label[parm.method()] << "..." << std::endl; }
    entropy = func[parm.method()](parm.stream(), parm.num_k(), parm.prefixlength());

    if(kVerbose) { std::cout << parm.num_k() << "-th empirical entropy of file " << parm.filename() << " is "; }
    std::cout << entropy << std::endl;

    return 0;
}
