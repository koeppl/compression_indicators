#include <numeric>
#include <utility>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include "common.hpp"
#include "dcheck.hpp"

double second_entropy(std::istream& is, const size_t maxlength) {
    constexpr size_t counts_length = std::numeric_limits<uint8_t>::max()+1;
    constexpr size_t number_keys = std::numeric_limits<uint16_t>::max()+1;
    std::vector<size_t> counts[number_keys];
    for(size_t i = 0; i < number_keys; ++i) {
	counts[i].resize(counts_length, 0);
    }
    size_t length = 0;
    if(is.eof()) return 0;
    uint16_t key = is.get();
    ++length;
    if(is.eof()) return 0;
    key <<=8;
    key |= is.get();
    ++length;

    while(!is.eof()) {
	const uint8_t read_char = is.get();
	if(is.eof()) { break; }
	DCHECK_LT(read_char, counts_length);
	++counts[key][read_char];
	++length;
	display_read_process(length, maxlength);
	if(length == maxlength) { break; }
	key <<= 8;
	key |= read_char;
    }
    double ret = 0;
    for(size_t i = 0; i < number_keys; ++i) {
	ret += entropy_in_container(counts[i]);
    }
    return ret/length;
}


int main(const int argc, char** argv) {
    common_parameters parm("computes the second entropy");
    parm.getopts(argc, argv, "vp:");

    if(kVerbose) { std::cout << "second empirical entropy of file " << parm.filename() << " is "; }
    std::cout << second_entropy(parm.stream(), parm.prefixlength()) << std::endl;

	return 0;
}
