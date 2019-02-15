#include <numeric>
#include <utility>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include "common.hpp"
#include "dcheck.hpp"

double first_entropy(std::istream& is, const size_t maxlength) {
    constexpr size_t counts_length = std::numeric_limits<uint8_t>::max()+1;
    std::vector<size_t> counts[counts_length];
    for(size_t i = 0; i < counts_length; ++i) {
	counts[i].resize(counts_length, 0);
    }
    size_t length = 0;
    if(is.eof()) return 0;
    uint8_t prev_char = is.get();
    ++length;
    while(!is.eof()) {
	const uint8_t read_char = is.get();
	if(is.eof()) { break; }
	DCHECK_LT(read_char, counts_length);
	++counts[prev_char][read_char];
	++length;
	display_read_process(length, maxlength);
	if(length == maxlength) { break; }
	prev_char = read_char;
    }
    // DCHECK_EQ(std::accumulate(counts,counts+counts_length,0ULL), length);
    double ret = 0;
    for(size_t i = 0; i < counts_length; ++i) {
	ret += entropy_in_container(counts[i]);
    }
    return ret/length;
}




int main(const int argc, char** argv) {
    common_parameters parm("computes the first entropy");
    parm.getopts(argc, argv, "vp:");

    if(kVerbose) { std::cout << "first empirical entropy of file " << parm.filename() << " is "; }
    std::cout << first_entropy(parm.stream(), parm.prefixlength()) << std::endl;

	return 0;
}
