#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cmath>
#include <numeric>
#include "dcheck.hpp"
#include "common.hpp"

using namespace std;

size_t effective_alphabet(std::istream& is, const size_t maxlength) {
	constexpr size_t counts_length = std::numeric_limits<uint8_t>::max()+1;
	uint8_t counts[counts_length];
	std::fill(counts, counts+counts_length, 0);
	size_t length = 0;
	while(!is.eof()) {
		const uint8_t read_char = is.get();
		if(is.eof()) break;
		++length;
		display_read_process(length, maxlength);
		DCHECK_LT(read_char, counts_length);
		counts[ read_char ] = 1;
		if(length == maxlength) { break; }
	}
	return std::accumulate(counts, counts+counts_length, 0ULL);
}


int main(const int argc, char** argv) {
    common_parameters parm("computes the effective alphabet size");
    parm.getopts(argc, argv, "vp:");

    if(kVerbose) { std::cout << "effective alphabet size of file " << parm.filename() << " is "; }
    std::cout << effective_alphabet(parm.stream(), parm.prefixlength()) << std::endl;

	return 0;
}
