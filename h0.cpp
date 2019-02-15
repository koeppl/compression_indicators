#include <numeric>
#include <utility>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include "common.hpp"
#include "dcheck.hpp"

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

int main(const int argc, char** argv) {
    common_parameters parm("computes the zeroth entropy");
    parm.getopts(argc, argv, "vp:");

    if(kVerbose) { std::cout << "zeroth empirical entropy of file " << parm.filename() << " is "; }
    std::cout << zeroth_entropy(parm.stream(), parm.prefixlength()) << std::endl;
    return 0;
}
