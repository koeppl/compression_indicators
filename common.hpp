#pragma once

#include "dcheck.hpp"
#include <cstdint>
#include <iostream>

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

template<class container_t>
double entropy_in_container(const container_t& counts) {
	const size_t length = std::accumulate(counts.begin(), counts.end(), 0ULL);
	DCHECK_EQ(counts.size(), std::numeric_limits<uint8_t>::max()+1);
	double ret = 0;
	for(size_t i = 0; i < std::numeric_limits<uint8_t>::max()+1; ++i) {
		if(counts[i] == 0) continue;
		ret += counts[i] * std::log2(double_div(length, static_cast<size_t>(counts[i])));
	}

	return ret;
}


#if __GNUC__ <= 7
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif//__GNUC__

size_t filesize(const std::string& filename) {
	try {
	    return fs::file_size(filename); 
	} catch(fs::filesystem_error& e) {
	    std::cerr << "file " << filename << " is not readable." << std::endl;
	    std::cerr << e.what() << std::endl;
		throw e;
	}
}

#include <unistd.h>

struct common_parameters {
    size_t m_prefixlength = 0;
    std::string m_filename;
    const char*const m_program_description;
    std::istream* m_stream = nullptr;
    protected:

    virtual void fallback_option(const int, const int, char** argv) {
	printUsage(argv);
	abort();
    }
    
    public:
    size_t prefixlength() const { return m_prefixlength; }

    const char* filename() const {
	return m_filename.empty() ? "stdin" : m_filename.c_str();
    }
    std::istream& stream() {
	return *m_stream;
    }
    ~common_parameters() {
	if(m_stream != nullptr && m_stream != &std::cin) {
	    delete m_stream;
	    m_stream = nullptr;
	}
    }

    common_parameters(const char*const program_description)
	: m_program_description(program_description)
    {
    }


    virtual void printUsage(char** argv) {
	std::cerr << "Usage: " << argv[0] << " [-p prefix] [-v] [filename]\n";
	std::cerr << m_program_description;
	std::cerr << 
"\noptions: \n"
"-p prefix \t for reading only up to  `prefix` characters. A length of 0 means to read everything.\n"
"-v \t for verbose output\n"
"filename \t if no filename is given, the input is assumed to be streamed from stdin."
;
    }
    
    void getopts(const int argc, char** argv, const char*const optionstring) {
	int c;
	while((c = getopt (argc, argv, optionstring)) != -1) {
	    switch(c) {
		case 'v':
		    kVerbose = 1;
		    break;
		case 'p':
		    m_prefixlength = strtoul(optarg, NULL, 10);
		    break;
		default:
		    fallback_option(c, argc, argv);
		    break;
	    }
	}
	if(optind < argc) {
	    m_filename = argv[optind];
	}
	if(kVerbose) {
	    std::cout << "Parameters: \n - filename = " << (m_filename.empty() ? "[STDIN]" : m_filename) << "\n - prefix = " << m_prefixlength << std::endl;
	}

	if(m_filename.empty()) {
	    m_stream = &std::cin;
	} else {
	    if(m_prefixlength == 0) { m_prefixlength = filesize(m_filename); }

	    m_stream = new std::ifstream(m_filename);
	    if(!m_stream->good()) {
		std::cerr << "file " << m_filename << " is not readable." << std::endl;
		abort();
	    }
	}
    }

};

