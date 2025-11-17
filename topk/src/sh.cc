/**
	SH: Substring HeavyKeeper
	Copyright (C) 2024 Roberto Grossi and Veronica Guerrini

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.#include <iostream>
**/
#include "fp/rk61.hpp"

#include <string>
#include <chrono>
#include <ctime>
#include <vector>
#include <cmath>
#include <limits>
#include <omp.h>
#include <sys/time.h>
#include <tuple>
#include "defs.h"
#include "krfp.h"
#include "utils.h"
#include "unordered_dense.h"

#include "heavykeeperplus.hpp"

using namespace std;

template<typename Out>
void output_top_k(unsigned char * sequence, vector<pair<int64_t,uint32_t>> &topK, Out& out) {
    fp::RabinKarp61 rk(257);

    for (auto &u : topK) {
        uint64_t fp = 0;
		for(uint32_t j=0; j<u.second;j++)
            fp = rk.push(fp, sequence[u.first+j]);

		out << fp << "," << u.second << "," << 0 << endl;
    }
}

int main(int argc, char* argv[])
{
	if(argc<3)
	{	
		cerr << "Usage: " << argv[0] << " text-file " << " K " << " [outfile] " << " [prefix] " << endl;
		return -1;
	}
	
	INT K = atoll(argv[2]);
	int64_t L = 1;
    
    INT max_n = (argc > 4) ? atoll(argv[4]) : std::numeric_limits<INT>::max();
	
	ifstream seq(argv[1], ios::in | ios::binary);
    	unsigned char * input_seq_char = NULL;
	unsigned char * sequence;
	INT n = 0;
    	char c;
    	while (n < max_n && seq.get(c))
    	{
       		if(n == 0 || n % ALLOC_SIZE)	input_seq_char = ( unsigned char * ) realloc ( input_seq_char,   ( n + ALLOC_SIZE ) * sizeof ( unsigned char ) );
        	input_seq_char[n] = c;
       		n++;
    	}
    	seq.close();
	
	if( n == 0 || n % ALLOC_SIZE)    input_seq_char = ( unsigned char * ) realloc ( input_seq_char,   ( n + ALLOC_SIZE ) * sizeof ( unsigned char ) );
	if( c == '\n')	input_seq_char[--n] = 255;	
	else		input_seq_char[n] = 255;	
	//Temporary
	sequence = input_seq_char;
    	sequence[++n]='\0';
	
	cerr << "Text is of length n = " << n - 1 << "." << endl;
	cerr << "K = " << K << endl;

	std::chrono::steady_clock::time_point  stream_begin = std::chrono::steady_clock::now();

	vector<pair<int64_t,uint32_t>> topK; 
	auto cms_hk_ncols =(256*K<268435456?256*K:268435456); // temporary: min between 256*K and 2^28
	HeavyKeeperPlus HeavySketch(cms_hk_ncols, K);
	
	HeavySketch.extractTopK(sequence, n, K, L, topK);

	std::chrono::steady_clock::time_point  stream_end = std::chrono::steady_clock::now();
	cerr<<"Top-K computation time: "<< std::chrono::duration_cast<std::chrono::milliseconds>(stream_end - stream_begin).count() << "[ms]." << std::endl;	
    
    if(argc > 3) {
        ofstream of(argv[3]);
        output_top_k(sequence,topK, of);
        of.close();
    } else {
        output_top_k(sequence,topK, cout);
    }
	
	free (sequence);
	return 0;
}
