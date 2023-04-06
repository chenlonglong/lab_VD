#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <tuple>
#include <array>
#include <queue>
#include <stdlib.h>
#include <time.h>

using namespace std;

int main(int argc, char** argv)
{
	uint32_t 	read_len 	= 64;
    uint32_t 	hap_len 	= 128;
	uint32_t 	data_num	= 100;

    for(int i=0;i<argc;i++) {
        std::string argv_str(argv[i]);

		// Turn the argv into uint32_t (stoul) and int (stoi)		
		if(argv_str == "-read_len") read_len = std::stoul(argv[i+1]);
		if(argv_str == "-hap_len") 	hap_len = std::stoul(argv[i+1]);
		if(argv_str == "-data_num") data_num = std::stoul(argv[i+1]);
    }
	
	std::cerr << "\n";
	std::cerr << "========= Start Generating dataset =========" << "\n";
    std::cerr << "Read Length: " << read_len << "\n";
	std::cerr << "Hap Length: " << hap_len << "\n";
	std::cerr << "Data num (set of hap-read): " << data_num << "\n";
	std::cerr << "========================================" << "\n\n";

	if (read_len <= 20) {
		std::cerr << "Read length should be > 20! \n";
		return 1;
	}
    if (read_len >= hap_len) {
        std::cerr << "Read length should be < Haplotype length\n";
            return 1;
    }
	if ((read_len % 4) != 0) {
		std::cerr << "Read length % 4 should be 0! \n";
		return 1;
	}
	if ((hap_len % 4) != 0) {
		std::cerr << "Haplotype length % 4 should be 0! \n";
		return 1;
	}

	std::ofstream fout_data("./test_data/random_pattern.txt");
	std::ofstream fout_bin("./test_data/random_pattern.bin", ios::out | ios::binary);
	if(!fout_bin) {
      cout << "Cannot open binary file!" << endl;
      return 1;
   	}

	srand (2);		//random seed
	int random_num;	// random number
	std::string haplotype_seq, read_seq, read_BQ;

	char* hap_buffer = new char[32];
	char* read_buffer = new char[32];
	char* read_BQ_buffer = new char[32];
	int temp_byte = 0;
	
	for (int i = 0; i < data_num; i++){
        // generate haplotype
        haplotype_seq = "";
        read_seq = "";
        for(int k = 0; k < hap_len; k++) {
            random_num = rand();

            if      ((random_num % 4) == 0) haplotype_seq.append("A");
            else if ((random_num % 4) == 1)	haplotype_seq.append("C");
            else if ((random_num % 4) == 2)	haplotype_seq.append("G");
            else if ((random_num % 4) == 3)	haplotype_seq.append("T");
        }
        
        // generate read by copying haplotype substring
        uint32_t tmp_pos = rand()%(hap_len - read_len);
        read_seq = haplotype_seq.substr(tmp_pos, read_len);

		if (i < data_num/5){ // exact match
            // do nothing
        }
		else if (i < data_num*2/5){ // 1~10 mismatch
			for (int m = 0; m < i; m++){
				int offset = rand()%read_len;	// position in read seq
                if (rand()%2 == 0) {
                    if (haplotype_seq[tmp_pos + offset] == 'A') read_seq.replace(offset, 1, "C");
                    else								     read_seq.replace(offset, 1, "A");
                }
			}
		}
		else if (i < data_num*3/5){	// insertion 3-7 bp
			int offset = rand()%read_len;	// position in read seq
            int insertion_len = rand()%5 + 3;

			for (int m = 0; m < insertion_len; m++){
                int base_seed = rand()%4;
                if      (base_seed==0) read_seq.insert(offset + m, "A");
                else if (base_seed==1) read_seq.insert(offset + m, "C");
                else if (base_seed==2) read_seq.insert(offset + m, "G");
                else if (base_seed==3) read_seq.insert(offset + m, "T");
			}
            read_seq = read_seq.substr(0, read_len);
		}
        else if (i < data_num*4/5) {// deletion: 3-7 bp, erase substring of read
            int deletion_len = rand()%5 + 3;
			int offset = rand()%(read_len - deletion_len);	// position in read seq

            for (int m = 0; m < deletion_len; m++) {
                read_seq.erase(offset, 1);
            }
			for (int m = 0; m < deletion_len; m++){
                int base_seed = rand()%4;
                if      (base_seed==0) read_seq.append("A");
                else if (base_seed==1) read_seq.append("C");
                else if (base_seed==2) read_seq.append("G");
                else if (base_seed==3) read_seq.append("T");
			}
        }
		else {	// generate random reads
            read_seq.clear();

			random_num = rand();

			for(int k = 0; k < read_len; k++) {
				random_num = rand()%4;

				if      (random_num == 0) read_seq.append("A");
				else if (random_num == 1) read_seq.append("C");
				else if (random_num == 2) read_seq.append("G");
				else if (random_num == 3) read_seq.append("T");
			}
		}

        for (int k = 0; k < read_len; k++) {
            uint32_t BQ_seed = rand()%100;
            if      (BQ_seed < 90) read_BQ.append("3");
            else if (BQ_seed < 95) read_BQ.append("2");
            else if (BQ_seed < 97) read_BQ.append("1");
            else                   read_BQ.append("0");
        }

		// write binary file (Hap)
		for (int k = 0; k < hap_len; k++){
			if (haplotype_seq[k] == 'A')		temp_byte += 0;
			else if (haplotype_seq[k] == 'C')	temp_byte += 1;
			else if (haplotype_seq[k] == 'G')	temp_byte += 2;
			else if (haplotype_seq[k] == 'T')	temp_byte += 3;

			if (k % 4 == 3) {
				hap_buffer[k/4] = temp_byte;
				temp_byte = 0;
			}
			else temp_byte = temp_byte * 4;
		}
		fout_bin.write (hap_buffer, (hap_len/4));

		// write binary file (Read)
		temp_byte = 0;
		for (int k=0; k<read_len; k++){
			if (read_seq[k] == 'A')			temp_byte += 0;
			else if (read_seq[k] == 'C')	temp_byte += 1;
			else if (read_seq[k] == 'G')	temp_byte += 2;
			else if (read_seq[k] == 'T')	temp_byte += 3;

			if (k % 4 == 3) {
				read_buffer[k/4] = temp_byte;
				temp_byte = 0;
			}
			else temp_byte = temp_byte * 4;
		}
		fout_bin.write (read_buffer, (read_len/4));
		
		// write binary file (Read BQ)
		temp_byte = 0;
		for (int k=0; k<read_len; k++){
			if      (read_BQ[k] == '0') temp_byte += 0;
			else if (read_BQ[k] == '1')	temp_byte += 1;
			else if (read_BQ[k] == '2')	temp_byte += 2;
			else if (read_BQ[k] == '3')	temp_byte += 3;

			if (k % 4 == 3) {
				read_buffer[k/4] = temp_byte;
				temp_byte = 0;
			}
			else temp_byte = temp_byte * 4;
		}
		fout_bin.write (read_BQ_buffer, (read_len/4));


		// Write data
		fout_data << haplotype_seq << "\n";
		fout_data << read_seq << "\n";
		fout_data << read_BQ << "\n\n";

		// reset read seq
		read_seq = "";
        read_BQ = "";
	}
	
	fout_data.close();
	fout_bin.close();

	delete[] read_buffer;
	delete[] read_BQ_buffer;
	delete[] hap_buffer;

    return 0;
}
