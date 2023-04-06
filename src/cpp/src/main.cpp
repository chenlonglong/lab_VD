#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <tuple>
#include <array>
#include <queue>
#include <cmath>

#include "VD.cpp"
using namespace std;

// For VC parameters
#define CONST_MATCH_SCORE		1
#define CONST_MISMATCH_SCORE	-4
#define CONST_GAP_OPEN			-6
#define CONST_GAP_EXTEND		-1

int main(int argc, char** argv)
{
	string	data_path 	= "";

    for(int i=0;i<argc;i++) {
        string argv_str(argv[i]);

		// Turn the argv into uint32_t (stoul) and int (stoi)
        if(argv_str == "-data_path")    data_path   = argv[i+1];
    }


    if (argc < 2) {
        cerr << "You should enter data path!!!\n\n";
        return -1;
    }

	// Score parameters
	int match_score 		= CONST_MATCH_SCORE;
	int mismatch_score 		= CONST_MISMATCH_SCORE;
	int gap_open_score		= CONST_GAP_OPEN;
	int gap_extend_score 	= CONST_GAP_EXTEND;
	
	cout << "========= Start Smith-Waterman =========" << "\n";
    cout << "Data path: "          << data_path << "\n";
	cout << "Match score: "        << match_score << "\n";
	cout << "Mismatch score: "     << mismatch_score << "\n";
	cout << "Gap open score: "     << gap_open_score << "\n";
	cout << "Gap extend score: "   << gap_extend_score << "\n";
	cout << "========================================" << "\n\n";

	// Read in File
    ifstream fin_data(data_path);
    if(fin_data) cerr << "Testdata file ok!\n";
    else {
        cerr << "Testdata file err!\n";
        return 1;
    }

	// Write File
	ofstream fout_golden("./test_data/pattern_ans.txt");
	if(!fout_golden) {
        cerr << "Cannot open file!" << endl;
        return 1;
   	}else {
        cerr << "Open output file successfully\n";
    }

	string haplotype_seq;
	string read_seq;
	string read_BQ;


    cout << "hi\n\n";

    float gatk_m2m = 0.998; //0.9998; //0.998;
    float gatk_i2m = 0.9;
    float gatk_m2i = 0.001; //0.0001; //0.001;
    float gatk_i2i = 0.1;

    float BQ0_mismatch_prior = 0.631; 
    float BQ1_mismatch_prior = 0.0794;
    float BQ2_mismatch_prior = 0.0032;
    float BQ3_mismatch_prior = 0.0002;
    float BQ0_match_prior    = 1.0 - BQ0_mismatch_prior;
    float BQ1_match_prior    = 1.0 - BQ1_mismatch_prior;
    float BQ2_match_prior    = 1.0 - BQ2_mismatch_prior;
    float BQ3_match_prior    = 1.0 - BQ3_mismatch_prior;
    BQ0_mismatch_prior /= 3;
    BQ1_mismatch_prior /= 3;
    BQ2_mismatch_prior /= 3;
    BQ3_mismatch_prior /= 3;

    int total_bits             = 16;
    int decimal_place_bits     = 10;
    int the_bias_score         = pow(2, total_bits - 1) - 1;
    int log_m2m                = round(log10(gatk_m2m)*pow(2, decimal_place_bits));
    int log_i2m                = round(log10(gatk_i2m)*pow(2, decimal_place_bits));
    int log_m2i                = round(log10(gatk_m2i)*pow(2, decimal_place_bits));
    int log_i2i                = round(log10(gatk_i2i)*pow(2, decimal_place_bits));
    int log_BQ3_match_prior    = round(log10(BQ3_match_prior   )*pow(2, decimal_place_bits));        //0;
    int log_BQ3_mismatch_prior = round(log10(BQ3_mismatch_prior)*pow(2, decimal_place_bits));    //-3968;
    int log_BQ2_match_prior    = round(log10(BQ2_match_prior   )*pow(2, decimal_place_bits));       //-2;
    int log_BQ2_mismatch_prior = round(log10(BQ2_mismatch_prior)*pow(2, decimal_place_bits));    //-2944;
    int log_BQ1_match_prior    = round(log10(BQ1_match_prior   )*pow(2, decimal_place_bits));      //-18;
    int log_BQ1_mismatch_prior = round(log10(BQ1_mismatch_prior)*pow(2, decimal_place_bits));    //-1920;
    int log_BQ0_match_prior    = round(log10(BQ0_match_prior   )*pow(2, decimal_place_bits));     //-227;
    int log_BQ0_mismatch_prior = round(log10(BQ0_mismatch_prior)*pow(2, decimal_place_bits));     //-896;

    int highest_score;
    int data_count = 0;
    while (!fin_data.eof()) {
		fin_data >> haplotype_seq;
		fin_data >> read_seq;
		fin_data >> read_BQ;
        data_count ++;
		cout << "Set " << dec << data_count << "\n";
		cout << "Haplotype seq: \n"     << haplotype_seq << "\n";
		cout << "Read seq: \n"          << read_seq << "\n";
		cout << "Read base qiality: \n" << read_BQ << "\n";

		// Execute Viterbi Algorithm for Pair-HMM
        highest_score = viterbi_decoding(
            haplotype_seq,
            read_seq,
            read_BQ,
            the_bias_score,
            log_m2m,
            log_i2m,
            log_m2i,
            log_i2i,
            log_BQ0_match_prior,
            log_BQ0_mismatch_prior, 
            log_BQ1_match_prior,
            log_BQ1_mismatch_prior, 
            log_BQ2_match_prior,
            log_BQ2_mismatch_prior, 
            log_BQ3_match_prior,
            log_BQ3_mismatch_prior);


		cout << "Highest score: " << dec << highest_score << "\n";
        cout << "=============================\n";

		fout_golden << dec << highest_score << "\n";
	}


	fin_data.close();
	fout_golden.close();

    return 0;
}
