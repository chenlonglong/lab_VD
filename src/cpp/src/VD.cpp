
int viterbi_decoding(
    std::string haplotype, std::string read, std::string read_BQ,
    int bias_score,
    int m2m, int i2m, int m2i, int i2i,
    int BQ0_match, int BQ0_mismatch, 
    int BQ1_match, int BQ1_mismatch, 
    int BQ2_match, int BQ2_mismatch, 
    int BQ3_match, int BQ3_mismatch) 
{
    int align_scores[read.length()+1][haplotype.length()+1]  = {};
    int insert_scores[read.length()+1][haplotype.length()+1] = {};
    int delete_scores[read.length()+1][haplotype.length()+1] = {};
    
    for (int j = 0; j < haplotype.length() + 1; ++j) {
        align_scores [0][j] = -pow(2, 15);
        insert_scores[0][j] = -pow(2, 15);
        delete_scores[0][j] = bias_score;
    }
    for (int i = 1; i < read.length() + 1; ++i) {
        align_scores [i][0] = -pow(2, 15);
        insert_scores[i][0] = -pow(2, 15);
        delete_scores[i][0] = -pow(2, 15);
    }


    int scores[read.length()][haplotype.length()] = {0};

    int weight = 0;
    int temp_align_score = 0;
    int temp_insert_score = 0;
    int temp_delete_score = 0;
    int highest_score = -pow(2, 15);


    for (int i=0; i < read.length(); ++i) {
        for (int j=0; j < haplotype.length(); ++j) {
            // align scores operation
            char read_base = read[i];
            char haplotype_base = haplotype[j];
            if (read_BQ[i]=='0') {
                if (read_base == haplotype_base)
                    weight = BQ0_match;
                else
                    weight = BQ0_mismatch;
            } else if (read_BQ[i]=='1') {
                if (read_base == haplotype_base)
                    weight = BQ1_match;
                else
                    weight = BQ1_mismatch;
            } else if (read_BQ[i]=='2') {
                if (read_base == haplotype_base)
                    weight = BQ2_match;
                else
                    weight = BQ2_mismatch;
            } else if (read_BQ[i]=='3') {
                if (read_base == haplotype_base)
                    weight = BQ3_match;
                else
                    weight = BQ3_mismatch;
            } 


            temp_align_score  = align_scores[i][j]  + weight + m2m;
            temp_insert_score = insert_scores[i][j] + weight + i2m;
            temp_delete_score = delete_scores[i][j] + weight + i2m;
            // align_scores[i+1][j+1] = log10sumpow10(temp_align_score, temp_insert_score, temp_delete_score);
            if (temp_delete_score > temp_align_score && temp_delete_score > temp_insert_score)
                align_scores[i+1][j+1] = temp_delete_score;
            else if (temp_insert_score > temp_align_score)
                align_scores[i+1][j+1] = temp_insert_score;
            else
                align_scores[i+1][j+1] = temp_align_score;
            
            // bound align score to most negative
            if (align_scores[i+1][j+1] < -pow(2, 15))
                align_scores[i+1][j+1] = -pow(2, 15);

            // insert scores operation
            temp_align_score  = align_scores[i][j+1]  + m2i;
            temp_insert_score = insert_scores[i][j+1] + i2i;
            // insert_scores[i+1][j+1] = log10sumpow10(temp_align_score, temp_insert_score);
            if (temp_insert_score > temp_align_score)
                insert_scores[i+1][j+1] = temp_insert_score;
            else
                insert_scores[i+1][j+1] = temp_align_score;
            
            if (insert_scores[i+1][j+1] < -pow(2, 15))
                insert_scores[i+1][j+1] = -pow(2, 15);
            

            // delete scores operation
            temp_align_score  = align_scores[i+1][j]  + m2i;
            temp_delete_score = delete_scores[i+1][j] + i2i;
            // delete_scores[i+1][j+1] = log10sumpow10(temp_align_score, temp_delete_score);
            if (temp_delete_score > temp_align_score)
                delete_scores[i+1][j+1] = temp_delete_score;
            else
                delete_scores[i+1][j+1] = temp_align_score;
            
            if (delete_scores[i+1][j+1] < -pow(2, 15))
                delete_scores[i+1][j+1] = -pow(2, 15);

            // scores operation
            temp_align_score  = align_scores[i+1][j+1];
            temp_insert_score = insert_scores[i+1][j+1];
            temp_delete_score = delete_scores[i+1][j+1];
            // scores[i][j] = log10sumpow10(temp_align_score, temp_insert_score);
            if (temp_insert_score > temp_align_score)
                scores[i][j] = temp_insert_score;
            else
                scores[i][j] = temp_align_score;

            if (i == read.length()-1) {
                if (scores[i][j] > highest_score)
                    highest_score = scores[i][j];
            }
        }
    } // end of dynamic programming


    // return highest_score;
    return highest_score;
}

