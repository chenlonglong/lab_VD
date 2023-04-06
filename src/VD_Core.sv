`timescale 1ns/100ps



module VD_Core( 
    input                                               clk,
    input                                               rst_n,
    input refresh,
    output                                              o_ready,
    input                                               i_valid,
    input         [2*`HAP_MAX_LENGTH-1:0]               i_sequence_A,
    input         [2*`READ_MAX_LENGTH-1:0]              i_sequence_B,
    input         [2*`READ_MAX_LENGTH-1:0]              i_sequence_B_qualities,
    input         [$clog2(`HAP_MAX_LENGTH):0]           i_seq_A_length,
    input         [$clog2(`READ_MAX_LENGTH):0]          i_seq_B_length,
    input         [$clog2(`GENO_SRAM_WORD_AMOUNT)-1:0]  i_geno_address_ID,

    input                                               i_ready,
    output                                              o_valid,
    output signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]      o_alignment_score,
    output        [$clog2(`GENO_SRAM_WORD_AMOUNT)-1:0]  o_geno_address_ID
);
    integer i, j, k, l;
    genvar gv;
    localparam  S_idle                  = 4'd0,
                S_calculate             = 4'd1,
                S_select_highest        = 4'd2,
                S_output_score          = 4'd3;

    localparam MOST_NEGATIVE = {1'b1, {(`DP_PAIRHMM_SCORE_BITWIDTH-1){1'b0}}};
    localparam MOST_POSITIVE = {1'b0, {(`DP_PAIRHMM_SCORE_BITWIDTH-1){1'b1}}};




    ///////////////////////////// main registers ////////////////////////////////
    reg [3:0]                                               state, state_n;

    reg [$clog2(`HAP_MAX_LENGTH+`READ_MAX_LENGTH):0]        end_count, end_count_n;
    reg [2*`HAP_MAX_LENGTH-1:0]                             hap_seq, hap_seq_n;
    reg [2*`READ_MAX_LENGTH-1:0]                            read_seq, read_seq_n;
    reg [2*`READ_MAX_LENGTH-1:0]                            read_base_qualities, read_base_qualities_n;
    reg [$clog2(`HAP_MAX_LENGTH)-1:0]                       hap_length, hap_length_n;
    reg [$clog2(`READ_MAX_LENGTH)-1:0]                      read_length, read_length_n;
    reg [$clog2(`GENO_SRAM_WORD_AMOUNT)-1:0]                geno_address_ID, geno_address_ID_n;
    
    reg [$clog2(`HAP_MAX_LENGTH+`READ_MAX_LENGTH):0]        counter,     counter_n;
    reg [$clog2(`HAP_MAX_LENGTH+`READ_MAX_LENGTH):0]        counter_d,   counter_d_n;
    reg [$clog2(`HAP_MAX_LENGTH+`READ_MAX_LENGTH):0]        counter_dd,  counter_dd_n;
    reg [$clog2(`HAP_MAX_LENGTH+`READ_MAX_LENGTH):0]        counter_ddd, counter_ddd_n;
    reg [8-1:0]                                             shift_counter, shift_counter_n;

    reg first_pe_en, first_pe_en_n;
    
    reg [`READ_MAX_LENGTH*`DP_PAIRHMM_SCORE_BITWIDTH-1:0]   row_max_scores, row_max_scores_n;
    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]             row_max_score[0:`READ_MAX_LENGTH-1];
    always@(*) begin 
        for (i=0;i<`READ_MAX_LENGTH;i=i+1) row_max_score[i] = row_max_scores[i*`DP_PAIRHMM_SCORE_BITWIDTH+:`DP_PAIRHMM_SCORE_BITWIDTH];
    end


    //----------------------------------------------------------------------------------------
    wire PE_refresh = (refresh | state==S_idle);
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  PE_o_INDEL_dd[0:`READ_MAX_LENGTH-1];
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  PE_o_I_d[0:`READ_MAX_LENGTH-1];
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  PE_o_A_dd[0:`READ_MAX_LENGTH-1];
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  PE_o_A_d_add_M2I[0:`READ_MAX_LENGTH-1];
    wire                                          PE_o_valid_d[0:`READ_MAX_LENGTH-1];
    wire                                          PE_o_valid_dd[0:`READ_MAX_LENGTH-1];
    wire                                          PE_o_valid_ddd[0:`READ_MAX_LENGTH-1];
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  PE_o_Max[0:`READ_MAX_LENGTH-1];
    wire [1:0]                                    PE_o_hap_base[0:`READ_MAX_LENGTH-1];

    

    always@(*) begin
        state_n = state;
        case(state)
            S_idle:             state_n = (i_valid) ? S_calculate : state;
            S_calculate:        state_n = (counter_ddd == end_count) ? S_select_highest : state;
            S_select_highest:   state_n = (|shift_counter) ? state : S_output_score;
            S_output_score:     state_n = (i_ready) ? S_idle : state;
        endcase
    end

    always@(*) begin
        end_count_n             = end_count;
        hap_seq_n               = hap_seq;
        read_seq_n              = read_seq;
        read_base_qualities_n   = read_base_qualities;
        hap_length_n            = hap_length;
        read_length_n           = read_length;
        geno_address_ID_n       = geno_address_ID;

        counter_n               = counter;
        counter_d_n             = counter_d;
        counter_dd_n            = counter_dd;
        counter_ddd_n           = counter_ddd;
        shift_counter_n         = shift_counter;
        first_pe_en_n           = first_pe_en;
        row_max_scores_n        = row_max_scores;

        if (refresh) begin
            end_count_n = 0;
            hap_seq_n =  {(2*`HAP_MAX_LENGTH){1'b0}};
            read_seq_n = {(2*`READ_MAX_LENGTH){1'b0}};
            read_base_qualities_n = {(2*`READ_MAX_LENGTH){1'b0}};
            hap_length_n = 0;
            read_length_n = 0;
            geno_address_ID_n = 0;

            counter_n = 0;
            counter_d_n   = 0;
            counter_dd_n  = 0;
            counter_ddd_n = 0;
            first_pe_en_n = i_valid;
            row_max_scores_n = {`READ_MAX_LENGTH{MOST_NEGATIVE}};
            shift_counter_n = 0;
            
        end else begin
            case(state)
            S_idle: begin
            // ** TODO
            end
            S_calculate: begin
            // ** TODO
            end
            S_select_highest: begin
            // ** TODO
            end
            endcase
        end
    end



    generate
        for (gv=0;gv<`READ_MAX_LENGTH;gv=gv+1) begin: PEs
            if (gv==0) begin
                VD_PE u_VD_PE(
                    .clk                    (clk),
                    .rst_n                  (rst_n),
                    .refresh                (PE_refresh),
                    .i_en                   (first_pe_en),
                    .i_hap_base             (hap_seq[2*`HAP_MAX_LENGTH-1-:2]),
                    .i_read_base            (read_seq[2*`READ_MAX_LENGTH-1-2*gv-:2]),
                    .i_read_base_quality    (read_base_qualities[2*`READ_MAX_LENGTH-1-2*gv-:2]),
                    .i_A_top_add_M2I        (MOST_NEGATIVE),
                    .i_I_top                (MOST_NEGATIVE),
                    .i_A_diag               (MOST_NEGATIVE),
                    .i_INDEL_diag           (MOST_POSITIVE),
                    .o_INDEL_dd             (PE_o_INDEL_dd[gv]),
                    .o_I_d                  (PE_o_I_d[gv]),
                    .o_A_dd                 (PE_o_A_dd[gv]),
                    .o_A_d_add_M2I          (PE_o_A_d_add_M2I[gv]),
                    .o_valid_d              (PE_o_valid_d[gv]),
                    .o_valid_dd             (PE_o_valid_dd[gv]),
                    .o_valid_ddd            (PE_o_valid_ddd[gv]),
                    .o_Max                  (PE_o_Max[gv]),
                    .o_hap_base             (PE_o_hap_base[gv])
                );
            end else begin
                VD_PE u_VD_PE(
                    .clk                    (clk),
                    .rst_n                  (rst_n),
                    .refresh                (PE_refresh),
                    .i_en                   (PE_o_valid_d[gv-1]),
                    .i_hap_base             (PE_o_hap_base[gv-1]),
                    .i_read_base            (read_seq[2*`READ_MAX_LENGTH-1-2*gv-:2]),
                    .i_read_base_quality    (read_base_qualities[2*`READ_MAX_LENGTH-1-2*gv-:2]),
                    .i_A_top_add_M2I        (PE_o_A_d_add_M2I[gv-1]),
                    .i_I_top                (PE_o_I_d[gv-1]),
                    .i_A_diag               (PE_o_A_dd[gv-1]),
                    .i_INDEL_diag           (PE_o_INDEL_dd[gv-1]),
                    .o_INDEL_dd             (PE_o_INDEL_dd[gv]),
                    .o_I_d                  (PE_o_I_d[gv]),
                    .o_A_dd                 (PE_o_A_dd[gv]),
                    .o_A_d_add_M2I          (PE_o_A_d_add_M2I[gv]),
                    .o_valid_d              (PE_o_valid_d[gv]),
                    .o_valid_dd             (PE_o_valid_dd[gv]),
                    .o_valid_ddd            (PE_o_valid_ddd[gv]),
                    .o_Max                  (PE_o_Max[gv]),
                    .o_hap_base             (PE_o_hap_base[gv]) 
                );
            end
        end
    endgenerate


    assign o_ready                  = (state==S_idle);
    assign o_valid                  = (state==S_output_score);
    assign o_alignment_score        = row_max_scores[`READ_MAX_LENGTH*`DP_PAIRHMM_SCORE_BITWIDTH-1-:`DP_PAIRHMM_SCORE_BITWIDTH];
    assign o_geno_address_ID        = geno_address_ID;


    always@(posedge clk) begin
        if (!rst_n) begin
            state               <= S_idle;
            end_count           <= 0;
            hap_seq             <= {(2*`HAP_MAX_LENGTH){1'b0}};
            read_seq            <= {(2*`READ_MAX_LENGTH){1'b0}};
            read_base_qualities <= {(2*`READ_MAX_LENGTH){1'b0}};
            hap_length          <= 0;
            read_length         <= 0;
            geno_address_ID     <= 0;

            counter             <= 0;
            counter_d           <= 0;
            counter_dd          <= 0;
            counter_ddd         <= 0;
            shift_counter       <= 0;
            first_pe_en         <= 0;
            row_max_scores      <= {`READ_MAX_LENGTH{MOST_NEGATIVE}};
        end else begin
            state               <= state_n;
            end_count           <= end_count_n;
            hap_seq             <= hap_seq_n;
            read_seq            <= read_seq_n;
            read_base_qualities <= read_base_qualities_n;
            hap_length          <= hap_length_n;
            read_length         <= read_length_n;
            geno_address_ID     <= geno_address_ID_n;

            counter             <= counter_n;
            counter_d           <= counter_d_n;
            counter_dd          <= counter_dd_n;
            counter_ddd         <= counter_ddd_n;
            shift_counter       <= shift_counter_n;
            first_pe_en         <= first_pe_en_n;
            row_max_scores      <= row_max_scores_n;
        end
    end

endmodule


//----------------
// NTU DCS
// Yen-Lung Chen
//----------------
// insertion and deletion are both dased on query sequency

module VD_PE( 

    ///////////////////////////////////// I/Os //////////////////////////////////////
    input clk,
    input rst_n,
    input refresh,

    input                                           i_en,
    input [1:0]                                     i_hap_base,          // haplotype base
    input [1:0]                                     i_read_base,          // read base
    input [1:0]                                     i_read_base_quality,  // read base quality

    input signed  [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  i_A_top_add_M2I,
    input signed  [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  i_I_top,

    input signed  [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  i_A_diag,
    input signed  [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  i_INDEL_diag,

    output signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  o_INDEL_dd,
    output signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  o_I_d,
    output signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  o_A_dd,
    output signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  o_A_d_add_M2I,

    output                                          o_valid_d,
    output                                          o_valid_dd,
    output                                          o_valid_ddd,
    output signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]  o_Max,
    output [1:0]                                    o_hap_base
);

    localparam signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] MOST_NEGATIVE = {1'b1, {(`DP_PAIRHMM_SCORE_BITWIDTH-1){1'b0}}};
    localparam signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] MOST_POSITIVE = {1'b0, {(`DP_PAIRHMM_SCORE_BITWIDTH-1){1'b1}}};

    // DFF
    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] D_d, D_d_n;
    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] I_d, I_d_n;
    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] A_d, A_d_n;
    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] A_dd, A_dd_n;

    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] INDEL_dd, INDEL_dd_n;
    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] V_dd, V_dd_n;
    reg signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] Max,  Max_n;

    reg output_valid_d, output_valid_d_n;
    reg output_valid_dd, output_valid_dd_n;
    reg output_valid_ddd, output_valid_ddd_n;

    reg [1:0] hap_base, hap_base_n;

    always@(*) begin
        output_valid_d_n   = output_valid_d;
        output_valid_dd_n  = output_valid_dd;
        output_valid_ddd_n = output_valid_ddd;
        hap_base_n = hap_base;
        if (refresh) begin
            output_valid_d_n   = 0;
            output_valid_dd_n  = 0;
            output_valid_ddd_n = 0;
            hap_base_n = 0;
        end else begin
            output_valid_d_n   = i_en;
            output_valid_dd_n  = output_valid_d;
            output_valid_ddd_n = output_valid_dd;
            hap_base_n = i_hap_base;
        end
    end    

    // deletion
    
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH:0] D_left_add_I2I = D_d + `CONST_I2I;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] D_left_add_I2I_bounded;
    // ** TODO
    assign D_left_add_I2I_bounded = ;
    always@(*) begin
        D_d_n = D_d;
        if (refresh) D_d_n = MOST_NEGATIVE;
        // ** TODO
        else if (i_en) D_d_n = ;
    end


    // insertion
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH:0] I_top_add_I2I;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] I_top_add_I2I_bounded;
    // ** TODO
    assign I_top_add_I2I = ;
    assign I_top_add_I2I_bounded = ;
    always@(*) begin
        I_d_n = I_d;
        if (refresh) I_d_n = MOST_NEGATIVE;
        // ** TODO
        else if (i_en) I_d_n = ;
    end

    // INDEL
    always@(*) begin
        INDEL_dd_n = INDEL_dd;
        if (refresh) INDEL_dd_n = MOST_NEGATIVE;
        // ** TODO
        else if (output_valid_d) INDEL_dd_n = ;
    end

    // alignment
    //* match prior
    reg signed [`CONST_MATCH_BITWIDTH-1:0] the_match_prior;
    always@(*) begin
        case(i_read_base_quality)
        2'd0: the_match_prior = `CONST_BQ0_MATCH_SCORE;
        2'd1: the_match_prior = `CONST_BQ1_MATCH_SCORE;
        2'd2: the_match_prior = `CONST_BQ2_MATCH_SCORE;
        2'd3: the_match_prior = `CONST_BQ3_MATCH_SCORE;
        endcase
    end
    //* mismatch prior
    reg signed [`CONST_MISMATCH_BITWIDTH-1:0] the_mismatch_prior;
    always@(*) begin
        case(i_read_base_quality)
        2'd0: the_mismatch_prior = `CONST_BQ0_MISMATCH_SCORE;
        2'd1: the_mismatch_prior = `CONST_BQ1_MISMATCH_SCORE;
        2'd2: the_mismatch_prior = `CONST_BQ2_MISMATCH_SCORE;
        2'd3: the_mismatch_prior = `CONST_BQ3_MISMATCH_SCORE;
        endcase
    end
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]    prior;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH:0]      A_diag_add_M2M;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH:0]      INDEL_add_I2M;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH:0]      greater_diag;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH:0]      diag_add_prior;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0]    diag_bounded;
    // ** TODO
    prior               = ; 
    A_diag_add_M2M      = ; 
    INDEL_add_I2M       = ; 
    greater_diag        = ; 
    diag_add_prior      = ; 
    assign diag_bounded = ;




    always@(*) begin
        A_d_n = A_d;
        if (refresh) A_d_n = MOST_NEGATIVE;
        // ** TODO
        else if (i_en) A_d_n = ;

        A_dd_n = A_dd;
        if (refresh) A_dd_n = MOST_NEGATIVE;
        // ** TODO
        else A_dd_n = ;
    end
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH:0] A_d_add_M2I;
    wire signed [`DP_PAIRHMM_SCORE_BITWIDTH-1:0] A_d_add_M2I_bounded;
    // ** TODO
    assign A_d_add_M2I = ;
    assign A_d_add_M2I_bounded = ;

    // V & Max
    always@(*) begin
        V_dd_n = V_dd;
        if (refresh) V_dd_n = MOST_NEGATIVE;
        // ** TODO
        else V_dd_n = ;

        Max_n = Max;
        if (refresh) Max_n = MOST_NEGATIVE;
        // ** TODO
        else Max_n = ;
    end



    // Sequential citcuit
    always@(posedge clk) begin
        if (!rst_n) begin
            D_d                 <= MOST_NEGATIVE;
            I_d                 <= MOST_NEGATIVE;
            A_d                 <= MOST_NEGATIVE;
            A_dd                <= MOST_NEGATIVE;
            INDEL_dd            <= MOST_NEGATIVE;
            V_dd                <= MOST_NEGATIVE;
            Max                 <= MOST_NEGATIVE;
            output_valid_d      <= 0;
            output_valid_dd     <= 0;
            output_valid_ddd    <= 0;
            hap_base            <= 0;
        end else begin
            D_d                 <= D_d_n;
            I_d                 <= I_d_n;
            A_d                 <= A_d_n;
            A_dd                <= A_dd_n;
            INDEL_dd            <= INDEL_dd_n;
            V_dd                <= V_dd_n;
            Max                 <= Max_n;
            output_valid_d      <= output_valid_d_n;
            output_valid_dd     <= output_valid_dd_n;
            output_valid_ddd    <= output_valid_ddd_n;
            hap_base            <= hap_base_n;
        end
    end


    assign o_INDEL_dd       = INDEL_dd;
    assign o_I_d            = I_d;
    assign o_A_dd           = A_dd;
    assign o_A_d_add_M2I    = A_d_add_M2I_bounded;
    assign o_valid_d        = output_valid_d;
    assign o_valid_dd       = output_valid_dd;
    assign o_valid_ddd      = output_valid_ddd;
    assign o_Max            = Max;
    assign o_hap_base       = hap_base;
endmodule
