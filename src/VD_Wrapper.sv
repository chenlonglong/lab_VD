

//* ----- Global parameters -----
`define HAP_MAX_LENGTH              400
`define READ_MAX_LENGTH             150

//* DP for Genotyping
`define DP_PAIRHMM_SCORE_BITWIDTH   16

`define CONST_M2M                   -1
`define CONST_M2I                   -3072
`define CONST_I2I                   -1024
`define CONST_I2M                   -47

`define CONST_MATCH_BITWIDTH        16
`define CONST_BQ0_MATCH_SCORE       -443//-227
`define CONST_BQ1_MATCH_SCORE       -37//-18
`define CONST_BQ2_MATCH_SCORE       -1//-2
`define CONST_BQ3_MATCH_SCORE       0

`define CONST_MISMATCH_BITWIDTH     16
`define CONST_BQ0_MISMATCH_SCORE    -693//-896
`define CONST_BQ1_MISMATCH_SCORE    -1615//-1920
`define CONST_BQ2_MISMATCH_SCORE    -3043//-2944
`define CONST_BQ3_MISMATCH_SCORE    -4276//-3968



module VD_Wrapper (
    input         avm_rst,
    input         avm_clk,
    output  [4:0] avm_address,
    output        avm_read,
    input  [31:0] avm_readdata,
    output        avm_write,
    output [31:0] avm_writedata,
    input         avm_waitrequest
);

localparam RX_BASE     = 0*4;
localparam TX_BASE     = 1*4;
localparam STATUS_BASE = 2*4;
localparam TX_OK_BIT   = 6;
localparam RX_OK_BIT   = 7;

// Feel free to design your own FSM!

// Remember to complete the port connection
VD_core vd_core(
    .clk                (avm_clk),
    .rst                (avm_rst),

    .o_ready            (),
    .i_valid            (),
    .i_sequence_ref     (),
    .i_sequence_read    (),
    .i_seq_ref_length   (),
    .i_seq_read_length  (),
    
    .i_ready            (),
    .o_valid            (),
    .o_alignment_score  (),
    .o_column           (),
    .o_row              ()
);


// ** TODO
always_comb begin
    
end

// ** TODO
always_ff @(posedge avm_clk or posedge avm_rst) begin
    if (avm_rst) begin
        

    end
    else begin
        

    end
end

endmodule
