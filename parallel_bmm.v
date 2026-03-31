module parallel_bmm #(parameter N = 256, parameter W = 32, parameter ALPHA = 7, parameter BETA = 8) (input clk, reset,start ,input [N-1:0]A,B,M, input [2*N-1:0]mu,output reg [N-1:0]Z);

reg [N-1:0]A_reg,B_reg,M_reg,Z_reg;
reg [2*N-1:0]mu_reg;
wire [N-1:0]S_H,S_L;
reg [N-1:0]S_H_reg,S_L_reg,T_final,T_reg;
reg [N/2+N-1:0]T = 0;
wire done1, done2;


wire [N/2-1:0] A_H,A_L;

assign A_H = A_reg[N-1:N/2];
assign A_L = A_reg[(N/2)-1:0];

barrett_bmm_general #(.N(N), .W(W), .ALPHA(ALPHA), .BETA(BETA)) inst1 (.A(A_H), .B(B_reg), .M(M_reg), .mu(mu_reg), .clk(clk), .reset(reset), .start(start), .done(done1), .Z(S_H));

barrett_bmm_general #(.N(N), .W(W), .ALPHA(ALPHA), .BETA(BETA)) inst2 (.A(A_L), .B(B_reg), .M(M_reg), .mu(mu_reg), .clk(clk), .reset(reset), .start(start), .done(done2), .Z(S_L));

always@(posedge clk) begin

    A_reg = A;
    B_reg = B;
    M_reg = M;
    mu_reg = mu;
    S_H_reg = S_H;
    S_L_reg = S_L;
    T_reg = T_final;
    Z = Z_reg;

end

always@(*) begin

//    T = (S_H_reg * (2 ** (N/2)));

   T[N/2+N-1:N/2] = S_H_reg ;
 
  // T_final = T[N-1:0] - ((T[N-1:0] >= M) ? M : 0);

    T_final = T % M;
   
   if((T_reg + S_L_reg) >= M)
    
        Z_reg = (T_reg + S_L_reg) - M ;
        
   else

        Z_reg = (T_reg + S_L_reg) ;
end
endmodule 


module barrett_bmm_general #(
    parameter N  = 256,    // modulus width
    parameter W  = 32,     // Comba word size
    parameter ALPHA = 7,  // General α value
    parameter BETA  = 8   // General β value
)(
    input               clk,
    input               reset,
    input               start,
    input  [N/2-1:0]      A,
    input  [N-1:0]      B,
    input  [N-1:0]      M,
    input  [2*N-1:0]    mu,     // mu = floor( 2^(2N+α) / M )

    output reg          done,
    output reg [N-1:0]  Z
);

    // ============================================================
    // STEP 1: P = A * B
    // ============================================================
    wire [N/2+N-1:0] P_full;
    wire done_P;
    reg start_P;

    pipelined_comba_mult #(
        .AW(N/2),
        .BW(N),
        .W(W)
    ) mult_AB (
        .clk(clk),
        .start(start_P),
        .A(A),
        .B(B),
        .done(done_P),
        .P(P_full)
    );

    // ============================================================
    // STEP 2: PS = P >> (N - β)
    // ============================================================
    localparam SHIFT_PS = N - BETA;

    wire [N/2 + W - 1:0] PS = P_full[ N/2 + N - 1 : SHIFT_PS ];

    // ============================================================
    // STEP 3: q1 = PS * μ
    // ============================================================
    wire [N/2 + 2*N + W - 1 : 0] q1_full;  // very wide
    wire done_q1;
    reg start_q1;

    pipelined_comba_mult #(
        .AW(N/2 + W),
        .BW(2*N),
        .W(W)
    ) mult_PS_mu (
        .clk(clk),
        .start(start_q1),
        .A(PS),
        .B(mu),
        .done(done_q1),
        .P(q1_full)
    );

    // ============================================================
    // STEP 4: q = q1 >> (2N + α + β)
    // ============================================================
    localparam SHIFT_Q = N + ALPHA + BETA;

    wire [N-1:0] q = q1_full[ N/2 + 2*N + W - 1 : SHIFT_Q ];

    // ============================================================
    // STEP 5: T = q * M
    // ============================================================
    wire [2*N-1:0] T_full;
    wire done_T;
    reg start_T;

    pipelined_comba_mult #(
        .AW(N),
        .BW(N),
        .W(W)
    ) mult_qM (
        .clk(clk),
        .start(start_T),
        .A(q),
        .B(M),
        .done(done_T),
        .P(T_full)
    );

    // ============================================================
    // STEP 6: Z_temp = P - T
    // ============================================================
    wire [N-1:0] Z_temp = P_full[N-1:0] - T_full[N-1:0];

    // ============================================================
    // STEP 7: correction Z ≥ M
    // ============================================================
    wire [N-1:0] Z_corr = (Z_temp >= M) ? (Z_temp - M) : Z_temp;

    // ============================================================
    // FSM CONTROLLER
    // ============================================================
    reg [2:0] state;
    localparam S_IDLE = 0,
               S_P    = 1,
               S_Q1   = 2,
               S_QM   = 3,
               S_OUT  = 4;

    always @(posedge clk) begin
        if (reset) begin
            start_P  <= 0;
            start_q1 <= 0;
            start_T  <= 0;
            Z        <= 0;
            done     <= 0;
            state    <= S_IDLE;
        end else begin
            done <= 0;

            case (state)

                S_IDLE: begin
                    start_P  <= 0;
                    start_q1 <= 0;
                    start_T  <= 0;

                    if (start) begin
                        start_P <= 1;
                        state   <= S_P;
                    end
                end

                S_P: begin
                    
                    if (done_P) begin
                        start_P <= 0;
                        start_q1 <= 1;
                        state    <= S_Q1;
                    end
                end

                S_Q1: begin
                    
                    if (done_q1) begin
                        start_q1 <= 0;
                        start_T <= 1;
                        state   <= S_QM;
                    end
                end

                S_QM: begin
                    
                    if (done_T) begin
                        start_T <= 0;
                        Z    <= Z_corr;
                        done <= 1;
                        state <= S_OUT;
                    end
                end

                S_OUT: begin
                    state <= S_IDLE;
                end
            endcase
        end
    end

endmodule



 // PIPELINED COMBA MULTIPLIER

module pipelined_comba_mult #(
    parameter AW = 128,    // bitwidth of A
    parameter BW = 256,    // bitwidth of B
    parameter W  = 32      // word size
)(
    input  wire             clk,
    input  wire             start,
    input  wire [AW-1:0]    A,
    input  wire [BW-1:0]    B,
    output reg              done,
    output reg [AW+BW-1:0]  P
);

    // ----------------------------
    // Number of words
    // ----------------------------
    localparam NA = (AW + W - 1) / W;
    localparam NB = (BW + W - 1) / W;
    localparam NW = NA + NB;

    // ----------------------------
    // Split A and B into words
    // ----------------------------
    reg [W-1:0] Aword [0:NA-1];
    reg [W-1:0] Bword [0:NB-1];
    
    integer i;



    // ----------------------------
    // Pipeline registers
    // ----------------------------
    reg [15:0] k;  // column index
    reg busy = 0;

    reg [2*W+15:0] carry;
    reg [W-1:0] C [0:NW-1];   // output words

    integer j,m;
    reg [2*W+15:0] temp;
    
    reg [W+1:0] booth_reg;
    reg [2*W -1:0] partial_sum, T;

    // ----------------------------
    // MAIN FSM
    // ----------------------------
    always @(posedge clk) begin
        done <= 0;

        if (start && !busy) begin
        
            for (i = 0; i < NA; i = i + 1)
                 Aword[i] <= A[(i*W) +: W];

            for (i = 0; i < NB; i = i + 1)
                 Bword[i] <= B[(i*W) +: W];
                
            busy  <= 1;
            k     <= 0;
            carry <= 0;
        end

        else if (busy) begin
            // Compute one column per cycle
            temp = carry;

            for (j = 0; j < NA; j = j + 1)
                if ((k-j) >= 0 && (k-j) < NB) begin
                    
                    
                    booth_reg =   {Bword[k-j], 1'b0};
                    partial_sum = 0;
                    
                    for (m = 0; m <= W / 2; m = m + 1) begin
                        case (booth_reg[2:0])
                            3'b000, 3'b111: partial_sum = partial_sum;            
                            3'b001, 3'b010: partial_sum = partial_sum + (Aword[j] << ( m<<1)); 
                            3'b011:          partial_sum = partial_sum + (Aword[j] << ((m<<1) + 1)); 
                            3'b100:          partial_sum = partial_sum - (Aword[j] << ((m<<1) + 1)); 
                            3'b101, 3'b110:  partial_sum = partial_sum - (Aword[j] << ((m<<1))); 
                        endcase
                        booth_reg = booth_reg >> 2; 
                    end
                    
                    T = partial_sum;
                    
                    temp = temp + T;
                    
                    end

            // Extract output word
            C[k] <= temp[W-1:0];

            // Next carry
            carry <= temp >> W;

            // Advance to next column
            k <= k + 1;

            // Finished all columns?
            if (k == NW-1) begin
                busy <= 0;

                // Assemble output product`

                done <= 1;
            end
        end
    end
    
    always@(*) begin
        for (j = NW; j > 0; j = j - 1)
                    P[(j*W)-1 -: W] = C[j-1];
    end

endmodule







