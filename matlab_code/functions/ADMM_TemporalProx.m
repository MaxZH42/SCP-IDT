temp_diff = SP1_k - SP_T_prev;                          % Forward temporal gradient

V_t = temp_diff + y3_k;                                 % Proximal operator argument

DSP_T_k = sign(V_t) .* max(abs(V_t) - tau_T/rho_T, 0);  % Soft-thresholding for L1-norm temporal TV prior

y3_k = y3_k + temp_diff - DSP_T_k;                      % Dual variable update
