temporalSwitch = temporalTV && ~markBegin && (temporal_iter <= 1);

for iter = 1:max_iter

    eval ADMM_SP1Optimize;

    eval ADMM_HessianProx;

    y1_k = y1_k + delta;

    eval ADMM_ProjProx;
    y2_k = y2_k + (SP1_k - SP2_k);
    SP_curr = SP2_k;

    if temporalSwitch
        eval ADMM_TemporalProx;
    end

    if (~temporalSwitch) && temporalTV && (~markBegin) && (iter == temporal_iter)
        DSP_T_k = SP1_k - SP_T_prev;
        temporalSwitch = true;
    end

    if iter >= min_iter
        rel_change = norm(SP_curr(:) - SP_last(:)) / (norm(SP_last(:)) + eps);
        if rel_change < 5e-3
            fprintf(': ADMM converged at iter %d', iter);
            break;
        end
    end

    SP_last = SP_curr;

end

v_re = SP_curr(:,:,:,1);
v_im = SP_curr(:,:,:,2);
SP_T_prev = SP_curr;

n_re = sqrt(((N_medium.^2 + v_re) + sqrt((N_medium^2 + v_re).^2 + v_im.^2)) / 2);
n_im = v_im ./ n_re ./ 2;

RI=n_re+1i*n_im;
