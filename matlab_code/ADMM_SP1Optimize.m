SP1_k = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 2, 'like', conj_PTF_Iten);

AHI1_k1 = conj_PTF_Iten + rho_S * (...
    conj(Dxx).*fft2(DSP_k(:,:,:,1) + y1_k(:,:,:,1)) + ...
    conj(Dyy).*fft2(DSP_k(:,:,:,2) + y1_k(:,:,:,2)) + ...
    conj(Dxy).*fft2(DSP_k(:,:,:,3) + y1_k(:,:,:,3)));

AHI2_k1 = conj_ATF_Iten + rho_S * (...
    conj(Dxx).*fft2(DSP_k(:,:,:,4) + y1_k(:,:,:,4)) + ...
    conj(Dyy).*fft2(DSP_k(:,:,:,5) + y1_k(:,:,:,5)) + ...
    conj(Dxy).*fft2(DSP_k(:,:,:,6) + y1_k(:,:,:,6)));

AHI1_k1 = AHI1_k1 + rho_S * fft2(SP2_k(:,:,:,1) - y2_k(:,:,:,1));
AHI2_k1 = AHI2_k1 + rho_S * fft2(SP2_k(:,:,:,2) - y2_k(:,:,:,2));

if ~temporalSwitch
    SP1_k(:,:,:,1) = real(ifft2((AHI1_k1.*sum_ATF - AHI2_k1.*conj_term1) ./ Normalized_term));
    SP1_k(:,:,:,2) = real(ifft2((AHI2_k1.*sum_PTF - AHI1_k1.*conj_term2) ./ Normalized_term));
else
    AHI1_k1 = AHI1_k1 + rho_T * fft2(SP_T_prev(:,:,:,1) + DSP_T_k(:,:,:,1) - y3_k(:,:,:,1));
    AHI2_k1 = AHI2_k1 + rho_T * fft2(SP_T_prev(:,:,:,2) + DSP_T_k(:,:,:,2) - y3_k(:,:,:,2));
    SP1_k(:,:,:,1) = real(ifft2((AHI1_k1.*sum_ATF_wT - AHI2_k1.*conj_term1) ./ Normalized_term_wT));
    SP1_k(:,:,:,2) = real(ifft2((AHI2_k1.*sum_PTF_wT - AHI1_k1.*conj_term2) ./ Normalized_term_wT));
end