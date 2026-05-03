clc;

modified        = false;
co02            = 0;
Pupil_AO        = System_Aperture;
markEnd         = false;
markBegin       = true;

if timeStop == 1
    slidingWindow = false;
end

if ~slidingWindow
    inter       = LED_Num;
end

frameTotal      = LED_Num / inter;

sum_PTF         = 0;
sum_ATF         = 0;
conj_PTF_Iten   = 0;
conj_ATF_Iten   = 0;
conj_term1      = 0;
conj_term2      = 0;

for i = 1:LED_Num

    P = PTF_4D_Used(:,:,:,i);
    A = ATF_4D_Used(:,:,:,i);

    sum_PTF = sum_PTF + real(P).^2 + imag(P).^2;
    sum_ATF = sum_ATF + real(A).^2 + imag(A).^2;

    conj_term1 = conj_term1 + conj(P) .* A;

end

conj_term2 = conj(conj_term1);

SP_T_prev = [];

clear P A