% Credit to: https://github.com/Waller-Lab/3DQuantitativeDPC

direction = [1,-1];
SP_re = SP1_k(:,:,:,1) + y2_k(:,:,:,1);
SP_im = SP1_k(:,:,:,2) + y2_k(:,:,:,2);

if direction(1) == 1
    % Positivity constraint on the real part of scattering potential
    SP_re(SP_re < 0) = 0;
else
    % Negativity constraint on the real part of scattering potential
    SP_re(SP_re > 0) = 0;
end

if direction(2) == 1
    % Positivity constraint on the imaginary part of scattering potential
    SP_im(SP_im < 0) = 0;
else
    % Negativity constraint on the imaginary part of scattering potential
    SP_im(SP_im > 0) = 0;
end

SP2_k = cat(4, SP_re, SP_im);