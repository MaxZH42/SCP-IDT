disp('Calculating Optical Transfer Functions');

PixelShift_col = zeros(1, LED_Num);
PixelShift_row = zeros(1, LED_Num);

for i = 1:LED_Num
    PixelShift_col(i) = round(LED_NAx(i) * System_Pixelsize * Pic_Nx);
    PixelShift_row(i) = round(LED_NAy(i) * System_Pixelsize * Pic_Ny);
end

% Calculate 4D optical transfer function

PTF_4D = single(zeros(Pic_Ny, Pic_Nx, Pic_Nz, LED_Num));
ATF_4D = single(zeros(Pic_Ny, Pic_Nx, Pic_Nz, LED_Num));
PTF_3D = single(zeros(Pic_Ny, Pic_Nx, Pic_Nz));
ATF_3D = single(zeros(Pic_Ny, Pic_Nx, Pic_Nz));

for i = 1:LED_Num

    LED_pos = i;

    LED_kx = LED_NAx(LED_pos);
    LED_ky = LED_NAy(LED_pos);

    G  = real(1 ./ sqrt((N_medium/LED_Lambda).^2 - ((fx2D-LED_kx).^2 + (fy2D-LED_ky).^2)));
    Gf = real(1 ./ sqrt((N_medium/LED_Lambda).^2 - ((fx2D+LED_kx).^2 + (fy2D+LED_ky).^2)));

    uv_vector1 = real(sqrt((N_medium/LED_Lambda).^2 - ((fx2D-LED_kx).^2 + (fy2D-LED_ky).^2)));
    uv_vector2 = real(sqrt((N_medium/LED_Lambda).^2 - ((fx2D+LED_kx).^2 + (fy2D+LED_ky).^2)));

    Pupil  = circshift(System_Aperture,  [PixelShift_row(LED_pos),  PixelShift_col(LED_pos)]);
    Pupilf = circshift(System_Aperture, -[PixelShift_row(LED_pos),  PixelShift_col(LED_pos)]);

    LED_kxy = sqrt((N_medium/LED_Lambda).^2 - (LED_kx^2 + LED_ky^2));

    for j = 1:Pic_Nz
        PTF_3D(:,:,j) = ...
            (Pupil  .* sin(2*pi.*Depth_Set(j).*(uv_vector1 - LED_kxy)) .* G + ...
             Pupilf .* sin(2*pi.*Depth_Set(j).*(uv_vector2 - LED_kxy)) .* Gf) + ...
            1i * (Pupil  .* cos(2*pi.*Depth_Set(j).*(uv_vector1 - LED_kxy)) .* G - ...
                  Pupilf .* cos(2*pi.*Depth_Set(j).*(uv_vector2 - LED_kxy)) .* Gf);

        ATF_3D(:,:,j) = ...
            -(Pupil  .* cos(2*pi.*Depth_Set(j).*(uv_vector1 - LED_kxy)) .* G + ...
              Pupilf .* cos(2*pi.*Depth_Set(j).*(uv_vector2 - LED_kxy)) .* Gf) + ...
             1i * (Pupil  .* sin(2*pi.*Depth_Set(j).*(uv_vector1 - LED_kxy)) .* G - ...
                   Pupilf .* sin(2*pi.*Depth_Set(j).*(uv_vector2 - LED_kxy)) .* Gf);

        PTF_3D(:,:,j) = fftshift(PTF_3D(:,:,j));
        ATF_3D(:,:,j) = fftshift(ATF_3D(:,:,j));
    end

    PTF_4D(:,:,:,LED_pos) = single(0.5 * Z_step * k_G^2 .* PTF_3D);
    ATF_4D(:,:,:,LED_pos) = single(0.5 * Z_step * k_G^2 .* ATF_3D);
end

PTF_4D_Used = PTF_4D;
ATF_4D_Used = ATF_4D;

conj_PTF_4D_Used = conj(PTF_4D_Used);
conj_ATF_4D_Used = conj(ATF_4D_Used);

clear DC PTF_4D ATF_4D PTF_3D ATF_3D PixelShift_col PixelShift_row G Gf uv_vector1 uv_vector2 Pupil Pupilf;
clear LED_pos LED_kx LED_ky LED_kxy LED_NAx LED_NAy LED_NAz LED_FDL LED_XN LED_YN;
clear LED_Expect_NA LED_Radius LED_Radius_ALL LED_Theta LED_Theta_ALL;