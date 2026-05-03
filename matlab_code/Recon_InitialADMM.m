SP1_k   = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 2);
DSP_k   = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 6);
y1_k    = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 6);
SP2_k   = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 2);
y2_k    = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 2);
delta   = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 6);

if temporalTV
    DSP_T_k = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 2);
    y3_k    = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 2);
end

if isempty(SP_T_prev)
    SP_T_prev = zeros(Pic_Ny, Pic_Nx, Pic_Nz, 2);
end

temporalSwitch = false;