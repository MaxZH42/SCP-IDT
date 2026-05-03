% Title : Spatiotemporal-Continuity-Prior-Guided Intensity Diffraction Tomography
% Author: Zihao Zhou, Ning Zhou

clc; clear; close all;

addpath('.\functions');
addpath1 = '.\IDT';                     % Captured images
addpath2 = '.\BG';                      % Captured background images


%% 1. System Parameters

% LED parameters
LED_Lambda      = 0.517;                % Illumination wavelength (μm)
LED_Expect_NA   = 0.745;                % Target illumination numerical aperture (NA)
LED_Radius      = 97.09*1000;           % Radial distance of LED array from optical axis (μm)
LED_Num         = 28;                   % Total number of illumination angles (LEDs)

% Angular positions of individual LEDs
theta_need      = -0.58;
LED_Theta       = (50.00+theta_need)*pi/90 : 2*pi/LED_Num : (230.00+theta_need)*pi/90-2*pi/LED_Num;

% Axial distance from LED plane to sample plane
LED_FDL         = sqrt((1-LED_Expect_NA.^2).*LED_Radius.^2)./LED_Expect_NA;

% OBJ parameters
OBJ_NA          = 0.75;                 % Objective lens numerical aperture (NA)
OBJ_Mag         = 40;                   % Objective lens magnification

% Camera parameters
Cam_Nx          = 1200;                 % Sensor width in pixels
Cam_Ny          = 1200;                 % Sensor height in pixels
Cam_Pixelsize   = 6.5;                  % Pixel size of the camera (μm)

% Volume parameters
Pic_Nx          = 500;                  % Width for reconstruction
Pic_Ny          = 500;                  % Height for reconstruction
Pic_Nz          = 20;                   % Number of axial (depth) reconstruction planes
Z_step          = 0.2;                  % Axial step size (μm)

% Other parameters
N_medium        = 1.33;                 % Refractive index of the background/immersion medium
k_G             = 2*pi/LED_Lambda;      % Wave number of illumination light (rad/μm)
Lambda_n        = LED_Lambda/N_medium;  % Wavelength of light within the medium (μm)


%% 2. LED Illumination Spatial Frequencies

% Convert LED array polar coordinates to Cartesian positions in the illumination plane
LED_Radius_ALL   = repmat(LED_Radius,1,LED_Num);
LED_Theta_ALL    = LED_Theta;
[LED_XN, LED_YN] = pol2cart(LED_Theta_ALL, LED_Radius_ALL);

LED_NAx = zeros(1,LED_Num);
LED_NAy = zeros(1,LED_Num);
LED_NAz = zeros(1,LED_Num);

for i = 1:LED_Num
    LED_NAx(i) = LED_XN(i) ./ sqrt(LED_XN(i)^2 + LED_YN(i)^2 + LED_FDL^2) / LED_Lambda;
    LED_NAy(i) = LED_YN(i) ./ sqrt(LED_XN(i)^2 + LED_YN(i)^2 + LED_FDL^2) / LED_Lambda;
    LED_NAz(i) = real(sqrt((N_medium/LED_Lambda)^2 - LED_NAx(i)^2 - LED_NAy(i)^2));
end


%% 3. Optical Transfer Function (OTF)

% Effective pixel size at the sample plane
System_Pixelsize = Cam_Pixelsize / OBJ_Mag;

% Spatial coordinates at the object plane
X = (-fix(Pic_Nx/2):1:fix((Pic_Nx-1)/2)) * System_Pixelsize;
Y = (-fix(Pic_Ny/2):1:fix((Pic_Ny-1)/2)) * System_Pixelsize;
Z = (-fix(Pic_Nz/2):1:fix((Pic_Nz-1)/2)) * Z_step;

% Frequency domain sampling intervals
[x3D, y3D, z3D] = meshgrid(X, Y, Z);

delta_x = 1/(System_Pixelsize * Pic_Nx);
delta_y = 1/(System_Pixelsize * Pic_Ny);
delta_z = 1/(Z_step * Pic_Nz);

Fx = (-fix(Pic_Nx/2):1:fix((Pic_Nx-1)/2)) * delta_x;
Fy = (-fix(Pic_Ny/2):1:fix((Pic_Ny-1)/2)) * delta_y;
Fz = (-fix(Pic_Nz/2):1:fix((Pic_Nz-1)/2)) * delta_z;

[fx2D , fy2D]       = meshgrid(Fx, Fy);
[Theta, Fre_radius] = cart2pol(fx2D, fy2D);

fz2D = real(sqrt((N_medium/LED_Lambda).^2 - fx2D.^2 - fy2D.^2));

% System aperture function
Fre_cutoff      = OBJ_NA / LED_Lambda;
System_Aperture = double(~(Fre_radius > Fre_cutoff));

PixelShift_col_AO = zeros(1, LED_Num);
PixelShift_row_AO = zeros(1, LED_Num);

for i = 1:LED_Num
    PixelShift_col_AO(i) = round(LED_NAx(i) * System_Pixelsize * Pic_Nx);
    PixelShift_row_AO(i) = round(LED_NAy(i) * System_Pixelsize * Pic_Ny);
end

% Axial depth profile for OTF generation
Depth_Set = Z;

eval OTF_Generate;


%% 4. Reconstruction Preparation

% Define Region of Interest (ROI)
Select_RangeX       = Pic_Nx;
Select_RangeY       = Pic_Ny;
Select_CenterPointX = Cam_Nx / 2;
Select_CenterPointY = Cam_Ny / 2;

colStart = round(Select_CenterPointX - Select_RangeX/2);
rowStart = round(Select_CenterPointY - Select_RangeY/2);

roiCols = colStart : (colStart + Select_RangeX - 1);
roiRows = rowStart : (rowStart + Select_RangeY - 1);

assert(min(roiCols) >= 1 && max(roiCols) <= Cam_Nx, 'ROI column indices exceed camera width.');
assert(min(roiRows) >= 1 && max(roiRows) <= Cam_Ny, 'ROI row indices exceed camera height.');


%% 5. Reconstruction

% Method parameters
timeStop         = 1;                   % Total number of timestamps to process in the sequence
slidingWindow    = true;                % Sliding window processing
inter            = 4;                   % Reconstruction interval
temporalTV       = true;                % Temporal Total Variation (TV) regularization

% Regularization parameters
% => Hessian
rho_S            = 1e1;                 % ADMM augmented Lagrangian penalty parameter for spatial subproblem
tau_S            = 0.005 * rho_S;       % Proximal step size for Hessian prior
max_iter         = 80;
min_iter         = 30;

% => Temporal TV
rho_T            = 8e-3;                % ADMM augmented Lagrangian penalty parameter for temporal subproblem
tau_T            = 0.10 * rho_T;        % Proximal step size for temporal TV prior
temporal_iter    = 15;

eval Recon_InitialVariables;

% Regularization initialization

fprintf('Spatial Strength: tau_S/rho_S =  %d\n', tau_S/rho_S);
if temporalTV
    fprintf('Temporal Strength: tau_T/rho_T =  %d\n', tau_T/rho_T);
end

% 2D finite difference kernels for Hessian operator
Dxx = zeros(Pic_Ny, Pic_Nx); Dxx(1,1)=-2; Dxx(1,2)=1; Dxx(1,end)=1;
Dxx = fft2(Dxx);
Dyy = zeros(Pic_Ny, Pic_Nx); Dyy(1,1)=-2; Dyy(2,1)=1; Dyy(end,1)=1;
Dyy = fft2(Dyy);
Dxy = zeros(Pic_Ny, Pic_Nx); Dxy(2,2)=1/4; Dxy(2,end)=-1/4; Dxy(end,2)=-1/4; Dxy(end,end)=1/4;
Dxy = fft2(Dxy);

if ~modified

    reg_H = abs(Dxx).^2 + abs(Dyy).^2 + abs(Dxy).^2;

    % Augment system matrices with spatial regularization penalty
    sum_PTF = sum_PTF + rho_S * (reg_H + 1);
    sum_ATF = sum_ATF + rho_S * (reg_H + 1);
    Normalized_term = sum_PTF .* sum_ATF - conj_term1 .* conj_term2;

    % Precompute terms with temporal regularization penalty
    sum_PTF_wT = sum_PTF + rho_T;
    sum_ATF_wT = sum_ATF + rho_T;
    Normalized_term_wT = sum_PTF_wT .* sum_ATF_wT - conj_term1 .* conj_term2;
    
    modified = true;

end

pathPrefix = './IDT_Hessian';
if slidingWindow
    pathPrefix = [pathPrefix, '_SW/'];
else
    pathPrefix = [pathPrefix, '/'];
end

% Main Reconstruction Loop
for timeStamp = 1:timeStop

    fprintf('\nReconstructing Timestamp %d', timeStamp);

    for frame = 1:frameTotal

        if frameTotal ~= 1
            fprintf('\n-> Frame %d', frame);
        end

        % => Read files and pre-calculate
        eval Recon_ReadFiles;

        % => Initialize ADMM variables
        eval Recon_InitialADMM;

        % => Execute core reconstruction solver
        eval Recon_Solver;

        % => Output
        eval Recon_WriteResults;

        if markEnd
            markEnd = false;
            disp([newline, 'Reconstruction Complete']);
            break;
        elseif markBegin
            markBegin = false;
        end

    end
end
