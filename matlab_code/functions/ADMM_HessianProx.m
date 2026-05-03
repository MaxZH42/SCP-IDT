H_k = zeros(Pic_Ny,Pic_Nx,Pic_Nz, 6, 'like', SP1_k);

% Real part Hessian components (Channels 1-3)
H_k(:,:,:,1) = circshift(SP1_k(:,:,:,1),[0,1,0]) - 2*SP1_k(:,:,:,1) + circshift(SP1_k(:,:,:,1),[0,-1,0]); % xx
H_k(:,:,:,2) = circshift(SP1_k(:,:,:,1),[1,0,0]) - 2*SP1_k(:,:,:,1) + circshift(SP1_k(:,:,:,1),[-1,0,0]); % yy
H_k(:,:,:,3) = ( circshift(SP1_k(:,:,:,1),[1,1,0]) - circshift(SP1_k(:,:,:,1),[1,-1,0]) - circshift(SP1_k(:,:,:,1),[-1,1,0]) + circshift(SP1_k(:,:,:,1),[-1,-1,0]) ) / 4;

% Imaginary part Hessian components (Channels 4-6)
H_k(:,:,:,4) = circshift(SP1_k(:,:,:,2),[0,1,0]) - 2*SP1_k(:,:,:,2) + circshift(SP1_k(:,:,:,2),[0,-1,0]); % xx
H_k(:,:,:,5) = circshift(SP1_k(:,:,:,2),[1,0,0]) - 2*SP1_k(:,:,:,2) + circshift(SP1_k(:,:,:,2),[-1,0,0]); % yy
H_k(:,:,:,6) = ( circshift(SP1_k(:,:,:,2),[1,1,0]) - circshift(SP1_k(:,:,:,2),[1,-1,0]) - circshift(SP1_k(:,:,:,2),[-1,1,0]) + circshift(SP1_k(:,:,:,2),[-1,-1,0]) ) / 4;

% L1 Proximal step (Soft-thresholding)

V = H_k - y1_k;
DSP_k = sign(V) .* max(abs(V) - tau_S/rho_S, 0);

delta = DSP_k - H_k;

% L1 & L2 Proximal step (Soft-thresholding)
%
% V = H_k - y1_k;
% DSP_k = zeros(size(V), 'like', V);
% 
% Vr = V(:,:,:,1:3);
% nr = sqrt(abs(V(:,:,:,1)).^2 + abs(V(:,:,:,2)).^2 + abs(V(:,:,:,3)).^2);
% scale_r = max(0, 1 - (tau_S/rho_S) ./ (nr + eps));
% DSP_k(:,:,:,1) = scale_r .* V(:,:,:,1);
% DSP_k(:,:,:,2) = scale_r .* V(:,:,:,2);
% DSP_k(:,:,:,3) = scale_r .* V(:,:,:,3);
% 
% ni = sqrt(abs(V(:,:,:,4)).^2 + abs(V(:,:,:,5)).^2 + abs(V(:,:,:,6)).^2);
% scale_i = max(0, 1 - (tau_S/rho_S) ./ (ni + eps));
% DSP_k(:,:,:,4) = scale_i .* V(:,:,:,4);
% DSP_k(:,:,:,5) = scale_i .* V(:,:,:,5);
% DSP_k(:,:,:,6) = scale_i .* V(:,:,:,6);
% 
% delta = DSP_k - H_k;
