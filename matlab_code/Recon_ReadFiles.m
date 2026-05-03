globalSlide = (timeStamp-1) * frameTotal + (frame-1);

if globalSlide == 0
    % Initialize intensity buffer
    I_buf = zeros(Select_RangeY, Select_RangeX, LED_Num, 'double');

    for i = 1:LED_Num
        str   = [addpath1,'\1_',num2str(i,'%05d'),'.tif'];
        strBG = [addpath2,'\1_',num2str(i,'%05d'),'.tif'];

        Iraw  = im2double(imread(str));
        BGraw = im2double(imread(strBG));

        Icrop  = Iraw(roiRows, roiCols);
        BGcrop = BGraw(roiRows, roiCols);

        assert(isequal(size(Icrop),  [Select_RangeY, Select_RangeX]), 'ROI crop size mismatch for Iraw.');
        assert(isequal(size(BGcrop), [Select_RangeY, Select_RangeX]), 'ROI crop size mismatch for BGraw.');

        I_buf(:,:,i) = Icrop ./ BGcrop;
        % I_buf(:,:,i) = Icrop;
    end
else
    % Sliding window update
    repStart = mod(((globalSlide-1) * inter), LED_Num);

    for k = 1:inter
        gIdx  = globalSlide * inter + LED_Num - inter + k;
        bgIdx = mod(gIdx-1, LED_Num) + 1;

        str   = [addpath1,'\1_',num2str(gIdx,'%05d'),'.tif'];
        strBG = [addpath2,'\1_',num2str(bgIdx,'%05d'),'.tif'];

        Iraw  = im2double(imread(str));
        BGraw = im2double(imread(strBG));

        Icrop  = Iraw(roiRows, roiCols);
        BGcrop = BGraw(roiRows, roiCols);

        assert(isequal(size(Icrop),  [Select_RangeY, Select_RangeX]), 'ROI crop size mismatch for Iraw.');
        assert(isequal(size(BGcrop), [Select_RangeY, Select_RangeX]), 'ROI crop size mismatch for BGraw.');

        I_buf(:,:,repStart+k) = Icrop ./ BGcrop;
    end

    if timeStamp == timeStop
        markEnd = true;
    end
end

Intensity     = zeros(Select_RangeY, Select_RangeX, LED_Num);
conj_PTF_Iten = zeros(Pic_Ny, Pic_Nx, Pic_Nz);
conj_ATF_Iten = zeros(Pic_Ny, Pic_Nx, Pic_Nz);

for i = 1:LED_Num
    IntTemp = I_buf(:,:,i);
    % IntTemp = IntTemp ./ mean(IntTemp(:));
    Intensity(:,:,i) = IntTemp - 1;

    F2_I = fft2(Intensity(:,:,i));
    conj_PTF_Iten = conj_PTF_Iten + conj_PTF_4D_Used(:,:,:,i) .* F2_I;
    conj_ATF_Iten = conj_ATF_Iten + conj_ATF_4D_Used(:,:,:,i) .* F2_I;
end