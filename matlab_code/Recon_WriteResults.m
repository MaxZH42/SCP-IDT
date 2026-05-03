outputRI = real(RI);

output = uint16((outputRI - 1.3200) ./ (1.4200 - 1.3200) * 65535);

if ~exist(pathPrefix, 'dir')
    mkdir(pathPrefix);
end

for writeIndex = 1:1:Pic_Nz
    if ~slidingWindow
        writePath = [pathPrefix,'Z',sprintf('%02d',writeIndex),'T',sprintf('%06d',timeStamp),'.tif'];
    else
        writePath = [pathPrefix,'Z',sprintf('%02d',writeIndex),'T',sprintf('%06d',globalSlide+1),'.tif'];
    end
    imwrite(output(:,:,writeIndex), writePath);
end