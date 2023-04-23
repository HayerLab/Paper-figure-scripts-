function [im_ch1, im_ch2, im_ch3,im_ch4, im_ch5, im_ch6] = nd2read_hayer_6chan(filename,finfo,counter,varargin)
tic
%finfo = nd2finfo(filename);
%disp(['analyzing file structure used ', sprintf('%0.2f', toc), ' seconds'])

c = numel(num2str(counter)); 

im_ch1 = zeros(finfo.img_width, finfo.img_height, 'uint16');
im_ch2 = zeros(finfo.img_width, finfo.img_height, 'uint16');
im_ch3 = zeros(finfo.img_width, finfo.img_height, 'uint16');
if finfo.ch_count == 6
  im_ch4 = zeros(finfo.img_width, finfo.img_height, 'uint16');
   im_ch5 = zeros(finfo.img_width, finfo.img_height, 'uint16');
    im_ch6 = zeros(finfo.img_width, finfo.img_height, 'uint16');
end

fid = fopen(filename, 'r');
fseek(fid, finfo.file_structure(strncmp(strcat('ImageDataSeq|',num2str(counter),'!'), ...
  {finfo.file_structure(:).nameAttribute}, 14+c)).dataStartPos, 'bof');

tic
% Image extracted from ND2 has image width defined by its first dimension.
if finfo.padding_style == 1
  if finfo.ch_count == 6
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
       im_ch3(:, ii) = temp(1, :);
        im_ch1(:, ii) = temp(2, :);
        im_ch2(:, ii) = temp(3, :);
        im_ch4(:, ii) = temp(4, :);
         im_ch5(:, ii) = temp(5, :);
          im_ch6(:, ii) = temp(6, :);
        fseek(fid, 2, 'cof');
    end
  else
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(1, :);
        im_ch2(:, ii) = temp(2, :);
       im_ch3(:, ii) = temp(3, :);
        fseek(fid, 2, 'cof');
    end
  end
else
  if finfo.ch_count == 6
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(1, :);
        im_ch2(:, ii) = temp(2, :);
        im_ch3(:, ii) = temp(3, :);
        im_ch4(:, ii) = temp(4, :);
         im_ch5(:, ii) = temp(5, :);
          im_ch6(:, ii) = temp(6, :);
    end
  else
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(1, :);
        im_ch2(:, ii) = temp(2, :);
        im_ch3(:, ii) = temp(3, :);
    end 
  end
end

fclose(fid);

im_ch1 = permute(im_ch1, [2 1]);
im_ch2 = permute(im_ch2, [2 1]);
im_ch3 = permute(im_ch3, [2 1]);
if finfo.ch_count == 6
    im_ch4 = permute(im_ch4, [2 1]);
    im_ch5 = permute(im_ch5, [2 1]);
    im_ch6 = permute(im_ch6, [2 1]);
end
if any(strcmpi(varargin, 'use_ch4'))
  im_ch3 = im_ch4;
end
  

%disp(['reading complete image data used ', sprintf('%0.2f', toc), ' seconds'])
end