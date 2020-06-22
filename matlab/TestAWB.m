% Test AWB
clear; clc;
close all;

I = imread('image4.bmp');
% J = imread('TestComparison.jpg');

% subplot(2,2,1), imshow(I);
% subplot(2,2,2), imshow(I(:,:,1));
% subplot(2,2,3), imshow(I(:,:,2));
% subplot(2,2,4), imshow(I(:,:,3));

% Calculate Illumination estimation 
% p = 93:0.5:94;

% for i=1:length(p)
% profile on;
    O = PerformAWB(I, 90);
% profile viewer;
    
%     subplot(1,2,1), imshow(I)
%     subplot(1,2,2), imshow(O);
  figure, montage({I, O},'Size',[1 2]);
    
%     imwrite(O, strcat('AWB_', num2str(p(i)), '.png'));
       
% end




