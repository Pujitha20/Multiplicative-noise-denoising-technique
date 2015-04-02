clc;
clear all;

im = imread('C:\Users\AK PUJITHA\Desktop\10805256_1714587932099999_1523826375_n.jpg');
grayImage = rgb2gray(im);
mask = grayImage >300; % Mask is bright stuff.
% Fill in the body
mask = imfill(mask, 'holes'); % Mask is whole solid body.
% OR it in with the zeros
mask = mask | (grayImage == 0); % Mask now includes pure zeros.
% Extract pixels that are not masked
darkNonZeroOutsidePixels = grayImage(~mask);
