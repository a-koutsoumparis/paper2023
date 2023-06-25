clc;
close all;
clear;

% Import images from folder
imagefiles = dir('*.tif');

% Put all images into one variable "images"
for i=1:size(imagefiles,1)
currentfilename = imagefiles(i).name;
currentimage = imread(currentfilename);
a{i,1} = currentimage;
end

% Find the maximum value and the index
for i=1:size(a,1)
[maxval(i),idx(i)]= max(a{i}(:));
end

% Convert the index into x-y coordinates
for i=1:size(idx,2)
[row,col]= ind2sub(size(a{i}), idx);
end

% Take the 12x12 (RIS size) values around the maximum
for i=1:size(row,2)
x(i) = (col(i)-5);
y(i) = (row(i)-5);
s{i} = [x(i) y(i) 11 11];
end

% Isolate the signal
for i=1:size(a,1)
I{i} = imcrop(a{i},s{i});
end

% Background is the pixel with the least intensity
for i=1:size(I,2)
minval{i} = min(min(I{i}));
end

% Subtract background
for i=1:size(I,2)
ImSub{i} = I{i} - minval{i};
end

% Average RIS intensity
for i=1:size(ImSub,2)
m(i) = mean2(ImSub{i});
end

% Calculate speed for 5 second time interval
for i=2:size(row,2)
speed(i) = (sqrt(((row(i) - row(i-1)).^2)+((col(i)-col(i-1)).^2)))./5;
end