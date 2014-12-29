function A = transcolor(image)
%TRANSCOLOR Summary of this function goes here
%   Detailed explanation goes here

[row, col, ~] = size(image);
image = im2double(image);

L = [1; 0; 0];
R = eye(3);
X = reshape(image(:), row * col, 3);

A = reshape(X * R * L, row, col);

subplot(1, 2, 1);
imshow(A)
subplot(1, 2, 2);
imshow(rgb2gray(image));



end

