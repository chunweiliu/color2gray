function tones = decolor(image)
%DECOLOR_NONLINEAR2009 (implementation of a paper in Siggraph Asia 2009)
% g(x, y) = L + f(theta) * C (where L, f, c are three channels of LCH color
% space)
% input
%   - image (n x m x 3 double)
% output
%   - tones (n x m double)

%% 0. Preprocessing
image = im2double(image);
image = image ./ 255;
[row, col, ~] = size(image);

% Define user parameters
lambda = row * col;
alpha = 1e0;

% Convert image to other spaces
addpath('colorspace');
imlab = rgb2lab(image);
imluv = rgb2luv(image);
imlch = rgb2lch(image);

%% 1. Compute color differences G
% Nayatani's chromatic lightness model
% adaptlum = 20.;
% adaptlump = adaptlum .^ 0.4495; %  0.8147419482
% whiteu = 0.20917;
% whitev = 0.48810;
% kbr = 0.8147419482;
theta = atan2((imluv(:, :, 3) - 0.48810), (imluv(:, :, 2) - 0.20917));  % atan2 but not atan!
qtheta = - 0.01585 - 0.03017 * cos(theta) - 0.04556 * cos(2 * theta)...
         - 0.02677 * cos(3 * theta) - 0.00295 * cos(4 * theta)...
         + 0.14592 * sin(theta) + 0.05084 * sin(2 * theta)...
         - 0.01900 * sin(3 * theta) - 0.00764 * sin(4 * theta);
suv = 13. * ((imluv(:, :, 2) - 0.20917) .^ 2 + (imluv(:, :, 3) - 0.48810) .^ 2) .^ 0.5;
LHK = imluv(:, :, 1) + (-0.1340 .* qtheta + 0.0872 * 0.8147419482) .* suv .* imluv(:, :, 1);

% Compute difference on x and y direction
% First, compute the sign of g.
deltalhk_x = circshift(LHK, [0 -1]) - circshift(LHK, [0 1]);
deltalhk_y = circshift(LHK, [-1 0]) - circshift(LHK, [1 0]);
deltalab_x = circshift(imlab, [0 -1]) - circshift(imlab, [0 1]);  % c(x+1, y) - c(x-1, y)
deltalab_y = circshift(imlab, [-1 0]) - circshift(imlab, [1 0]);  % c(x, y+1) - c(x, y-1)

deltal_x = deltalab_x(:, :, 1);
deltal_y = deltalab_y(:, :, 1);
deltalab3_x = deltalab_x(:, :, 1) .^ 3 + deltalab_x(:, :, 2) .^ 3 + deltalab_x(:, :, 3) .^ 3;
deltalab3_y = deltalab_y(:, :, 1) .^ 3 + deltalab_y(:, :, 2) .^ 3 + deltalab_y(:, :, 3) .^ 3;
            
% If sign deltaLHK = 0, then try deltaL, otherwise deltaLab3
signg_x = sign(deltalhk_x);
idx = signg_x == 0;
signg_x(idx) = sign(deltal_x(idx));
idx = signg_x == 0;
signg_x(idx) = sign(deltalab3_x(idx));
signg_y = sign(deltalhk_y);
idx = signg_y == 0;
signg_y(idx) = sign(deltal_y(idx));
idx = signg_y == 0;
signg_y(idx) = sign(deltalab3_y(idx));

% Second, compute the g
% r = 3.59210244843;  % 2.54 * 1.41421356237;
gx = (deltalab_x(:, :, 1) .^ 2 + (alpha .* ((deltalab_x(:, :, 2) .^ 2 + deltalab_x(:, :, 3) .^ 2) .^ 0.5) ./ 3.59210244843) .^ 2) .^ 0.5;
gy = (deltalab_y(:, :, 1) .^ 2 + (alpha .* ((deltalab_y(:, :, 2) .^ 2 + deltalab_y(:, :, 3) .^ 2) .^ 0.5) ./ 3.59210244843) .^ 2) .^ 0.5;

% To sum up, G become:
Gx = signg_x .* gx;
Gy = signg_y .* gy;

%% 2. Optimization
% Compute M_s, b_s
L = imlch(:, :, 1);
C = imlch(:, :, 2);
H = imlch(:, :, 3);
T = zeros(row, col, 9);
for n = 1:4
    T(:, :, n) = C .* cos(n .* H);
    T(:, :, n + 4) = C .* sin(n .* H);
end
T(:, :, 9) = C;

[U, V] = gradient(T);
[Lx, Ly] = gradient(L);

% Resterilization
p = (Gx(:) - Lx(:))';
q = (Gy(:) - Ly(:))';
u = reshape(shiftdim(U, 2), 9, []);  % Move the 9 dim in the first dim and resterilize the matrix
v = reshape(shiftdim(V, 2), 9, []);

b_s = sum(bsxfun(@times, u, p) + bsxfun(@times, v, q), 2);
M_s = u * u' + v * v';

% Solve the energy function based on M_s and b_s
x = (M_s + lambda .* eye(9)) \ b_s;

%% 3. Render gray-scale image based on L, C and x (f(theta))
ftheta = reshape(x' * reshape(shiftdim(T, 2), 9, []), row, col);
tones = mat2gray(L + ftheta .* C);

%% 4. Add saliency
addpath('simpsal');
map = imresize(simpsal(image), [row, col]);
tones = mat2gray(tones + map);

end