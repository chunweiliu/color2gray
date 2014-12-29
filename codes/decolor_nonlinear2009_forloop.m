function tones = decolor_nonlinear2009_forloop(image)
%DECOLOR_NONLINEAR2009 (implementation of a paper in Siggraph Asia 2009)
% g(x, y) = L + f(theta) * C (where L, f, c are three channels of LCH color
% space)
% input
%   - image (n x m x 3 double)
% output
%   - tones (n x m double)

% Predefine arguments

% Processing
% 0. rgb2lab
imlab = rgb2lab(image);
imluv = rgb2luv(image);
imlch = rgb2lch(image);

% 1. Compute color differences G
% For loop version
[row, col, ~] = size(image);
Gx = zeros(row, col);
Gy = zeros(row, col);
for n = 2:row - 1
    for m = 2:col - 1        
        Gx(n, m) = colordiff(imlab(n + 1, m, :), imlab(n - 1, m, :), imluv(n + 1, m, :), imluv(n - 1, m, :));
        Gy(n, m) = colordiff(imlab(n, m + 1, :), imlab(n, m - 1, :), imluv(n, m + 1, :), imluv(n, m - 1, :));
    end
end

% 2. Optimization
% Compute M_s, b_s
L = imlch(:, :, 1);
C = imlch(:, :, 2);
H = imlch(:, :, 3);
T = zeros(row, col, 9);
for n = 1:row
    for m = 1:col        
        T(n, m, :) = C(n, m) .* ...
                     [cos(H(n, m)), cos(2 * H(n, m)), cos(3 * H(n, m)), cos(4 * H(n, m)), ...
                      sin(H(n, m)), sin(2 * H(n, m)), sin(3 * H(n, m)), sin(4 * H(n, m)), 1];
    end
end
[U, V] = gradient(T);
[Lx, Ly] = gradient(L);

M_s = zeros(9, 9);
b_s = zeros(9, 1);
for n = 2:row - 1
    for m = 2:col - 1
        % Compute u, v, p, q
        p = Gx(n, m) - Lx(n, m);
        q = Gy(n, m) - Ly(n, m);
        u = reshape(U(n, m, :), 9, []);
        v = reshape(V(n, m, :), 9, []);
        
        b_s = b_s + (p * u + q * v);
        M_s = M_s + (u * u' + v * v');        
    end
end

% Solve the energy function based on M_s and b_s
x = (M_s + row .* col .* eye(9)) \ b_s;

% 3. Render gray-scale image based on L, C and x (f(theta))
tones = zeros(row, col);
for n = 1:row
    for m = 1:col
        f = reshape(T(n, m, :), 1, 9) * x;
        tones(n, m) = L(n, m) + f * C(n, m);
    end
end
tones = mat2gray(tones);


end

function g = colordiff(labi, labj, luvi, luvj, alpha)
%COLORDIFF color difference of two vector in Lab color space
% input
%   labi, labj: 1 x 3 double
%   laui, lauj: 1 x 3 double
% output
%   g: 1 x 1 double, color difference 
if nargin < 5
    alpha = 0.1;
end
deltal = labi(1) - labj(1);
deltal2 = (labi(1) - labj(1)) .^ 2;
deltaa2 = (labi(2) - labj(2)) .^ 2;
deltab2 = (labi(3) - labj(3)) .^ 2;
r = 3.59210244843;  % 2.54 * 1.41421356237;
g = (deltal2 + (alpha .* ((deltaa2 + deltab2) .^ 0.5) ./ r) .^ 2) .^ 0.5;

lhki = l2lhk(luvi(1), luvi(2), luvi(3));
lhkj = l2lhk(luvj(1), luvj(2), luvj(3));
deltalhk = lhki - lhkj;
if deltalhk == 0
    if deltal == 0
        signg = sign(deltal2 .^ 1.5 + deltaa2 .^ 1.5 + deltab2 .^ 1.5);
    else
        signg = sign(deltal);
    end
else
    signg = sign(deltalhk);
end

g = signg * g;

end

function lhk = l2lhk(l, u, v)
% adaptlum = 20.;
% adaptlump = adaptlum .^ 0.4495;
% whiteu = 0.20917;
% whitev = 0.48810;

% Adapting luminance dependency
kbr = 0.2717 .* (6.469 + 6.362 * 20 .^ 0.4495) / (6.469 + 20 .^ 0.4495);
% kbr = 0.8147419482;

% Saturation
suv = 13. * ((u - 0.20917) .^ 2 + (v - 0.48810) .^ 2) .^ 0.5;

% Theta
theta = atan((u - 0.20917) ./ (v - 0.44810));

qtheta = - 0.01585 - 0.03017 * cos(theta) - 0.04556 * cos(2 * theta)...
         - 0.02677 * cos(3 * theta) - 0.00295 * cos(4 * theta)...
         + 0.14592 * sin(theta) + 0.05084 * sin(2 * theta)...
         - 0.01900 * sin(3 * theta) - 0.00764 * sin(4 * theta);

%lhk = l + (-0.1340 * qtheta + 0.0872 * 0.8147419482) * suv * l;
lhk = l + (-0.1340 * qtheta + 0.0872 * kbr) * suv * l;

end