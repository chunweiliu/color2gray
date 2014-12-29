function tones = decolor_sparse(image)
    %sparse_dr(im,seg,sali,pdata,psparse,colorT, sAccumulateT)
%SPARSE_DR Dimensionality reduction using sparse linear model
% input:
%   im          row x col x 3   uint8   | input color image
%   seg         row x col       int32   | EGBIS output segmentation
%   sali        row x col       double  | saliency map
%   pdata       1 x 1           double  | fit data
%   psparse     1 x 1           double  | fit sparse
%   colorT      1 x 1           double  | threshold of color disimilar
%   sAccumlate  1 x 1           double  | threshold of accumlate saliency
%
% output:
%   out         row x col       uint8   | output gray image
%   L           1 x 3           double  | output linear projection
%   R           3 x 3           double  | output rotation matrix
% others for visualization
%   D           3 x k           double
%   

%% 0. Preprocessing
%image = im2double(image);
%image = image ./ 255;
[row, col, ~] = size(image);
%imlab = rgb2lab(image);

% Define user parameters
% Segmentation
sigma = 0.5;
k = 500;
minpixel = 20;
% Optimization of L
pdata = 1;              % data fitting for sparse coding [tau](1)
psparse = 0.05^0.5;     % regulization for sparse coding [sigma](0.05^0.5)
% Color difference
tol_color = 45;

%% 1. Compute Saliency
addpath('saliency/simpsal');
salmap = imresize(simpsal(image), [row, col]);

%% 2. Compute Segmentation
addpath('segmentation');
[segmap, nseg] = EGBIS(image, sigma, k, int32(minpixel));

%% 3. Find the dictionary D
% Find the saliency and mean color
imlab = rgb2lab(image);
L = imlab(:, :, 1);
A = imlab(:, :, 2);
B = imlab(:, :, 3);
C = zeros(3, nseg);
s = zeros(1, nseg);
for n = 1:nseg
    idx = segmap == n;
    
    C(1, n) = mean(L(idx));
    C(2, n) = mean(A(idx));
    C(3, n) = mean(B(idx));
    
    s(n) = sum(salmap(idx));
end

% Merge the dictionary by color (unique with torrence 20)
addpath('uniquetol');
[C_merged, C_idx, ~] = uniquetol(C', tol_color, 'rows');
C_merged = C_merged';
s_merged = s(C_idx);

% Sort the dictionary by saliency
[s_sorted, s_idx] = sort(s_merged, 'descend');
D_sorted = C_merged(:, s_idx);

D = D_sorted;
s = s_sorted;

%% 4. Find the linear projection L
% find a linear projection L = L*R* preserve good structure in high dimension
%
% find L*
% Gkioulekas and Zickler NIPS 2011
[evector, evalue] = eigs(D * D');
evector_m = evector(:, 1); % columns are the corresponding eigenvectors.
evalue_m = evalue(1, 1);

L = diag(foo(evalue_m, pdata, psparse)) * evector_m';

%% 5. Find the linear projection R
% find R*
% fit a rotation R for the linear projection L
%s = sMerge(:,1:size(D,2));

k       = size(D,2);
delta_x = zeros(3,k^2);
delta_s = zeros(1,k^2);
for n = 1:k
    for m = 1:k
        delta_x(:,(n-1)*k +m) = D(:,n) - D(:,m); %D have been sorted        
        delta_s(:,(n-1)*k +m) = s(n);
    end
end

x1 = delta_x;
d1 = delta_s;

%
% Generate nMax 3d random matrices and see which is the best
%
nMax   =  10000;
vpower =  2;
fMin   =  Inf;

%rng('default');
rand('seed',0);
r  = rand(4, 1);
R  = rotationmat3D(r(1),[r(2),r(3),r(4)]);

v  = L*R*x1;
%fMax = (v.^vpower)*d1';
fMax = ((v.^vpower)*d1')/size(d1,2);

% Reset the random number generator to default
randMatrix = rand(nMax, 4, 1);

for iter = 1:nMax
    r = randMatrix(iter,:,:);
    RTmp = rotationmat3D(r(1),[r(2),r(3),r(4)]);
    v = L*RTmp*x1;    
    %f = (v.^vpower)*d1';
    f = ((v.^vpower)*d1')/size(d1,2);
    
    if f < fMin
        fMin = f;
    end
    
    if f > fMax
        fMax = f;        
        R = RTmp;
    end        
end

%% 6.
% dimensionality reduction
x = reshape(shiftdim(double(image) ./ 255., 2), 3, []);
% x = zeros(3,row*col);
% for n = 1:3
%     x(n,:) = channels{n}(:);
% end

%L = [1, 0, 0];
%R = eye(3);
y = L*R*x;
%y = min(max(y,0),1);

% a = x;
% aMean = mean(a,2);
% aScatter = zeros(size(a,1),size(a,1));
% for i = 1:size(a,2)
%     aScatter = aScatter + (a(:,i)-aMean)*(a(:,i)-aMean)';
% end
% [evector,evalue] = eigs(aScatter);
% %evalue
% fprintf('largest evalue divided by total in x = %3.2f\n', ...
%     evalue(1)/sum(diag(evalue)));

% map the output to gray scale display range, and correct domain
out = mat2gray(y);
tones = reshape(out,row,col);

end

function out = foo(n, tau, sigma)
out = ((4 * tau^4 * n)/...
       (sigma^4 + 4*tau^2*sigma^2*n + 4*tau^4*n^2))^0.5;
end

function R= rotationmat3D(r,Axis)
%function R= rotationmat3D(radians,Axis)
%
% creates a rotation matrix such that R * x 
% operates on x by rotating x around the origin r radians around line
% connecting the origin to the point "Axis"
%
% example:
% rotate around a random direction a random amount and then back
% the result should be an Identity matrix
%
%r = rand(4,1);
%rotationmat3D(r(1),[r(2),r(3),r(4)]) * rotationmat3D(-r(1),[r(2),r(3),r(4)])
%
% example2: 
% rotate around z axis 45 degrees
% Rtest = rotationmat3D(pi/4,[0 0 1])
%
%Bileschi 2009

if nargin == 1
   if(length(rotX) == 3)
      rotY = rotX(2);
      rotZ = rotZ(3);
      rotX = rotX(1);
   end
end

% useful intermediates
L = norm(Axis);
if (L < eps)
   error('axis direction must be non-zero vector');
end
Axis = Axis / L;
L = 1;
u = Axis(1);
v = Axis(2);
w = Axis(3);
u2 = u^2;
v2 = v^2;
w2 = w^2;
c = cos(r);
s = sin(r);
%storage
R = nan(3);
%fill
R(1,1) =  u2 + (v2 + w2)*c;
R(1,2) = u*v*(1-c) - w*s;
R(1,3) = u*w*(1-c) + v*s;
R(2,1) = u*v*(1-c) + w*s;
R(2,2) = v2 + (u2+w2)*c;
R(2,3) = v*w*(1-c) - u*s;
R(3,1) = u*w*(1-c) - v*s;
R(3,2) = v*w*(1-c)+u*s;
R(3,3) = w2 + (u2+v2)*c;
end