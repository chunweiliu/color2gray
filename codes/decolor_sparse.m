function [out] = decolor_sparse(im)
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
[row, col, ~] = size(im);
% parameters setting
sigma       = 0.5;          % smooth input image (0.5)
k           = 50;           % threshold function (500)
minp        = 20;           % min pixel (20)
pdata       = 1;            % data fitting for sparse coding [tau](1)
psparse     = 0.05^0.5;     % regulization for sparse coding [sigma](0.05^0.5)
colorT      = 0.25;         % color threshold in RGB domain (0.3)
sAccumulateT= 0.6;          % accumulated color saliency (0.6)

%% 1. Compute Saliency
addpath('saliency/simpsal');
sali = imresize(simpsal(image), [row, col]);

%% 2. Compute Segmentation
addpath('segmentation');
[seg, nseg] = EGBIS(uint8(image), sigma, k, int32(minp));

%%

% convert data to [0,1]
if isinteger(im)
    im = im2double(im);
end

if isinteger(sali)
    sali = im2double(sali);
    if size(sali,3) == 3 % if loading a reference iamge directly
        sali = rgb2gray(sali); 
    end   
end

% decide number of segments and image dimensions
%nseg            = max(max(seg));
[row,col,nch]   = size(im);
channels        = cell(nch,1);
for n = 1:nch
    channels{n} = im(:,:,n); 
end

% compute saliency and gray information for each segment
D0          = zeros(nch,nseg);  % dictionary of superpixels
s0          = zeros(1,nseg);    % super saliency 
segSize0    = zeros(1,nseg);    % super size

for n = 1:nseg
    idx = seg == n;
    
    s0(n) = sum(sali(idx));
    segSize0(n) = sum(sum(idx));
        
    D0(1,n) = mean(channels{1}(idx));
    D0(2,n) = mean(channels{2}(idx));
    D0(3,n) = mean(channels{3}(idx));    
end

% sort saliency over segments and the corresponding gray values
[sSort, sIDX]   = sort(s0,'descend');
DSort           = D0(:,sIDX);
segSizeSort     = segSize0(sIDX);

% merge saliency over segments
DMerge          = DSort(:,1);
sMerge          = sSort(1);
segSizeMerge    = segSizeSort(1);
mIDX            = sIDX(1);

curser = 1;
sAccumulate = 0;
sTotal = sum(sSort);
%while curser < size(D0,2) && size(DMerge,2) < nColorTheme
while curser < size(D0,2) && sAccumulate < sAccumulateT
    curser = curser + 1;
    this_color = DSort(:,curser);
    
    too_similar = false;
    for n = 1:size(DMerge,2)
        mean_color_difference = norm(this_color-DMerge(:,n));
        if mean_color_difference < colorT
            too_similar = true;
            break;
        end
    end
    
    if ~too_similar
        sAccumulate = sAccumulate + sSort(curser)/sTotal;
        
        mIDX = [mIDX sIDX(curser)];
        DMerge = [DMerge DSort(:,curser)];
        sMerge = [sMerge sSort(curser)];
        segSizeMerge = [segSizeMerge segSizeSort(curser)];           
    end
end

%
% choose atoms for dictionary
%
% prepare small size image pixels
sm = imresize(im,0.1);
x = zeros(3,size(sm,1)*size(sm,2));
for n = 1:nch
    x(n,:) = reshape(sm(:,:,n),1,size(sm,1)*size(sm,2));
end

D = DMerge;

rankD           = -1;
currentRankD    = rank(D);
while currentRankD ~= rankD && currentRankD >= 3
    
    rankD = currentRankD;
    a = zeros(size(D,2),size(x,2));
    
    % project x to a
    for n = 1:size(x,2)
        sols = SolveLasso(D, x(:,n), size(D,2));       
        if sum(isnan(sols)) ~= 0, sols = zeros(size(D,2),1); end
        a(:,n) = sols;
    end
    
    % scatter matrix 
    aMean = mean(a,2);
    aScatter = zeros(size(a,1),size(a,1));
    for i = 1:size(a,2)
        aScatter = aScatter + (a(:,i)-aMean)*(a(:,i)-aMean)';
    end
    [evector,evalue] = eigs(aScatter);
    %evalue
    fprintf('largest evalue divided by total in a = %3.2f\n', ...
        evalue(1)/sum(diag(evalue)));
    
    % update dictionary
    currentRankD = rank(aScatter);
    D = DMerge(:, 1:currentRankD);
    
end

%
% find a linear projection L = L*R* preserve good structure in high dimension
%
% find L*
% Gkioulekas and Zickler NIPS 2011
[evector,evalue]    = eigs(D*D');
evector_m           = evector(:,1); % columns are the corresponding eigenvectors.
evalue_m            = evalue(1,1);

L = diag(foo(evalue_m,pdata,psparse))*evector_m';

% find R*
% fit a rotation R for the linear projection L
s = sMerge(:,1:size(D,2));

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

% dimensionality reduction
x = zeros(nch,row*col);
for n = 1:nch
    x(n,:) = channels{n}(:);
end
y  = L*R*x;

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
out = reshape(out,row,col);

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