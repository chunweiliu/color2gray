function [tones,recolor] = decolor_linear2006(picture,effect,scale,noise)
% Created: 4 June 2005
% Version: 8 September 2005
% Copyright 2005, Mark Grundland. All rights resaved. 
% Licensed for noncommercial, academic use only.

% Developed by Mark Grundland.
% Computer Laboratory, University of Cambridge, UK.
% Technical Report UCAM-CL-TR-649, University of Cambridge.
% Web page: http://www.eyemaginary.com/Portfolio/Publications.html
% Any questions? Contact mark @ eyemaginary.com .

% Usage:   Coverts RGB picture to grayscale tones,
%             while endeavoring to express color contrasts as tone changes.
%          Can also recolor the original picture 
%             to incorporate the enhanced grayscale tones.
% Options: effect specifies how much the picture's achromatic content 
%             should be altered to accommodate the chromatic contrasts
%             (default effect = 0.5 )
%          scale in pixels is the typical size of 
%             relevant color contrast features
%             (default scale = sqrt(2*min(size(picture.dimensions)) )
%          noise quantile indicates the amount of noise in the picture
%             enabling the dynamic range of the tones to be appropriately scaled
%             (default noise = 0.001)
% Assumes: valid RGB 0 <= picture <= 1
%          positive 0 <= effect <= 1  
%          positive 1 <= scale << min(picture.dimensions)
%          small quantile 0 <= noise << 0.5 
% Applies: Standard NTSC color to grayscale conversion
%             is used as the default achromatic color channel.
% Example: tones=decolorize(picture,0.5,25,0.001)


% Examine inputs
frame=[size(picture,1), size(picture,2)];
pixels=frame(1)*frame(2);
if nargin<2 || isempty(effect)
    effect=0.5;
end;
if nargin<3 || isempty(scale)
    scale=sqrt(2*min(frame));
end;
if nargin<4 || isempty(noise)
    noise=0.001;
end;

% Reset the random number generator
randn('state',0);
tolerance=100*eps;

% Define the YPQ color space
colorconvert=[0.2989360212937753847527155, 0.5870430744511212909351327,  0.1140209042551033243121518;
              0.5,                         0.5,                         -1;
              1,                          -1,                            0]';
colorrevert= [1,                           0.1140209042551033243121518,  0.6440535265786729530912086;
              1,                           0.1140209042551033243121518, -0.3559464734213270469087914;
              1,                          -0.8859790957448966756878482,  0.1440535265786729530912086]';
colorspan=[ 0,  1;
           -1,  1;
           -1,  1];          
maxluminance=1;
scaleluminance=0.66856793424088827189;
maxsaturation=1.1180339887498948482; 
alter=effect*(maxluminance/maxsaturation);

% Covert picture to the YPQ color space
picture=reshape(picture,[pixels,3]);
image=picture*colorconvert;
original=image;
chroma=sqrt(image(:,2).*image(:,2)+image(:,3).*image(:,3));

% Pair each pixel with a randomly chosen sample site
mesh=reshape(cat(3,repmat((1:frame(1))',[1,frame(2)]),repmat((1:frame(2)),[frame(1),1])),[pixels,2]);
displace=(scale*sqrt(2/pi))*randn(pixels,2);
look=round(mesh+displace);
redo=find((look(:,1)<1));
look(redo,1)=2-rem(look(redo,1),frame(1)-1);
redo=find((look(:,2)<2));
look(redo,2)=2-rem(look(redo,2),frame(2)-1);
redo=find((look(:,1)>frame(1)));
look(redo,1)=frame(1)-1-rem(look(redo,1)-2,frame(1)-1);
redo=find((look(:,2)>frame(2)));
look(redo,2)=frame(2)-1-rem(look(redo,2)-2,frame(2)-1);
look=look(:,1)+frame(1)*(look(:,2)-1);

% Calculate the color differences between the paired pixels
delta=image-image(look,:);
contrastchange=abs(delta(:,1));
contrastdirection=sign(delta(:,1));
colordifference=picture-picture(look,:);
colordifference=sqrt(sum(colordifference.*colordifference,2))+eps;

% Derive a chromatic axis from the weighted sum of chromatic differences between paired pixels
weight=1-((contrastchange/scaleluminance)./colordifference);
weight(find(colordifference<tolerance))=0;
axis=weight.*contrastdirection;
axis=delta(:,2:3).*[axis, axis];
axis=sum(axis,1);

% Project the chromatic content of the picture onto the chromatic axis
projection=image(:,2)*axis(1)+image(:,3)*axis(2);
projection=projection/(quantiles(abs(projection),1-noise)+tolerance);

% Combine the achromatic tones with the projected chromatic colors and adjust the dynamic range
image(:,1)=image(:,1)+effect*projection;
imagerange=quantiles(image(:,1),[noise; 1-noise]);
image(:,1)=(image(:,1)-imagerange(1))/(imagerange(2)-imagerange(1)+tolerance);
targetrange=effect*[0; maxluminance]+(1-effect)*quantiles(original(:,1),[noise; 1-noise]);
image(:,1)=targetrange(1)+(image(:,1)*(targetrange(2)-targetrange(1)+tolerance));
image(:,1)=min(max(image(:,1),original(:,1)-alter.*chroma),original(:,1)+alter.*chroma);
image(:,1)=min(max(image(:,1),0),maxluminance);

% Return the results
tones=image(:,1)/maxluminance;
tones=reshape(tones,frame);
if nargout>1
    recolor=image*colorrevert;
    recolor=cat(3,reshape(recolor(:,1),frame),reshape(recolor(:,2),frame),reshape(recolor(:,3),frame));
    recolor=min(max(recolor,0),1);
end;








function r = quantiles(x,q);
% Usage:   Finds quantiles q of data x
% Assumes: quantiles 0 <= q <= 1
% Example: r=quantiles(x,[0; 0.5; 1]);
%           xmin=r(1); xmedian=r(2); xmax=r(3)

tolerance=100*eps;
q=q(:);
x=sort(x);
n=size(x,1);
k=size(x,2);
e=1/(2*n);
q=max(e,min(1-e,q));
q=n*q+0.5;
p=min(n-1,floor(q));
q=q-p;
q(find(tolerance>q))=0;
q(find(q>(1-tolerance)))=1; 
q=repmat(q,[1,k]);
r=(1-q).*x(p,:)+q.*x(p+1,:);
