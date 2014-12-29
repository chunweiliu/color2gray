function LCH = rgb2lch(image)
%RGB2LCH Summary of this function goes here
%   LCH = (lightness, chroma, hue angle)

image = double(image) ./ 255.;
LCH = colorspace('LCH', image);

end

