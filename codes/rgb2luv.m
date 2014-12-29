function Luv_back = rgb2luv(Im_rgb)
% takes an image in Rgb coords and returns either a vector or an image
% in Luv coords, depending on Type (Image or vector).

Vector = shiftdim(Im_rgb, 2);
Vector = double(reshape(Vector, [3, size(Im_rgb, 1)*size(Im_rgb, 2)]));
Vector_nz = 4096*(sum(Vector, 1) > 0);
Vector_z = not(Vector_nz);
Vector = max(Vector, [Vector_z; Vector_z; Vector_z]);	% to avoid all zero pixels

Matrix = [...
      0.490 0.310 0.200; ...
      0.177 0.812 0.011; ...
      0.000 0.010 0.990];

XYZ = Matrix*Vector;		% now get a vector in XYZ coords

Y0 = 1;		% is it true?
X0 = 1;
Z0 = 1;

u0 = 4*X0/(X0+15*Y0+3*Z0);
v0 = 9*Y0/(X0+15*Y0+3*Z0);

Luv_back = zeros(size(Vector));
Luv_back(1, :) = 25*(100*XYZ(2, :)/Y0).^(0.333)-16;
Zeros = zeros(size(Luv_back(1, :)));
Luv_back(1, :) = max(Luv_back(1, :), Zeros);
u = 4*XYZ(1, :) ./ (XYZ(1, :)+15*XYZ(2, :)+3*XYZ(3, :));
v = 9*XYZ(2, :) ./ (XYZ(1, :)+15*XYZ(2, :)+3*XYZ(3, :));
Luv_back(2, :) = 13*Luv_back(1, :).*(u-u0);
Luv_back(2, :) = max(Luv_back(2, :), Zeros);
Luv_back(3, :) = 13*Luv_back(1, :).*(v-v0);
Luv_back(3, :) = max(Luv_back(3, :), Zeros);
Luv_back = min(Luv_back, [Vector_nz; Vector_nz; Vector_nz]);	% zeroing out
Luv_back = Luv_back/(718.3176/256);
Luv_back = reshape(Luv_back, [3, size(Im_rgb, 1), size(Im_rgb, 2)]);
Luv_back = shiftdim(Luv_back, 1);

