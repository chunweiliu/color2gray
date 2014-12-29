clear all

image_dir = '../images/';
image_name = '*g';

images = dir([image_dir, image_name]);

for n = 1:numel(images)
    
    % For every images in the folder
    image = [image_dir, images(n).name];

    % Preprocessing image to desired format
    if ischar(image), image = imread(image); end
    if isinteger(image), image = im2double(image); end

    % Comparison
    S.mat2gray = @mat2gray;  % Color image
    S.rgb2gray = @rgb2gray;
    S.decolor_linear2006 = @decolor_linear2006;
    S.decolor_nonlinear2009 = @decolor_nonlinear2009;
    S.decolor_sparse = @decolor_sparse;
    out = structfun(@(f) f(image), S, 'UniformOutput', false);

    % % Display output
    fn = fieldnames(S);
    for m = 1:numel(fn)
        subplot(numel(images), numel(fn), (n - 1) * numel(fn) + m)
        imshow(out.(fn{m}))
        
        name = images(n).name;
        out_name = sprintf('../out/%s_%d.png', name(1:end-4), m);       
        imwrite(out.(fn{m}), out_name, 'png');
    end

end