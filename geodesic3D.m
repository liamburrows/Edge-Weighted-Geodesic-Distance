function gd = geodesic3D(image)


[n,m,sl] = size(image);

slice = round(sl/2);
mask0 = double(roipoly(image(:,:,slice)));
mask = zeros(size(image));
mask(:,:,slice) = mask0;

[gx,gy,gz]=gradient(image);
grad = gx.^2 + gy.^2 + gz.^2;

f = 1e-3 + 1000*grad;

gd = fast_sweeping3D_mex(f,mask);

gd = gd - min(gd(:));
gd = gd./max(gd(:));



end