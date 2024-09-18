function gd = geodesic2D(image)


[n,m] = size(image);

mask = double(roipoly(image));

[gx,gy]=gradient(image);
grad = gx.^2 + gy.^2;

f = 1e-3 + 1000*grad;

gd = fast_sweeping_mex(f,mask);

gd = gd - min(gd(:));
gd = gd./max(gd(:));



end