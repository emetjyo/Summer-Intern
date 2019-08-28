% given input image size output uv map
function [lat_d,lon_d,old_w,old_h] = getOldLatLon(eq_w,eq_h)
% get old lat lon
factor_r2d = 180/pi;
old_w = repmat([1:eq_w],[eq_h,1]);
old_h = repmat([1:eq_h]',[1,eq_w]);
u_map = (old_w+0.5)/eq_w; % acually unused
v_map = (old_h+0.5)/eq_h;
lon_r = (u_map-0.5)*2*pi; lon_d = factor_r2d*lon_r;
lat_r = (0.5-v_map)*pi; lat_d = factor_r2d*lat_r; %(89~-90)
end
% new_h = (0.5-lat/pi)*eq_h-0.5
% new_w = (lon/2/pi+0.5)*eq_w-0.5