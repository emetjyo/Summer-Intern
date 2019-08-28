function out_img = real_method2(old_lat_d,shift_phi,eq_img,old_w,eq_h)
factor_r2d = 180/pi;
%% 1. gen a new warped eq
%%%%% warp function

boundary_d = asin(0.5*sin(shift_phi))*factor_r2d;
up_t = old_lat_d > (boundary_d*ones(size(old_lat_d)));
theta_map = 90 - old_lat_d;
theta_map(~up_t) = 90 + old_lat_d(~up_t);
%     x_map = asin(shift_phi*sin(lat_r));%theta_map/factor_r2d)); bug#1
%     x_map = asin(shift_phi*sin(theta_map/factor_r2d)); % if use R*phi , large scale obtained
x_map = asin(sin(shift_phi)*sin(theta_map/factor_r2d));

A_map = theta_map+x_map;
B_map = theta_map-x_map;
new_old_lat_d_A = 90-A_map;
new_old_lat_d_B = B_map-90;
% get mapping latitude (longitude unchanged)
old_lat_d(up_t) = new_old_lat_d_A(up_t);
old_lat_d(~up_t) = new_old_lat_d_B(~up_t);
% get he new h locaions (uv-> pixel locations)
lat_r = old_lat_d/factor_r2d;
new_h = (0.5-lat_r/pi)*eq_h-0.5;

out_img = remap_bilinear(eq_img,new_h,old_w);

end