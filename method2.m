% method 2 : move he upper view downward and generae a virual lower view
% generate a warped virtual eq
clear all;
% %% measure d1
datapath = './video_capture/';
nameu = strcat(datapath,'cube/py_cubeu_0.jpg');
named = strcat(datapath,'cube/py_cubed_0.jpg');
name_eq = strcat(datapath,'./capture_up.jpg');
per_w = 1080;per_h = 720; fov_w = 90;radius = 128;
th = 0; phi =0;
per1 = im2double(imread(nameu));
per2 = im2double(imread(named));
%%% calibrate longitude shift (assumme fl fr is known)
% fl = im2double(imread(strcat(datapath,'capture_down.jpg')));fr = im2double(imread(strcat(datapath,'capture_up.jpg')));
% [opt_d,eq_cali] = shift_longitude(fr,fl,90,0,0,720,1080);
% [per_cali,~,~] = eq2perspective(eq_cali,fov_w,th,phi,per_h,per_w,radius);
% per2 = per_cali;
% imwrite(per_cali,strcat(datapath,'cali_d.jpg'));
%%% measure dy
% dy = measure_dy(per1,per2);


%%

% D : fovw = 90 , phi = 10
%
% name_eq = strcat(datapath,'cube/c90_warp90.jpg');

factor_r2d = 180/pi;w_fov_r = fov_w/factor_r2d;
eq_img = im2double(imread(name_eq));
% eq_img = im2double(imread('./street/l_D1.jpg'));
[eq_h,eq_w,~] = size(eq_img);
%%
h_fov_r = 2*atan((per_h/per_w)*tan(w_fov_r/2));
virual_per_h = 2*radius*tan(h_fov_r/2);
dh = virual_per_h/per_h ;
dy = 5; % measure from eq and perspective

apple = atan(dy/360)*180/pi;
R_per = 360/tan(33.69/180*pi);
dy = 2;
orange = atan(dy/R_per)*180/pi;
%%

a_phi = [5.7 6.3 7.1 8.1 8.7 9.5 10.3 11.3 12.5 14 15.9 16.6 18.4 21.8 26 30 33.7 38 42 45 50 53 57  60 65 70 75 80 85]/factor_r2d;

% a_phi = [ 50 53 57 60 65 70 75 80 85 90]/factor_r2d;
a_phi = [ 20 30 40]/factor_r2d;
diff_cell = cell(1+numel(a_phi),3);
diff_cell{1,1} = abs(per1-per2);
diff_cell{1,2} = 0;
diff_cell{1,3} = per1;
ic = 2;
R=1;
% warp812 = eq_img;

for a = 1:numel(a_phi)
    
    new_phi = atan(dy/R);
    new_phi = a_phi(a);
    [lat_d,lon_d,old_w,old_h] = getOldLatLon(eq_w,eq_h);shift_phi = new_phi;
    
    disp(['R = ',num2str(R),', shift_phi(deg.) = ',num2str(shift_phi*factor_r2d)])
    
    deg = round(shift_phi*factor_r2d*10); % for output deg display only
    up_t = lat_d> shift_phi*factor_r2d;
    % compute mapping old lat lon -> new lat lon
    theta_map = 90-lat_d;
    theta_map(~up_t) = 90+lat_d(~up_t); % max = 129 ?
    % x_map = asin(shift_phi*sin(lat_r));%theta_map/factor_r2d)); bug#1
%     x_map = asin(shift_phi*sin(theta_map/factor_r2d)); % if use R*phi , large scale obtained
    x_map = asin(sin(shift_phi)*sin(theta_map/factor_r2d))...
        *factor_r2d;
    
    A_map812 = atan(sin(theta_map/factor_r2d)./(sin(shift_phi)+cos(theta_map/factor_r2d)));
    %% theta map is deg. x_map is rad
    % if replace x_map with x_map*factor -> distortion severe?
    A_map = theta_map+x_map; % i think it's right
    B_map = theta_map-x_map; %
    new_lat_d_A = 90-A_map;
    new_lat_d_B = B_map-90; % -(90-B_map)
%     new_lat_d_A = 90 - A_map812*factor_r2d;
%     new_lat_d_B = A_map812*factor_r2d-90;
    
    % get the latitude on new sphere
    lat_d(up_t) = new_lat_d_A(up_t);
    lat_d(~up_t) = new_lat_d_B(~up_t);
    % ge he new h locaions
    lat_r = lat_d/factor_r2d;
    new_h = (0.5-lat_r/pi)*eq_h-0.5;
%     warp812 = remap_bilinear(warp812,new_h,old_w);
vec_lat = lat_r(:,1); % size = eq_h,1
vec_theta = theta_map(:,1);
%% 0820 upsample eq
dt_map = ones(eq_h,1);
dt_map = 
%%
    out_img = remap_bilinear(eq_img,new_h,old_w);
    % imwrite(out_img,warpn);
    
    %% gen a new perspective
    %% cubemap
    cube_r = 1080;
    [new_per,lat2,lon2] = eq2perspective(out_img,90,0,0,cube_r,cube_r,radius);
    %%
%     [new_per,lat2,lon2] = eq2perspective(out_img,fov_w,th,phi,per_h,per_w,radius);
    diff_per = abs(new_per-per2);
    diff_cell{ic,1} = diff_per;
    diff_cell{ic,2} = deg/10;
    diff_cell{ic,3} = new_per;
    ic = ic+1;
    imwrite(new_per,strcat('./video_capture/0820debug/per_'...
        ,int2str(round(a_phi(a)*factor_r2d)),'.jpg'));
%     imwrite(diff_per,'./video_capture/cube/c90_diff90.jpg')
    imwrite(out_img,strcat('./video_capture/0820debug/eq_'...
        ,int2str(round(a_phi(a)*factor_r2d)),'.jpg'));
end
figure(7);imshow(new_per);
figure(8);imshow(out_img);
% save('./video_capture/cube/double_cube90_difmap_phi.mat','diff_cell');
