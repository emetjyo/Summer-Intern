
%%% theta: along y axis(longitude) , pi: along z-axis(latitude)
function [out_img,lat,lon] = eq2perspective(img,w_fov,th,phi,out_h,out_w,radius)
if ~exist('radius','var')
    radius = 1;
end

out_c_w = ((out_w+1)/2); out_c_h = ((out_h+1)/2);
factor_r2d = 180/pi;
factor_d2r = pi/180;
w_fov_r = w_fov*factor_d2r;
[eq_h,eq_w,~] = size(img);
eq_c_w = ((eq_w+1)/2); eq_c_h = ((eq_h+1)/2);

h_fov_r = 2*atan((out_h/out_w)*tan(w_fov_r/2));

virtual_per_w = 2*radius*tan(w_fov_r/2);
virtual_per_h = 2*radius*tan(h_fov_r/2);
dw = virtual_per_w/(out_w-1) ; % [1,...,out_w]*dw = [~,...,virtual_per_w]
dh = virtual_per_h/(out_h-1) ;
%%% construct D ratio 
% x: point to perspective img 
% y: point to longitude/width 
% z: point to latitude/height
D_x = zeros(out_h,out_w)+radius;
D_y = repmat(([1:out_w]-out_c_w)*dw,[out_h,1]); % [1:n] = 1Xn 
D_z = -repmat(([1:out_h]'-out_c_h)*dh,[1,out_w]);
% distance to perspective plane
D_map = sqrt(D_x.^2 + D_y.^2 + D_z.^2);
DR = radius./D_map; Dr_x = DR.*D_x; Dr_y = (DR.*D_y); Dr_z = (DR.*D_z);
xyz(:,:,1) = Dr_x; xyz(:,:,2) = Dr_y; xyz(:,:,3) = Dr_z;
xyz = reshape(xyz,[],3)';
R1 = rotz(th);
R22 = roty(-phi);

% R2 = rotationVectorToMatrix(R1_*(phi)*factor_d2r);% cannot use -phi here
% R2*R1(phi after theta) == R1*R22 (theta after phi,theta is less) 
xyz = R1*R22*xyz; xyz = xyz';xyz = reshape(xyz,out_h,out_w,[]);
% if non rotation => go to rotation
Dr_z = xyz(:,:,3); Dr_y = xyz(:,:,2); Dr_x = xyz(:,:,1);
lat = asin(Dr_z/radius); % phi radians

theta = atan(Dr_y./Dr_x);
lon = theta;
idx1 = Dr_x > 0;
idx2 = Dr_y > 0;

idx3 = ((~idx1) & (idx2));
idx4 = ((~idx1) & (~idx2));


lon(idx3) = theta(idx3) + pi;
lon(idx4) = theta(idx4) - pi;


lon = lon*factor_r2d;
lat = -lat*factor_r2d;
lon = lon / 180 * eq_c_w + eq_c_w; %% u = atan2(x,z)/2pi+0.5
lat = lat / 90 * eq_c_h + eq_c_h; %% v = y*0.5+0.5


out_img = zeros(out_h,out_w,3);
% heuristic cv2.remap 
lon(lon<1) = lon(lon<1)+eq_w; lon(lon>eq_w) = lon(lon>eq_w)-eq_w;
lat(lat<1) = lat(lat<1)+eq_h; lat(lat>eq_h) = lat(lat>eq_h)-eq_h;


%%% accelerate ? 
out_acc = zeros(size(out_img));
for ch = 1:size(img,3)
    imgch = img(:,:,ch);
    fl_omap = floor(lon); fl_amap = floor(lat);
    ce_omap = ceil(lon); ce_amap = ceil(lat);
    a_map = lon-fl_omap; b_map = lat - fl_amap;
    fl_omap(fl_omap==0) = eq_w; fl_amap(fl_amap==0) = eq_h;
    ce_omap(ce_omap>eq_w) = 1; ce_amap(ce_amap>eq_h) = 1;
    %
    p1_map = imgch(fl_amap+fl_omap*eq_h-eq_h); p2_map = imgch(fl_amap+ce_omap*eq_h-eq_h);
    p3_map = imgch(ce_amap+fl_omap*eq_h-eq_h); p4_map = imgch(ce_amap+ce_omap*eq_h-eq_h);
    %
    out_acc(:,:,ch) = (a_map.*p2_map+(1-a_map).*p1_map).*(1-b_map)+...
        (a_map.*p4_map+(1-a_map).*p3_map).*b_map; 
end

out_img = out_acc;
end


