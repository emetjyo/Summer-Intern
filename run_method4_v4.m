%% method 4 ; directly interpolation
data_path = './video_capture/';
fu = 'capture_up.jpg';
fd = 'capture_d_cali.jpg';
save_path = strcat(data_path,'method4/');
cube_map_n = 'out_cubes_interp_v5.jpg';
max_disp = 20;
max_dy = 0;
fov = 90; cube_r = 512;
profile on;
%%
eq_img_u = im2double(imread(strcat(data_path,fu)));
eq_img_d = im2double(imread(strcat(data_path,fd)));
[eq_h,eq_w,~] = size(eq_img_u);
d2r = pi/180; r2d = 180/pi;
dy = zeros(1,4);
%% generate cube face : up
%% measure dy

%% measure lower region(1/4) dy
dy_b = zeros(1,4);
dy_u = dy_b;
lower_mask = zeros(cube_r);
lower_mask(end-round(cube_r/4):end,:) = 1;
upper_mask = rot90(lower_mask,2);


%% front right back left
%% look up table
sample_x = [0 2 6 8 11 13 15]; % dy
sample_y = [0 15 30 45 60 75 90];
luta = interp1(sample_x,sample_y,0:15)/180*pi;
%% parameter re-used
d_per = cube_r/2;
r = -d_per+1 : 1 : d_per; % a horizontal vector
theta_a = abs(atand(r/d_per)); % 1*512
theta_a_map = repmat(theta_a , [cube_r 1]);% <direction> 512*512
r_map = repmat(r , [cube_r 1]); % 512*512
dy_frbl = zeros(cube_r,cube_r,4);

warp_cube = cell(1,4);
%% initialized and refine lower-upper 1/4
for i = 1:4
    %% measure dy
    [face_u,~,~] = eq2perspective(eq_img_u,fov,90*(i-1),0,cube_r,cube_r);
    [face_d,~,~] = eq2perspective(eq_img_d,fov,90*(i-1),0,cube_r,cube_r);
    warp_cube{1,i} = face_u;
    dy(1,i) = measure_dy_m4(face_u,face_d,max_disp);
    dy_b(1,i) = measure_dy_m4(face_u,face_d,max_disp,lower_mask);
    dy_u(1,i) = measure_dy_m4(face_u,face_d,max_disp,upper_mask);
    %% initialized dy
    dy_frbl(:,:,i) = dy(1,i);
    dy_frbl(end-round(cube_r/4):end,:,i) = dy_b(1,i);
    dy_frbl(1:round(cube_r/4),:,i) = dy_u(1,i);
    %% interpolation 1/2 - 1/4 bottom and top
    vec_dy = reshape(dy_frbl(:,end,i),1,[]); % reshape tp horizonal
    theta_b = atand(3/4)-theta_a;
    dya = dy(1,i); 
    lower = theta_a + theta_b;
    up_a = theta_b./lower;
    up_b = theta_a./lower;
    %% interpolate bottom region
    dyb =  dy_b(1,i);
    temp = round(dya*up_a + dyb*up_b);
    t = (r/d_per < 3/4).*(r/d_per>=0);
    vec_dy(t==1) = temp(t==1);
    %% interpolate upper region
    dyb =  dy_u(1,i);
    temp = round(dya*up_a + dyb*up_b);
    t = (r/d_per > -3/4).*(r/d_per < 0);
    vec_dy(t==1) = temp(t==1);
    %% store back
    dy_frbl(:,:,i) = repmat(vec_dy',[1 cube_r]);
    if max(vec_dy(:)) > max_dy
        max_dy = max(vec_dy(:));
    end
   
end

%% initialize - faces actually can be optimized

[lat_d,lon_d,old_w,old_h] = getOldLatLon(eq_w,eq_h);

% for i = 1:numel(dy)
%     face = 90*(i-1);
%     rmp_img = real_method2(lat_d,luta(dy(1,i)),eq_img_u,old_w,eq_h);
%     [warp_cube{1,i},~,~] = eq2perspective(rmp_img,fov,face,0,cube_r,cube_r);
%     % [lut_per,~,~] = method2_2perspective_2in1(lat_d,luta(i),eq,old_w,eq_h, 90,face,0,cube_r,cube_r,128);
% %     imwrite(u_per,strcat(save_path,'lut_per_',int2str(face),'_',int2str(dy(1,i)),'.jpg'));
% %     warp_cube{1,i} = im2double(imread(strcat(save_path,'lut_per_',int2str(face),'_',int2str(dy(1,i)),'.jpg')));
% end

%% refine boundary of front - right - left - bottom dy_frbl->dy_map
%% for theta in each cube face > 22.5
dy_map = cell(1,4);
for i = 1:4
    a = i;
    theta_b = 90 - theta_a; % 1*512
    lower = theta_a + theta_b;
    up_a = theta_b ./ lower;
    up_b = theta_a ./ lower;
    
    dya = dy_frbl(:,end,a); % 512*1
    %% positive sign part
    if i>3
        b= 1;
    else
        b = i+1;
    end
    dyb = dy_frbl(:,end,b); % 512*1
    dy_map{1,i} = round(dya*up_a + dyb*up_b);  % 512*512
    temp = dy_frbl(:,:,a);%% input
    dy_map{1,i}(theta_a_map < fov/2/2) = temp(theta_a_map < fov/2/2); %
    %% negative sign
    if i<2
        b= 4;
    else
        b = i-1;
    end
    dyb = dy_frbl(:,end,b); % 512*1
    t = round(dya*up_a + dyb*up_b); % 512*512
    dy_map{1,i}(atand(r_map/d_per) < -fov/2/2) = t(atand(r_map/d_per) < -fov/2/2);
    %% keep max_dy: save computation
    if max_dy < max(dy_map{1,i}(:))
        max_dy = max(dy_map{1,i}(:));
    end
end
dy_map_3d = cell(1,4);

%% measure top bottom (5.a)
%% ver.4
dy_top = round(sum(dy(:))/4);
dy_bot = round(sum(dy(:))/4);
%% ver.5
dy_top = round(sum(dy_u(:))/4);
dy_bot = round(sum(dy_b(:))/4);
%% initialize top bottom
rmp_imgt = real_method2(lat_d,luta(dy_top),eq_img_u,old_w,eq_h);
rmp_imgb = real_method2(lat_d,luta(dy_bot),eq_img_u,old_w,eq_h);

[face_top,~,~] = eq2perspective(rmp_imgt,fov,0,90,cube_r,cube_r);
[face_bot,~,~] = eq2perspective(rmp_imgb,fov,0,-90,cube_r,cube_r);



%% refine top-bottom face
%% split up down to 4 region
r_row = repmat( -d_per+1:d_per  ,[cube_r 1]);
r_col = repmat([-d_per+1:d_per]',[1 cube_r]);
% maybe only one t is needed
mask = cell(1,4);
mask{1,1} = (r_col > r_row) .* (r_col > -r_row); % up-front region
mask{1,2} = (r_col < r_row) .* (-r_col < r_row ); %= rot90(t1)
mask{1,3} = (r_col < r_row) .* (r_col <-r_row);
mask{1,4} = (r_col > r_row) .* (-r_row > r_col);
%% interpolate top face

%  top parameters
dy_top_map = zeros(cube_r);
rot_n = [3 0 1 2]; % how many to rot90
%  bottom parameters
dy_bot_map = zeros(cube_r);
rotb_n = [1 0 3 2 ];
t2b_i = [3 2 1 4];% t1234{top2bot_i(i)} t1234(i)
%% interpolate top bot from [front right back left]
for i = 1:4 
    % dyb : 1*512'  = 512*1
    % dya : scalar
    %% top : interpolation
    dyb = dy_map{1,i}(1,:)'; % on the +-45 deg boundary
    %       theta varies along horizontal direction
    theta_b = 45-theta_a;
    lower = theta_a + theta_b;
    up_a = theta_b./lower;
    up_b = theta_a./lower;
    % rotate for the mask and assign
    tmp_dy = rot90(round( dy_top*up_a + dyb*up_b),rot_n(i)); % 512*512
    dy_top_map (mask{1,i}==1) = tmp_dy(mask{1,i}==1);
    %% bottom interpolation
    dyb = dy_map{1,i}(end,:)';
    theta_b = 45-theta_a;
    lower = theta_a + theta_b;
    up_a = theta_b./lower;
    up_b = theta_a./lower;
    % rotate for the mask and assign
    tmp_dy = rot90(round(dy_bot*up_a + dyb*up_b),rotb_n(i)); % 512*512
    dy_bot_map (mask{1, t2b_i(i) }==1) = tmp_dy(mask{1, t2b_i(i) }==1);
    %% update max_dy
    if max_dy < max(dy_top_map(:)) 
        max_dy = max(dy_top_map(:));
    end
    if max_dy < max(dy_bot_map(:))
        max_dy = max(dy_bot_map(:));
    end
    
    
    %% save dy_map
    dy_map_3d{1,i} = repmat(dy_map{1,i},[1 1 3]); 
end

%% paste back refined top bot
for d = 1:max_dy
    rmp_img = real_method2(lat_d,luta(d),eq_img_u,old_w,eq_h);
    %% lrbf
    for f = 1:4
        face = 90*(f-1);
        [u_per,~,~] = eq2perspective(rmp_img,fov,face,0,cube_r,cube_r);
        warp_cube{1,f}(dy_map_3d{1,f}==d) = u_per(dy_map_3d{1,f}==d);
    end
    %% top bottom
    [temp_top,~,~] = eq2perspective(rmp_img,fov,0,90,cube_r,cube_r);
    [temp_bot,~,~] = eq2perspective(rmp_img,fov,0,-90,cube_r,cube_r);
    face_top(dy_top_map ==d) = temp_top(dy_top_map==d);
    face_bot(dy_bot_map ==d) = temp_bot(dy_bot_map==d);
end

%% cubemap construction
cube_map = zeros(3*cube_r,4*cube_r,3);
face = [4 1 2 3];
for i = 1:4
    cube_map(cube_r+1:2*cube_r,cube_r*(i-1)+1:cube_r*i,:) = warp_cube{1,face(i)};
    figure(1);imshow(cube_map)
end
cube_map(1:cube_r,cube_r+1:cube_r*2,:) = face_top;
cube_map(2*cube_r+1:3*cube_r,cube_r+1:2*cube_r,:) = face_bot;
figure(1);imshow(cube_map)
imwrite(cube_map,strcat(save_path,cube_map_n));


profile viewer;
