%% method 4 ; directly interpolation
data_path = './video_capture/';
fu = 'capture_up.jpg';
fd = 'capture_d_cali.jpg';
save_path = strcat(data_path,'method4/');
cube_map_n = 'out_cubes_interp_v3.jpg';
max_disp = 20;
max_dy = 0;
fov = 90; cube_r = 512;
re_interpolate = true;
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
%% loop for measurement
for i = 1:4
    [face_u,~,~] = eq2perspective(eq_img_u,fov,90*(i-1),0,cube_r,cube_r);
    [face_d,~,~] = eq2perspective(eq_img_d,fov,90*(i-1),0,cube_r,cube_r);
    dy(1,i) = measure_dy_m4(face_u,face_d,max_disp);
    
    
end


%% 2.a warp each face with a dy
%% look up table
sample_x = [0 2 6 8 11 13 15]; % dy
sample_y = [0 15 30 45 60 75 90];
luta = interp1(sample_x,sample_y,0:15)/180*pi;
%% parameter re-used
d = cube_r/2;
r = -d+1 : 1 : d; % a horizontal vector
theta_a = abs(atand(r/d));
     


% look up table - faces
warp_cube = cell(1,4);
[lat_d,lon_d,old_w,old_h] = getOldLatLon(eq_w,eq_h);

for i = 1:numel(dy)
    face = 90*(i-1);
    rmp_img = real_method2(lat_d,luta(dy(1,i)),eq_img_u,old_w,eq_h);
    [u_per,~,~] = eq2perspective(rmp_img,fov,face,0,cube_r,cube_r);
    % [lut_per,~,~] = method2_2perspective_2in1(lat_d,luta(i),eq,old_w,eq_h, 90,face,0,cube_r,cube_r,128);
    imwrite(u_per,strcat(save_path,'lut_per_',int2str(face),'_',int2str(dy(1,i)),'.jpg'));
    warp_cube{1,i} = im2double(imread(strcat(save_path,'lut_per_',int2str(face),'_',int2str(dy(1,i)),'.jpg')));
end
%% refine boundary of front - right - left - bottom dy->dy_
if re_interpolate
    %% for theta in each cube face > 22.5
    dy_ = cell(1,4);
    for i = 1:4
        a = i;
        theta_b = 90 - theta_a;
        lower = theta_a + theta_b;
        up_a = theta_b ./ lower;
        up_b = theta_a ./ lower;
        
        dya = dy(1,a); % dya = dy_col(1,a)'
        %% positive sign part
        if i>3
            b= 1;
        else
            b = i+1;
        end
        dyb = dy(1,b); % dyb = dy_col(1,b)'
        dy_{1,i} = round(dya*up_a + dyb*up_b);
        dy_{1,i}(theta_a < fov/2/2) = dy(i); % 
        %% negative sign
        if i<2
            b= 4;
        else
            b = i-1;
        end
        dyb = dy(1,b);
        t = round(dya*up_a + dyb*up_b);
        dy_{1,i}(atand(r/d) < -fov/2/2) = t(atand(r/d) < -fov/2/2);
        %% keep max_dy: save computation
        if max_dy < max(dy_{1,i}(:))
            max_dy = max(dy_{1,i}(:));
        end
    end
    %% warping according to dy map
    % for each disp warp eq
    %%% for each face , paste corresponding pixel
    
    for d = 1:max_dy
        rmp_img = real_method2(lat_d,luta(d),eq_img_u,old_w,eq_h);
        for f = 1:4
            face = 90*(f-1);
            if size(dy_{1,f},1)==1
                dy_map = repmat(dy_{1,f},[cube_r 1 3]);
            else
                dy_map = dy_{1,f};
            end
            [u_per,~,~] = eq2perspective(rmp_img,fov,face,0,cube_r,cube_r);
            
            warp_cube{1,f}(dy_map==d) = u_per(dy_map==d);
        end
        
    end
    %% save img and reload
    for f = 1:4
        face = 90*(f-1);
        imwrite(warp_cube{1,f},strcat(save_path,'per_',int2str(face),'_interp.jpg'));
        warp_cube{1,f} = im2double(imread(strcat(save_path,'per_',int2str(face),'_interp.jpg')));
    end
end
%% measure top bottom (5.a)
face_topbot = cell(1,2);
dy_tb = round(sum(dy(:))/4);
%% warping raw top bottom
rmp_img = real_method2(lat_d,luta(dy_tb),eq_img_u,old_w,eq_h);

[face_topbot{1,1},~,~] = eq2perspective(rmp_img,fov,0,90,cube_r,cube_r);
[face_topbot{1,2},~,~] = eq2perspective(rmp_img,fov,0,-90,cube_r,cube_r);



%% refine top-bottom face
if re_interpolate
    %% split up down to 4 region
    d = cube_r/2;
    r = -d+1:d;
    r_row = repmat( -d+1:d  ,[cube_r 1]);
    r_col = repmat([-d+1:d]',[1 cube_r]);
    % maybe only one t is needed
    t1234 = cell(1,4);
    t1234{1,1} = (r_col > r_row) .* (r_col > -r_row); % up-front region
    t1234{1,2} = (r_col < r_row) .* (-r_col < r_row ); %= rot90(t1)
    t1234{1,3} = (r_col < r_row) .* (r_col <-r_row);
    t1234{1,4} = (r_col > r_row) .* (-r_row > r_col);
    %% interpolate top face
    dy_max_top = 0;
    dy_max_bot = 0;
    dya = dy_tb;
    %  top parameters
    theta_a = abs(atand(r/d));
    dy_top_map = zeros(cube_r);
    % dy_top = dy_ % todo : for differnt top bottom boundary
    dy_tb_ = cell(1,4); % can be cut
    rot_n = [3 0 1 2]; % how many to rot90
    %  bottom parameters
    dy_bot_map = zeros(cube_r);
    rotb_n = [1 0 3 2 ];
    dy_bot_  = dy_ ; %% if use different measurement in frlb -dy
    top2bot_i = [3 2 1 4];% t1234{top2bot_i(i)} t1234(i)
    
    for i = 1:4 % front right back left
        % size
        % dyb : 1*512'  = 512*1
        % dya : scalar
        %% top
        dyb = dy_{1,i}'; % on the +-45 deg boundary
        %       theta varies along horizontal direction
        theta_b = 45-theta_a;
        lower = theta_a + theta_b;
        up_a = theta_b./lower;
        up_b = theta_a./lower;
        
        dy_tb_{1,i} = round(dya*up_a + dyb*up_b); % 512*512
        % rotate for the mask and assign
        dy_tb_{1,i} = rot90(dy_tb_{1,i},rot_n(i));
        dy_top_map (t1234{1,i}==1) = dy_tb_{1,i}(t1234{1,i}==1);
        %% bottom
        %% interpolate bottom map
        dyb = dy_bot_{1,i}';
        theta_b = 45-theta_a;
        lower = theta_a + theta_b;
        up_a = theta_b./lower;
        up_b = theta_a./lower;
        
        dy_tb_{1,i} = round(dya*up_a + dyb*up_b); % 512*512
        % rotate for the mask and assign
        dy_tb_{1,i} = rot90(dy_tb_{1,i},rotb_n(i));
        dy_bot_map (t1234{1, top2bot_i(i) }==1) = ...
            dy_tb_{1,i}(t1234{1, top2bot_i(i) }==1);
        if dy_max_top < max(dy_top_map(:))
            dy_max_top = max(dy_top_map(:));
        end
        if dy_max_bot < max(dy_bot_map(:))
            dy_max_bot = max(dy_bot_map(:));
        end
    end
    %% to  paste back orignal dy ? v5-v6
    
    %% paste back refined top bot
    for d = 1:max(dy_max_bot,dy_max_top)
        rmp_img = real_method2(lat_d,luta(d),eq_img_u,old_w,eq_h);
        
        [temp_top,~,~] = eq2perspective(rmp_img,fov,0,90,cube_r,cube_r);
        [temp_bot,~,~] = eq2perspective(rmp_img,fov,0,-90,cube_r,cube_r);
        face_topbot{1,1}(dy_top_map ==d) = temp_top(dy_top_map==d);
        face_topbot{1,2}(dy_bot_map ==d) = temp_bot(dy_bot_map==d);
    end
    
end

%% cubemap construction
cube_map = zeros(3*cube_r,4*cube_r,3);
face = [4 1 2 3];
for i = 1:4
    cube_map(cube_r+1:2*cube_r,cube_r*(i-1)+1:cube_r*i,:) = warp_cube{1,face(i)};
    figure(1);imshow(cube_map)
end
cube_map(1:cube_r,cube_r+1:cube_r*2,:) = face_topbot{1,1};
cube_map(2*cube_r+1:3*cube_r,cube_r+1:2*cube_r,:) = face_topbot{1,2};
figure(1);imshow(cube_map)
imwrite(cube_map,strcat(save_path,cube_map_n));

% final_face = cell(1,6);
% final_face{1,1:4} = warp_cube{1,1:4};
% final_face{1,5} = face_topbot{1,1};
% final_face{1,6} = face_topbot{1,2};


