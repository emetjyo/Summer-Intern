%% use look up table and disparity to gen a warp cubeface
clear all; close all;

data_path = './video_capture/';
eq = im2double(imread('./video_capture/capture_up.jpg'));

face = 0;
disp = imread(strcat('./video_capture/cube/disp/refine_lab',int2str(face),'.jpg'));
refper = im2double(imread(strcat('./video_capture/cube/py_cubed_',int2str(face),'.jpg')));
% refine_paste or paste
out_n = strcat('./video_capture/method3_result/refine_paste_',int2str(face),'.jpg');
opp_n = strcat(data_path,'cube/result',int2str(face),'_gau.jpg');
%% for 90
if face ==90
    disp =  imread(strcat('./video_capture/cube/disp/refine_lab90.jpg'));
    eq = im2double(imread('./video_capture/c90_warp90.jpg'));
    opp_n = strcat(data_path,'cube/','result90sec_uni.jpg');
end
%%
[eq_h,eq_w,~] = size(eq);

cube_r = 1080; radius = 128;
max_disp = 15; % the limit warping can handle
dq = 255/max_disp;
labels = disp/dq;

%% look up table
sample_x = [0 2 6 8 11 13 15];
sample_y = [0 15 30 45 60 75 90];

luta = interp1(sample_x,sample_y,0:15)/180*pi;
phi_map = luta(labels+1);
factor_d2r = pi/180;
% look up table - faces
[lat_d,lon_d,old_w,old_h] = getOldLatLon(eq_w,eq_h);
% ref75 = im2double(imread('./video_capture/lut_per/debug_deg75.jpg'));
for i = 1:numel(luta)
    rmp_img = real_method2(lat_d,luta(i),eq,old_w,eq_h);
% [lut_per,~,~] = method2_2perspective_2in1(lat_d,luta(i),eq,old_w,eq_h, 90,face,0,cube_r,cube_r,128);
    [lut_per,~,~] = eq2perspective(rmp_img,90,face,0,cube_r,cube_r,128);
    imwrite(lut_per,strcat(data_path,'cube/lut_per/','lut_per_',int2str(face),'_',int2str(i-1),'.jpg'));
%     figure(2);imshow(abs(ref75-lut_per));
end
%% read lut
lut = cell(16,1);
for i = 1:numel(lut)
    im = im2double(imread(strcat(data_path,'cube/lut_per/','lut_per_',int2str(face),'_',int2str(i-1),'.jpg')));
    lut{i,1} = im;
end
%%
paste = double(zeros(size(lut{1,1})));

for i = 0:max(labels(:))
    
    t = repmat(labels==i,[1,1,3]);
    paste(t) = lut{i+1,1}(t);
end
figure;imshow(paste)
% imwrite(paste,out_n);
figure(3);subplot(1,2,1);imshow(abs(lut{1,1}-refper));title('input diff')
subplot(1,2,2);imshow(abs(paste-refper));title('output diff')

% opp = im2double(imread(opp_n));
% figure(4);subplot(1,2,1);imshow(abs(opp-refper));title('opp vs ref')
% subplot(1,2,2);imshow(abs(paste-refper));title('new vs ref')




