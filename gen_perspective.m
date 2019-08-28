clear all; close all;
path = '../theta_img/';
path = './video_capture/';
% f12 = 'l_D'; id = '4';
% name1 = strcat(path,f12,'3.jpg'); name2 = strcat(path,f12,id,'.jpg');
% img1 = im2double(imread(name1)); img2 = im2double(imread(name2));
% adj_ang = 0;
% 
% out_path = strcat('./street/');
% outname = strcat(out_path,f12,'_per_',id,'.jpg');
%%% compute longitude shift
% fr = img2;fl = img1;
% opt_d = shift_longitude(fr,fl);(fr,fl,fovw_d,phi,theta,ph,pw)
% temp = zeros(size(fl)); temp(:,1+opt_d:end,:) = fl(:,1:end-opt_d,:);
% temp(:,1:opt_d,:) = fl(:,end-opt_d+1:end,:);
% img1 = temp;
% imwrite(img1,'../theta_img/B_d_lon.jpg');

%%%
img1 = im2double(imread('./video_capture/py_u.jpg'));
img1 = im2double(imread('./video_capture/0820debug/eq_1.jpg'));
% img2 = im2double(imread(name2));
%%% input parameters , angle all in degree (fov , theta , phi)
radius = 128;
w_fov = 90;
% output img size
out_w = 1080; out_h = 720;out_h = 1080;
th = 0; 
phi = 0; phi = 90;
% % opt_d = shift_longitude(img2,img1,w_fov,phi,th,out_h,out_w);
% % fl = img1;fr = img2;
% % for move-right
% temp = zeros(size(fl)); temp(:,1+opt_d:end,:) = fl(:,1:end-opt_d,:);
% temp(:,1:opt_d,:) = fl(:,end-opt_d+1:end,:);
% img1 = temp;
% % imwrite(img1,'../theta_img/D2_cali2D1.jpg');
% % for move-left
% 
% temp = zeros(size(fr)); temp(:,1:end-opt_d,:) = fr(:,1+opt_d:end,:);
% temp(:,end-opt_d+1:end,:) = fr(:,1:opt_d,:);
% img2 = temp;
% imwrite(img2,'../theta_img/l_D2.jpg')
% 
% img2 = im2double(imread('../theta_img/D1.jpg'));
% step 1 : show perspectice
[out_img1,lat1,lon1] = eq2perspective(img1,w_fov,th,phi,out_h,out_w,radius);
% [out_img2,lat2,lon2] = eq2perspective(img2,w_fov,th,phi+adj_ang,out_h,out_w,radius);
% figure(1);subplot(2,1,1);imshow(out_img1);title('lower view'); % subplot(2,1,2);imshow(out_img2);title('upper view');

% imwrite(out_img1,strcat(out_path,of1));
% imwrite(out_img2,strcat(out_path,of2));
imwrite(out_img1,'./video_capture/0820debug/face1.jpg')
% imwrite(out_img2,outname);





