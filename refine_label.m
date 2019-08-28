clear all;close all;
%% function refine_img_ =  
addpath filtercode
addpath L0smoothing
disp1 = double(imread('./video_capture/cube/disp/pycube180.jpg'));
left_img = im2double(imread('./video_capture/cube/py_cubeu_180_b.jpg'));
save_n = strcat('video_capture/cube/disp/','refine_lab180_.jpg');
% %% mask region for 0 deg
% mh1 = 334; mh2 = 390; mw1 = 528; mw2 = 607;
%% mask region for 270 deg
% mh1 = 970; mh2 = 1079; mw1 = 239; mw2 = 395;
%% mask region for 90 deg (right)
% mh1 = 405; mh2 = 635; mw1 = 647; mw2 = 851;
%% mask region for 180 deg
p1 = double(imread('./video_capture/cube/py_cubeu_180_b.jpg'));
p2 = double(imread('./video_capture/cube/py_cubed_180_b.jpg'));
mask_region = (sum(p1,3)==0) ;
figure;imshow(mask_region);
face = 180;

%%
q = 15;
dq = 255/q;
lab1 = round(disp1/dq);

inlab1 = lab1;
bad_disp = lab1==0;
if exist('mh1','var')
bad_disp(mh1:mh2,mw1:mw2,:) = 1;
end
if exist('mask_region','var')
    bad_disp = (lab1==0)|mask_region ;
end

figure(1); subplot(2,2,1);title('before LR check')
imshow(lab1,[]); subplot(2,2,2);
imshow(bad_disp);


%% weighted median filter
gamma_c = 0.1;          % \sigma_c in eq. (6)
gamma_d = 9;            % \sigma_s in eq. (6)
winsize = 29;          % filter kernel of weighted median
numDisp = 16;
lab1(bad_disp) = -1;
% refine_img = weightedMedianMatlab(rot90(left_img,3),rot90(lab1,3),winsize,gamma_c,gamma_d);
% lab1(lab1==0) = refine_img(lab1==0);
refine_img = fillPixelsReference(rot90(left_img), rot90(lab1), gamma_c, gamma_d, winsize, numDisp);
subplot(2,2,3);
imshow(rot90(refine_img,3),[]); title('after LR check');
subplot(2,2,4);
imshow((rot90(refine_img,3)~=inlab1).*bad_disp);

%%
if exist('face','var') && face==180
    refine_img(rot90(mask_region)) = -1;
    refine_img = rot90(fillPixelsReference(rot90(left_img,2), rot90(refine_img,2), gamma_c, gamma_d, winsize, numDisp),2);
    refine_img(rot90(mask_region)) = 1;
end
refine_img2 = weightedMedianMatlab(rot90(left_img),refine_img,winsize,gamma_c,gamma_d);
refine_img2 = rot90(refine_img2,3);
%%
t1 = inlab1~=0; t2 = refine_img2==0;
too_over = t1 & t2;
refine_img2(too_over) = inlab1(too_over);
%%
figure(2);
subplot(2,2,1); imshow(rot90(refine_img,3),[]); title('before WMF2');
subplot(2,2,2); imshow(refine_img2,[]);
% refine_img3 = weightedMedianMatlab(rot90(left_img,3),refine_img2,winsize,gamma_c,gamma_d);
%%
refine_img3 = refine_img2;
t1 = inlab1~=0; t2 = refine_img3==0;
too_over = t1 & t2;
refine_img3(too_over) = inlab1(too_over);
%%
subplot(2,2,3); imshow(refine_img2,[]);
subplot(2,2,4); imshow(refine_img3,[]);


figure(7);
refine_img_ = medfilt2(refine_img3, [19 19]);
subplot(1,3,1);imshow(refine_img_,[]);title('MF 1st');
refine_img_ = medfilt2(refine_img_, [9 9]);
subplot(1,3,2);imshow(refine_img_,[]);title('MF 2nd');
refine_img_ = medfilt2(refine_img_, [5 5]);
subplot(1,3,3);imshow(refine_img_,[]);title('MF 3rd');


figure(10);subplot(1,2,1);imshow(refine_img_,[]);title('over');
t1 = inlab1~=0; t2 = refine_img_ ==0;
too_over = t1 & t2;
refine_img_(too_over) = inlab1(too_over);


subplot(1,2,2);imshow(refine_img_,[]);title('not too over');
refine_img_ = min(max(refine_img_*dq,0),255)/255.0;
% imwrite(refine_img_,save_n_notover);
mask_region = (sum(p1,3)==0) |(sum(p2,3)==0);
refine_img_ (mask_region) = 0;
refine_img_ = medfilt2(refine_img_, [19 19]);
imwrite(refine_img_,strcat('video_capture/cube/disp/','refine_lab180_mask.jpg'));

% %% paste face onto cube map
% img_left = im2double(imread('./video_capture/cube/result270_uni.jpg'));
% img_front = im2double(imread('./video_capture/cube/result0uni.jpg'));
% img_right = im2double(imread('./video_capture/cube/result90sec_uni.jpg'));
% img_back = im2double(imread('./video_capture/cube/py_cubeu_180.jpg'));
% 
% img_top = im2double(imread('./video_capture/cube/py_cubeu_top.jpg'));
% img_bot = im2double(imread('./video_capture/cube/py_cubeu_bot.jpg'));
% 
% [h,~,~] = size(img_top);
% 
% cube = zeros(3*h,4*h,3);
% 
% cube(h+1:2*h,1:h,:) = img_left;
% cube(h+1:2*h,h+1:2*h,:) = img_front;
% cube(h+1:2*h,2*h+1:3*h,:) = img_right;
% cube(h+1:2*h,3*h+1:4*h,:) = img_back;
% cube(1:h,h+1:2*h,:) = img_top;
% cube(2*h+1:end,h+1:2*h,:) = img_bot;
% imwrite(cube,'../py360convert-master/temp_result_cmuni2.jpg')
