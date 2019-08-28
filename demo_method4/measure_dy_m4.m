% measure dy between 2 perspective image , 
% input with unspecify lower/upper
% so pad per1 with 2*disp_max row
function [dy] = measure_dy_m4(per1,per2,max_disp,region)
[h,w,~] = size(per1);

if ~exist('max_disp','var')
    max_disp = 25;
end
if ~exist('region','var')
    %% gaussian kernal : not used yet
    x = 1:w; y =1:h;
    [X,Y] = meshgrid(x,y);
    X = X-w/2; Y = Y-h/2;
    sigx = w/6; sigy = h/6;
    region = 2*exp(- (X.^2/(2*sigx^2)+Y.^2/(2*sigy^2)) );
end
dou_disp = max_disp*2;

pad1 = zeros(h+dou_disp,w,3);
% pad max_disp/2 row-zeros in upper region of per1
pad1(max_disp+1:max_disp+h,:,:) = per1;
temp1 = zeros(size(per1));

err = zeros(dou_disp+1,1);
region = repmat(region,[1 1 3]);
%%
% d = 1 :map to per1 up shift max_disp
% d = 26 :map to original per1
% d = dou_disp :map to per1 down shift max_disp
for d = 1:dou_disp+1 
    temp1(:,:,:) = pad1(d:d+h-1,:,:) ;
    er_map = (temp1 - per2).*region;
    err(d) = sum(sqrt(er_map(:).^2))/(h*w);

end

%% display
[opt_er, opt_disp] = sort(err);

dy = opt_disp(1)-round(max_disp+1)
if dy>0
    disp('img1 is upper view')
else 
    disp('img1 is lower view')
end
% dy<0 : per1 is lower view , otherwise per2 is lower view
end