% measure dy between 2 perspective image , input with unspecify lower/upper
% tag
function [dy] = measure_dy(per1,per2)
max_disp = 50;
[h,w,~] = size(per1);
pad1 = zeros(h+max_disp,w,3);
% pad 25-row 0 in upper region of per1
pad1(max_disp/2+1:max_disp/2+h,:,:) = per1;
temp1 = zeros(size(per1));
er_ed = round(h/2);
err = zeros(max_disp+1,1);
% d = 1 map to per1 shift -25
% d = 26 map to original per1
for d = 1:max_disp+1 
    temp1(:,:,:) = pad1(d:d+h-1,:,:) ;
    er_map = temp1(1:er_ed,:,:) - per2(1:er_ed,:,:);
    err(d) = sum(sqrt(er_map(:).^2))/(er_ed*w);

end
[opt_er, opt_disp] = sort(err)
top_5_energy = opt_er(1:5)
dy = opt_disp(1)-26
if dy>0
    disp('img1 is upper view')
else 
    disp('img1 is lower view')
end
% dy<0 : per1 is lower view , otherwise per2 is lower view
end