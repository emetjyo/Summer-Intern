function new_map = remap_bilinear(eq_img,new_i_map_h,new_i_map_w,new_map)
[eq_h,eq_w,~] = size(eq_img);
% h_map
fl_h_map = floor(new_i_map_h); ce_h_map = ceil(new_i_map_h);
amap = new_i_map_h-fl_h_map;
ce_h_map(ce_h_map<1) = eq_h;  ce_h_map(ce_h_map>eq_h) = 1;
fl_h_map(fl_h_map<1) = 1; fl_h_map(fl_h_map>eq_h) = 1;
% w_map
fl_w_map = floor(new_i_map_w); ce_w_map = ceil(new_i_map_w);
bmap = new_i_map_w-fl_w_map;
fl_w_map(fl_w_map<1) = eq_w;%eq_w;
t = fl_w_map>eq_w; fl_w_map(t) = eq_w;%1;
t = ce_w_map>eq_w; ce_w_map(t) = eq_w;%1;
ce_w_map(ce_w_map<1) = eq_w;
if ~exist('new_map','var')
    new_map = zeros(size(bmap,1),size(bmap,2),3);
end
if size(size(eq_img),2)==3
    ch_n = 3;
else
    ch_n = 1;
end
for ch = 1:ch_n
    if size(size(eq_img),2)==3
        imgch = eq_img(:,:,ch);
    else
        imgch = eq_img(:,:);
    end
    
    p1 = imgch((fl_w_map-1)*eq_h+fl_h_map);
    p2 = imgch((ce_w_map-1)*eq_h+fl_h_map);
    
    p3 = imgch((fl_w_map-1)*eq_h+ce_h_map);
    p4 = imgch((ce_w_map-1)*eq_h+ce_h_map);
    
    new_map(:,:,ch) = (1-bmap).*((amap.*p3)+(1-amap).*p1)+...
        bmap.*(amap.*p4+(1-amap).*p2);
    
end

end