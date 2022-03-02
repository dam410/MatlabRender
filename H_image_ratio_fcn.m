function [H_image] = H_image_ratio_fcn(h,l)
        H_image = [0,-(l-1)/2,l/2;-(l-1)/2,0,h/2;0,0,1];
end
