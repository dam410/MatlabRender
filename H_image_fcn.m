function [H_image] = H_image_fcn(h,l)
        H_image = [0,-(h-1)/2,h/2+0.5;-(l-1)/2,0,l/2+0.5;0,0,1];
end
