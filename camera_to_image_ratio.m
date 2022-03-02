function [U] = camera_to_image_ratio(P,K,h,l)
        H_image = H_image_ratio_fcn(h,l);
        U = H_image*K*P;
        U = U./U([3,3,3],:);
end

