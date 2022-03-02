function [U] = camera_to_image(P,K,h,l)
        H_image = H_image_fcn(h,l);
        U = H_image*K*P;
        U = U./U([3,3,3],:);
end

