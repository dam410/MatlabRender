function [P] = image_to_camera(U,K,h,l)
        inv_K = inv(K);
        inv_H_image = inv(H_image_fcn(h,l));
        if size(U,1) == 2
                U = [U;ones(1,size(U,2))];
        end
        H_image = H_image_fcn(h,l);
        P = inv_K*inv_H_image*U;
        P = P./P([3,3,3],:);
end

