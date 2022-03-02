function [P] = image_to_camera_direct(U,K,h,l)
        inv_K = inv(K);
        if size(U,1) == 2
                U = [U;ones(1,size(U,2))];
        end
        P = inv_K*U;
        P = P./P([3,3,3],:);
end

