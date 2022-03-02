function [X_0,Y_0,H_image,lignes,colonnes] = image_coordinate_transform(h,l)
        %% Calculate the 3D coordinates of all the points
        lignes = -1:2/(l-1):1;
        colonnes = -1:2/(h-1):1;
        [grid_X,grid_Y] = meshgrid(lignes,colonnes);
        X_0 = grid_X(:);
        Y_0 = grid_Y(:);
        H_image = H_image_fcn(h,l);
end

