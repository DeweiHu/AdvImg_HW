%clc, clear, close all;

%load TestHw8_I4L1_5.mat

% Resize D to get the same shape
dynamic = imresize(dynamic,size(static));
% Normalize intensity
dynamic = rescale(dynamic,0,255);
static = rescale(static,0,255);

%% Main
global alpha, alpha = 0.5;
global iter_num, iter_num = 150;
global ini_shape, ini_shape = size(static);

% Initialization
tx = zeros(ini_shape);
ty = zeros(ini_shape);

for i = 1:numlevel
    
    figure
    colormap gray
    subplot(2,2,1),imagesc(static),title('static')
    subplot(2,2,2),imagesc(dynamic),title('dynamic')
    
    % Update the shape 
    L = ini_shape(1)/2^(numlevel-i);
    shape = [L,L];
    tx = imresize(tx,shape);
    ty = imresize(ty,shape);
    D = imresize(dynamic,shape);
    S = imresize(static,shape);
    
    % Update the optical flow
    [tx,ty] = Demon(D,S,std,2*tx,2*ty,L);  
end

function [tx,ty] = Demon(dynamic,static,std,tx,ty,size)
    global iter_num, global alpha
    D_new = dynamic;
    [Gx,Gy] = gradient(static);
    G_mag = Gx.^2+Gy.^2;
    for i = 1:iter_num
        diff = (static-D_new);
        Vx = (diff.*Gx)./(G_mag+(alpha*diff.^2));
        Vy = (diff.*Gy)./(G_mag+(alpha*diff.^2));
        Vx(isnan(Vx))=0; 
        Vy(isnan(Vy))=0;
        
        Vx_smooth = imgaussfilt(Vx,std);
        Vy_smooth = imgaussfilt(Vy,std);

        tx = Vx_smooth+tx;
        ty = Vy_smooth+ty;
        
        [X,Y] = meshgrid(1:size);
        
        move_X = X+tx;
        move_Y = Y+ty;

        D_new = interp2(X,Y,dynamic,move_X,move_Y);
        D_new(isnan(D_new))=0;

        subplot(2,2,3),imagesc(D_new),title('Evolving Process')
        subplot(2,2,4),imagesc(abs(static-D_new)),title('Difference')
        drawnow
    end
end

