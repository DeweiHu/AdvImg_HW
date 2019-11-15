clc;clear;close all;

load TestHw8_I4L3_15.mat

% Resize D to get the same shape
dynamic = imresize(dynamic,size(static));
% Normalize intensity
dynamic = rescale(dynamic,0,255);
static = rescale(static,0,255);
%%
figure,
colormap gray
subplot(2,2,1),imagesc(static),title('static')
subplot(2,2,2),imagesc(dynamic),title('dynamic')

global alpha, alpha = 1;
global iter_num, iter_num = 150;

if numlevel == 1
    tx = zeros(256,256);
    ty = zeros(256,256);
    dynamic1 = dynamic;
    static1 = static;
    [Tx,Ty] = Demon(dynamic1,static1,std,tx,ty,256);
end

if numlevel == 3
    tx = zeros(64,64);
    ty = zeros(64,64);
    dynamic1 = imresize(dynamic,[64,64]);
    static1 = imresize(static,[64,64]);
    
    [Tx1,Ty1]=Demon(dynamic1,static1,std,tx,ty,size(dynamic1,2));
    
    Tx1_resize = imresize(Tx1,[128,128]);
    Ty1_resize = imresize(Ty1,[128,128]);
    
    figure,
    colormap gray
    subplot(2,2,1),imagesc(static),title('static')
    subplot(2,2,2),imagesc(dynamic),title('dynamic')
    
    
    dynamic2 = imresize(dynamic,[128,128]);
    static2 = imresize(static,[128,128]);
    
    [Tx2,Ty2] = Demon(dynamic2,static2,std,Tx1_resize,...
        Ty1_resize,size(dynamic2,2));
    
    Tx2_resize = imresize(Tx2,[256,256]);
    Ty2_resize = imresize(Ty2,[256,256]);
    
    figure,
    colormap gray
    subplot(2,2,1),imagesc(static),title('static')
    subplot(2,2,2),imagesc(dynamic),title('dynamic')
    
    dynamic3 = imresize(dynamic,[256,256]);
    static3 = imresize(static,[256,256]);
    
    [Tx3,Ty3]=Demon(dynamic3,static3,std,Tx2_resize,...
        Ty2_resize, size(dynamic3,2));
    
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