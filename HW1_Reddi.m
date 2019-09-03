clc,clear,close all

% Image Generator
size_b = 256;
size_f1 = 150;
size_f2 = 100;
img_b = uint8(zeros(size_b,size_b));
img_f1 = uint8(zeros(size_f1,size_f1));
img_f2 = uint8(zeros(size_f2,size_f2));

bg = imnoise(img_b,'gaussian',0.1,0.01);
fg1 = imnoise(img_f1,'gaussian',0.5,0.01);
fg2 = imnoise(img_f2,'gaussian',1,0.01);

pos = randi([1,size_b-size_f1],1,4);

img = bg;
img(pos(1):pos(1)+size_f1-1,pos(2):pos(2)+size_f1-1) = fg1;
img(pos(3):pos(3)+size_f2-1,pos(4):pos(4)+size_f2-1) = fg2;

%imwrite(mat2gray(img),'testimage.jpg')
[k1,k2] = Reddi(img)
%% Reddi's Method
function [k1,k2] = Reddi(img)
    [count,binLoc] = imhist(img);
    siz = size(img);
    output = uint8(zeros(siz(1),siz(2)));
    % initial
    k1 = 1./3*255;
    k2 = 2./3*255;
    e1 = 10;
    e2 = 10;

    figure
    subplot(2,2,1);imshow(img);title('Input image')
    subplot(2,2,2);imhist(img);xlim([0,256]);title('optimize threshold')
    hold on

    while abs(e1)>1 || abs(e2)>1
        e1 = .5*(RegionMean(1,k1,img)+RegionMean(k1+1,k2,img))-k1;
        e2 = .5*(RegionMean(k1+1,k2,img)+RegionMean(k2+1,256,img))-k2;
        % update
        k1 = k1+round(e1);
        k2 = k2+round(e2);
        % visualization
        xline(k1,'Color',rand(1,3),'LineWidth',1.5)
        hold on
        xline(k2,'Color',rand(1,3),'LineWidth',1.5)
        hold on
        pause(1)
    end

    for i = 1:siz(1)
        for j = 1:siz(2)
            if img(i,j) >= k2
                output(i,j) = 255;
            elseif img(i,j) < k2 && img(i,j) >= k1
                output(i,j) = 128;
            end
        end
    end

    subplot(2,2,3);imshow(output);title('Output image')
end
    
    
function m = RegionMean(a,b,img)
    [count,binLoc] = imhist(img);
    mul = count.*binLoc;
    m = sum(mul(a:b),1)./sum(count(a:b),1);
end
