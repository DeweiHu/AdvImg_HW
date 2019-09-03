clc, clear, close all

% Image Generator
size_b = 256;
size_f = 80;
img_b = uint8(zeros(size_b,size_b));
img_f = uint8(zeros(size_f,size_f));

bg = imnoise(img_b,'gaussian',0.2,0.01);
fg = imnoise(img_f,'gaussian',0.8,0.03);

pos = randi([1,size_b-size_f],1,2);

img = bg;
img(pos(1):pos(1)+size_f-1,pos(2):pos(2)+size_f-1) = fg;

k = Otsu(img)


%% Otsu's Method
function th = Otsu(img)
    [count,binLoc] = imhist(img);
    mul = count.*binLoc;
    total = sum(count,1);
    mean = sum(sum(img))/total;

    sig = zeros(256,1);

    for k = 1:256
       p_1 = sum(count(1:k),1)/total;
       p_2 = 1-p_1;
       mean_1 = sum(mul(1:k),1)/sum(count(1:k),1);
       mean_2 = sum(mul(k+1:256),1)/sum(count(k+1:256),1);
       sig(k,:) = p_1*(mean_1-mean).^2+p_2*(mean_2-mean).^2;
    end

    [~,th] = max(sig);

    output = sign(img-th)*255;

    % plot the inter-class variance
    figure
    subplot(2,2,1); imshow(img),title('Input image')
    subplot(2,2,2); imhist(img),title('Histogram')
    subplot(2,2,3);
    plot(binLoc,sig,'LineWidth',2),title('Inter-class variance')
    hold on
    xline(th,'Color','r','LineWidth',1.5)
    xlim([0,256])
    subplot(2,2,4); imshow(output),title('Output image')
end