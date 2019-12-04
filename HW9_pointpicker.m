clc, clear, close all

hand1 = imread("hand1.tif");
hand2 = imread("hand2.tif");
hand3 = imread("hand3.tif");
hand4 = imread("hand4.tif");
hand5 = imread("hand5.tif");
hand6 = imread("hand6.tif");
hand7 = imread("hand7.tif");
hand8 = imread("hand8.tif");
hand9 = imread("hand9.TIF");
hand10 = imread("hand10.tif");

n =58;
pointset10 = zeros(n,2);
figure
imshow(hand10)
hold on
for i = 1:n
    [x,y] = ginput(1);
    pointset10(i,:) = [x,y];
    plot(pointset10(i,1),pointset10(i,2),'r.','MarkerSize',25);
end
hold off
%%
