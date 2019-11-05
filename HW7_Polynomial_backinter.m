%clc, clear, close all

%load HW7_practice2.mat

%% Least-Square Solution
n = max(size(xc1));
Y_back = [xc1',yc1'];
Y_forward = [xc2',yc2'];

% Form A matrix by the 2nd order polynomial
A_back = [ones(n,1),xc2',yc2',xc2'.*yc2',xc2'.*xc2',yc2'.*yc2'];
A_forward = [ones(n,1),xc1',yc1',xc1'.*yc1',xc1'.*xc1',yc1'.*yc1'];

% Coefficients of the polynomial transformation T
Transform_back = pinv(A_back)*Y_back;
Transform_forward = pinv(A_forward)*Y_forward;

%% Apply Transform_back and do Interpolation for image
global num, num = 256*256;

% Form matrix of integer coordinates(scatter coordinate)
Int = ones(num,6);
for col = 1:256
    y_vec = [1:256]';
    x_vec = col*ones(256,1);
    Int((col-1)*256+1:col*256,2:end) = [x_vec,y_vec,x_vec.*y_vec,...
                       x_vec.*x_vec,y_vec.*y_vec];
end

% Get matrix of non-integer coordinates(scatter coordinate)
Non_Int = Int*Transform_back;

% Filter out the negative coordinate values  
for i = 1:num
    if Non_Int(i,1)<0 || Non_Int(i,2)<0
        Non_Int(i,:) = [0,0];
    end
end
% Interpolate
[X,Y] = meshgrid(1:1:256);
vec_rg = interp2(X,Y,source,Non_Int(:,1),Non_Int(:,2));
img_rg = reshape(vec_rg,size(source));

% Apply Transform_forward for points
points_rg = A_forward*Transform_forward;

%% Grids to show deformation field
Grid = zeros(256,256);
for i = 1:256
    if mod(i,8)==0
        Grid(i,:) = 255;
        Grid(:,i) = 255;
    end
end
[X,Y] = meshgrid(1:1:256);
vec_field = interp2(X,Y,Grid,Non_Int(:,1),Non_Int(:,2));
img_field = reshape(vec_field,size(source));

%% Visualization
figure(1)
colormap gray
subplot(2,2,1),imagesc(source),hold on,scatter(xc1,yc1,500,'g.')
hold off,title('Source')
subplot(2,2,2),imagesc(target),hold on,scatter(xc2,yc2,500,'r.')
hold off,title('Target')
subplot(2,2,3),imagesc(img_rg),hold on,scatter(points_rg(:,1),points_rg(:,2),500,'y.')
hold on,scatter(xc2,yc2,500,'r.')
hold off,title('Source to Target')
subplot(2,2,4),imagesc(img_field),title('Deformation Field')


