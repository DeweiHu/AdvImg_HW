%clc, clear, close all 

%load HW7_practice6.mat

%% Least-Square Solution
global n, n = max(size(xc1));

% Form matrix L_forward, L_back
K_forward = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i == j
            K_forward(i,j) = 0;
        else
            r = (xc1(i)-xc1(j))^2+(yc1(i)-yc1(j))^2;
            K_forward(i,j) = r*log(r);
        end
    end
end

K_back = zeros(n,n);
for i = 1:n
    for j = 1:n
        if i == j
            K_back(i,j) = 0;
        else
            r = (xc2(i)-xc2(j))^2+(yc2(i)-yc2(j))^2;
            K_back(i,j) = r*log(r);
        end
    end
end

P_forward = [ones(n,1),xc1',yc1'];
P_back = [ones(n,1),xc2',yc2'];

O = zeros(3,3);
L_forward = [P_forward,K_forward;O,P_forward'];
L_back = [P_back,K_back;O,P_back'];

Y_forward = zeros(n+3,2);
Y_forward(1:n,:) = [xc2',yc2'];
Y_back = zeros(n+3,2);
Y_back(1:n,:) = [xc1',yc1'];

% Transform_forward
Transform_forward = pinv(L_forward)*Y_forward;
Transform_back = pinv(L_back)*Y_back;

%% Apply Transform_back and do interpolation for image
global num, num = 256*256;

PK_back = ones(num,n+3);

for col = 1:256
    y_vec = [1:256]';
    x_vec = col*ones(256,1);
    PK_back((col-1)*256+1:col*256,2:3) = [x_vec,y_vec];
end

for i = 1:num
    dist_list = zeros(1,n);
    for j = 1:n
        r = (PK_back(i,2)-xc2(j))^2+(PK_back(i,3)-yc2(j))^2;
        dist_list(j) = r*log(r);
    end
    PK_back(i,4:end) = dist_list;
end

Non_Int = PK_back*Transform_back;

for i = 1:num
    if Non_Int(i,1)<0 || Non_Int(i,2)<0
        Non_Int(i,:) = [0,0];
    end
end

[X,Y] = meshgrid(1:1:256);
vec_rg = interp2(X,Y,source,Non_Int(:,1),Non_Int(:,2));
img_rg = reshape(vec_rg,size(source));

% Apply Transform_forward for points
points_rg = L_forward*Transform_forward;
points_rg = points_rg(1:n,:);

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
subplot(2,2,3),imagesc(img_rg),hold on,scatter(points_rg(:,1),points_rg(:,2),700,'y.')
hold on,scatter(xc2,yc2,300,'r.')
hold off,title('Source to Target')
subplot(2,2,4),imagesc(img_field),title('Deformation Field')