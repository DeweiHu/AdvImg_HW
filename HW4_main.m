%clc,clear,close all

%load testhw4_5.mat

%% Step 1: Interpolate on initial contour
[px,py] = EvenSpace(initcontourx,initcontoury,nsample);

figure(1)
imagesc(binim),title('Initialization')
hold on 
plot(initcontourx,initcontoury,'-yo','MarkerSize',20,'LineWidth',2)
hold on
plot(px,py,'-.r*','MarkerSize',10,'LineWidth',2)
hold off

pause(2)
close

%% Step 2: Form A matrix
% input: alpha, beta
% output: A

unit = [beta,-4*beta-alpha,6*beta+2*alpha,-4*beta-alpha,beta];
A = zeros(nsample,nsample);
A(1,1:3) = unit(3:5);
A(1,nsample-1:nsample) = unit(1:2);

% from 1st row, get the pentadiagonal matrix A
for i = 2:nsample
    for j = 1:nsample
        if j == 1
            A(i,j) = A(i-1,nsample);
        else
            A(i,j) = A(i-1,j-1);
        end
    end
end

%% Step 3: Specify the External Force
% input: forcetype, binim, std, support
% output: Fext_x, Fext_y

%Get E_ext by convolve with a first order derivative of Gaussian filter
gauss = fspecial('gaussian',2*support,std);

if forcetype == 1
    [Fext_x, Fext_y] = Gradient(binim,gauss);
end
if forcetype == 2
    [Fext_x, Fext_y] = Dist_map(binim);
end
if forcetype == 3
    [Fext_x,Fext_y] = GVF(binim,itergvf);
end

%% Step 4: Iteration
I = eye(nsample);
iter = 1;

while iter <= Niter
    % [OPEN THE RING!] and make it a vector
    buffer_x = px(1:nsample)';
    buffer_y = py(1:nsample)';
    
    % 2d interpolation to find the F_ext at points
    F_ext = Interpolate(buffer_x,buffer_y,Fext_x,Fext_y,nsample);
    
    % Balloon force
    [bal_x,bal_y] = Balloon(buffer_x,buffer_y);
    
    % Update the position of points on the contour
    buffer_x = (I+gamma*A)\(buffer_x+extcoef*gamma*F_ext(:,1))+balcoef*bal_x;
    buffer_y = (I+gamma*A)\(buffer_y+extcoef*gamma*F_ext(:,2))+balcoef*bal_y;
    
    % Resample to make points evenly space [CLOSE THE RING!]
    Re_x = [buffer_x;buffer_x(1)]';
    Re_y = [buffer_y;buffer_y(1)]';
    [px,py] = EvenSpace(Re_x,Re_y,nsample);
    
    % Visualization
    imagesc(binim),title('Evolving Process')
    hold on
    plot(px,py,'-rx','LineWidth',2,'MarkerSize',10)
    drawnow 
    
    iter = iter+1;

end


function [Gx,Gy] = Gradient(binim,gauss)
    [dx,dy] = gradient(gauss);

    Gx = conv2(binim,dx,'same');
    Gy = conv2(binim,dy,'same');
end
    
function [Fext_x, Fext_y] = Dist_map(binim)
    dmap = bwdist(binim);
    [grad_x,grad_y]=gradient(dmap);
    Fext_x = -double(grad_x);
    Fext_y = -double(grad_y);
end

function [u,v] = GVF(binim,itergvf)
    dim = size(binim);
    
    % compute b(x,y),c1(x,y),c2(x,y)
    [fx,fy] = gradient(binim);
    b = fx.^2+fy.^2;
    c1 = b.*fx;
    c2 = b.*fy;
    r = 0.05;
    kernel = r*[0,1,0;1,-4,1;0,1,0];

    % initialize u(x,y) & v(x,y)
    u = zeros(dim);
    v = zeros(dim);
    iter = 1;

    while iter <= itergvf
        % Update u(x,y) & v(x,y)
        u = conv2(u,kernel,'same')+u-b.*u+c1;
        v = conv2(v,kernel,'same')+v-b.*v+c2;
        iter = iter+1;
    end
end

function [px,py] = EvenSpace(initcontourx,initcontoury,nsample)
    Np = max(size(initcontourx));
    dist_vec = zeros(Np,1);

    for i = 2:Np
        dist = Dist([initcontourx(i),initcontoury(i)],...
            [initcontourx(i-1),initcontoury(i-1)]);
        dist_vec(i) = dist+dist_vec(i-1);
    end

    delta = dist_vec(end)/nsample;

    px = interp1(dist_vec,initcontourx,0:delta:dist_vec(end)); 
    py = interp1(dist_vec,initcontoury,0:delta:dist_vec(end));
    
    function dist = Dist(p1,p2)
    dist = sqrt((p1(1,1)-p2(1,1))^2+(p1(1,2)-p2(1,2))^2);
    end

end

function F_ext = Interpolate(new_x,new_y,Fext_x,Fext_y,nsample)
    % Get F_ext for points on the contour by 2d interpolation
    F_ext = zeros(nsample,2);
    for i = 1:nsample
        % Define the location of 4 corners
        lx = floor(new_x(i));
        hx = ceil(new_x(i));
        ly = floor(new_y(i));
        hy = ceil(new_y(i));
        [X,Y] = meshgrid(ly:hy,lx:hx);

        % value in Fext_x and Fext_y
        Vx = [Fext_x(ly,lx),Fext_x(ly,hx);Fext_x(hy,lx),Fext_x(hy,hx)];
        Vy = [Fext_y(ly,lx),Fext_y(ly,hx);Fext_y(hy,lx),Fext_y(hy,hx)];

        % position to be interpolated
        Xq = new_y(i);
        Yq = new_x(i);

        % interpolation
        F_ext(i,1) = interp2(X,Y,Vx,Xq,Yq);
        F_ext(i,2) = interp2(X,Y,Vy,Xq,Yq);
    end
end

function [opt_x,opt_y] = Balloon(buffer_x,buffer_y)
    kernel = [1;0;-1];

    ipt_x = [buffer_x(end);buffer_x;buffer_x(1)];
    ipt_y = [buffer_y(end);buffer_y;buffer_y(1)];

    nx = -conv(ipt_x,kernel,'valid');
    ny = conv(ipt_y,kernel,'valid');

    coef = zeros(size(nx));
    for i = 1:max(size(nx))
        coef(i) = 1/sqrt(nx(i)^2+ny(i)^2);
    end

    opt_x = coef.*ny;
    opt_y = coef.*nx;
    
end