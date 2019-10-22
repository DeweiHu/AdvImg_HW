%clc,clear,close all

%load testhw5_3.mat

%% Step 1: Visualize the initial state
global dist_map, global sample, global num
sample = [pointx;pointy];
num = max(size(pointx));

I = zeros(256,256);
area = double(roipoly(I,xc,yc));
binim = edge(area,'sobel');
dist_map = double(bwdist(binim));

figure(1)
imagesc(dist_map),title("Initial State")
hold on
plot(pointx,pointy,'-y*','MarkerSize',10,'LineWidth',2)
hold off

pause(2)
close

%% Step 2: Optimization

dist = Interpolate(pointx,pointy);

% Define the vector of parameters to optimize
theta = 0;
tx = 0;
ty = 0;
v_0 = [theta,tx,ty];

fun = @Error;

options = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8,...
            'MaxFunEvals',1000,'OutputFcn',@Optimplot);
v = fminsearch(fun,v_0,options);

function stop = Optimplot(v,~,state)
    global dist_map
    stop = false;    
    reg = Transform(v*1000);
    imagesc(dist_map),title('Evolving Process')
    hold on
    plot(reg(1,:),reg(2,:),'-y*','MarkerSize',10,'LineWidth',2)
    drawnow
    
    if state == 'done'
        stop = true;
    end
end

function error = Error(v)
    moved = Transform(1000*v);
    dist = Interpolate(moved(1,:),moved(2,:));
    error = dot(dist,dist);
end

function moved = Transform(v)
    global sample, global num
    moved = zeros(size(sample));
    R = [cos(v(1)),-sin(v(1));sin(v(1)),cos(v(1))];
    T = [v(2);v(3)];
    for i = 1:num
        moved(:,i) = R*sample(:,i)+T;
    end
end

function dist = Interpolate(pointx,pointy)
    global num, global dist_map
    dist = zeros(num,1); 
    
    % uniform the index (plot to image)
    X = pointy;
    Y = pointx;

    for i = 1:num
        % Define the location of 4 corners
        lx = floor(X(i));
        hx = ceil(X(i));
        ly = floor(Y(i));
        hy = ceil(Y(i));
        
        t = (X(i)-lx)/(hx-lx);
        u = (Y(i)-ly)/(hy-ly);
        % bi-linear interpolation
        dist(i) = (1-t)*(1-u)*dist_map(lx,ly)+t*(1-u)*dist_map(hx,ly)+...
                  t*u*dist_map(hx,hy)+(1-t)*u*dist_map(lx,hy);
    end
end
