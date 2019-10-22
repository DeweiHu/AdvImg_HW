%clc,clear,close all

%load testhw5_6.mat

%% Visualize the Initial State
global sample, global num
sample = [pointx;pointy];
num = max(size(pointx));

I = zeros(256,256);
area = double(roipoly(I,xc,yc));
binim = edge(area,'sobel');

figure(1)
imagesc(binim),title('Initial State')
hold on
plot(pointx,pointy,'-ro','MarkerSize',10,'LineWidth',2)
hold on
scatter(xc,yc,100,'yo')
hold off

pause(2)
close
%%
P = [pointx;pointy];
X = P;
iter = 1;

while iter < 50
    % Step 1: Find the pseudo-correspondence: closest points Y
    Y = Link(X(1,:),X(2,:),xc,yc);

    % Step 2: Pairwise Rigid Registration
    [R,t] = Register(X,Y);

    % Step3: Update
    X = R*X+t;
    
    % plot process
    figure(2)
    imagesc(binim),title('Evolving Process')
    hold on
    plot(X(1,:),X(2,:),'-r*','MarkerSize',10,'LineWidth',2)
    drawnow

    iter = iter+1;
end

function [R,t] = Register(X,Y)
    Xc = X-repmat(mean(X,2),1,size(X,2));
    Yc = Y-repmat(mean(Y,2),1,size(Y,2));
    [U,~,V] = svd(Xc*Yc');
    R = V*diag([1,det(V*U)])*U';
    t = mean(Y,2)-R*mean(X,2);
end

function Y = Link(pointx,pointy,xc,yc)
    global num
    Y = zeros(2,num);
    % Find the pseudo-corresponding point for each sample points
    for k = 1:num
        p = [pointx(k);pointy(k)];
        y = zeros(2,5);
        % Compute the distance of the point to each line element
        for l = 1:5
            X1 = [xc(l);yc(l)];
            X2 = [xc(l+1);yc(l+1)];
            nm = (p(1)-X1(1))*(X2(1)-X1(1))+(p(2)-X1(2))*(X2(2)-X1(2));
            dm = dot(X2-X1,X2-X1);
            u = nm/dm;
            if u >= 1
                y(:,l) = X2;
            elseif u <= 0
                y(:,l) = X1;
            else
                y(:,l) = X1+u*(X2-X1);
            end
        end
        % Find the minimum distance among these 5 distances
        for i = 1:5
            dist = sum((y-p).^2,1);
            idx = find(dist==min(dist));
            Y(:,k) = y(:,idx(1));
        end
    end
end