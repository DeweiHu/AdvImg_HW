%clc, clear, close all

%load TestActiveShape_2019_2.mat
load Mean.mat
load Eigenvectors.mat

%% ICP
global Num, Num = 58;
global l_num, l_num = max(size(xy))-1;

b = zeros(116,1);
Y = xy';
Mean = reshape(mean,[2,Num]);
Eigen = zeros(116,116);
Eigen(:,1:4) = COEFF(:,1:4);

figure(1)
subplot(2,2,1)
plot(Y(1,:),Y(2,:),'-g*','MarkerSize',10,'LineWidth',2)
xlim([0,400])
ylim([0,400])
hold on
plot(Mean(1,:),Mean(2,:),'-r','LineWidth',2)
hold off
title('Initial')


subplot(2,2,2)
plot(Y(1,:),Y(2,:),'-g*','MarkerSize',10,'LineWidth',2),title('Progressing')
xlim([0,400])
ylim([0,400])
hold on
for i = 1:20
    X = reshape(mean+Eigen*b,[2,Num]);
    iterNum = 100;
    iter = 1;
    
    while iter <= iterNum
        Xc = Link(X(1,:),X(2,:),Y(1,:),Y(2,:));
        [R,t] = Register(X,Xc);
        X = R*X+t;
        iter = iter+1;
    end
    
    % Project back to mean space
    Xm = R\Xc-t;
    
    % Shape Adaptation
    b = COEFF'*(Xm(:)-mean);
    
    plot(X(1,:),X(2,:),'-b*','MarkerSize',10,'LineWidth',2)
    drawnow
end
hold off

subplot(2,2,3)
plot(Y(1,:),Y(2,:),'-g*','MarkerSize',10,'LineWidth',2)
xlim([0,400])
ylim([0,400])
hold on
plot(X(1,:),X(2,:),'-b','LineWidth',2)
hold off
title('Final Result')

pause(7)
close

function [R,t] = Register(X,Y)
    Xc = X-repmat(mean(X,2),1,size(X,2));
    Yc = Y-repmat(mean(Y,2),1,size(Y,2));
    [U,~,V] = svd(Xc*Yc');
    R = V*diag([1,det(V*U)])*U';
    t = mean(Y,2)-R*mean(X,2);
end

function Y = Link(pointx,pointy,xc,yc)
    global Num, global l_num
    Y = zeros(2,Num);
    % Find the pseudo-corresponding point for each sample points
    for k = 1:Num
        p = [pointx(k);pointy(k)];
        y = zeros(2,l_num);
        % Compute the distance of the point to each line element
        for l = 1:l_num
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
        % Find the minimum distance among these l_num distances
        for i = 1:l_num
            dist = sum((y-p).^2,1);
            idx = find(dist==min(dist));
            Y(:,k) = y(:,idx(1));
        end
    end
end