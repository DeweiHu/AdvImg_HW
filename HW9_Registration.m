%clc, clear, close all

filenames = {'pointset1.mat','pointset2.mat','pointset3.mat','pointset4.mat'...
    'pointset5.mat','pointset6.mat','pointset6.mat','pointset7.mat'...
    'pointset8.mat','pointset9.mat','pointset10.mat'};
for i = 1:numel(filenames)
    load(filenames{i})
end

global Num, Num = 10; 
data = {pointset1,pointset2,pointset3,pointset4,pointset5...
    pointset6,pointset7,pointset8,pointset9,pointset10};

global pNum, pNum = max(size(pointset1));

%% Rigid Point Registration

% Stage 1
fixed = data{1}';
update = cell(1,Num);
update{1} = fixed;

for i = 2:Num
    move = data{i}';
    [R,t] = Register(move,fixed);
    update{i} = R*move+t;
end

% Stage 2
iterNum = 100;
iter = 0;
while iter < iterNum
    averset = GetAverage(update);
    for i = 1:Num
        move = update{i};
        [R,t] = Register(move,averset);
        update{i} = R*move+t;
    end
    iter = iter+1;
end


figure(1)
subplot(2,2,1)
for i = 1:Num
    plot(data{i}(:,1),data{i}(:,2),'-b*','MarkerSize',10,'LineWidth',2)
    hold on
end
hold off
title('Input shapes')

subplot(2,2,2)
for i = 1:Num
    plot(update{i}(1,:),update{i}(2,:),'-b*','MarkerSize',10,'LineWidth',2)
    hold on
end
hold off
title('Realigned shapes')

subplot(2,2,3)
plot(averset(1,:),averset(2,:),'-r','LineWidth',2)
xlim([0 400])
ylim([0,400])
title('Mean shape')

pause(3)
close

%% Eigenshapes
mean = averset(:);
A = zeros(2*pNum,Num);
for i = 1:Num
    A(:,i) = update{i}(:)-mean;
end
cov = A*A'/(Num-1);

% PCA
[COEFF,latant,explained] = pcacov(cov);

% Eigenshapes
Eigenshape = cell(1,Num-1);
factor = -100;
step = 10;
while factor <= 100
    for i = 1:Num-1
        vector = mean+COEFF(:,i)*factor;
        Eigenshape{i} = reshape(vector,[2,pNum]);
    end
    factor = factor+step;

    figure(2)
    subplot(2,2,1)
    plot(Eigenshape{1}(1,:),Eigenshape{1}(2,:),'-b','LineWidth',2)
    xlim([0,400])
    ylim([0,400])
    title('First Eigenshape')
    subplot(2,2,2)
    plot(Eigenshape{2}(1,:),Eigenshape{2}(2,:),'-b','LineWidth',2)
    xlim([0,400])
    ylim([0,400])
    title('Second Eigenshape')
    subplot(2,2,3)
    plot(Eigenshape{3}(1,:),Eigenshape{3}(2,:),'-b','LineWidth',2)
    xlim([0,400])
    ylim([0,400])
    title('Third Eigenshape')
    subplot(2,2,4)
    plot(Eigenshape{4}(1,:),Eigenshape{4}(2,:),'-b','LineWidth',2)
    xlim([0,400])
    ylim([0,400])
    title('Fourth Eigenshape')
    
    drawnow
    pause(0.3)
end
pause(3)
close


function [R,t] = Register(X,Y)
    Xc = X-repmat(mean(X,2),1,size(X,2));
    Yc = Y-repmat(mean(Y,2),1,size(Y,2));
    [U,~,V] = svd(Xc*Yc');
    R = V*diag([1,det(V*U)])*U';
    t = mean(Y,2)-R*mean(X,2);
end

function averset = GetAverage(update)
    global Num;
    sumset = zeros(2,58);
    for i = 1:Num
        sumset = sumset+update{i};
    end
    averset = sumset/Num;
end

