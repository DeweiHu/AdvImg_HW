clc, clear, close all

load testhw2_3.mat

par = EM(testima,[60,100,150],[10,10,10],[.3,.3,.4])

function par = EM(ipt,mu,sigma,prior)

%ipt = testima;
%global ipt; 

global cnt; global loc;

[cnt,loc] = imhist(ipt);

% Initialization
% mu = [1/4,1/2,3/4]*double(max(max(ipt)));
% sigma = [10,10,10];
% prior = [.3,.3,.4];

% Updating
iter = 0;
max_iter = 500;
M = 256^2;

figure
subplot(2,2,1),imshow(ipt),title('Original image')
subplot(2,2,2),imhist(ipt),title('Histogram')
subplot(2,2,3)
x = 0:1:255;
y = prior(1)*sum(cnt)*normpdf(x,mu(1),sigma(1))+...
    prior(2)*sum(cnt)*normpdf(x,mu(2),sigma(2))+...
    prior(3)*sum(cnt)*normpdf(x,mu(3),sigma(3));
plot(loc,cnt,'Color','b','LineWidth',2),title('Initialization')
hold on
plot(x,y,'Color','r','LineWidth',2)
hold off

while iter < max_iter
    
    % Visualization
    x = 0:1:255;
    y = prior(1)*sum(cnt)*normpdf(x,mu(1),sigma(1))+...
        prior(2)*sum(cnt)*normpdf(x,mu(2),sigma(2))+...
        prior(3)*sum(cnt)*normpdf(x,mu(3),sigma(3));
    subplot(2,2,4)
    plot(loc,cnt,'Color','b','LineWidth',2),title('Fitting histogram')
    hold on
    plot(x,y,'Color','r','LineWidth',2)
    drawnow
    hold off
    
    % weighted(pixel # for each class and each intensity)
    w = PixelClass(mu,sigma,prior);
    
    % update the prior
    prior = sum(w)/M;
    % update the mean
    mu = sum(w.*loc)./(M*prior);
    % update the variance
    var = sum(w.*(loc-mu).^2)./(M*prior);
    sigma = var.^.5;
    % update the buffer
    par = [mu,sigma,prior];
    
    % Error value
    error = Err(par);
    %disp(error)
    
    iter = iter+1;
    disp([iter,max_iter])
end


% function to compute wji
function w = PixelClass(mu,sigma,prior)
    %global cnt
    % Gaussian value for each pixel and class
    gauss = zeros(256,3);
    for i = 0:255
        for j = 1:3
            gauss(i+1,j) = normpdf(i,mu(j),sigma(j));
        end
    end
    % Gaussian value * prior = numerator
    numerator = gauss.*prior;
    w = (numerator./sum(numerator,2)).*cnt;
end

% error function
function error = Err(par)
    error = 0;
    for i = 1:256
        error = error + (cnt(i)-par(7)*sum(cnt)*normpdf(loc(i),par(1),par(4))...
        -par(8)*sum(cnt)*normpdf(loc(i),par(2),par(5))-(sum(cnt)...
        -par(9))*sum(cnt)*normpdf(loc(i),par(3),par(6)))^2;
    end
end

end