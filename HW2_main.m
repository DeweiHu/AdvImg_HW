clc, clear, close all

button = 'y';

while button == 'y' || button == 'Y'
    % Type in dir to load image 
    prompt = ('Please type in the image directory:\n');
    dir = input(prompt,'s');
    struct = load(dir);
    img = cell2mat(struct2cell(struct));

    % Select function
    prompt = ('Which method do you want to try?\n');
    req = input(prompt,'s');
    
    if strcmp(req,'Optimization') || strcmp(req,'optimization')
        prompt = ('Input 3 initial values of mean:\n');
        mu = input(prompt);
        
        prompt = ('Input 3 initial values of standard deviation:\n');
        sigma = input(prompt);
        
        prompt = ('Input 2 initial values of prior probability:\n');
        prior = input(prompt);
        
        par = Optimize(img,mu,sigma,prior);
        
    elseif   strcmp(req,'EM') || strcmp(req,'em')
        prompt = ('Input 3 initial values of mean:\n');
        mu = input(prompt);
        
        prompt = ('Input 3 initial values of standard deviation:\n');
        sigma = input(prompt);
        
        prompt = ('Input 3 initial values of prior probability:\n');
        prior = input(prompt);
        
        par = EM(img,mu,sigma,prior);
        
    else
        fprintf('Request invalid\n');
    end
   
    prompt = 'Do you want to try more?[y/n]:\n';
    button = input(prompt,'s');
    
    if button == 'n' || button == 'N'
        fprintf('Thanks\n');
        break
    end
end

function v = Optimize(ipt,mu,sigma,prior)
global cnt,global loc
[cnt,loc] = imhist(ipt);

% show the image and histogram
figure
subplot(2,2,1),imshow(ipt);title('Original image');
subplot(2,2,2),imhist(ipt);title('Histogram')

% Define the vector to optimize
v_0 = [mu,sigma,prior*sum(cnt)];

x = 0:1:255;

y = v_0(7)*normpdf(x,v_0(1),v_0(4))...
    +v_0(8)*normpdf(x,v_0(2),v_0(5))...
    +(sum(cnt)-v_0(7)-v_0(8))*normpdf(x,v_0(3),v_0(6));

subplot(2,2,3)
plot(x,cnt,'LineWidth',2)
hold on
plot(x,y,'LineWidth',2),title('Initialization'),legend('data hist','initial')
hold off

% Optimization
fun = @Err;

tic

options = optimset('Display','iter','MaxIter',1e+4,'OutputFcn',@Optimplot);
[v,~,exitflag] = fminsearch(fun,v_0,options);

toc

means = v(1:3);
sd = v(4:6);
prior = [v(7:8)/sum(cnt),1-v(7)/sum(cnt)-v(8)/sum(cnt)];

fprintf('Mean of Gaussian Distributions:\n')
disp(means)
fprintf('Standard deviation of Gaussian Distributions:\n')
disp(sd)
fprintf('Prior of Gaussian Distributions:\n')
disp(prior)

function stop = Optimplot(v,~,state)
    stop = false;
    a = 0:1:255;
    b = v(7)*normpdf(a,v(1),v(4))...
        +v(8)*normpdf(a,v(2),v(5))+(sum(cnt)-v(7)-v(8))*normpdf(a,v(3),v(6));
    
    subplot(2,2,4)
    plot(a,cnt,'Color','b','LineWidth',2)
    hold on
    plot(a,b,'Color','r','LineWidth',2),title('Fitting Result')
    drawnow
    hold off
    
    if state == 'done'
        stop = true;
    end
end

function error = Err(v)
    error = 0;
    for i = 1:256
        error = error + (cnt(i)-v(7)*normpdf(loc(i),v(1),v(4))...
        -v(8)*normpdf(loc(i),v(2),v(5))-(sum(cnt)-v(7)-v(8))*normpdf(loc(i),v(3),v(6)))^2;
    end
end

end

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
max_iter = 300;
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
plot(x,y,'Color','r','LineWidth',2),legend('data hist','initial')
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
    
    iter = iter+1;
    disp([iter,max_iter])
end
% show final reaults
means = par(1:3);
sd = par(4:6);
prior = par(7:9);

fprintf('Mean of Gaussian Distributions:\n')
disp(means)
fprintf('Standard deviation of Gaussian Distributions:\n')
disp(sd)
fprintf('Prior of Gaussian Distributions:\n')
disp(prior)


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

end