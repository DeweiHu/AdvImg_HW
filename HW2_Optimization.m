clc, clear, close all

load testhw2_3.mat

v = Optimize(testima,[60,120,180],[10,10,10],[.3,.3])

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
plot(x,y,'LineWidth',2),title('Initialization')

% Optimization
fun = @Err;

tic

options = optimset('Display','iter','MaxIter',1e+4,'OutputFcn',@Optimplot);
[v,~,exitflag] = fminsearch(fun,v_0,options);

toc

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