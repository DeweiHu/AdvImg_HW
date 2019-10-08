clc, clear, close all

button = 'y';

while button == 'y' || button == 'Y'
    % Type in dir to load image 
    prompt = ('Please type in the image directory:\n');
    dir = input(prompt,'s');
    load(dir);
    
    b_field = bf_correction(ima_att,param);
    
    prompt = 'Do you want to try more?[y/n]:\n';
    button = input(prompt,'s');
    
    if button == 'n' || button == 'N'
        fprintf('Thanks\n');
        break
    end
end
    

function b_field = bf_correction(ima_att,param)

% Pre-process the input image and params
ipt_img = double(ima_att);
mu = [param(1),param(4)];
sigma = [param(2),param(5)];
prior = [param(3),param(6)];

% Visualiza the input information
[cnt,loc] = imhist(uint8(ipt_img)); 
figure(1)
subplot(2,2,1),plot(loc,cnt/sum(cnt),'LineWidth',2),title('Original Distribution')
subplot(2,2,2),mesh(ipt_img),title('Original Image 3D surf')
subplot(2,2,3),imshow(ind2rgb(uint8(ipt_img),jet)),title('Original image')
pause(2)
close
%% define the evolving plot
figure(2)
subplot(1,2,1),imshow(ind2rgb(uint8(ipt_img),jet)),title('input')

C_Total = zeros(6,1);

it = 10;
while it>0 
    ipt = uint8(ipt_img);    
    % Find the theta by optimization
    par = EM(ipt,mu,sigma,prior);
    theta = par;
    
    % Prepare the y matrix
    % reshape the ipt into vector
    ipt_vec = ipt_img(:);
    
    % compute the posterior probability
    pst = zeros(sum(cnt),2);
    for i = 1:sum(cnt)
        [pst(i,1),pst(i,2)] = Post(ipt_vec(i),theta);
    end
    
    % get the alpha
    var = [theta(3)^2,theta(4)^2];
    alpha = pst./var;
    
    % get g tilde
    mu = [theta(1),theta(2)];
    g_tilde = sum(alpha.*mu,2)./sum(alpha,2);
    
    % get y vector
    Y = ipt_vec - g_tilde;
    
    % Prepare the basis function matrix
    basis = zeros(sum(cnt),6);
    for y = 1:256
        for x = 1:256
            basis((y-1)*256+x,:) = [x^2,y^2,x*y,x,y,1];
        end
    end
    
    % get solution of C
    C = pinv(basis)*Y;
    
    % Do the correction
    opt_vec = ipt_vec-basis*C;
    opt_img = reshape(opt_vec,[256,256]);
    opt = uint8(opt_img);

    subplot(1,2,2),imshow(ind2rgb(opt,jet)),title('output of 10 iterations')
    drawnow
    pause(.4)
    
    % Update the image and the Total C
    C_Total = C_Total+C;
    ipt_img = opt_img;
    
    it = it-1;
end
close
%% Compute the bias field
b_field = reshape(basis*C_Total,[256,256]);

% Visualize 2 slices in the middle of the image
% Original
img1 = ima_att(:,128);
img2 = ima_att(128,:);
img1_c = opt_img(:,128);
img2_c = opt_img(128,:);
bf1 = b_field(:,128);
bf2 = b_field(128,:);

figure(3)
subplot(3,2,1),plot(img1,'LineWidth',1.3);title('Original Image (:,128)')
subplot(3,2,2),plot(img1_c,'LineWidth',1.3);title('Corrected Image (:,128)')
subplot(3,2,3),plot(img2,'LineWidth',1.3);title('Original Image (128,:)')
subplot(3,2,4),plot(img2_c,'LineWidth',1.3);title('Corrected Image (128,:)')
subplot(3,2,5),plot(bf1,'LineWidth',1.5);title('Correction Surface (:,128)')
subplot(3,2,6),plot(bf2,'LineWidth',1.5);title('Correction Surface (128,:)')

% Print the result of the algorithm
x = 0:1:255;
y = theta(5)*normpdf(x,theta(1),theta(3))+...
    theta(6)*normpdf(x,theta(2),theta(4));
[cnt,loc] = imhist(uint8(ima_att));

figure(4)
subplot(2,2,1),plot(loc,cnt/sum(cnt),'LineWidth',2)
hold on
plot(x,y,'r','LineWidth',2);title('Original and Updated Distribution')
hold off
subplot(2,2,2),mesh(b_field),colormap hsv;title('Bias Field')
subplot(2,2,3),imshow(ind2rgb(uint8(ima_att),jet));title('Original Image')
subplot(2,2,4),imshow(ind2rgb(opt,jet));title('Image After Correction')

function [post1,post2] = Post(g_i,theta)
    class1 = normpdf(g_i,theta(1),theta(3))*theta(5);
    class2 = normpdf(g_i,theta(2),theta(4))*theta(6);
    factor = class1+class2;
    post1 = class1/factor;
    post2 = class2/factor;
end

function par = EM(ipt,mu,sigma,prior)

    [cnt,loc] = imhist(ipt);

    % Updating
    iter = 0;
    max_iter = 100;
    M = 256^2;

    while iter < max_iter

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
    end
    fprintf('fitting done, updating...\n')


    % function to compute wji
    function w = PixelClass(mu,sigma,prior)
        %global cnt
        % Gaussian value for each pixel and class
        gauss = zeros(256,2);
        for i = 0:255
            for j = 1:2
                gauss(i+1,j) = normpdf(i,mu(j),sigma(j));
            end
        end
        % Gaussian value * prior = numerator
        numerator = gauss.*prior;
        w = (numerator./sum(numerator,2)).*cnt;
    end

end

end