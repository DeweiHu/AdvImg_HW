%clc,clear,close all

%load TestHw6_1.mat

global im_fix, global im_move, global lower, global upper, global bin
im_fix = mitest;
im_move = mitestrot;
lower = minval;
upper = maxval;
bin = nbins;

theta = 0;

fun = @Error;

options = optimset('Display','iter','TolFun',1e-8,'TolX',1e-8,...
            'MaxFunEvals',1000,'OutputFcn',@Optimplot);

angle = fminsearch(fun,theta,options)*1e5;

function stop = Optimplot(theta,~,state)
    global im_fix, global im_move
    stop = false;    
    
    J = imrotate(im_move,theta*1e5,'bilinear','crop');
    Hgram = Joint_Histogram(J);
    
    figure(1)
    colormap jet
    subplot(2,2,1),imagesc(im_fix),title('Original Image')
    subplot(2,2,2),imagesc(im_move),title('Rotated Image')
    subplot(2,2,3),imagesc(J),title('Evolving Process')
    subplot(2,2,4),imagesc(Hgram),title('Joint-Histogram')
    hold on
    
    drawnow
    
    pause(.5)
    
    if state == 'done'
        hold off
        stop = true;
    end
end

function error = Error(theta)
    global im_move
    rotated = imrotate(im_move,theta*1e5,'bilinear','crop');
    Hgram = Joint_Histogram(rotated);
    error = -Mutual_Information(Hgram);
end

function Hgram = Joint_Histogram(rotated)
    global im_fix, global lower, global upper, global bin
    
    % Fit in number of bins
    bin_size = (upper-lower)/bin;

    % Initialize the lower bound
    low = lower;
    im_x = zeros(size(im_fix));
    im_y = zeros(size(rotated));

    for i = 1:bin
        % Mark pixels of interest with 1~64 correspond to bins
        im_x(im_fix>=low & im_fix<low+bin_size) = i;
        im_y(rotated>=low & rotated<low+bin_size) = i;
        low = low+bin_size;
    end

    % Make Joint-Histogram
    Hgram = zeros([bin,bin]);
    for i = 1:bin

        % Get locations that in ith bin in im_x
        [row,col] = find(im_x==i);  
        len = max(size(row));

        % Get bin number of these locations in im_y
        for j = 1:len
            bin_num = im_y(row(j),col(j));
            if bin_num ~= 0
                Hgram(i,bin_num) = Hgram(i,bin_num)+1;
            end
        end

    end
end

function MI = Mutual_Information(Hgram)
    % Normalize to get prob and mariginal prob
    pxy = Hgram / sum(Hgram,'all');
    px = sum(pxy,2);  % marginal all y/columns
    py = sum(pxy,1);  % marginal all x/rows
    px_py = px*py;
    nzo = pxy > 0;
    MI = sum(pxy(nzo).*log(pxy(nzo)./px_py(nzo)));
end