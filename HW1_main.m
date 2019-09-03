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
    
    if strcmp(req,'Otsu') || strcmp(req,'otsu')
        th = Otsu(img);
    elseif   strcmp(req,'Reddi') || strcmp(req,'reddi')
        [k1,k2] = Reddi(img);
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

function th = Otsu(img)
    [count,binLoc] = imhist(img);
    mul = count.*binLoc;
    total = sum(count,1);
    mean = sum(sum(img))/total;

    sig = zeros(256,1);

    for k = 1:256
       p_1 = sum(count(1:k),1)/total;
       p_2 = 1-p_1;
       mean_1 = sum(mul(1:k),1)/sum(count(1:k),1);
       mean_2 = sum(mul(k+1:256),1)/sum(count(k+1:256),1);
       sig(k,:) = p_1*(mean_1-mean).^2+p_2*(mean_2-mean).^2;
    end

    [~,th] = max(sig);

    output = sign(img-th)*255;

    % plot the inter-class variance
    figure
    subplot(2,2,1); imshow(img),title('Input image')
    subplot(2,2,2); imhist(img),title('Histogram')
    subplot(2,2,3);
    plot(binLoc,sig,'LineWidth',2),title('Inter-class variance')
    hold on
    xline(th,'Color','r','LineWidth',1.5)
    xlim([0,256])
    subplot(2,2,4); imshow(output),title('Output image')
end

function [k1,k2] = Reddi(img)
    [count,binLoc] = imhist(img);
    siz = size(img);
    output = uint8(zeros(siz(1),siz(2)));
    % initial
    k1 = 1./3*255;
    k2 = 2./3*255;
    e1 = 10;
    e2 = 10;

    figure
    subplot(2,2,1);imshow(img);title('Input image')
    subplot(2,2,2);imhist(img);xlim([0,256]);title('optimize threshold')
    hold on

    while abs(e1)>1 || abs(e2)>1
        e1 = .5*(RegionMean(1,k1,img)+RegionMean(k1+1,k2,img))-k1;
        e2 = .5*(RegionMean(k1+1,k2,img)+RegionMean(k2+1,256,img))-k2;
        % update
        k1 = k1+round(e1);
        k2 = k2+round(e2);
        % visualization
        xline(k1,'Color',rand(1,3),'LineWidth',1.5)
        hold on
        xline(k2,'Color',rand(1,3),'LineWidth',1.5)
        hold on
        pause(1)
    end

    for i = 1:siz(1)
        for j = 1:siz(2)
            if img(i,j) >= k2
                output(i,j) = 255;
            elseif img(i,j) < k2 && img(i,j) >= k1
                output(i,j) = 128;
            end
        end
    end

    subplot(2,2,3);imshow(output);title('Output image')
end
        
function m = RegionMean(a,b,img)
    [count,binLoc] = imhist(img);
    mul = count.*binLoc;
    m = sum(mul(a:b),1)./sum(count(a:b),1);
end
