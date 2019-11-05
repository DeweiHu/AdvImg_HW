clc,clear,close all

button = 'y';

while button == 'y' || button == 'Y'
    % Type in dir to load image 
    prompt = ('Please type in the image directory:\n');
    dir = input(prompt,'s');
    load(dir);
    
    if method == 1
        run('/home/dewei/Desktop/HW2019F/AdvImg_HW/HW7_Polynomial_backinter.m')
    else
        run('/home/dewei/Desktop/HW2019F/AdvImg_HW/HW7_TPS_backinter.m')
    end
    
    prompt = 'Do you want to try more?[y/n]:\n';
    button = input(prompt,'s');
    
    if button == 'n' || button == 'N'
        fprintf('Thanks\n');
        break
    end
end