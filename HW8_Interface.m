clc,clear,close all

button = 'y';

while button == 'y' || button == 'Y'
    % Type in dir to load image 
    prompt = ('Please type in the image directory:\n');
    dir = input(prompt,'s');
    load(dir);
    
    run('/home/dewei/Desktop/HW2019F/AdvImg_HW/HW8_main.m')
    
    prompt = 'Do you want to try more?[y/n]:\n';
    button = input(prompt,'s');
    
    if button == 'n' || button == 'N'
        fprintf('Thanks\n');
        break
    end
end