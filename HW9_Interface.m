clc,clear,close all

button = 'y';

iter = 1;
while button == 'y' || button == 'Y'
    % Type in dir to load image 
    prompt = ('Please type in the pointset directory:\n');
    dir = input(prompt,'s');
    load(dir);
    
    if iter == 1
        run('/home/dewei/Desktop/HW2019F/AdvImg_HW/HW9_Registration.m')
        iter = iter+1;
    end
    run('/home/dewei/Desktop/HW2019F/AdvImg_HW/HW9_main.m')
    
    
    prompt = 'Do you want to try more?[y/n]:\n';
    button = input(prompt,'s');
    
    if button == 'n' || button == 'N'
        fprintf('Thanks\n');
        break
    end
end