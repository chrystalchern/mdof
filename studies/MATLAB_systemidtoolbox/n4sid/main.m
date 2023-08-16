clear, clc

nt = length(importdata("input_1.txt"));

inputs_raw =  [importdata("input_1.txt"), importdata("input_2.txt"), ...
              importdata("input_3.txt"), importdata("input_4.txt")];

outputs_raw = [importdata("output_1.txt"), importdata("output_2.txt"), ...
              importdata("output_3.txt"), importdata("output_4.txt")];

[nt, n_inputs] = size(inputs_raw);
[~, n_outputs] = size(outputs_raw);

for d = 1:16
    decimation = 1:d:nt;
    
    inputs = zeros(length(decimation),n_inputs);
    outputs = zeros(length(decimation),n_outputs); 
    
    for i = 1:length(input_channels)
        input = importdata("berke_"+input_channels(i)+".txt");
        input = input(decimation);
        inputs(:,i) = input;
    end
    
    for i = 1:length(output_channels)
        output = importdata("berke_"+output_channels(i)+".txt");
        output = output(decimation);
        outputs(:,i) = output;
    end
    
    Ts = 0.01;
    
    sys = n4sid(inputs,outputs,12,'Ts',Ts,'Form','modal','Feedthrough',1);
    
    [~,f] = modalfrf(sys);
    [fd, zeta] = modalfit(sys,f,6)
    periods = 1./fd
end
