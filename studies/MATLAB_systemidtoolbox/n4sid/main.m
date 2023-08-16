clear, clc

inputs_raw =  [importdata("input_1.txt"), importdata("input_2.txt"), ...
              importdata("input_3.txt"), importdata("input_4.txt")];

outputs_raw = [importdata("output_1.txt"), importdata("output_2.txt"), ...
              importdata("output_3.txt"), importdata("output_4.txt")];

[nt, n_inputs] = size(inputs_raw);
[~, n_outputs] = size(outputs_raw);

Ts = 0.01;
Form = 'modal';
Feedthrough = 1;
order = 12;
nModes = order/2-2;

for d = 1:16
    decimation = 1:d:nt;
    
    inputs = inputs_raw(decimation,:);
    outputs = outputs_raw(decimation,:);

    sys = n4sid(inputs,outputs,order,'Ts',Ts,'Form',Form,'Feedthrough',Feedthrough);
    
    [~,f] = modalfrf(sys);
    [fd, damping] = modalfit(sys,f,nModes);
    periods = 1./fd;
    
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;

    result_filename = "decimation_"+d+".mat";
    save(result_filename, "A", "B", "C", "D", "periods", "damping")
end
