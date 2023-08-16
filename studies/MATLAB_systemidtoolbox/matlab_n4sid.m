input_channels = [25 2 7 18];
output_channels = [23 13 15 20];

d = 1;
decimation = 1:d:nt;

nt = length(importdata("berke_25.txt"));
inputs = zeros(length(decimation),length(input_channels));
outputs = zeros(length(decimation),length(input_channels));

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
[fd, zeta] = modalfit(sys,f,3)
periods = 1./fd
