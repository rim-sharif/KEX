
fname = "FAL_6_viterbiremapped.hdf5";
info = h5info(fname); % Retreive information from HDF5 container 
reads = info.Groups.Groups; % Retreive DNA reads

N = length(reads); % Number of reads
n = randi(N); % Pic a random read
read_name = reads(n).Name; % Obtain name of picked read

% Load data from read
Signal = double(h5read(fname,strcat(read_name,"/Dacs"))); %current signal
x = Signal; x = x - mean(x); x = x ./ std(x);       %normerar signalen (s.k. z-score normalization)
Ref_to_signal = double(h5read(fname,strcat(read_name,"/Ref_to_signal")));   %linjerna
Reference = double(h5read(fname,strcat(read_name,"/Reference")));   % bokstäverna

% Pick random segment from read
DNA = ['A','C','G','T']; % Neucleotide letters for reference
M = length(Reference); % Number of bases in reference
L = 20; % Number of bases to plot
m1 = randi(M-L);        %startpunkt
m2 = m1 + L;            %slutpunkt

% Plot signal and reference
figure(1); clf; hold on;
xb = x(Ref_to_signal(m1)+1:Ref_to_signal(m2+1));        %varför +1?
plot(xb,'b');

% Plot borders and reference
hold on;
t1 = Ref_to_signal(m1);
for m = m1:m2
    plot((Ref_to_signal(m) - t1 + 1)*[1 1],[-5 5],'r'); % Plot border line
    text((Ref_to_signal(m)+Ref_to_signal(m+1))/2-t1, 4, DNA(Reference(m)+1)); % Plot reference base
end

% Picture formatting
box on;
axis([1 (Ref_to_signal(m2+1)-t1) -5 5]);
