
% ---------------------Code supervisors------------------------------
fname = "FAL_6_viterbiremapped.hdf5";   %Container of data
info = h5info(fname); % Retreive information from HDF5 container 
reads = info.Groups.Groups; % Retreive DNA reads

N = length(reads); % Number of reads
n = randi(N); % Pic a random read
read_name = reads(n).Name; % Obtain name of picked read

% Load data from DNA read
Signal = double(h5read(fname,strcat(read_name,"/Dacs"))); %Ström signal
x = Signal; x = x - mean(x); x = x ./ std(x);
Ref_to_signal = double(h5read(fname,strcat(read_name,"/Ref_to_signal"))); %streck
Reference = double(h5read(fname,strcat(read_name,"/Reference")));   %Bokstäverna

% Pick random segment from read
DNA = ['A','C','G','T']; % Neucleotide letters for reference
M = length(Reference); % Number of bases in reference
L = 20; % NUmber of bases to plot
m1 = randi(M-L);
m2 = m1 + L;

% Plot signal and reference
figure(1); clf; hold on;
xb = x(Ref_to_signal(m1)+1:Ref_to_signal(m2+1));
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
%----------------------
% 
% hold on
% t1 = Ref_to_signal(m1);
% vec = [];
% for m = m1+1:m2
%     sampel = Ref_to_signal(m) - t1;
%     vec = [vec, sampel];
%     t1 = Ref_to_signal(m);
% end
% val = unique(vec)
% cnt = histc(vec,val)
% 
% figure(2)
% plot(val,cnt);
% axis([-1 30 -1 10]);

% -------------------------Egen kod-------------------------
 % t1 = Ref_to_signal(m1);
 % vec = [];
 % for m = m1+1:m2
 %     antal_sampel = Ref_to_signal(m) - t1;
 %     vec = [vec, antal_sampel];
 %     t1 = Ref_to_signal(m);      %update start point to next bar
 % end

 vec = [];
 M = length(Reference); % Number of bases in reference
 L = 500; % Number of bases to plot
 s1 = randi(M-L);        %startpunkt
 s2 = s1 + L;            %slutpunkt
 t1 = Ref_to_signal(s1);

 for s = s1+1:s2
     antal_sampel = Ref_to_signal(s) - t1;
     vec = [vec, antal_sampel]; %radvektor
     t1 = Ref_to_signal(s);      %update start point to next bar
 end

 figure(2); clf; grid on; hold on;
 uniques = unique(vec); %returnerar värden utan repetetition
 N = hist(vec, length(uniques));    %liknande graf som joakims
 plot(N)
 
 
 figure(3); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'gamma');
 hold off;
 
 figure(4); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'gev');
 hold off;
 
 figure(5); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'birnbaumsaunders'); %funkar inte ibland pga zeros
 hold off;
 
 figure(6); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'inverse gaussian');
 hold off;
 
 figure(7); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'nbin');
 hold off;
 
 figure(8); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'loglogistic');
 hold off;
 
 figure(9); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'lognormal');
 hold off;
 
 figure(10); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'rician');
 hold off;
 
 figure(11); clf; grid on; hold on;
 N = histfit(vec, length(uniques), 'wbl');
 hold off;
%  distributionFitter
 
%  
%  [T,edges] = histcounts(vec,length(uniques));
%  %plot(uniques,T);
%  
%  T1 = T'
%  pd  = fitdist(T1, 'gamma');
%  y = pdf(pd, length(uniques));
%  figure(4); clf; grid on; hold on
%  plot(length(uniques), y);
 
 
 
%  figure(4); clf; grid on; hold on;
 
 
 %Använda fitdist i samma plot se vilekn som passar bäst.
 %distributionFitter(N)
 
 %prblemet är att veta fördelningen utan att veta fördelningen. och hur rita
 %utan att ha med alla bars. Finns säkert någon funktion i matlab men tar
 %för långt att googla fram.
