function vec = NumberOfSamplesPerBase(Ref_to_signal)

vec = [];
t1 = Ref_to_signal(1);

for s = 2:length(Ref_to_signal)
    antal_sampel = Ref_to_signal(s) - t1;
    vec = [vec, antal_sampel];  
    t1 = Ref_to_signal(s);      %update start point to "next" bar
end

end