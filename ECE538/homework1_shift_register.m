clear all
clf
register=[1 0 0 0 0 0 0];
for ri=1:128,
m15(ri)=register(1,7);
register(2:7)=register(1:6);
register(1,1)=rem((register(1,1)+m15(1,ri)),2);
end
m15=2*m15-1;
stem((0:127),m15)