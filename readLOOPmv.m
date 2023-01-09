% A 1D Tissue Model of the Cardiac Purkinje-Ventricular System
% by Mengya Yuan, Heqiang Lian, and Pan Li*, 
% as in Yuan M, Lian H, Li P (2023) Spatiotemporal patterns of early afterdepolarizations underlying abnormal T-wave morphologies in a tissue model of the Purkinje-ventricular system. PLoS ONE 18(1): e0280267. https://doi.org/10.1371/journal.pone.0280267
% *Email: pan.li@nih.gov

% To plot use matlab:
%  [t,S,X,T,Si] = readLOOPmv('*.txt',201,1000,16);

function [t,S,X,T,Si] = untitled2(infile, num_ints, S2, run)
fid = fopen(infile, 'rb');
format = [repmat('%f ', 1, num_ints), '%s'];
S = textscan(fid, format, run*(S2+500)/0.4);
S = [S{1:num_ints}];
fclose(fid);
t = S(:,1);
a = 1:1:num_ints;
[X,T]=meshgrid(a,t);
S = S(:,2:num_ints);
X = X(:,2:num_ints);
T = T(:,2:num_ints);
Si=[];
for i = 1:run
  Si = [Si S((((S2+500)/0.4)*(i-1)+1):((S2+500)/0.4)*i,1:(num_ints-1))];
surf(transpose(Si));
shading interp;
colormap(jet);
view(2);
end
