% A 1D Tissue Model of the Cardiac Purkinje-Ventricular System
% by Mengya Yuan, Heqiang Lian, and Pan Li*, 
% as in Spatiotemporal Patterns of Early Afterdepolarizations Underlying Abnormal T-wave Morphologies in a Tissue Model of the Purkinje-Ventricular System, PLoS ONE, Jan 2023
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
