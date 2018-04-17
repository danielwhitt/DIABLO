% This file creates a grid in the y-direction for input to DIABLO
% Written by John Taylor, 10/23/2005
% The grid definitions are produced as follows:
%   1 First, define the G grid based on a specified function with
%     loops over j=0,N+1 where the 0 and N+1 correspond to the bottom
%     and top wall values respectively.  Since matlab is unable to store
%     elements in the zero index, all arrays are indexed starting at 1.
%   2 Then, define the fractional grid GF halfway between neighboring
%     G grid locations
%   3 Stretch both GF and G so that GF(0) and GF(N) now correspond
%     to the upper and lower walls.
clear S0 zarray
load gridgyf321.mat
gyf321 = gyf(1:end-1)';
gyf321(1)=0;
load gridgyf161.mat
gyf = gyf(1:end-1)';
gyf(1)=0;
h5fn1 = '/glade/scratch/dwhitt/diablo_runoct17/W1024_2x_WT17/W1024_18_cool_rednoise_161_2x/start.h5'
h5fn2 = '/glade/scratch/dwhitt/diablo_runoct17/NW256_WT17/W1152_2/start.h5'
%nx=1024;
nx=1152;
%nx=128;
%nz=1024;
%nz=128;
nz=1152;
% 1024x1024
%chksize=[1024 1 8];
% 1024x128
%chksize=[1024 1 4];
chksize=[1152 1 8];
%chksize=[128 1 8];
nyLES=321;
nyb = 161;
clear TH1 U V W P
TH1 = repmat(interp1(gyf',nanmean(reshape(permute(h5read(h5fn1,'/Timestep/TH1'),[2 1 3]),[161 1024^2]),2),gyf321','pchip')',[nx 1 nz])+1e-5.*rand(nx,nyLES,nz);
U = repmat(interp1(gyf',nanmean(reshape(permute(h5read(h5fn1,'/Timestep/U'),[2 1 3]),[161 1024^2]),2),gyf321','pchip')',[nx 1 nz])+1e-5.*rand(nx,nyLES,nz);
V = repmat(interp1(gyf',nanmean(reshape(permute(h5read(h5fn1,'/Timestep/V'),[2 1 3]),[161 1024^2]),2),gyf321','pchip')',[nx 1 nz]);
W = repmat(interp1(gyf',nanmean(reshape(permute(h5read(h5fn1,'/Timestep/W'),[2 1 3]),[161 1024^2]),2),gyf321','pchip')',[nx 1 nz])+1e-5.*rand(nx,nyLES,nz);
P = repmat(interp1(gyf',nanmean(reshape(permute(h5read(h5fn1,'/Timestep/P'),[2 1 3]),[161 1024^2]),2),gyf321','pchip')',[nx 1 nz]);
if sum(isnan(U(:)))>0 || sum(isnan(TH1(:)))>0
display('NAN warning')
end
%figure
%plot(gyf,Ninterp(1,:,1),gyf,Pinterp(1,:,1),gyf,Zinterp(1,:,1),gyf,Dinterp(1,:,1));
display('write to hdf?')
pause
% Set the dimensions of the grid 


%h5create(h5fn,Nvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);
%h5create(h5fn,Pvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);
%h5create(h5fn,Zvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);
%h5create(h5fn,Dvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);

% Now, write the grid to file
h5write(h5fn2,'/Timestep/TH1',TH1)
h5write(h5fn2,'/Timestep/U',U)
h5write(h5fn2,'/Timestep/V',V)
h5write(h5fn2,'/Timestep/W',W)
h5write(h5fn2,'/Timestep/P',P)

% % For testing the grid stretching
% for j=2:N-1
%   r(j)=(GF(j+1)-GF(j))/(GF(j)-GF(j-1));
% end
% disp('The maximum grid-stretching ratio is:'),max(r)
% disp(' ');
% disp('grid.h5 has been written to the current directory');
% 



