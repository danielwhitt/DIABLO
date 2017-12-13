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
load gridgyf161.mat
gyf = gyf(1:end-1)';
gyf(1)=0;
h5fn = '/glade/scratch/dwhitt/diablo_run/NW256_2_cool_rednoise_161_2x_nofrontstart3d_rerunaug52017/start.h5'
%matfn = 'LEStest_comparebioaugust19.mat'
matfn = 'LESini_dz05H80_workspace.mat'
load(matfn)
zarray = fliplr(zarray);
%nx=1024;
nx=256;
%nx=128;
%nz=1024;
%nz=128;
nz=256;
% 1024x1024
%chksize=[1024 1 8];
% 1024x128
%chksize=[1024 1 4];
chksize=[256 1 16];
%chksize=[128 1 8];
nyb = 161;
clear Ninterp Pinterp Zinterp Dinterp
Ninterp = repmat(interp1(zarray,S0(1:nyb),gyf),[nx 1 nz]);
Pinterp = repmat(interp1(zarray,S0((nyb+1):(2*nyb)),gyf),[nx 1 nz]);
Zinterp = repmat(interp1(zarray,S0((2*nyb+1):3*nyb),gyf),[nx 1 nz]);
Dinterp = repmat(interp1(zarray,S0((3*nyb+1):4*nyb),gyf),[nx 1 nz]);
Ninterp(:,1,:)=Ninterp(:,2,:);
Pinterp(:,1,:)=Pinterp(:,2,:);
Zinterp(:,1,:)=Zinterp(:,2,:);
Dinterp(:,1,:)=Dinterp(:,2,:);
Ninterp(:,end,:)=Ninterp(:,end-2,:);
Pinterp(:,end,:)=Pinterp(:,end-2,:);
Zinterp(:,end,:)=Zinterp(:,end-2,:);
Dinterp(:,end,:)=Dinterp(:,end-2,:);
Ninterp(:,end-1,:)=Ninterp(:,end-2,:);
Pinterp(:,end-1,:)=Pinterp(:,end-2,:);
Zinterp(:,end-1,:)=Zinterp(:,end-2,:);
Dinterp(:,end-1,:)=Dinterp(:,end-2,:);
if sum(isnan(Dinterp(:)))>0 || sum(isnan(Ninterp(:)))>0
display('NAN warning')
end
%figure
%plot(gyf,Ninterp(1,:,1),gyf,Pinterp(1,:,1),gyf,Zinterp(1,:,1),gyf,Dinterp(1,:,1));
display('write to hdf?')
pause
Nvarname=['/Timestep/TH2']
Pvarname=['/Timestep/TH3']
Zvarname=['/Timestep/TH4']
Dvarname=['/Timestep/TH5']
% Set the dimensions of the grid 


%h5create(h5fn,Nvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);
%h5create(h5fn,Pvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);
%h5create(h5fn,Zvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);
%h5create(h5fn,Dvarname,[nx 161 nz],'ChunkSize',chksize,'FillValue',0.0);

% Now, write the grid to file
h5write(h5fn,Nvarname,Ninterp)
h5write(h5fn,Pvarname,Pinterp)
h5write(h5fn,Zvarname,Zinterp)
h5write(h5fn,Dvarname,Dinterp)

% % For testing the grid stretching
% for j=2:N-1
%   r(j)=(GF(j+1)-GF(j))/(GF(j)-GF(j-1));
% end
% disp('The maximum grid-stretching ratio is:'),max(r)
% disp(' ');
% disp('grid.h5 has been written to the current directory');
% 



