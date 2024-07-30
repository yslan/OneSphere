warning off;clear all; close all; format compact; profile off; diary off; restoredefaultpath;warning on;
pause(.1);hdr;

verbose = 1;
fldr_out = 'outputs_box';
nbx = 10;

mkdir(fldr_out);

% cube sphere at radius = 1
tag = 'cubeSph';
R0 = 1.0;
[X,Quad,Qfront] = gen_fbox_qsph(nbx); X = X * R0; nX = size(X,1);
%ifig = 1; plot_quad(ifig,X,Quad);
draw_quad_vtk(X,Quad,repmat(1:6,nbx*nbx,1),fldr_out,tag,-1);

% boundary layer on sphere
tag = 'sph1';
R1 = 1.2;
nlayers = 2;
ratio = 1.4;
[X,Hexes,Hfront,Hbc,Hcurve] = fbox_extrude_quad2sph(X,Quad,Qfront,nbx,R0,R1,nlayers,ratio);
nX=size(X,1); bc_set = chk_bcid([],Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>1);
draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);

% bounding box of the sph
tag = 'box1';
Lbx = R1*sqrt(2); % face center to origin
nlayers = 2;
ratio = 1.0;
geom1.type=2;
geom1.coef=R1;
geom2.type=1;
geom2.coef=Lbx;
[X,Hexes,Hfront,Hbc,Hcurve] = fbox_extrude_general(X,Hexes,Hfront,Hbc,Hcurve,nbx,geom1,geom2,nlayers,ratio);
nX = size(X,1); bc_set = chk_bcid(bc_set,Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>1);
draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);

% outer box
% bot
tag = 'boxb';
box1.nelz = -5;   box1.zz = [-5.0,-Lbx,1/1.1];
box1.nely = -nbx; box1.yy = [-Lbx,Lbx,1.0];
box1.nelx = -nbx; box1.xx = [-Lbx,Lbx,1.0];
[Xnew,Hnew,Hbc_new,bout1] = genbox(box1,verbose);
Hexes = [Hexes;Hnew+nX]; X = [X; Xnew]; nX=size(X,1); 
Hcurve= cat(3,Hcurve,zeros(6,12,size(Hnew,1)));
Hbc = update_Hbc(Hbc,Hbc_new); bc_set = chk_bcid(bc_set,Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>1);

draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);

% top
tag = 'boxt';
box2.nelz = -7;   box2.zz = [Lbx,10.0,1.3];
box2.nely = -nbx; box2.yy = [-Lbx,Lbx,1.0];
box2.nelx = -nbx; box2.xx = [-Lbx,Lbx,1.0];
[Xnew,Hnew,Hbc_new,bout2] = genbox(box2,verbose);
Hexes = [Hexes;Hnew+nX]; X = [X; Xnew]; nX=size(X,1); 
Hcurve= cat(3,Hcurve,zeros(6,12,size(Hnew,1)));
Hbc = update_Hbc(Hbc,Hbc_new); bc_set = chk_bcid(bc_set,Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>1);
draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);

% surrounding
tag = 'boxs';
box3.nelz = nbx+abs(box1.nelz)+abs(box2.nelz);
box3.zz = uniquetol([linspace(-Lbx,Lbx,nbx+1)';bout1.zz;bout2.zz],1e-6);

eetmp = [-6,-nbx,-6];
xxtmp = [-3,-Lbx,Lbx,3];
ratio = 1.3;
rrtmp = [ratio, 1.0, 1./ratio];
for iy=1:3
   box3.nely = eetmp(iy);
   box3.yy = [xxtmp(iy), xxtmp(iy+1), rrtmp(iy)];
for ix=1:3
   box3.nelx = eetmp(ix);
   box3.xx = [xxtmp(ix), xxtmp(ix+1), rrtmp(ix)];
   if (ix~=2 || iy~=2)
      [Xnew,Hnew,Hbc_new,~] = genbox(box3,verbose);
      Hexes = [Hexes;Hnew+nX]; X = [X; Xnew]; nX=size(X,1);
      Hcurve= cat(3,Hcurve,zeros(6,12,size(Hnew,1)));
      Hbc = update_Hbc(Hbc,Hbc_new);
   end
end
end
bc_set = chk_bcid(bc_set,Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>0);
draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);

% set CBC
bc_tmp{1} = [1,2,3,4,5,6];                   % sphere
bc_tmp{2} = [35,55,65,75,85,95,105,115,125]; % inlet
bc_tmp{3} = [46,56,66,76,86,96,106,116,126]; % oulet
bc_tmp{4} = [54,84,104];                     % side, xmin
bc_tmp{5} = [72,92,122];                     % side, xmax
bc_tmp{6} = [51,61,71];                      % side, ymin
bc_tmp{7} = [103,113,123];                   % side, ymax
CBC = set_cbc(Hbc,bc_tmp);
draw_Hexes_vtk(X,Hexes,CBC,fldr_out,[tag '_cbc'],-2);

% dump re2
fout='out';
fname=[fldr_out '/' fout];
%dump_nek_rea(fname,X,Hexes,CBC,0,[],verbose);
%dump_nek_re2(fname,X,Hexes,CBC,0,[],verbose);
dump_nek_re2(fname,X,Hexes,CBC,1,Hcurve,verbose);
