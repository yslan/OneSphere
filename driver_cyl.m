warning off;clear all; close all; format compact; profile off; diary off; restoredefaultpath;warning on;
pause(.1);hdr;

verbose = 1;
fldr_out = 'outputs_cyl';
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

% cylinder layer
tag = 'cyl1';
Rcyl = 4.4;
nlayers = 7;
ratio = 1/1.3;
geom1 = geom2;
geom2.type=3;
geom2.coef=Rcyl;
ratio = [0.25,0.45,0.6,0.75,0.87,0.95,1.0];
[X,Hexes,Hfront,Hbc,Hcurve] = fbox_extrude_general(X,Hexes,Hfront,Hbc,Hcurve,nbx,geom1,geom2,nlayers,ratio);
nX = size(X,1); bc_set = chk_bcid(bc_set,Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>1);
draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);

% bot / top cylinder via extrusion
func_cfill = @(X,c) c*ones(size(X,1),1);
tag = 'cbot';
bid_set=[25,35];
zmin = -5;
nlayers = 5;
ratio = 1/1.1;

func_geom_pj = @(X) [X(:,1:2),func_cfill(X,zmin)];
[X,Hexes,Hbc] = extrude_onto_geom(X,Hexes,Hbc,bid_set,func_geom_pj,nlayers,ratio);
bc_set = chk_bcid(bc_set,Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>1);
draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);


tag = 'ctop';
bid_set=[26,36];
zmax = 10;
nlayers = 7;
ratio = 1.3;
func_geom_pj = @(X) [X(:,1:2),func_cfill(X,zmax)];
[X,Hexes,Hbc] = extrude_onto_geom(X,Hexes,Hbc,bid_set,func_geom_pj,nlayers,ratio);
bc_set = chk_bcid(bc_set,Hbc,tag,1);
chk_hex(X,Hexes,tag,verbose); chk_hex_metric(X,Hexes,tag,verbose>1);
draw_Hexes_vtk(X,Hexes,0,fldr_out,tag,-4);
draw_Hexes_vtk(X,Hexes,Hbc,fldr_out,[tag '_f'],-2);


% set CBC
bc_tmp{1} = [1,2,3,4,5,6]; % sphere
bc_tmp{2} = [46];          % inlet
bc_tmp{3} = [56];          % oulet
bc_tmp{4} = [31,32,33,34]; % side, xmin
CBC = set_cbc(Hbc,bc_tmp);
Hcurve = fix_cyl_curve_via_BC(X,Hexes,Hbc,Hcurve,bc_tmp{4});
draw_Hexes_vtk(X,Hexes,CBC,fldr_out,[tag '_cbc'],-2);

% dump re2
fout='out';
fname=[fldr_out '/' fout];
%dump_nek_rea(fname,X,Hexes,CBC,0,[],verbose);
dump_nek_re2(fname,X,Hexes,CBC,1,Hcurve,verbose);
