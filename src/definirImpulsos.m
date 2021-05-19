function [kgrid, medium, mask, c_indice] = definirImpulsos(N)
% c diferente, d igual , sem absorção e linear

c_meio = 1584;
c_coracao = 1400;
%     c_osso = c_meio+2*I*step;
%     c_pulmao = c_meio-I*step;

d_meio = 1060;

Nx = N;           % 24; number of grid points in the x (row) direction
Ny = N;           % number of grid points in the y (column) direction
dx = 0.4/Nx;        % grid point spacing in the x direction [m]
dy = 0.4/Ny;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = c_meio*ones(Nx, Ny);         % [m/s]
medium.density = d_meio*ones(Nx, Ny);             % [kg/m^3]

arquivo=sprintf('torso4.png');
torso=loadImage(arquivo); % load distribution from an image e escala [0,1].
torso=resize(torso,[Nx, Ny]);
torso=fix(torso);
mask.torsoMask = boolean(torso);


coracao = zeros(Nx,Ny);
coracao(Nx*0.6:Ny*0.6+3,Ny*0.6:Ny*0.6+3) = ones(4,4);
mask.coracaoMask=boolean(coracao);

medium.sound_speed = medium.sound_speed + (c_coracao - c_meio)*coracao ;    % [m/s]

c_indice.c_meio = c_meio;
c_indice.c_coracao = c_coracao;

return;


