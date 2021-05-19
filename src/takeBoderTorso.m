function [bordaTorso,torso] = takeBoderTorso(kgrid)
% Dado o kgrid, retorna uma matriz de marcara com os pontos 
% que pertencem ao torso e a borda externa do torso
% Com hardcoding pode ser a borda interna.

arquivo=sprintf('torso4.png');
torso=loadImage(arquivo); % load distribution from an image e escala [0,1].
torso=resize(torso,[kgrid.Nx, kgrid.Ny]);
torso=fix(torso);
torso = im2bw(torso,0.5);

%  borda externa do torso
bordaTorso = imdilate(torso,strel('diamond', 1));
bordaTorso = bordaTorso-torso;

% % borda interna do torso
% erode = imerode(torso,strel('diamond', 1));
% bordaTorso = torso - erode;
% torso = erode;

    
end