function [kgrid, medium, mask, c_indice] = poucoHeterogeoC(N,I)
% Cria uma objeto com uma pequena variação em heterogeniedade em C
% Entrada:
% N = Tamanho do grid de saida do objeto
% I = Contraste entre os indices de velocidade do som no objeto

% Saida:
% kgrid = kgrid gerado
% medium = medium gerada
% mask = struct que contem as mascaras de cada uma das 4 regioes (torso, pulmão, coração, coluna e meio)
% c_indice = struct que contem os a velocidad do som em cada uma das regioes


% Define os indices de velocidade
step = 20;
c_meio = 1584;
c_coracao = c_meio - step;
c_osso = c_meio+2*I*step;
c_pulmao = c_meio-I*step;

% Define o indice de densidade
d_meio = 1060;

% Cria o kgrid
Nx = N;           % 24; number of grid points in the x (row) direction
Ny = N;           % number of grid points in the y (column) direction
dx = 0.4/Nx;        % grid point spacing in the x direction [m]
dy = 0.4/Ny;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = c_meio*ones(Nx, Ny);         % [m/s]
medium.density = d_meio*ones(Nx, Ny);             % [kg/m^3]

% load distribution from an image e escala [0,1].
arquivo=sprintf('torso4.png');
torso=loadImage(arquivo); 
torso=resize(torso,[Nx, Ny]);
torso=fix(torso);
mask.torsoMask = boolean(torso);

% Define a regiao do pulmão e altera a velocidade do som nessa regiao
arquivo=sprintf('pulmoes_anat4.png');
pulmoes=loadImage(arquivo); 
pulmoes=resize(pulmoes,[Nx, Ny]); 
pulmoes=fix(pulmoes);
mask.pulmoesMask=boolean(pulmoes);

medium.sound_speed = medium.sound_speed + (c_pulmao - c_meio)*pulmoes ;    % [m/s]

% Define a regiao da coluna e altera a velocidade do som nessa regiao
arquivo=sprintf('coluna_anat4.png');
coluna=loadImage(arquivo); % load distribution from an image e escala [0,1]
coluna=resize(coluna,[Nx, Ny]); % resize the image to match the size of the computational grid and assign
coluna=fix(coluna);
mask.colunaMask=boolean(coluna);

medium.sound_speed = medium.sound_speed + (c_osso - c_meio)*coluna ;    % [m/s]

% Define a regiao do coração e altera a velocidade do som nessa regiao
arquivo=sprintf('coracao_anat4.png');
coracao=loadImage(arquivo); % load distribution from an image e escala [0,1]
coracao=resize(coracao,[Nx, Ny]); % resize the image to match the size of the computational grid and assign
coracao=fix(coracao);
mask.coracaoMask=boolean(coracao);

medium.sound_speed = medium.sound_speed + (c_coracao - c_meio)*coracao ;    % [m/s]

% Seta os indices de velocidades na saida da função
c_indice.c_meio = c_meio;
c_indice.c_coracao = c_coracao;
c_indice.c_osso = c_osso;
c_indice.c_pulmao = c_pulmao;



return;
