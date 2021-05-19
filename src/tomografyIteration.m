function [deltaMedium, erro, xcf,lags] = tomografyIteration(kgrid, medium, source, sensorRealData, sensorTorso, sensorInterno, conversorInterno, ignorarPercentualInicial)
% Realiza uma iteação do algoritmo de reconstrução tomografico.
% Uma interação é caracterizada pela integral do laplaicano da propragação
% direta vezes a retro-propagação do valor residual.

% Entradas:
% kgrid: kgrid da imagem tomografica
% medium: medium da imagem tomografica descober até a iteração atual
% source: Fontes da propagação direta (apenas uma das fontes será diferente de zero)
% sensorRealData: Dados dos receptores aquisitados durante a fase de ensaio
% do 
% sensorTorso: sensores localizados ao redor do torso (quadrado envolta do torso ou na borda do mesmo)
% sensorInterno: sensores que estão dentro do torso
% conversorInterno: indice linear dos sensores internos
% ignorarPercentualIncial: percentual incial dos dados aquisitados que será
% ignorado durante o calculo do deltaMedium

% Saidas
%deltaMedium = matriz do tamanha do medium com o calculo da integral do 
% laplaicano da propragação direta vezes a retro-propagação do residual
% calculados ponto a ponto (zero nos pontos que não estão no interior do torso)
% erro(N sensores Torso, M stepTime) = Erro entre os dados reais
% aquisitados nos sensores do Torso e a simulação da propagação direta.
% xcf(N sensores Torso, K lag) = Correlação entre os dados reais
% aquisitados nos sensores do Torso e a simulação da propagação direta para
% K diferentes lags(deslocamento). K varia de 1 a 41
% lags = traduz os valores de lags de 0 a para -20 a +20


%% PROPAGATION

% Junta os sensores no torso, com os sensores interno
sensores.mask = sensorTorso.mask  | sensorInterno.mask;
sensores.record = {'p'};

% Define PML de acordo com o tamanho do grid
if(kgrid.Nx >= 256)
    PMLSize = [20,20];
elseif (kgrid.Nx >= 96)
    PMLSize = [20,10];
else
    PMLSize = [10,5];
end

% Executa a simulação da propagação
% Executado usando precisão unica
sensores_data = kspaceFirstOrder2D(kgrid, medium, source, sensores,'PMLSize',PMLSize,'DataCast', 'single','PlotSim' ,false);
% sensores_data = kspaceFirstOrder2D(kgrid, medium, source, sensores,'PMLSize',PMLSize,'DataCast', 'single','PlotSim' ,true);

% divide as informações dos sensores em duas partes.
% A primeira parte referente aos sensores do torso
% A segunda referente aos sensores interno
sensorTorso_data.p = zeros(sum(sum(sensorTorso.mask)),length(kgrid.t_array));
sensorInterno_data.p = zeros(sum(sum(sensorInterno.mask)),length(kgrid.t_array));
for t=1:length(kgrid.t_array)
    unmasked_sensores_data = unmaskSensorData(kgrid, sensores, sensores_data.p(:,t));
    sensorTorso_data.p(:,t) = unmasked_sensores_data(sensorTorso.mask~=0);
    sensorInterno_data.p(:,t) = unmasked_sensores_data(sensorInterno.mask~=0);
end


%% BACKPROPAGATION

% Copia as informações coletas nos sensores da borda do torso, para os
% novos sensores que também estão na borda do torso.
% Necessário para a  propagação retrograda
sourceBack.p_mask = sensorTorso.mask;
sourceBack.p = sensorTorso_data.p - sensorRealData.p; % seleciona apenas os sensores da bordo do torso

% Calcula o erro e a correlação entre o sinal na borda do torso real e o simulado
erro = sourceBack.p;
 for i=1:size(sensorTorso_data.p,1)
[xcf(i,:),lags] = crosscorr(sensorRealData.p(i,:),sensorTorso_data.p(i,:));
end

% inverte a ordem temporal dos sensores
sourceBack.p = flip(sourceBack.p, 2);
sourceBack.p_mode = 'dirichlet';

% simulação backpropagation
sensorInternoBack_data = kspaceFirstOrder2D(kgrid, medium, sourceBack, sensorInterno,'PMLSize', PMLSize, 'DataCast', 'single', 'PlotSim',false);
% sensorInternoBack_data = kspaceFirstOrder2D(kgrid, medium, sourceBack, sensores,'PMLSize', PMLSize, 'DataCast', 'single', 'PlotSim',true);

% inverte a ordem temporal dos sensores Interno
sensorInternoBack_data.p  = flip(sensorInternoBack_data.p, 2);


%% DELTA CALCULATION

%Calcula laplaciano de uk * zk
Mult = zeros(length(conversorInterno), size(sensorInterno_data.p,2));
pInterno = zeros(kgrid.Nx, kgrid.Ny);
for t=floor(1+size(sensorInterno_data.p,2)*ignorarPercentualInicial):size(sensorInterno_data.p,2) %para cada instante de tempo, depois de ignorar uma porcentagem incial dos dados
    
    % Cria uk em forma matricial
    pInterno(conversorInterno) = sensorInterno_data.p(:,t);
    
    %Calcula laplaciano de uk
    Lapla = 4 * del2(pInterno, kgrid.dx ,kgrid.dy); %Fator 4 nessário, pois del2 calcula 0.25*laplaciano(u)
  
    %Calcula laplaciano de uk * zk 
    Mult(:,t) = Lapla(conversorInterno) .* sensorInternoBack_data.p(:,t);
end


%Calcula integral de laplaciano de uk * zk
deltaMedium = zeros(size(medium.sound_speed));
deltaMedium(conversorInterno) = trapz(Mult,2) * kgrid.dt;

    
end