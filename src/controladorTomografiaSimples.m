function [resp,residual,corr,delays,MSE] = controladorTomografiaSimples(ensaio, N, varreduras, cinicial)
% Computa o algoritmo de reconstru��o tomografico simples.
% Exite hardcoding do codigo para varios parametros

% Entrada:
% ensaio: Ensaio do objeto a ser tomografado

% N: Tamanho de grid de reconstru��o tomografia

% varreduras: Numero de varreduras do algoritmo. Ser�o executados
% N varreduras * K fontes itera��es.

% cinicial: Imagem tomografica inicial.  Caso seja 0 a imagem tomografica
% inicial � definida apartir de uma deforma��o do imagem real do objeto. Se
% 1, define um meio inicial neutro. Caso seja uma matriz, essa matriz ser� 
% o mediumInicial.

% imagem tomografica inicial � um meio neutro, isto �, todos os valores s�o
% setados para 1584

% Saida:
% resp(N,N,iteracoes/5): imagem tomografica obtida a cada 5 itera��es do algoritmo

% residual(iteracoes): Soma dos valores residuais calculado entre os dados dos sensores 
% localizados bordo do torso coletados durante a fase de ensaio com os
% dados obtidos durante a simula��o da propaga��o de cada itera��o. A soma 
% compreende (M sensores * t step-time.

% Maxima corr(iteracoes): Maxima correla��o normalizada entre os dados dos sensores coletados no
% ensaio e na propaga��o durante a itera��o N. Calculamos primeiro a correla��o
% dos M sensores, depois fazemos a m�dia dessas M correla��es e ent�o
% encontramos e delay que maxima a m�dia.
% Delay: delay que maxima a correla��o normalizada, tem que trazer o valor
% para o centro, isto �, (delayReal = delay - lenght(delayReal)/2)

% MSE(iteracoes,8): Calcula o MSE e o MSE normalizado de cada uma das se��es do 
% objetos, isto �, Pulm�o, cora��o, ossso, meio. Mairoes detalhes do que
% cada coluna representa ver a fun��o calculaMSE.m

%% Kgrid 

%PArametros do kgrid para realizar as simulacaoes de propaga��o e backpropaga��o
Nx = N;
Ny = N;
dx = 0.4/Nx;
dy = 0.4/Ny;
CFL = 0.3;
cmim = 500;
cmax = 4000;
Tf =  0.4 * 10^-3;
nSource = size(ensaio.source.p,1);

kgridTomo = makeGrid(Nx, dx, Ny, dy);
[kgridTomo.t_array, dt] = makeTime(kgridTomo, [cmim cmax], CFL, Tf);
ignorarPercentualInicial = 0.01;

%% mediumInicial

% Analisa o caso da mediumInicial

if(cinicial == 0)
    % Cria uma tomografia inicial baseado na imagem do objeto a ser tomografado
    % A imagagem inicial pode sofre una defora��o em forma ou em intensidade.
    iterationForma = 0;%40 %Deforma as fronterias das regioes do objetos
    deslocamento = 0; %4 %Deforma o centro das regioes do objetos
    cIntensity = 0.0;%0.05  %Deforma o centro das regioes do objetos
    [mediumInicial, ~,~] = aproximacaoInicialPoucoHeterogeneoC(...
        ensaio.medium, ensaio.mask, iterationForma, deslocamento, ensaio.c_indice,cIntensity); 
    
    % Resize estimativa inicial para ser coerente com o tamanho do grid de simula��o
    mediumInicial.sound_speed = resize(mediumInicial.sound_speed, [kgridTomo.Nx kgridTomo.Ny]);  % Resize estimativa inicial para ser coerente com o tamanho do grid de simula��o
    mediumInicial.density = resize(mediumInicial.density, [kgridTomo.Nx kgridTomo.Ny]);
    c0 = 1584*ones(size(mediumInicial.sound_speed));
    
else
    % Definimos incicialmente o mediumInicial igual ao ensaio.medium para
    % depois altera-lo segundo o cinicial selecionado
    mediumInicial = ensaio.medium;   
    
    % Resize estimativa inicial para ser coerente com o tamanho do grid de simula��o
    mediumInicial.sound_speed = resize(mediumInicial.sound_speed, [kgridTomo.Nx kgridTomo.Ny]); 
    mediumInicial.density = resize(mediumInicial.density, [kgridTomo.Nx kgridTomo.Ny]);
    c0 = 1584*ones(size(mediumInicial.sound_speed));
    
    if(cinicial == 1)
        mediumInicial = c0;
    else
        mediumInicial.sound_speed = cinicial;
    end
end

% Define a variavel f do algortimo baseado no mediun incial
f =( mediumInicial.sound_speed ./ c0).*( mediumInicial.sound_speed ./ c0)-1;
    
%% Sources

%Interpola no espa�o. Ajusta o index localiza��o das fontes de um grid
%maior(ensaio) para um grid(menor) tomografia. A interpola��o � feita
%atrav�s do nearest point da coordenada cartesiano.
%Aten�ao que a ordem das fontes pode ser alterada, isto �, a fonte que
%antes era a n-�sima fonte em  ensaio.source, pode n�o ser mais a n-�sima 
%fonte  em kgridTomo . Reorder_index_sourceTomo traz a rela��o de ordem.
[cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.source.p_mask);
[sourceTomo.p_mask,~,reorder_index_sourceTomo] = cart2grid(kgridTomo, cart_mask);

%Interpola no tempo. Faz uma interpola�ao, alterando a fun��o de excita��o
%dos sources que antes estava com tempo inicial e final e time-step de
%ensaio para uma fun��o com 'dominio' em kgridTome. 
F = griddedInterpolant(0:Tf/length(ensaio.func):Tf-Tf/length(ensaio.func),ensaio.func, 'cubic');
func = F(kgridTomo.t_array);
sourceTomo.p =  repmat(func,nSource,1);

%% Sensor

% define os sensores na borda do torso e no torso
% A borda do torso � usado para fazer a compara��o entre os dados reais e
% simualados 
% O torso � usado para calcular o incremento da imagem tomografica

% Define os sensores Internos no torse
[~,torso] = takeBoderTorso(kgridTomo); %Est� fun��o est� sendo usada apenas para pegar os pontos do torso
sensorInterno.mask = torso;
sensorInterno.record = {'p'};
conversorInterno = find(sensorInterno.mask == 1); %Indice linear dos pontos que fazem parte do torso interno
torsoMenor = imerode(torso,strel('diamond', 1)); %matriz que contem os pontos do torso sem a borda

% Define os sensores na borda do torso fazendo uma interpola��o espacial
% entre o ensaio e o KgridTomo.
% Apesar da ordem dos sensores serem diferentes, da mesma forma que ocorre
% com as fontes, isso n�o altera o funcionamento do algortimo, pois os dados
% reais tambem v�o sofrer essa mesma reordena��o
[cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.sensorTorso.mask);
[sensorTomo.mask, ~, ~] = cart2grid(kgridTomo, cart_mask);
sensorTomo.record = {'p'};


% Codigo que define os sensores internos em todo o dominio interno aos
% sensores e n�o apenas os pontos dentro do torso.
% sensorInterno.mask = ones(N,N);
% gY = repmat(1:96,96,1);
% gX = repmat((1:96)',1,96);
% sensorInterno.mask (gY<=12/96*N) = 0;
% sensorInterno.mask (gY>=87/96*N) = 0;
% sensorInterno.mask (gX<=22/96*N) = 0;
% sensorInterno.mask (gX>=77/96*N) = 0;
% sensorInterno.record = {'p'};
% conversorInterno = find(sensorInterno.mask == 1);


%% SensorData

% Ajustar os sensor_data para o cenario da simula��o
%Interpola no tempo
for i=1:length(ensaio.sensor_data) %Para cada fonte emissora
    for j=1:size(ensaio.sensor_data(i).p,1) %Para cada sensor
        F = griddedInterpolant(ensaio.kgrid.t_array,ensaio.sensor_data(i).p(j,:), 'linear');
        sensorRealData(i).p(j,:) = F(kgridTomo.t_array)*ensaio.kgrid.Nx/N; % Fator de escala, por que ele � necessario?
                                                                           % Por que a amplitude do sinal simulado em um 
                                                                           % grid de tamanha N1 e em tamanho N2 s�o diferentes?
    end
end

%Interpola no espa�o
[cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.sensorTorso.mask);
for i=1:length(ensaio.sensor_data) %Para cada fonte emissora
    sensorRealData(i).p = interpCartData(kgridTomo, sensorRealData(i).p, cart_mask, sensorTomo.mask);
end

%% Algoritmo

% calcula o mse entre o ensaio e o meio incial
[MSE(1,:)] = calculaMSE(ensaio, mediumInicial.sound_speed);

% Prealoca a resposta da imagem tomografica
resp = zeros([size(mediumInicial.sound_speed) varreduras*nSource/5+1]);
resp(:,:,1) = mediumInicial.sound_speed;
imshow(mediumInicial.sound_speed,[]);

% Defini��o dos parametros iniciais
indIteracao = 0;
indResposta = 0;
    
%Roda as itera��es
for i=1:varreduras
    P = randperm(nSource); %Ordena aleatoriamente as fontes emissoras.
    for j=1:nSource
        indIteracao = indIteracao+1;
        
        % Define as fontes para simula��ode propaga��o e retropropaga��o.
        % Em cada itera��o apenas uma das fontes devem emitir sinais, as
        % outras devem ser nulas.
        % Observe que temos que usar o reorder_index_sourceTomo, pois
        % temos que ter a rela��o da fonte que queremos que excite o
        % sistema com os dadosReais gerados por elas. 
        % P(j) indica o numero fonte que queremos que esteja ativa de
        % acordo com a numera�� do kgrid ensaio
        % reorder_index_sourceTomo(P(j)) � a numero da fonte de acordo com 
        % kgridTome
        % sensorRealData(P(j)) s�o os dados reais obtidos pela ativa��o da
        % P(j) fonte numerada de acordo com  kgrid ensaio.
        sourceAux.p_mask = sourceTomo.p_mask;
        sourceAux.p = zeros(size(sourceTomo.p));
        sourceAux.p(reorder_index_sourceTomo(P(j)),:) = sourceTomo.p(P(j),:);
        
        %  Itera��o do algoritmo
        %  Calcula o a integral do laplaicano da propraga��o direta vezes a
        % retro-propaga��o do valor residual.
        % erro(nSensor,nTimeStep,itera��o) - Erro entre os dados reais do sensores e os simulados
        % xcf(nSensor,nLag,itera��o) - Correla��o entre os dados reais do sensores e os simulados com diferentes lags
        [deltaMedium, erro(:,:,indIteracao), xcf(:,:,indIteracao),lag] = tomografyIteration(kgridTomo, mediumInicial, sourceAux, sensorRealData(P(j)), sensorTomo, sensorInterno, conversorInterno, ignorarPercentualInicial);
        
        % Encontra o deslocamento temporal (lag) que maxima a correla��o
        soma = sum(xcf(:,:,indIteracao),1);
        [~, maxindice] = max(soma);
        lag = lag(maxindice);       

        % Calcula o deltaMedium de modo que o mesmo tenha media nula quando lag = 0;
        % Ou uma media maior ou menor que zero, por�m pequena, quando lag ~= 0
        media = mean( deltaMedium(torsoMenor));
        alteracaoMedia = 1+lag*0.025;
        deltaMedium = deltaMedium - media*alteracaoMedia;
        
        % Atualiza o mediumIncial
        f(torsoMenor) = f(torsoMenor) + 3 * deltaMedium(torsoMenor)./sum(sum(torsoMenor));
        mediumInicial.sound_speed = c0.*sqrt(1+f);
        
        % Imagem tomografica para as itera��es multiplas de 5
        if(mod(j,5)==1)
            figure;
            indResposta = indResposta+1;
            resp(:,:,indResposta+1) = mediumInicial.sound_speed;
            imshow(mediumInicial.sound_speed,[]);
        end
        
        % MSE das regioes para todas as itera��es
        [MSE(indIteracao+1,:)] = calculaMSE(ensaio, mediumInicial.sound_speed);
        
    end
end



%% Feedback

% Calcula a maxima correla��o da m�dia e a somatoria dos valores residuais 
% E identifica o delauy maximo.
corr = zeros(indIteracao,1);
delays = zeros(indIteracao,1);
residual = zeros(indIteracao,1);
for i = 1:indIteracao
    soma = sum(xcf(:,:,i),1);
    [corr(i), delays(i)] = max(soma);
    corr(i) = corr(i)./size(xcf,1);
    residual(i) = sum(sum(erro(:,:,i).*erro(:,:,i)));
end

        

end