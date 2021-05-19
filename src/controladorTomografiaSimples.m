function [resp,residual,corr,delays,MSE] = controladorTomografiaSimples(ensaio, N, varreduras, cinicial)
% Computa o algoritmo de reconstrução tomografico simples.
% Exite hardcoding do codigo para varios parametros

% Entrada:
% ensaio: Ensaio do objeto a ser tomografado

% N: Tamanho de grid de reconstrução tomografia

% varreduras: Numero de varreduras do algoritmo. Serão executados
% N varreduras * K fontes iterações.

% cinicial: Imagem tomografica inicial.  Caso seja 0 a imagem tomografica
% inicial é definida apartir de uma deformação do imagem real do objeto. Se
% 1, define um meio inicial neutro. Caso seja uma matriz, essa matriz será 
% o mediumInicial.

% imagem tomografica inicial é um meio neutro, isto é, todos os valores são
% setados para 1584

% Saida:
% resp(N,N,iteracoes/5): imagem tomografica obtida a cada 5 iterações do algoritmo

% residual(iteracoes): Soma dos valores residuais calculado entre os dados dos sensores 
% localizados bordo do torso coletados durante a fase de ensaio com os
% dados obtidos durante a simulação da propagação de cada iteração. A soma 
% compreende (M sensores * t step-time.

% Maxima corr(iteracoes): Maxima correlação normalizada entre os dados dos sensores coletados no
% ensaio e na propagação durante a iteração N. Calculamos primeiro a correlação
% dos M sensores, depois fazemos a média dessas M correlações e então
% encontramos e delay que maxima a média.
% Delay: delay que maxima a correlação normalizada, tem que trazer o valor
% para o centro, isto é, (delayReal = delay - lenght(delayReal)/2)

% MSE(iteracoes,8): Calcula o MSE e o MSE normalizado de cada uma das seções do 
% objetos, isto é, Pulmão, coração, ossso, meio. Mairoes detalhes do que
% cada coluna representa ver a função calculaMSE.m

%% Kgrid 

%PArametros do kgrid para realizar as simulacaoes de propagação e backpropagação
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
    % A imagagem inicial pode sofre una deforação em forma ou em intensidade.
    iterationForma = 0;%40 %Deforma as fronterias das regioes do objetos
    deslocamento = 0; %4 %Deforma o centro das regioes do objetos
    cIntensity = 0.0;%0.05  %Deforma o centro das regioes do objetos
    [mediumInicial, ~,~] = aproximacaoInicialPoucoHeterogeneoC(...
        ensaio.medium, ensaio.mask, iterationForma, deslocamento, ensaio.c_indice,cIntensity); 
    
    % Resize estimativa inicial para ser coerente com o tamanho do grid de simulação
    mediumInicial.sound_speed = resize(mediumInicial.sound_speed, [kgridTomo.Nx kgridTomo.Ny]);  % Resize estimativa inicial para ser coerente com o tamanho do grid de simulação
    mediumInicial.density = resize(mediumInicial.density, [kgridTomo.Nx kgridTomo.Ny]);
    c0 = 1584*ones(size(mediumInicial.sound_speed));
    
else
    % Definimos incicialmente o mediumInicial igual ao ensaio.medium para
    % depois altera-lo segundo o cinicial selecionado
    mediumInicial = ensaio.medium;   
    
    % Resize estimativa inicial para ser coerente com o tamanho do grid de simulação
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

%Interpola no espaço. Ajusta o index localização das fontes de um grid
%maior(ensaio) para um grid(menor) tomografia. A interpolação é feita
%através do nearest point da coordenada cartesiano.
%Atençao que a ordem das fontes pode ser alterada, isto é, a fonte que
%antes era a n-ésima fonte em  ensaio.source, pode não ser mais a n-ésima 
%fonte  em kgridTomo . Reorder_index_sourceTomo traz a relação de ordem.
[cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.source.p_mask);
[sourceTomo.p_mask,~,reorder_index_sourceTomo] = cart2grid(kgridTomo, cart_mask);

%Interpola no tempo. Faz uma interpolaçao, alterando a função de excitação
%dos sources que antes estava com tempo inicial e final e time-step de
%ensaio para uma função com 'dominio' em kgridTome. 
F = griddedInterpolant(0:Tf/length(ensaio.func):Tf-Tf/length(ensaio.func),ensaio.func, 'cubic');
func = F(kgridTomo.t_array);
sourceTomo.p =  repmat(func,nSource,1);

%% Sensor

% define os sensores na borda do torso e no torso
% A borda do torso é usado para fazer a comparação entre os dados reais e
% simualados 
% O torso é usado para calcular o incremento da imagem tomografica

% Define os sensores Internos no torse
[~,torso] = takeBoderTorso(kgridTomo); %Está função está sendo usada apenas para pegar os pontos do torso
sensorInterno.mask = torso;
sensorInterno.record = {'p'};
conversorInterno = find(sensorInterno.mask == 1); %Indice linear dos pontos que fazem parte do torso interno
torsoMenor = imerode(torso,strel('diamond', 1)); %matriz que contem os pontos do torso sem a borda

% Define os sensores na borda do torso fazendo uma interpolação espacial
% entre o ensaio e o KgridTomo.
% Apesar da ordem dos sensores serem diferentes, da mesma forma que ocorre
% com as fontes, isso não altera o funcionamento do algortimo, pois os dados
% reais tambem vão sofrer essa mesma reordenação
[cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.sensorTorso.mask);
[sensorTomo.mask, ~, ~] = cart2grid(kgridTomo, cart_mask);
sensorTomo.record = {'p'};


% Codigo que define os sensores internos em todo o dominio interno aos
% sensores e não apenas os pontos dentro do torso.
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

% Ajustar os sensor_data para o cenario da simulação
%Interpola no tempo
for i=1:length(ensaio.sensor_data) %Para cada fonte emissora
    for j=1:size(ensaio.sensor_data(i).p,1) %Para cada sensor
        F = griddedInterpolant(ensaio.kgrid.t_array,ensaio.sensor_data(i).p(j,:), 'linear');
        sensorRealData(i).p(j,:) = F(kgridTomo.t_array)*ensaio.kgrid.Nx/N; % Fator de escala, por que ele é necessario?
                                                                           % Por que a amplitude do sinal simulado em um 
                                                                           % grid de tamanha N1 e em tamanho N2 são diferentes?
    end
end

%Interpola no espaço
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

% Definição dos parametros iniciais
indIteracao = 0;
indResposta = 0;
    
%Roda as iterações
for i=1:varreduras
    P = randperm(nSource); %Ordena aleatoriamente as fontes emissoras.
    for j=1:nSource
        indIteracao = indIteracao+1;
        
        % Define as fontes para simulaçãode propagação e retropropagação.
        % Em cada iteração apenas uma das fontes devem emitir sinais, as
        % outras devem ser nulas.
        % Observe que temos que usar o reorder_index_sourceTomo, pois
        % temos que ter a relação da fonte que queremos que excite o
        % sistema com os dadosReais gerados por elas. 
        % P(j) indica o numero fonte que queremos que esteja ativa de
        % acordo com a numeraçã do kgrid ensaio
        % reorder_index_sourceTomo(P(j)) é a numero da fonte de acordo com 
        % kgridTome
        % sensorRealData(P(j)) são os dados reais obtidos pela ativação da
        % P(j) fonte numerada de acordo com  kgrid ensaio.
        sourceAux.p_mask = sourceTomo.p_mask;
        sourceAux.p = zeros(size(sourceTomo.p));
        sourceAux.p(reorder_index_sourceTomo(P(j)),:) = sourceTomo.p(P(j),:);
        
        %  Iteração do algoritmo
        %  Calcula o a integral do laplaicano da propragação direta vezes a
        % retro-propagação do valor residual.
        % erro(nSensor,nTimeStep,iteração) - Erro entre os dados reais do sensores e os simulados
        % xcf(nSensor,nLag,iteração) - Correlação entre os dados reais do sensores e os simulados com diferentes lags
        [deltaMedium, erro(:,:,indIteracao), xcf(:,:,indIteracao),lag] = tomografyIteration(kgridTomo, mediumInicial, sourceAux, sensorRealData(P(j)), sensorTomo, sensorInterno, conversorInterno, ignorarPercentualInicial);
        
        % Encontra o deslocamento temporal (lag) que maxima a correlação
        soma = sum(xcf(:,:,indIteracao),1);
        [~, maxindice] = max(soma);
        lag = lag(maxindice);       

        % Calcula o deltaMedium de modo que o mesmo tenha media nula quando lag = 0;
        % Ou uma media maior ou menor que zero, porém pequena, quando lag ~= 0
        media = mean( deltaMedium(torsoMenor));
        alteracaoMedia = 1+lag*0.025;
        deltaMedium = deltaMedium - media*alteracaoMedia;
        
        % Atualiza o mediumIncial
        f(torsoMenor) = f(torsoMenor) + 3 * deltaMedium(torsoMenor)./sum(sum(torsoMenor));
        mediumInicial.sound_speed = c0.*sqrt(1+f);
        
        % Imagem tomografica para as iterações multiplas de 5
        if(mod(j,5)==1)
            figure;
            indResposta = indResposta+1;
            resp(:,:,indResposta+1) = mediumInicial.sound_speed;
            imshow(mediumInicial.sound_speed,[]);
        end
        
        % MSE das regioes para todas as iterações
        [MSE(indIteracao+1,:)] = calculaMSE(ensaio, mediumInicial.sound_speed);
        
    end
end



%% Feedback

% Calcula a maxima correlação da média e a somatoria dos valores residuais 
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