function [resp,residual,corr,delays,MSE,erroInterpolacao] = controladorTomografiaInterpolacao(ensaio, N, iteration,varargin)

for input_index = 1:2:length(varargin)
    switch varargin{input_index}
        case 'cInicial'
            cInicial =  varargin{input_index+1};      
        % Verifica numero de Source
        case 'nSource'
            if(~exist('nSouce','var'))
                if(mod(ensaio.source.p,varargin{input_index+1}) == 0)
                    nSource = varargin{input_index+1};
                    distSource =  ensaio.source.p/varargin{input_index+1};
                else
                    error('nSource deve ser um divisor do numero de souce do ensaio');
                end
            else
                error('Numero de souce já foi definido: Use nSource ou distSource ');
            end
        case 'distSource'
            if(~exist('nSouce','var'))
                if(mod(nSouce,varargin{input_index+1}) == 0)
                    nSouce = ensaio.source.p/varargin{input_index+1};
                    distSource = varargin{input_index+1};
                else
                    error('distSource deve um divisor do numero de souce do ensaio');
                end
            else
                error('Numero de souce já foi definido: Use nSource ou distSource ');
            end
            
            
            % Verifica numero de sensores
        case 'nSensor'                
            nSensor = varargin{input_index+1};
%         case 'distSensor'
%             if(~exist('sensorIntercalado','var'))
%                 [sensorIntercalado] = defineSensor(ensaio.kgrid, ensaio.sensorTorso.mask, 'distSensor', varargin{input_index+1});
%                 if(varargin{input_index+1} ~= 1)
%                     precisaInterpolar = 1;
%                 else
%                       precisaInterpolar = 0;
%                  end
%             else
%                 error('sensorIntercalado já foi definido: Use nSensor ou distSensor ');
%             end
    end
end

if(~exist('nSource','var'))
    nSource = size(ensaio.source.p,1);
    distSource = 1;
end
if(~exist('nSensor','var'))
%     [sensorIntercalado] = defineSensor(ensaio.kgrid, ensaio.sensorTorso.mask, 'nsensor', sum(sum(ensaio.sensorTorso.mask))); %numero de sensores
     nSensor = 0;
end

    % Kgrid para realizar as simulacaoes de propagação e backpropagação
    Nx = N;           
    Ny = N;   
    dx = 0.4/Nx;   
    dy = 0.4/Ny;
    CFL = 0.3;
    cmim = 500;
    cmax = 4000;
    Tf =  0.4 * 10^-3;
	
    kgridTomo = makeGrid(Nx, dx, Ny, dy);    
    [kgridTomo.t_array, dt] = makeTime(kgridTomo, [cmim cmax], CFL, Tf);

    % Cria uma tomografia inicial
    iterationForma = 0;%40
    deslocamento = 0; %4
    cIntensity = 0.0;%0.05
    [mediumInicial, maskInicial,c_indiceInicial] = aproximacaoInicialPoucoHeterogeneoC(...
            ensaio.medium, ensaio.mask, iterationForma, deslocamento, ensaio.c_indice,cIntensity);
        
    % Resize estimativa inicial para ser coerente com o tamanho do grid de
    % simulação
    mediumInicial.density = resize(mediumInicial.density, [kgridTomo.Nx kgridTomo.Ny]);
    mediumInicial.sound_speed = resize(mediumInicial.sound_speed, [kgridTomo.Nx kgridTomo.Ny]);
    c0 = 1584*ones(size(mediumInicial.sound_speed));
    mediumInicial.sound_speed = c0;
    

    if(exist('cInicial','var'))
        mediumInicial.sound_speed = cInicial;
    end
    
    f =( mediumInicial.sound_speed ./ c0).*( mediumInicial.sound_speed ./ c0)-1;
    [MSE(1,:)] = calculaMSE(ensaio, mediumInicial.sound_speed);
    
    
    
    %Sources
    %Interpola no espaço
    [cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.source.p_mask);
    [sourceTomo.p_mask,~,reorder_index_sourceTomo] = cart2grid(kgridTomo, cart_mask);
    
    %Interpola no tempo
    F = griddedInterpolant(0:Tf/length(ensaio.func):Tf-Tf/length(ensaio.func),ensaio.func, 'cubic');
    func = F(kgridTomo.t_array);
    sourceTomo.p =  repmat(func,nSource,1);
    
    
    % define os sensores na borda e no torso
    %Obtem torsoMenor
    [~,torso] = takeBoderTorso(kgridTomo);
    torsoMenor = imerode(torso,strel('diamond', 1));

     %Obtem sensorTomo de acordo com sensorTorso (ensaio)
    [cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.sensorTorso.mask);
    [sensorTomo.mask, ~, ~] = cart2grid(kgridTomo, cart_mask);
    sensorTomo.record = {'p'};
    
%     sensorInterno.mask = torso;
%     sensorInterno.record = {'p'};
%     conversorInterno = find(sensorInterno.mask == 1);
        sensorInterno.mask = ones(N,N);
        gY = repmat(1:96,96,1);
        gX = repmat((1:96)',1,96);
        sensorInterno.mask (gY<=12/96*N) = 0;
        sensorInterno.mask (gY>=87/96*N) = 0;
        sensorInterno.mask (gX<=22/96*N) = 0;
        sensorInterno.mask (gX>=77/96*N) = 0;    
        sensorInterno.record = {'p'};
        conversorInterno = find(sensorInterno.mask == 1);
    
        
    % Ajustar os sensor_data para o cenario da simulação    
    %Interpola no tempo
    for i=1:length(ensaio.sensor_data) %Para cada fonte emissora
        for j=1:size(ensaio.sensor_data(i).p,1) %Para cada sensor
            F = griddedInterpolant(ensaio.kgrid.t_array,ensaio.sensor_data(i).p(j,:), 'linear');
            sensorRealData(i).p(j,:) = F(kgridTomo.t_array);
        end
    end
    
    %Interpola no espaço
    [cart_mask, ~] = grid2cart(ensaio.kgrid, ensaio.sensorTorso.mask);
    for i=1:length(ensaio.sensor_data) %Para cada fonte emissora        
        sensorRealData(i).p = interpCartData(kgridTomo, sensorRealData(i).p, cart_mask, sensorTomo.mask);
    end
    
    % define quais serão os sensores efetivamentes usados, chamadas de
    % intercalados
    [sensorTomoIntercalado] = defineSensor(kgridTomo, sensorTomo.mask, 'nsensor', nSensor); 
    [cart_sensorTomoIntercalado, ~] = grid2cart(kgridTomo, sensorTomoIntercalado.mask);
    [cart_sensorTomo, ~] = grid2cart(kgridTomo, sensorTomo.mask);
    
    % Relaciona os sensores intecalados com todos os sensores do torso.
    % A relação é mostrada pelo vetor  isIntercalado que tem tamanho igual
    % ao numero de sensores do torso. Cada posição mostra qual é o 
    % numero do sensor intercalado que corresponde com cada um
    % dos sensores do torso, caso não exista relação o numero é zero.
    isIntercalado = zeros(1,size(cart_sensorTomo,2));    
    for i=1:size(cart_sensorTomo,2)
        for j=1:size(cart_sensorTomoIntercalado,2)
            if(cart_sensorTomo(:,i) == cart_sensorTomoIntercalado(:,j))
                isIntercalado(i) = j;
            end
        end
    end
    
    %Identifica os 2 sensores intercalados que estão mais proxima de cada 
    %um dos sensores do torso. Caso o sensor intercalado seja um sensor do
    %torso a identificação é zero.
    % A identificação se da por meio de dois vetores, um que mostra a
    % coordanada e outro que mostra o numero do vizinho. O significado
    % de cada dimensão são:
    % coordenadasVizinhos = (X ou Y, numero do sensor,Vizinho 1 ou 2  EleMesmo, cat ou indice)
    % numeroVizinhos = (numero do vizinho 1 ou 2, numero do sensor)
    nSensorTomo = sum(sum(sensorTomo.mask));
    distanciaSensor = zeros(nSensorTomo, nSensor);    
    coordenadasVizinhos = zeros(2,nSensorTomo,3,2);
    numeroVizinhos = zeros(2,nSensorTomo);   
    for i=1:nSensorTomo
        for j=1:nSensor 
            distanciaSensor(i,j) = (cart_sensorTomo(1,i)-cart_sensorTomoIntercalado(1,j)).*...
                                   (cart_sensorTomo(1,i)-cart_sensorTomoIntercalado(1,j)) +...
                                   (cart_sensorTomo(2,i)-cart_sensorTomoIntercalado(2,j)).*...
                                   (cart_sensorTomo(2,i)-cart_sensorTomoIntercalado(2,j));      
        end
        [valorMim1,indiceMim1] = min(distanciaSensor(i,:));
        if (valorMim1~=0)
            distanciaSensor(i,indiceMim1) = inf;            
            numeroVizinhos(1,i) = find(isIntercalado == indiceMim1);            
            coordenadasVizinhos(:,i,1,1) = cart_sensorTomoIntercalado(:,indiceMim1);
            
            [~,indiceMim2] = min(distanciaSensor(i,:));
            numeroVizinhos(2,i) = find(isIntercalado == indiceMim2);
            coordenadasVizinhos(:,i,2,1) = cart_sensorTomoIntercalado(:, indiceMim2);
            
            coordenadasVizinhos(:,i,3,1) = cart_sensorTomo(:, i);
        end 
    end

    [coordenadaSource(:,:,1), ~] = grid2cart(kgridTomo, sourceTomo.p_mask);    
    [coordenadaSource(:,:,2)] = transformCat2Ind(coordenadaSource(:,:,1),kgridTomo);
    [coordenadasVizinhos(:,:,1,2)] = transformCat2Ind(coordenadasVizinhos(:,:,1,1),kgridTomo);
    [coordenadasVizinhos(:,:,2,2)] = transformCat2Ind(coordenadasVizinhos(:,:,2,1),kgridTomo);
    [coordenadasVizinhos(:,:,3,2)] = transformCat2Ind(coordenadasVizinhos(:,:,3,1),kgridTomo);
    
    resp = zeros([size(mediumInicial.sound_speed) iteration]);
    indExterno = 0;
    indInterno = 0;
    
    %Roda as iterações
    for i=1:iteration
        P = randperm(nSource);
        for j=1:nSource
            indExterno = indExterno+1;
            sourceAux.p_mask = sourceTomo.p_mask;
            sourceAux.p = zeros(size(sourceTomo.p));
            sourceAux.p(reorder_index_sourceTomo(P(j)),:) = sourceTomo.p(P(j),:);
            
            [sensorRealDataInterpolado, erroInterpolacao(:,indExterno)] = interpolateSensorData(sensorRealData(P(j)),mediumInicial, kgridTomo, coordenadasVizinhos, numeroVizinhos, coordenadaSource(:,reorder_index_sourceTomo(P(j)),:));
            
            [deltaMedium, erro(:,:,indExterno), xcf(:,:,indExterno),lags] = tomografyIteration(kgridTomo, mediumInicial, sourceAux, sensorRealDataInterpolado, sensorTomo, sensorInterno, conversorInterno);
            
            soma = sum(xcf(:,:,indExterno),1);
            [~, maxindice] = max(soma);
            lag = lags(maxindice);
            %                 if(lag >= 2)
            %                     mediaDelta = 1.05;
            %                 elseif (lag <= -2)
            %                     mediaDelta = 0.95;
            %                 end
            
            mediaDelta = 1+lag*0.025;
%             mediaDelta = 1;
            
            %                 [dis] = prctile(deltaMedium(torsoMenor),[15 70]);
            %                 dMaior = dis(2);
            %                 dMenor = dis(1);
            %                 dif(indGlobal) = dMaior-dMenor;
            %                 deltaMedium(deltaMedium > dMaior) = dMaior;
            %                 deltaMedium(deltaMedium< dMenor) = dMenor;
            med = mean( deltaMedium(torsoMenor));
            deltaMedium = deltaMedium - med*mediaDelta;
            %                 deltaMedium = deltaMedium - med;
            
            %                 f(torsoMenor) = f(torsoMenor) + 0.5 * deltaMedium(torsoMenor)./c0(torsoMenor);
            f(torsoMenor) = f(torsoMenor) + 3 * deltaMedium(torsoMenor)./c0(torsoMenor);
            
            mediumInicial.sound_speed = c0.*sqrt(1+f);
            % mediumInicial.sound_speed(torsoMenor) = mediumInicial.sound_speed(torsoMenor) + 0.5 * deltaMedium(torsoMenor);
            
            %                 % figure;
            %                 % imshow(deltaMedium,[]);
            if(mod(j,5)==1)
                figure;
                indInterno = indInterno+1;
                %                     mediumInicial.sound_speed = c0.*sqrt(1+f);
                resp(:,:,indInterno) = mediumInicial.sound_speed;
                imshow(mediumInicial.sound_speed,[]);
            end
            [MSE(indExterno+1,:)] = calculaMSE(ensaio, mediumInicial.sound_speed);
            
        end
        %         resp(:,:,i) = mediumInicial.sound_speed;
    end
    
    
    for i = 1:indExterno
        soma = sum(xcf(:,:,i),1);
        [corr(i), delays(i)] = max(soma);
        corr(i) = corr(i)./size(xcf,1);
        residual(i) = sum(sum(erro(:,:,i).*erro(:,:,i)));
    end
    
    

end