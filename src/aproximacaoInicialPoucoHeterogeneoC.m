function [medium,mask,c_indice] = aproximacaoInicialPoucoHeterogeneoC(medium, mask, iterationForma,deslocamento,c_indice,cIntensity)
% Cria uma medium aproximad do objeto tomografico.
% Essa aproximação ocorre alterando-se a forma das regioes do objeto, isto
% é, alterando-se o dominio do pulmão, coração, coluna e meio.
% Ou alterando-se o velocidade de propagação do som nessas regiões

% Entrada:
% medium: medium do objeto tomografado que deseja-se criar uma aproximação
% mask: Marcara das regiões (pulmão, coração, coluna e meio) que forma o objeto
% iterationForma: quantas deformações de forma serão executadas
% deslocamento: Amplitude maxima do deslocamento da forma
% c_indice: Velocidade do som das regiões que formam o objeto
% cIntensity: intensidade da variação do c_indice que o medium aproximado
% terá com relação ao medium do objeto tomografado

% Saida:
% medium: Medium aproximado
% mask: Mascará da regioes do objeto aproximado
% c_indice: Velocidade do som das regioes do objeto aproximado

%% Forma

% Deforma e desloca as mascaras da coluna, coração e pulmao
mask.colunaMask =  deforma(mask.colunaMask,0.5,iterationForma,deslocamento);
mask.coracaoMask =  deforma(mask.coracaoMask,0.5,iterationForma,deslocamento);

% separa o pulmao em duas partes e entao deforma cada pulmão separadamene,
% depois junta cada um dos pulmões
p = regionprops(mask.pulmoesMask,'PixelIdxList');
p1 = boolean(zeros(size(mask.pulmoesMask)));
p2 = p1;
p1(p(1).PixelIdxList) = 1;
p2(p(2).PixelIdxList) = 1;
p1 = deforma(p1,0.5,iterationForma,deslocamento);
p2 = deforma(p2,0.5,iterationForma,deslocamento);

mask.pulmoesMask = p1 | p2;

%% Intensidade

% Altera a velocidade de som no meio atráves de uma função gaussiana
% c_meio permanece igual
c_indice.c_coracao = c_indice.c_coracao*(1+0.1*cIntensity*randn); 
c_indice.c_osso = c_indice.c_osso*(1+2*cIntensity*randn);
c_indice.c_pulmao = c_indice.c_pulmao*(1+cIntensity*randn);

%% Criação do medium

% Altera a velocidade do som do  medium de entrada de acordo com as novas
% mask e os novos indices
medium.sound_speed(mask.torsoMask) =  c_indice.c_meio;
medium.sound_speed(mask.colunaMask) =  c_indice.c_osso;
medium.sound_speed(mask.pulmoesMask) =  c_indice.c_pulmao;
medium.sound_speed(mask.coracaoMask) =  c_indice.c_coracao;

end