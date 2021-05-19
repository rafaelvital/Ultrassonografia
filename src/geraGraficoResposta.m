% Gera graficos de MSE, MSE normalizado, correlação e residual

figure;
subplot(2,1,1)
hold on;
for i=1+4:4+4
    plot(MSE(:,i));
end
hold off;
xlabel('iteration');
ylabel('NMSE');
legend('pulmão','coracão','coluna','meio');

subplot(2,1,2)
hold on;
for i=1:4
     plot(MSE(:,i));
end
hold off;
xlabel('iteration');
ylabel('MSE');
legend('pulmão','coracão','coluna','meio');
saveas(gcf, ['C:\Users\vital\Desktop\resultadoTomo\' figu 'b'], 'fig');
saveas(gcf, ['C:\Users\vital\Desktop\resultadoTomo\' figu 'b'], 'png');

figure;
subplot(2,1,1);
plot(corr);
xlabel('iteration');
ylabel('normalized correlation');
subplot(2,1,2);
plot(residual);
xlabel('iteration');
ylabel('Residual');

saveas(gcf, ['C:\Users\vital\Desktop\resultadoTomo\' figu 'c'], 'fig');
saveas(gcf, ['C:\Users\vital\Desktop\resultadoTomo\' figu 'c'], 'png');

