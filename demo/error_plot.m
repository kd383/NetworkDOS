% error = zeros(20,1);
% for k = 1:5
%     for i = 50:50:1000
%         error(i/50) = error(i/50) + demof('',[],sign(randn(n,20)),i)/n;
%     end
% end
% error = error/5;

figure('outerposition',[0 0 900 900]);
hold on
plot(50:50:1000,error_nofilt,'r--','LineWidth',5)
plot(50:50:1000,error_onefilt,'g','LineWidth',5)
plot(50:50:1000,error_filt,'b','LineWidth',5)
plot(50:50:1000,error_nofilt,'rd','MarkerSize',20,'MarkerFaceColor','r')
plot(50:50:1000,error_onefilt,'g*','MarkerSize',20,'LineWidth',20,'MarkerFaceColor','m')
plot(50:50:1000,error_filt,'bo','MarkerSize',20,'MarkerFaceColor','b')
xlabel('#moments')
ylabel('rel. error')
xlim([50 1000])
grid on
set(gca,'xtick',200:200:1000,'FontSize',40,'FontWeight','bold');
legend('No Filter','Zero Filter', 'All Filter')
box on
