figure(1);

load('-ascii', "upper.txt");
load('-ascii', "lower.txt");
load('-ascii', "LCDM.txt");
x = LCDM(:, 1);
yupper = upper(:, 2);
ylower = lower(:, 2);
yLCDM = LCDM(:, 2);
yupper = yupper ./ yLCDM;
ylower = ylower ./ yLCDM;
yLCDM = yLCDM ./ yLCDM;
hold on;
plot(x(4000:end), yupper(4000:end));
plot(x(4000:end), ylower(4000:end));
plot(x(4000:end), yLCDM(4000:end));

xlim([0 3])
ylim([0.975 1.025])
xlabel('z')
title('7CPL')
legend({'Line 1', 'Line 2', 'Line 3'}, 'Location', 'northeast')
print("image", "-dpng")
