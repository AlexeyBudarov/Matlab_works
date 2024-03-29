axis equal;
grid on;
hold on;
axis([-5 5 -5 5]);

% Параметризованные функции
x1 = @(t) sqrt(2 * abs(cos(2 * t))) .* cos(t);
y1 = @(t) sqrt(2 * abs(cos(2 * t))) .* sin(t);

x2 = @(t) -sqrt(33/4) * cosh(t) - sqrt(33/8) * sinh(t) - 0.5;
y2 = @(t) sqrt(33/8) * sinh(t) + 1;

x3 = @(t) sqrt(33/4) * cosh(t) - sqrt(33/8) * sinh(t) - 0.5;
y3 = @(t) sqrt(33/8) * sinh(t) + 1;

lim_f1 = [0,2*pi];
lim_f2 = [-2,2];

A = [-1.6, 2];


Distance = @(v) sqrt((v(1,1) - v(2,1)).^2 + (v(1,2) - v(2,2)).^2);

% Поиск начального приближения д
dt = 0.1;
min_dist_3 = Inf;
min_dist_2 = Inf;
min_dist_1 = Inf;

for t_1 = lim_f1(1):dt:lim_f1(2)
    dist_1 = Distance([[x1(t_1),y1(t_1)]; A]);
    dist_2 = Distance([[x2(t_1),y2(t_1)]; A]);

    if dist_1 < min_dist_1
        min_dist_1 = dist_1;
        param1 = t_1;
    end
    if dist_2 < min_dist_2
        min_dist_2 = dist_2;
        param2 = t_1;
    end

    for t_2 = lim_f2(1):dt:lim_f2(2)
        dist_3 = Distance([[x1(t_1),y1(t_1)]; [x2(t_2),y2(t_2)]]);
        if dist_3 < min_dist_3
            min_dist_3 = dist_3;
            param3 = [t_1, t_2];
        end
    end
end

disp(['Approximate distance_1: ', num2str(min_dist_1)])
disp(['Approximate distance_2: ', num2str(min_dist_2)])
disp(['Approximate distance_3: ', num2str(min_dist_3)])

t = lim_f1(1):0.01:lim_f1(2);
plot(x1(t), y1(t), 'Color', 'b', 'LineWidth', 1.5)
t_temp = lim_f2(1):0.01:lim_f2(2);
plot(x2(t_temp), y2(t_temp), 'Color', 'g', 'LineWidth', 1.5)
plot(x3(t_temp), y3(t_temp), 'Color', 'g', 'LineWidth', 1.5)
plot(A(1), A(2), '*', 'Color','r', 'LineWidth', 2);


optFunction_1 = @(t) Distance([A;[x1(t), y1(t)]]);
[xval1, distance_1] = fminsearch(optFunction_1, param1);
plot(x1(xval1), y1(xval1),'*', 'color', 'r', 'LineWidth', 2);
plot([A(1), x1(xval1)], [A(2), y1(xval1)], 'LineWidth', 2, 'Color', 'k');
text((A(1) + x1(xval1))/2, (A(2) + y1(xval1))/2 + 0.1, sprintf('%.5f', distance_1), 'FontSize', 10);


optFunction_2 = @(t) Distance([A;[x2(t), y2(t)]]);
[xval2, distance_2] = fminsearch(optFunction_2, param2);
plot(x2(xval2), y2(xval2), '*', 'LineWidth', 2, 'Color', 'r');
plot([A(1), x2(xval2)], [A(2), y2(xval2)], 'LineWidth', 2, 'Color', 'k');
text((A(1) + x2(xval2))/2 - 0.5, (A(2) + y2(xval2))/2 + 0.2, sprintf('%.5f', distance_2), 'FontSize', 10);



optFunction_3 = @(t) Distance([[x1(t(1)),y1(t(1))];[x2(t(2)), y2(t(2))]]);
[xval3, distance_3]= fminsearch(optFunction_3, param3);
plot(x1(xval3(1)), y1(xval3(1)), '*', "LineWidth",2, 'Color', 'r');
plot(x2(xval3(2)), y2(xval3(2)), '*', "LineWidth",2, 'Color', 'r');
plot([x1(xval3(1)), x2(xval3(2))], [y1(xval3(1)), y2(xval3(2))], "LineWidth",2, 'Color', 'k');
text((x1(xval3(1)) + x2(xval3(2))) / 2 -0.08, (y1(xval3(1)) + y2(xval3(2))) / 2 -0.08, sprintf('%.4f', distance_3), 'FontSize', 10);



Line1 = [A(1) - x1(xval1), A(2) - y1(xval1)];
Line2 = [A(1) - x2(xval2), A(2) - y2(xval2)];
phi = rad2deg(acos((Line1(1) * Line2(1) + Line1(2) * Line2(2))/(distance_1 * distance_2)));
text(A(1) - 0.2, A(2) - 0.25, sprintf('%.1f°', phi), 'FontSize', 10);




