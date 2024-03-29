x_1 = @(tetta, alpha) 3.*sin(tetta).*cos(alpha);
y_1 = @(tetta, alpha) 3.*sin(tetta).*sin(alpha);
z_1 = @(tetta) 3.*cos(tetta);


x_2 = @(tetta, alpha) sqrt(3).*sin(tetta).*cos(alpha);
y_2 = @(tetta, alpha) sqrt((30./(-sqrt(29) + 11))).*sin(tetta).*sin(alpha);
z_2 = @(tetta) sqrt((30./(sqrt(29) + 11))).*cos(tetta);

distanceFunc = @(t) sqrt((x_1(t(1), t(2)) - x_2(t(3), t(4))).^2 + (y_1(t(1), t(2)) - y_2(t(3), t(4))).^2 + (z_1(t(1)) - z_2(t(3))).^2);

xlmt = [-1.5,1.5];
ylmt = [-2,2];


d_tetta = -pi/2:0.1:pi/2;
d_alpha = 0:0.1:2*pi;

min_distance = Inf;
for tetta= d_tetta
    for alpha = d_alpha
        x1 = x_1(tetta,alpha);
        y1 = y_1(tetta,alpha);
        z1 = z_1(tetta);
        for tetta2 = d_tetta
            for alpha2 = d_alpha
                x2 = x_2(tetta2, alpha2);
                y2 = y_2(tetta2, alpha2);
                z2 = z_2(tetta2);
                if (x2>=xlmt(1))&&(x2<=xlmt(2))&&(y2>=ylmt(1))&&(y2<=ylmt(2))&&(x1>=xlmt(1))&&(x1<=xlmt(2))&&(y1>=ylmt(1))&&(y1<=ylmt(2))
                    distance = distanceFunc([tetta, alpha, tetta2, alpha2]);
                    if distance < min_distance
                        min_distance = distance;
                        min_param = [tetta, alpha, tetta2, alpha2];
                    end
                end
            end
        end
    end
end

[params_, fval] = fmincon(distanceFunc, min_param, [], [], [], [], [], [], @nonlcon);
[paramss_, f_val] = fminsearch(distanceFunc, min_param);
disp(f_val)


A_temp = [x_1(paramss_(1), paramss_(2)),y_1(paramss_(1), paramss_(2)),z_1(paramss_(1))];
B_temp = [x_2(paramss_(3), paramss_(4)),y_2(paramss_(3), paramss_(4)),z_2(paramss_(3))];

A = [x_1(params_(1), params_(2)),y_1(params_(1), params_(2)),z_1(params_(1))];
B = [x_2(params_(3), params_(4)),y_2(params_(3), params_(4)),z_2(params_(3))];

disp(['First point A : ', num2str(A)]);
disp(['Second point B : ', num2str(B)]);
disp(['Min distance: ', num2str(fval)]);

f_1 = @(x,y,z) z.^2 + y.^2 + x.^2 - 9;               
% f_2 = @(x,y,z) z.^2 + (1.3).*x.^2 + (0.9).*y.^2 - x.*y - 3;
f_2 = @(x,y,z) z.^2./(30./(sqrt(29) + 11)) + y.^2./(30./(-sqrt(29) + 11)) + x.^2./3 - 1;
% interval = [-1.5 1.5 -2 2, -7 7];
interval = [-3 1.5 -2 2, -7 7];

% % Построение нормалей
syms x y z

grad_1 = [diff(f_1(x,y,z),x),diff(f_1(x,y,z),y),diff(f_1(x,y,z),z)];
grad_2 = [diff(f_2(x,y,z),x),diff(f_2(x,y,z),y),diff(f_2(x,y,z),z)];
grad_fun_1 = matlabFunction(grad_1);
grad_fun_2 = matlabFunction(grad_2);


n_A = -grad_fun_1(A(1), A(2), A(3))/norm(grad_fun_1(A(1), A(2), A(3)));
n_B = -grad_fun_2(B(1), B(2), B(3))/norm(grad_fun_2(B(1), B(2), B(3)));


A_normal = [A; A + n_A];
B_normal = [B; B + n_B];

 
% Вектор AB и его единичный вектор
AB = A - B;
AB_unit = AB/norm(AB);

% Вычисление косинуса угла между вектором AB и нормалью в точке A и в точке B
cos_angle_A = dot(AB_unit, n_A);
cos_angle_B = dot(AB_unit, n_B);

% Вычисление угла между AB и нормалью в точке A и в точке B
angle_A = acos(cos_angle_A) * 180 / pi;
angle_B = 180 - (acos(cos_angle_B) * 180 / pi);

disp(['The angle between AB and n_A: ', num2str(angle_A)]);
disp(['The angle between AB and n_B: ', num2str(angle_B)]);


z_temp = [A(3), B(3)];
x_temp = [A(1), B(1)];
y_temp = [A(2), B(2)];

%Построение 
fimplicit3(f_1, interval);
hold on;
fimplicit3(f_2, interval);

plot3(A(1), A(2), A(3), '*', 'LineWidth',5);
plot3(B(1), B(2), B(3), '*', 'LineWidth', 5);

plot3(A_temp(1), A_temp(2), A_temp(3), '*', 'LineWidth',5);
plot3(B_temp(1), B_temp(2), B_temp(3), '*', 'LineWidth', 5);

plot3(x_temp,y_temp , z_temp, 'r', 'LineWidth',3)
plot3(A_normal(:,1), A_normal(:,2), A_normal(:,3), 'g', 'LineWidth', 3)
plot3(B_normal(:,1), B_normal(:,2), B_normal(:,3), 'g', 'LineWidth', 3)

patch([A(1) B(1) A_normal(2,1)],[A(2) B(2) A_normal(2,2)],[A(3) B(3) A_normal(2,3)],'r','FaceAlpha',0.5)
patch([A(1) B(1) B_normal(2,1)],[A(2) B(2) B_normal(2,2)],[A(3) B(3) B_normal(2,3)],'y','FaceAlpha',0.5)

text(A(1),A(2),A(3),'A','FontSize', 20, 'FontWeight','bold')
text(B(1),B(2),B(3),'B','FontSize', 20, 'FontWeight','bold')
legend(['angle_A: ', num2str(angle_A)], ['angle_B:', num2str(angle_B)])
axis equal


function [c, ceq] = nonlcon(angle)
    ceq = [];
    c(1) = -1.5 - 3*sin(angle(1))*cos(angle(2));
    c(2) = 3*cos(angle(1))*cos(angle(2))-1.5;
    c(3) = -1.5 - sqrt(3).*sin(angle(3)).*cos(angle(4));
    c(4) = sqrt(3).*sin(angle(3)).*cos(angle(4)) - 1.5;
    c(5) = -2 - sqrt((30./(-sqrt(29) + 11))).*sin(angle(3)).*sin(angle(4));
    c(6) = sqrt((30./(-sqrt(29) + 11))).*sin(angle(3)).*sin(angle(4)) - 2;
    c(5) = -2 - 3.*sin(angle(1)).*sin(angle(2));
    c(6) = 3.*sin(angle(1)).*sin(angle(2)) -2;
end






