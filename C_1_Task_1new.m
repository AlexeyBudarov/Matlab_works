grid on;
f_1 = @(x,y) x.^2 - x.*sin(2.*y)./2+y.^2 - 3.*x -3;
f_2 = @(x,y) -x.^2+3.*x.*y+5.*y.^2+5.*x - 7;
ezplot(f_1)
hold on;
ezplot(f_2)

while true
    [x, y, button] = ginput(1);
    if button == 1
        temp_function = @(x) Fsistem(f_1,f_2,x);
        root = fsolve(temp_function, [x, y]);
        plot(root(1), root(2), '*r');
        hold on;
        text(root(1)+0.2, root(2)+0.2, sprintf('(%.5f, %.5f)', root(1), root(2)))
        disp([root(1), root(2)])
    elseif button == 3
        break;
    end
end



function f = Fsistem(f_1, f_2, x)
    f(1) = f_1(x(1),x(2));
    f(2) = f_2(x(1),x(2));
end



