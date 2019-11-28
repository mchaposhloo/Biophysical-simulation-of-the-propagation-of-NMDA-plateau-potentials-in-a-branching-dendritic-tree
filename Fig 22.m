delta = 0.03;
x_nullcline = 40*(-0.3:delta:1.1);
[x,y] = meshgrid(x_nullcline,40*(-0.15:(delta/4):0.2));
a = 0.2;
b = 0.00204;
c = 0.01;
v = x.*(a-x/40).*(x/40-1) - y;
w = b*x - c*y;

figure
quiver(x,y,v,w,2)
hold on
a = plot(x_nullcline, x_nullcline.*(a-x_nullcline/40).*(x_nullcline/40-1));
set(a,'LineWidth',1)
set(a,'Color','black')
hold on 
c = plot(x_nullcline, b/c*(x_nullcline));
set(c,'LineWidth',1)
hold on
felan1 = 7:1:8;
f = 8.5;
f2 = 8.531665;
felan3 = 8.54:1:10.54;
startx =   [felan1 f f2 felan3];
starty =  zeros(size(startx));
b = streamline(x,y,v,w,startx,starty);
set(b,'LineWidth',1)
hold on
felan2 = 8.5316596411;
h = streamline(x,y,v,w,felan2,0);
set(h,'Color','red')
set(h,'LineWidth',1.5)
grid on
xlim([-10 40])
ylim([-1 5])
xlabel('V')
ylabel('w')
legend('vector field','V-nullcline','w-nullcline', 'trajectories')