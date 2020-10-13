%% Rukaiya Batliwala
%Assignment4

%% Dirichlet k = 1
clc
clear all
close all



L = 1; Uo= 1;v = 1;A = 1;N = 200;N2 = 400;k =1 ; f =1; f2 =1 ;
h = L/(N+1);x = [0:h:L];h2 = L/(N2+1);x2 = [0:h2:L];

U = [];U2 = [];
  
for i = 1:1:length(x);
    U(i) =(((sinh(k*L - k*x(i)) + sinh(k*x(i)))/(sinh(k*L))) - 1)*(A/k^2) + Uo*(sinh(k*L  - k*x(i))/(sinh(k*L)));    
end

for i = 1:1:length(x2);
    U2(i) =(((sinh(k*L - k*x2(i)) + sinh(k*x2(i)))/(sinh(k*L))) - 1)*(A/k^2) + Uo*(sinh(k*L  - k*x2(i))/(sinh(k*L)));    
end

for i = 2:1:length(U)-1
    u(i) = (U(i+1)+ U(i-1)-(h^2)*(1+(k^2)*U(i)))/2;
end
u(1) = Uo;
u(end+1) = 0;

for i = 2:1:length(x2)-1
    u2(i) = (U2(i+1)+ U2(i-1)-(h2^2)*(1+(k^2)*U2(i)))/2;
end
u2(1)= Uo;
u2(end+1) = 1;

a = zeros(N);
for i = 1:N
    a(i,i) = -2 - h^(2)*k^(2);
end
for i = 1:1:N-1
    a(i+1,i) = 1;
    a(i,i+1) = 1;
end

 F = zeros(N,1);
 for i = 1:1:length(F)
 F(i) = (i+A)*h^(2);
 end
 F(1) = F(1)-Uo;
 
 uf= a\F;

a2 = zeros(N2);
for i = 1:N2
    a2(i,i) = -2-h2^(2)*k^(2);
end
for i = 1:1:N2-1
    a2(i+1,i) = 1;
    a2(i,i+1) = 1;
end
 
 F2 = zeros(N2,1);
 for i = 1:1:length(F2)
 F2(i) = (i+A)*h2^(2);
 end
F2(1) = F2(1)+2*h*v;
 
 uf2= a2\F2;

for i = 2:1:length(uf)-1
  ErrorN(i) = abs(u(i)-uf(i));
end

for i = 2:1:length(uf2)-1
 ErrorN2(i) = abs(u2(i) - uf2(i));
  end
 
 plot(x,u,'-b')
 hold on
 plot(x(2:end-1),uf,'-r')
legend('Exact Solution','Numerical Solution')
title('k = 1,Dirichlet solution')
xlabel('x value')
ylabel('U value')
grid on 
saveas(gcf,'Figure1.pdf')
 
 F_A = abs(log((max(ErrorN))/(max(ErrorN2)))/(log(2)));
 Dirichlet_1_Table = [u(end-20:end), uf(end-20:20)];
 writematrix(Dirichlet_1_Table,'K_1_Dirichlet.xls');

%%  Dirichlet k = 10
clear all





L = 1; Uo= 1;v = 1;A = 1;N = 90;N2 = 180;k =10 ; f =1; f2 =1 ;
h = L/(N+1);x = [0:h:L];h2 = L/(N2+1);x2 = [0:h2:L];

U = [];U2 = [];
  
for i = 1:1:length(x);
    U(i) =(((sinh(k*L - k*x(i)) + sinh(k*x(i)))/(sinh(k*L))) - 1)*(A/k^2) + Uo*(sinh(k*L  - k*x(i))/(sinh(k*L)));    
end

for i = 1:1:length(x2);
    U2(i) =(((sinh(k*L - k*x2(i)) + sinh(k*x2(i)))/(sinh(k*L))) - 1)*(A/k^2) + Uo*(sinh(k*L  - k*x2(i))/(sinh(k*L)));    
end

for i = 2:1:length(U)-1
    u(i) = (U(i+1)+ U(i-1)-(h^2)*(1+(k^2)*U(i)))/2;
end
u(1) = Uo;
u(end+1) = 0;

for i = 2:1:length(x2)-1
    u2(i) = (U2(i+1)+ U2(i-1)-(h2^2)*(1+(k^2)*U2(i)))/2;
end
u2(1)= Uo;
u2(end+1) = 1;

a = zeros(N);
for i = 1:N
    a(i,i) = -2 - h^(2)*k^(2);
end
for i = 1:1:N-1
    a(i+1,i) = 1;
    a(i,i+1) = 1;
end

 F = zeros(N,1);
 for i = 1:1:length(F)
 F(i) = (i+A)*h^(2);
 end
 F(1) = F(1)-Uo;
 
 uf= a\F;

a2 = zeros(N2);
for i = 1:N2
    a2(i,i) = -2-h2^(2)*k^(2);
end
for i = 1:1:N2-1
    a2(i+1,i) = 1;
    a2(i,i+1) = 1;
end
 
 F2 = zeros(N2,1);
 for i = 1:1:length(F2)
 F2(i) = (i+A)*h2^(2);
 end
F2(1) = F2(1)+2*h*v;
 
 uf2= a2\F2;

for i = 2:1:length(uf)-1
  ErrorN(i) = abs(u(i)-uf(i));
end

for i = 2:1:length(uf2)-1
 ErrorN2(i) = abs(u2(i) - uf2(i));
end
 
  figure()
 plot(x,u,'-b')
 hold on
 plot(x(2:end-1),uf,'-r')
 legend('Exact Solution','Numerical Solution')
title('k = 10,Dirichlet solution')
xlabel('x value')
ylabel('U value')
grid on 
saveas(gcf,'Figure2.pdf')
 
 F_A = abs(log((max(ErrorN))/(max(ErrorN2)))/(log(2)));
Dirichlet_2_Table = [u(end-20:end), uf(end-20:20)];
 writematrix(Dirichlet_2_Table,'K_10_Dirichlet.xls');


%% Neumann k = 1

L = 1;Uo= 1;v = 1;A = 1;N = 100;N2 = 200; k =1 ;

h = L/(N+1);
x = [0:h:L];
h2 = L/(N2+1);
x2 = [0:h2:L];

U = [];
dv_dx = v ;  %x = 0

for i = 1:1:length(x);
    U(i) = ((((cosh(k*x(i)))/(cosh(k*L)))-1)*(A/(k^2))) - ((v/k)*((sinh(k*L - k*x(i)))/(cosh(k*L))));
end
for i = 1:1:length(x2);
    U2(i) = ((((cosh(k*x2(i)))/(cosh(k*L)))-1)*(A/(k^2))) - ((v/k)*((sinh(k*L - k*x2(i)))/(cosh(k*L))));
end

u=[];
for i = 1:1:length(U)-2
    u(i) = -2*h*v + 4*U(i+2)-3*U(i+1);
end

for i = 1:1:length(U2)-2
     u2(i) = -2*h*v + 4*U2(i+2)-3*U2(i+1);
end
u(end+1) =  0;
u2(end+1) = 0;
a = zeros(N+1);
for i = 1:N+1
    a(i,i) = 2*h^(2)*k^(2);
end
for i = 1:1:N
    a(i+1,i) = 1;
    a(i,i+1) = 1;
end
a(1,2) = 2;

 F = zeros(N+1,1);
 for i = 1:1:length(F)
 F(i) = (i+A)*h^(2);
 end
 F(1) = F(1)- 2*h*v;
 
 uf= a\F;



x(end) = [];
x2(end) = [];
XL = length(x);
figure()

 hold on
 plot(x,u,'-b')
 hold on
 plot(x,uf,'-r')
 grid on
legend('Exact Solution','Numerical Solution')
title('k = 1,Neumann solution')
xlabel('x value')
ylabel('U value')
saveas(gcf,'Figure3.pdf')
Neumann_1_Table = [u(end-20:end), uf(end-20:20)];
 writematrix(Neumann_1_Table,'K_1_Neumann.xls');

 
%% Neumann k = 10
   
L = 1;Uo= 1;v = 1;A = 1;N = 200;N2 = 400; k =10 ;

h = L/(N+1);
x = [0:h:L];
h2 = L/(N2+1);
x2 = [0:h2:L];

U = [];
dv_dx = v ;  %x = 0

for i = 1:1:length(x);
    U(i) = ((((cosh(k*x(i)))/(cosh(k*L)))-1)*(A/(k^2))) - ((v/k)*((sinh(k*L - k*x(i)))/(cosh(k*L))));
end
for i = 1:1:length(x2);
    U2(i) = ((((cosh(k*x2(i)))/(cosh(k*L)))-1)*(A/(k^2))) - ((v/k)*((sinh(k*L - k*x2(i)))/(cosh(k*L))));
end

u=[];
for i = 2:1:length(x)-2
    u(i) = (-h^(2)*1 + U(i-1) + U(i+2))/2;
end

for i = 2:1:length(x2)-2
    u2(i) = (-h2^(2)*1 + U2(i-1) + U2(i+2))/2;
end
u(end+1) =  0;
u2(end+1) = 0;
a = zeros(N+1);
for i = 1:N+1
    a(i,i) = 2*h^(2)*k^(2);
end
for i = 1:1:N
    a(i+1,i) = 1;
    a(i,i+1) = 1;
end
a(1,2) = 2;

 F = zeros(N+1,1);
 for i = 1:1:length(F)
 F(i) = (i+A)*h^(2);
 end
 F(1) = F(1)- 2*h*v;
 
 uf= a\F;



x(end) = [];
x2(end) = [];
XL = length(x);
figure()

 hold on
 plot(x,u,'-b')
 hold on
 plot(x,uf,'-r')
 grid on
legend('Exact Solution','Numerical Solution')
title('k = 10,Neumann solution')
xlabel('x value')
ylabel('U value')
saveas(gcf,'Figure4.pdf')

Neumann_2_Table = [u(end-20:end), uf(end-20:20)];
 writematrix(Neumann_2_Table,'K.xls');
