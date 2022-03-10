function y=blackbody(T,Edistribution)
q = 1.602176565e-19;
h = 6.62606957e-34;
kbeV=8.6173324e-5;% technically kT/q, approx T=302K
c = 29979245800;
i=0;
% bb=0;
bb(:,1)=Edistribution';
for E=Edistribution%=0.001:0.001:4
    i=i+1;
    bb(i,2)=q*q*2*pi*1e3*(q^2/h^3/c^2)*power(E,2).*exp(-E./kbeV/T);%in unit of mA 
%     bb(i,2)=2*pi/(heV^3*c^2)*power(E,2)*exp(-E/kbeV/T);
end
y=bb;
end 