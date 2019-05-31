function z=Gaussian(x,y,mu,sigma)
x=x-mu(1);
y=y-mu(2);
        z=1/(2*pi*sigma^2) * exp(-1*(x.^2 + y.^2)/(2*sigma^2));
end