function new  = rk4(old,deriv,dt)
n = numel(deriv);

new = zeros(1,n);
k = zeros(4,n);

%first step (half step)
for i = 1:n
    k(1,i) = deriv{i}(old);
    new(i) = old(i) + k(1,i)*dt/2; 
end
    oldNew = new;

%2nd step (another half step using deriv(new))
for i = 1:n
    k(2,i) = deriv{i}(oldNew);
    new(i) = old(i) + k(2,i)*dt/2;
end
    oldNew = new;

%3rd step (full step)
for i = 1:n
    k(3,i) = deriv{i}(oldNew);
    new(i) = old(i) + k(3,i)*dt;
end
    oldNew =  new;
    
%get k4 (full step)
for i = 1:n
    k(4,i) = deriv{i}(oldNew);
end



%calculate phi
phi = zeros (1,n);
for i = 1:n
    phi(i) = (k(1,i)+2*k(2,i)+2*k(3,i)+k(4,i))/6;
end

%calculate final value;
for i = 1:n
new(i) = old(i) + phi(i)*dt;
end

end