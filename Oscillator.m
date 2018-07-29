classdef Oscillator
    properties(Access = public)
        mass = 5;
        cv = 4;%viscous damping constant
        cc = 0;%Coulomb damping constant
        k = 100;%spring constant
        isPendulum = false;
       pendulumLength = 1;
       

        tSpace;
        pSpace;
        vSpace;
        p0 = pi/2;
        v0 = 0;
        animationSpeed = 1;
        animationScale = 1; %zooms the animation in or out if its not a pendulum
        dt = 0.03;
        finalTime = 4;
        
         externalForce = @(X) 0;
         %{
         external force must take a 3 element vector, X
         X(1) : time 
         X(2) : position (angle if its a pendulum)
         X(3) : velocity (angular velocity if pendulum)
         %}
    end
    methods (Access = public)
        
        function this = Animate(this,fromData)

            if nargin == 1
            this = this.Calculate;
            end
            if nargin ==2 && fromData
            this = this.Extract;
            end
                
                close all;
                 figure(1);
                if this.isPendulum
                        w = this.pendulumLength;
                    for n = 1:length(this.tSpace)
                        drawline(0,0,-w*sin(this.pSpace(n)),-w*cos(this.pSpace(n)));
                        xlim([-w, w])
                        ylim([-w, w])
                        pause(this.dt);%10x speed
                    end
                
                end
                if ~this.isPendulum
                    m = this.mass;
                        w = this.animationScale;
                    ymax = 1.5;
                
                     for n = 1:length(this.tSpace)
                         clf;
                        rectangle('Position',[this.pSpace(n)-this.animationScale/3,-1/2,2*this.animationScale/3,1])
                        xlim([-w, w])
                        ylim([-ymax, ymax])
                        pause(this.dt/this.animationSpeed);
                        
                    end
                end
        end
        function Plot(this)
            plot(this.tSpace,this.pSpace)
        end
        function this = Calculate(this)
                    this.tSpace = 0:this.dt:this.finalTime;
                m = this.mass;
                if this.isPendulum
                w = this.pendulumLength;
                I = m*w^2/3;
                else
                    I = m;
                end
                numIterations = length(this.tSpace);
                X = zeros(3,numIterations);
                X(1,1) = this.tSpace(1);
                X(2,1) = this.p0;
                X(3,1) = this.v0;
                deriv = cell(1,3);
                deriv{1} = @(X) 1; % dt/dt
                deriv{2} = @(X) X(3); %theta dot
                deriv{3} = @(X) -(this.k*X(2)+this.cv*X(3)+this.cc*sign(X(3))-this.externalForce(X))/I; %theta dot dot
                    for i = (2:numIterations)
                        X(:,i) = this.rk4(X(:,i-1),deriv,this.dt);%rk4.m does 1 rk4 step
                    end
                this.tSpace = X(1,:);
                this.pSpace = X(2,:);
                this.vSpace = X(3,:);
        end
        function this = Extract(this)
        this.p0 = this.pSpace(1);
        this.v0 = this.vSpace(1);
        this.dt = this.tSpace(2)-this.tSpace(1);
        this.finalTime = this.tSpace(end);
        end
        function new = rk4(~,old,deriv,dt)
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
        function drawline(x1,y1,x2,y2)
            X = [x1,x2];
            Y = [y1,y2];
            plot (X,Y);
        end
        
    end
end