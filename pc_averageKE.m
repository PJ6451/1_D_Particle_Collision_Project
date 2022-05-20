function AK = pc_averageKE(time, number_of_particles, len, left_energy, right_energy)

%Description: Models a one dimensional chain with heat baths on either end
%providing energy to particles, which then colide with both each other and
%the baths. Outputs collision count of particles over the time limit

X =zeros(3,number_of_particles);                                            %plotting the differences between initial/terminal properties at given time instances might be easier with a matrix, using diff()

rng default;
reset(RandStream.getGlobalStream,sum(100*clock));

X(1,:) = linspace(len/(number_of_particles+1),len*(1-1/(number_of_particles+1)),number_of_particles);  %positions, evenly spaced
X(2,:) = rand(1,number_of_particles)+randi(10,1,number_of_particles);  % random masses between 1 and 2
X(3,:) = sqrt(-2./X(2,:)*1/((left_energy+right_energy)/2).*log(rand(1,number_of_particles)));  % To start, initialize with exponentially generated energies halfway between left and right.  Exponentially generated energies are more realisitic
                                                                                               % and I believe can be computed by -E*log(U) where U is a uniform random variable [might want to check this] and E is the energy

%CREATING A ROW VECTOR, ONLY CAN BE 0 OR 1 AS RANDI FROM 1 TO 2 MINUS 1 IS
%EITHER 1 OR 2
direction=randi(2,1,number_of_particles)-1;
direction(~direction)=-1;  %vector of +1 and -1 , to randomly assign directions, 0 BECOMES -1
X(3,:)=X(3,:).*direction;

z = 1;
time1 = 0;
Average_Kinetic = zeros(1,number_of_particles);
ttt = 0;

%Position 0 will be the cold bath, position "len" will be the hot bath
collisioncount = zeros(1,number_of_particles);
while time1 < time
    t = zeros(1,size(X,2));
    for i = 1:number_of_particles        
        t1 = inf;
        t2 = inf;
        j  = i+1;
        if i < number_of_particles
            if X(3,i) > 0                                                       %calculates collision times for paritlces going right (toward hot bath) 
                if X(3,i) > X(3,j)
                    t1 = (X(1,j)-X(1,i))/(X(3,i)-X(3,j));
                end
            elseif X(3,i) < 0                                                   %Calculates collision times for particles going left (toward cold bath)
                if i == 1
                    j = i+1;
                    t1 = -X(1,i)/X(3,i);
                    if X(3,i) > X(3,j)
                        t2 = (X(1,j)-X(1,i))/(X(3,i)-X(3,j));
                    end
                elseif X(3,i) > X(3,j)
                        t1 = (X(1,j)-X(1,i))/(X(3,i)-X(3,j));
                end
            end
            t(i) = min([t1 t2]);
        else
            if X(3,number_of_particles) > 0
                t(end) = (len-X(1,i))/X(3,i);
            else
                t(end) = inf;
            end
        end     
    end
   
    [t,I]  = min(t);
    
    if z == 1
        Average_Kinetic(z,:) = Average_Kinetic(z,:) + (t/time)*1/2*X(2,:).*(X(3,:).^2);
        ttt(z) = 0;
    else
        Average_Kinetic(z,:) = Average_Kinetic(z-1,:) + (t/time)*1/2*X(2,:).*(X(3,:).^2);
        ttt(z) = time1;
    end
    
    X(1,:) = X(1,:) + X(3,:)*t;
    if I==1  
        if X(1,I) <= 0  %% Collision with left wall
            X(1,I) = 0;
            X(2,I) = rand+randi(10,1);
            X(3,I) = sqrt(-2*left_energy*log(rand)/X(2,I));
        else  %% Collision with particle right of I
            J=I+1;
            collisioncount([I,J])=collisioncount([I,J])+1;
            v1 = X(3,I);
            v2 = X(3,J);
            X(3,I) = ((X(2,I)-X(2,J))*v1 + 2*X(2,J)*v2)/(X(2,I)+X(2,J));
            X(3,J) = (2*X(2,I)*v1 + (X(2,J)-X(2,I))*v2)/(X(2,I)+X(2,J));
        end
    elseif I == number_of_particles  %% Collision with right wall
        X(1,I) = len;
        X(2,I) = rand+randi(10,1);
        X(3,I) = -sqrt(-2*right_energy*log(rand)/X(2,I));
    else
        J = I+1;
        collisioncount([I,J])=collisioncount([I,J])+1;
        v1 = X(3,I);
        v2 = X(3,J);
        X(3,I) = ((X(2,I)-X(2,J))*v1 + 2*X(2,J)*v2)/(X(2,I)+X(2,J));
        X(3,J) = (2*X(2,I)*v1 + (X(2,J)-X(2,I))*v2)/(X(2,I)+X(2,J));
    end
    time1=time1+t;
    if sum(collisioncount) > time*10
        time1 = time;
    elseif t <= 0
        time1 = time;
    end
    z = z+1;
end

AK = [Average_Kinetic(end,:) X(2,:)];

end

