%ELEC 4700
%Assignment 1
%Part 3
%Tariq Aboushaer
%101064544

clear all
clearvars
close all
clc

    %% 
    global C X Y
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_n = 0.26*C.m_0;                 % effective mass of electrons
    
    %%
    %Starting with 100 particles to test then going to 1000
    Particles = 1000;
    
    %Assume Temp = 300 K.
    Temp = 300;
    
    %The nominal size of the region is 200 nm × 100 nm.
    Xregion = 2e-7;
    Yregion = 1e-7;
    
    % Typically the spacial step should be smaller than 1/100 of the region size
    Step_size = 1e-9;
    Timestep = 1000;
    
    %Calculating the Thermal Velocity
    Vth = sqrt(2*C.kb*(Temp/C.m_n));
    
    %Equation for the Change of velocity for each timestep
    Vchange = Step_size/Vth;
    
    % Mean time between collisions is 0.2e-12
    MeanTime = 2e-13;
    
    %Sigma can be calculated by
    Sigma = sqrt(C.kb*Temp/C.m_n)/4;
    
    
    %%
    %Creating 2 random numbers that are used for the angle and position.
    %The two numbers will be used for each particle at X coordinates.
    %The Y coordinates only require the position
    %Then the coordinates are multiplied by the region to make sure that
    %particles are within the region 
    X = rand(2,Particles);
    Xposition(1,:) = X(1,:)*Xregion;
    Y = rand(1,Particles);
    Yposition(1,:) = Y(1,:)*Yregion;
    
    %The angle of each particle
    angle(1,:) = X(2,:)*2*pi;
    
    %Initial velocity of each particle using vth
    Xvelocity = Vth*Vchange*cos(angle(1,:));
    Yvelocity = Vth*Vchange*sin(angle(1,:));  
    
    % Scattering probability
    Pscat = 1 - exp(-Vchange/MeanTime);
    MFPvector = zeros(1,Particles);
    
    %Now Particle are not allowed to be added in the rectangle
    %First we need to check if the particle is inside X boundary
    CheckXLeft = Xposition > 0.8e-7;
    CheckXRight = Xposition < 1.2e-7;
    CheckX = CheckXLeft & CheckXRight;
    
    % Check which box it's in by checking Y coordinates
    CheckTop = Yposition > 0.6e-7;
    CheckBottom = Yposition < 0.4e-7;
    
    %The box
    BoxTop = CheckTop & CheckX;
    BoxBottom = CheckBottom & CheckX;
    InsideBox = BoxTop | BoxBottom;
    
    %Randomizing particles inside the box
    while(sum(InsideBox) > 0)
        
        TX = rand(1,sum(InsideBox));
        TY = rand(1,sum(InsideBox));
        Xposition(InsideBox) = TX*Xregion;
        Yposition(InsideBox) = TY*Yregion;
        
        %Checking again
        CheckXLeft = Xposition > 0.8e-7;
        CheckXRight = Xposition < 1.2e-7;
        CheckX = CheckXLeft & CheckXRight;
        CheckTop = Yposition > 0.6e-7;
        CheckBottom = Yposition < 0.4e-7;
        BoxTop = CheckTop & CheckX;
        BoxBottom = CheckBottom & CheckX;
        InsideBox = BoxTop | BoxBottom;
    end
    
    %%
    %Starting the loop of with the time steps
    %In this loop the electrons are moving using the initail position and
    %initial velocity.
    for i = 1:Timestep
        
        % Scattering all Electrons
        % When s is less than s, the particle scatters.
        s = rand(1,Particles);
        Scattered = s < Pscat;
        
        % New random angle:
        angle(Scattered) = rand*2*pi;
        % New random velocity:
        MBDistribution = makedist('Normal',Vth,Sigma);
        Velocity = random(MBDistribution,1,Particles);
        Xvelocity(Scattered) = Vchange*Velocity(Scattered) .*cos(angle(Scattered));
        Yvelocity(Scattered) = Vchange*Velocity(Scattered) .*sin(angle(Scattered));
        
        % Keep track of how many particles scattered
        % to figure out mean free path
        MFPvector(~Scattered) = MFPvector(~Scattered)+Step_size;
        % Anything that scattered will be set to 0
        MFPvector(Scattered) = 0;
        
        %To keep all the particles within the region, once their velocity
        %and position goes over the boundary it will subtract region_x to
        %make it go to the opposite side.
        Xleft = Xposition + Xvelocity < 0;
        Xright = Xposition + Xvelocity > Xregion;
        Xposition(Xleft) = Xposition(Xleft) + Xvelocity(Xleft) + Xregion;
        Xposition(Xright) = Xposition(Xright) + Xvelocity(Xright) - Xregion;
        
        %For Y boundary we will check if the particle is over the limit and
        %reverse its velocity by setting it to negative.
        %Once the particle is back their position would be updated.
        Yout = (Yposition + Yvelocity > Yregion | Yposition + Yvelocity < 0);
        Yvelocity(Yout) = -1*Yvelocity(Yout);
        Yposition(1,:) = Yposition(1,:) + Yvelocity(1,:);
              

       
        % Check if particles are inside the X boundary
        CheckXLeft = (Xposition + Xvelocity) > 0.8e-7;
        CheckXRight = (Xposition + Xvelocity) < 1.2e-7;
        CheckX = CheckXLeft & CheckXRight;
        
        % Checking the Y boundary for the Bottom box
        CheckBottom = (Yposition + Yvelocity) < 0.4e-7;
        
        % Checking the bottom box
        BoxBottom = CheckBottom & CheckX;
        Xvelocity(BoxBottom) = -1*Xvelocity(BoxBottom);
        
        % Cheking if we need to multipy by -1
        CheckXLeft = (Xposition + Xvelocity) > 0.8e-7 + Step_size;
        CheckXRight = (Xposition + Xvelocity) < 1.2e-7 - Step_size;
        CheckX = CheckXLeft & CheckXRight;
        CheckYbottomBoxtop = Yposition < 0.4e-7 - Step_size;
        YBoxbottom = CheckX & CheckYbottomBoxtop;
        Yvelocity(YBoxbottom) = -1*Yvelocity(YBoxbottom);
        
        % Check if particles are inside X boundary
        CheckXLeft = (Xposition + Xvelocity) > 0.8e-7;
        CheckXRight = (Xposition + Xvelocity) < 1.2e-7;
        CheckX = CheckXLeft & CheckXRight;
        
        
        % Checking the Y boundary for the Top box
        CheckTop = (Yposition + Yvelocity) > 0.6e-7;
        
        % Checking the bottom box
        BoxTop = CheckTop & CheckX;
        Xvelocity(BoxTop) = -1*Xvelocity(BoxTop);
        
        % Cheking if we need to multipy by -1
        CheckXLeft = (Xposition + Xvelocity) > 0.8e-7 + Step_size;
        CheckXRight = (Xposition + Xvelocity) < 1.2e-7 - Step_size;
        CheckX = CheckXLeft & CheckXRight;
        CheckYtopBoxbottom = Yposition > 0.6e-7 + Step_size;
        YBoxtop = CheckX & CheckYtopBoxbottom;
        Yvelocity(YBoxtop) = -1*Yvelocity(YBoxtop);
        
         %D particles are the particles that are not in Xleft or Xright.
        D_particles = ~(Xleft | Xright | BoxTop | BoxBottom);
        Xposition(D_particles) = Xposition(D_particles) + Xvelocity(D_particles);
       
        %the two lines below are to save the timesteps for plotting
        Xplot(i,:) = Xposition(1,:);
        Yplot(i,:) = Yposition(1,:);

        % Calculating temperature
        Avelocity = sum(Velocity)/Particles;
        Currentt(1,i) =  Avelocity^2*C.m_n/(C.kb*2);
       
        % Calculating the Mean Free Path 
        MFP = sum(MFPvector)/Particles;
        avg_tbc = MFP/Avelocity;
        
    end
    
      %%

     figure(1)
     ElectronDensity = [Xposition',Yposition'];
     hist3(ElectronDensity,'CdataMode','auto')
     xlabel('X (m)')
     ylabel('Y (m)')
     title('Electron Density Map (TA 101064544)')
     colorbar
     view(2)
     
     
      Xcan = linspace(0,Xregion,10);
      Ycan = linspace(0,Yregion,10);
      Xbin = discretize(Xposition,Xcan);
      Ybin = discretize(Yposition,Ycan);
      Tbins = zeros(10,10);
    for x = 1:10
        for y=1:10
            LX = Xbin == x;
            LY = Ybin == y;
            L = LX & LY;
            
            Xsum = sum(Xvelocity(L)) / Vchange;
            Ysum = sum(Yvelocity(L)) / Vchange;
            
            AvgVelocity = sqrt((Xsum)^2 + (Ysum)^2);
            % Use average velocity to calculate temperature
            % and put it into the proper bin
            Tbins(x,y) = AvgVelocity^2*C.m_n/(C.kb*2);
        end
    end
    
        figure(2)
    surf(Tbins)
    title('Temperature Map (TA 101064544)')
    colorbar
    
    BoxPicX = [0.8e-7 0.8e-7  1.2e-7 1.2e-7];
    BoxTopY = [1e-7 0.6e-7 0.6e-7 1e-7];
    BoxBottomY = [0 0.4e-7 0.4e-7 0];
    
     for row = 1:Timestep
         
        figure(3)
        plot(Xplot(row,:),Yplot(row,:),'o')
        xlim([0 Xregion])
        ylim([0 Yregion])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('2-D Plot of Particle Trajectories (TA 101064544)')
        legend(['Number of Particles:' num2str(Particles)])
        hold on
        plot(BoxPicX,BoxTopY,'k')
        plot(BoxPicX,BoxBottomY,'k')
    end
