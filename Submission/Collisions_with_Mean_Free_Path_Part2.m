%ELEC 4700
%Assignment 1
%Part 2
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
    
    % This part is to set up the Histogram for the Maxwell-Boltzmann
    % Distribution. The average of Vth is used.
    
    MBDistribution = makedist('Normal',Vth,Sigma);
    Velocity = random(MBDistribution,1,Particles);
    figure(1)
    
    %Matlab ’hist' function
    hist(Velocity)
    title('Velocity Histogram (TA 101064544)')

    % Scattering probability
    Pscat = 1 - exp(-Vchange/MeanTime);
    MFPvector = zeros(1,Particles);
    
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
        
        %D particles are the particles that are not in Xleft or Xright.
        D_particles = ~(Xleft | Xright);
        Xposition(D_particles) = Xposition(D_particles) + Xvelocity(D_particles);
       
        %For Y boundary we will check if the particle is over the limit and
        %reverse its velocity by setting it to negative.
        %Once the particle is back their position would be updated.
        Yout = (Yposition + Yvelocity > Yregion | Yposition + Yvelocity < 0);
        Yvelocity(Yout) = -1*Yvelocity(Yout);
        Yposition(1,:) = Yposition(1,:) + Yvelocity(1,:);
              
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
    
      for paths = 1:Particles
        
        figure(2)
        plot(Xplot(:,paths),Yplot(:,paths),'-')
        xlim([0 Xregion])
        ylim([0 Yregion])
        xlabel('X (m)')
        ylabel('Y (m)')
        title('2-D Plot of Particle Trajectories (TA 101064544)')
        legend(['Number of Particles:' num2str(Particles)])
        hold on
    end
    for row = 1:Timestep
        
        figure(3)
        plot(row,Currentt(row),'.k');
        title('Temperature plot (TA 101064544)')
        xlabel('Time-step')
        ylabel('Temperature (K)')
        legend(['Current Temperature:' num2str(Currentt(row))])
        hold on
        pause(0.01)
    end