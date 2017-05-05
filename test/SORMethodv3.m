% This MATLAB m-file is to be used for the project for the Scientific
% Computing for Mechanical Engineers Course.
% The purpose of this m-file is to solve the 2D Helmholtz equation with
% various boundary conditions in a square domain in Cartesian Coordinates
% through the Successive Over-Relaxation Method
% This m-file was written by John Jimenez (PSID:1253921)
% Project Code: AHc2-3

clc,clearvars % Clears command window & workspace of unnecessary clutter from a previous run of another m-file.
disp('2D Helmholtz Equation Solver: SOR Method')
disp('Solves the 2D Helmholtz Equation within a square domain using the SOR Method')

%% Dimensions of the Square Domain and Some Important Parameters
%  Note: the domain of interest is the square:
% -pi<x<pi and -pi<y<pi

    % Start and End points for domain in x

        ax = -pi;    % Start point of domain in x
        bx =  pi;    % End point of domain in x

    % Start and End points for domain in y are the the same as that in x for Project Code.

        ay = -pi;   % Start point of domain in y
        by =  pi;   % End point of domain in y

    % A value for Lambda can be assigned directly but is left up to user input for greater versatility

        Lambda = input('Please input a value for Lambda:'); Lambda = -Lambda;
        G      = input('Please input a value for G. This will be used for SOR and is limited for values between 1 and 2:'); % Value needed for SOR.
    
            if G > 2
                disp('Please try again. Your selected value for G is incompatible for SOR.')
                disp('Restart process through typing in SORMethodv3.')
            return 
            end
  
    % Assigning an acceptable limit of error for convergence

        % Es = 10^-7;
        Es = input('Please input a value for the acceptable limit of error for convergence (e.g. 10^7):')

%% Boundary Conditions Setup and Evaluation of Relevant Functions

    % Defines the Interior nodes for both x and y through user input. However,these do not
    % include exterior boundary points.

        IntNodesX = input('Please input a value for the no. of nodes for the  x-coordinate:');
        IntNodesY = input('Please input a value for the no. of nodes for the  y-coordinate:');

    % Safeguard to ensure a consistent, square mesh.
        if IntNodesX ~= IntNodesY
                disp('Please try again. This m-file performs the task only for square matrices.')
                disp('Restart process through typing SORMethodv3.')
                return 
        end
        
    % Generating a range of linearly spaced vectors for x and y through the variables IntNodesX
    % and IntNodesY
    % IMPORTANT: THESE ARE THE INPUT VALUES FOR EVERY FUNCTION THAT NEEDS
    % THEM!

        xval = linspace(ax,bx,IntNodesX); % Start and End nodes for x can be redefined through changing ax and bx above.
        yval = linspace(ay,by,IntNodesY); % Start and End nodes for y can be redefined through changing ay and by above.
        [X,Y]= meshgrid(xval,yval);
        Y = flipud(Y);

    % Value of spacing for Each Node

        dx = (bx-ax)/IntNodesX;
        dy = dx;


    % Defining some coefficients that remain the same for all matrix
    % operations. Refer to report for details concerning the discretization

        A = Lambda - (2 * (dx^2)) - (2 * (dy^2));
        B = (dx^2) * (dy^2);
        U = zeros(IntNodesX,IntNodesY);

    % The following functions below are needed for the evaluation of the
    % boundary conditions as well as the function F(x,y)
    
        fb  =  Y*(by-Y).^2;
        gb  =  (by-Y).^2*(cos((pi*Y)/by));
        uw  =  ((ay*(by-ay).^2))+(((X-ax)/(bx-ax))*((cos((pi*ay)/by)*(by - ay).^2) - ((ay*(by-ay).^2))));
      
      % F = zeros(IntNodesX,IntNodesY); % This is for the Final Simulation
    
        F =  sin(pi*(X-ax)/(bx-ax)).*cos((0.5*pi)*(2*((Y-ay)/(by-ay))+1));
   
     % Dirichlet Boundary Conditions on x-axis
         U(2:IntNodesX-1,1) = Y(2:IntNodesX-1,1).*(by-Y(2:IntNodesX-1,1)).^2;
         U(2:IntNodesX-1,IntNodesX) = ((by-Y(2:IntNodesX-1,1)).^2).*cos((pi*Y(2:IntNodesX-1,1)/by));
     
     % Dirichlet Boundary Conditions on y-axis
        U(IntNodesY,:) = ((ay*(by-ay).^2))+(((X(IntNodesY,:)-ax)/(bx-ax))*((cos((pi*ay)/by)*(by - ay).^2) - ((ay*(by-ay).^2))));
     

    w = 2/(1+sin(pi*dx/(2*pi)));
    e = 1;

%% The Successive Over Relaxation (SOR) Method Execution
    while e > Es

    U1 = U; % Initial Guess
    for a = 2:IntNodesX-1
        for b = 2:IntNodesY-1
            W(a,b) = U(a,b); % Stores Results of Previous Iteration
            if a == 2
                U(1,b) = (F(1,b)*dx^2) - (2*U(2,b)) - U(1,b+1) - U(1,b-1); % Expression Derived from Discretized Equation. Refer to Report for details
                U(1,b) = (G * U(1,b)) + ((1-G) * W(1,b));
            else
                U(a,b) = (((F(a,b)*(dx^2)) - U(a-1,b) - U(a,b-1) - U(a+1,b) - U(a,b+1))/((Lambda*dx^2)-4)); 
                U(a,b) = (G * U(a,b)) + ((1-G) * W(a,b));
            end
        end
    end
    
    e = max(max(abs((U-U1)./U))); % Calculates Relative Error for comparison
    
end

Value = mean(mean(U(2:IntNodesX-1,2:IntNodesX-1).^2))

%% Plotting Results
    figure(1) % VECTOR FIELD U
    set(gcf,'units','normalized','position',[0.02 0.52 0.3 0.32]);
    % surf(X,Y,U);
    mesh(X,Y,U);
    rotate3d
    xlabel('x'); ylabel('y'); zlabel('U');
    title('Solution using SOR Method','fontweight','normal');