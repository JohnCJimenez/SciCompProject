% This MATLAB m-file is to be used for the project for the Scientific
% Computing for Mechanical Engineers Course. 
% The purpose of this m-file is to solve the 2D Helmholtz equation with
% various boundary conditions in a square domain in Cartesian Coordinates
% through Gaussian Elimination
% This m-file was written by John Jimenez (PSID:1253921)
    % Project Code: AHc2-3

clc,clearvars % Clears command window & workspace of unnecessary clutter from a previous run of another m-file.
disp('2D Helmholtz Equation Solver: Gaussian Elimination') 

%% Dimensions of the Square Domain and Some Important Parameters
    % Note: the domain of interest is the square: 
            % -pi<x<pi and -pi<y<pi 
    % However, in order to make the computation more efficient, 
    % the domain is then truncated to 
            %  0<x<pi and 0<y<pi
    % As such, the overall solution will reflect the changes through the
    % factor of 4 which will be introduced into the approximation. 

    % Start and End points for domain in x
        
        ax = 0;   % Start point of domain in x
        bx = pi;  % End point of domain in x

    % Start and End points for domain in y are the the same as that in x for Project Code. 
        
        ay = ax;  % Start point of domain in y 
        by = bx;  % End point of domain in y

    % A value for Lambda can be assigned directly but is left up to user input for greater versatility 
        
        Lambda = input('Please input a value for Lambda:'); 

    % Assigning an acceptable limit of error for convergence
        
        Es = 10^-7;
        
%% Assembling the Coefficient Matrix & the F(x,y) Vector such that: Au = F (A is a pentadiagonal matrix)

    % Defines the Interior nodes for both x and y through user input. However,these do not
    % include exterior boundary points.

        IntNodesX = input('Please input a value for the no. of nodes for the  x-coordinate:');
        IntNodesY = input('Please input a value for the no. of nodes for the  y-coordinate:');
    
    % Due to the Neumann boundary condition, a ghost node must be
    % introduced in the y-boundary condition. 
       
        IntNodesX1 = IntNodesX+2; 
        IntNodesY1 = IntNodesY+2;

    % Generating a range of linearly spaced vectors for x and y through the variables IntNodesX
    % and IntNodesY 
        % IMPORTANT: THESE ARE THE INPUT VALUES FOR EVERY FUNCTION THAT NEEDS
        % THEM!
        
        xval = linspace(ax,bx,IntNodesX1); % Start and End nodes for x can be redefined through changing ax and bx above. 
        yval = linspace(ay,by,IntNodesY1); % Start and End nodes for y can be redefined through changing ay and by above. 
    
    % Value of spacing for Each Node   
       
        dx = xval/length(xval);
        dy = yval/length(yval);
    
    % Defining some coefficients that remain the same for all matrix
    % operations. Refer to report for details concerning the discretization
    
        A = Lambda - (2 * (dx.^2)) - (2 * (dy.^2));
        B = (dx.^2) .* (dy.^2); 
        
    % The following functions below are needed for the evaluation of the 
    % boundary conditions as well as the function F(x,y)
        
        fb  =  @(y)   y*(by-y).^2;
        gb  =  @(y)   (by-y).^2*(cos((pi*y)/by));
        Fxy = @(x,y)  sin(pi*(x-ax)/(bx-ax))*cos((0.5*pi)*(2*((y-ay)/(by-ay))+1));
    
     
    
   





