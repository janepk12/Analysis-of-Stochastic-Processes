function x_filtered=smooth_MA(x,tau)

%number of points
num_pts=length(x);

%Radius of filter
r=(tau-1)/2;

%time stamps of filtered series
x_filtered=zeros(size(x));

%Moving average filter
for i=1:num_pts

    %----------------------------------------------------------------------   
    %effective radius at the boundary
    r_eff = min([i-1, num_pts-i, r]);
    
    %Fenstergrenzen:
    left=i-r_eff;
    right=i+r_eff;
    
    % Gewichte setzen (1 innerhalb, 0 außerhalb)
    weights=zeros(num_pts,1);
    weights(left:right) = 1;

    %----------------------------------------------------------------------
    %Calculate mean via least squares adjustment
    A=[ones(num_pts,1)];
    P=diag(weights);

    N=A'*P*A;
    n=A'*P*x;
    
    X_hat=N\n;

    x_filtered(i)=X_hat;
    
end