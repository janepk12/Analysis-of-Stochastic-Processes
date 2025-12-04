function x_filtered=smooth_MLS(t,x,tau,p)

%number of points
num_pts=length(x);

%Radius of filter
r=(tau-1)/2;

%time stamps of filtered series
x_filtered=zeros(size(x));

%gaussian function for the weights
gaussian=@(u) exp(-5*u.^2);

%Moving average filter
for i=1:num_pts

    %----------------------------------------------------------------------      
    %window limits
    left  = max(1, i-r);
    right = min(num_pts, i+r);
    
    %INdex of points within window
    idx = (left:right)';
    
    %points within window
    ti = t(idx);
    xi = x(idx);

    %mapping the t values of the moving window on [-1,1]
    ui = ti - t(i);
    ui = ui / max(abs(ui));
    
    %fitting polynomial inside window
    A=[];
    for j=0:p
        A=[A ui.^j];
    end

    %Gaussian weights
    weights=gaussian(ui);
    P=diag(weights);

    N=A'*P*A;
    x_hat=N\A'*P*xi;
    
    %get filtered value at u=0
    x_filtered(i)=x_hat(1);
    
end