function x_filtered=smooth_SG(x,tau,p)

%number of points
num_pts=length(x);

%Radius of filter
r=(tau-1)/2;

%time stamps of filtered series
x_filtered=zeros(size(x));

%mapping the t values of the moving window on [-1,1]
u=linspace(-1,1,tau)';

A=[];
for j=0:p
    A=[A u.^j];
end

%Get Kernel
C=(A'*A)\A';
kernel=C(1,:);

%Moving average filter
for i=1:num_pts

    %----------------------------------------------------------------------      
    %window limits
    left  = max(1, i-r);
    right = min(num_pts, i+r);
    
    %INdex of points within window
    idx = left:right;
    
    %x values of points within window
    xi = x(idx);

    %Clipping kernel for the limits
    kerneli = kernel((left:right)-i+r+1);
    kerneli = kerneli/sum(kerneli);  %work around

    x_filtered(i)=kerneli*xi;
    
end