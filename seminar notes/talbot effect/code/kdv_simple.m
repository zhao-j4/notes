%% Parameters 
    Nref=2^14;  % Size of spatial mesh
    N=2^14;      % Fourier modes
    Jmax=1;    % Refining the time steps

    Xvec=pi*(-Nref/2+1:Nref/2)'/Nref*2; % Spacial mesh
    indexvec=(-Nref/2+1:Nref/2)'; % Vector of Fourier indices
    indexvec(Nref/2)=1; % Correct zero index so can use to rescale Fourier coefficients
    T=pi; % Final time

    nonlinear_part=@(u) 1/2*dx(conv1(u,u)); % define function for KdV splitting

    
u0_hat=fftpi(heaviside(Xvec));

for l = 1:Jmax 
        M=10^(l + 3);
        h=T/M; % Time step

        % Schratz and Hofmanova
        
        w_hat=u0_hat;
        indexvec=(-N/2+1:N/2);
        indexvec(N/2)=1;

        for m=1:M

            wstep=w_hat;
            w_hat=kdv_resonance_based_first_order(wstep,h);

            if mod(m, 10^(l+1)) == 0
            w = ifftpi(w_hat);
            
            plot(Xvec, real(w));
            title(['Time: ', num2str(m/M) '*' num2str(T,3) ' with Time step:', num2str(h, 8)]);
            drawnow;
            end
        end
       
 end    

