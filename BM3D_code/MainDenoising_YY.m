

% Scale the image (BM3D processes inputs in [0,1] range)
scale_range = 0.7;

% Data = importdata('DATAFILE');
% Alpha_arr =  importdata('ALPHA_NOISE_PARAMETER_FILE');
% Sigma_arr =importdata('SIGMA_NOISE_PARAMETER_FILE');
Data = imageFinal((700-511):(701+511),(700-511):(701+511),:);
Alpha_arr = 0.0226*ones(size(imageFinal,3));
Sigma_arr = zeros(size(imageFinal,3));
Pjnum=size(imageFinal,3);
Dset = zeros(size(imageFinal((700-511):(701+511),(700-511):(701+511),:)));

%%
%DiffInd_hard = [14 34 37 55 57  63 64 65 66];
%DiffInd_med = [8 49 54];
ScaleFac_hard = 0.7;
ScaleFac_med = 0.85;
for Pnum=1:Pjnum
%for Pnum=1:68
    
    % Poisson-Gaussian noise parameters 
    curr_Alp = Alpha_arr(Pnum);     
    curr_Sigma = Sigma_arr(Pnum);
    curr_Mu = 0;
    
    z = Data(:,:,Pnum);
    

    
    %curr_Sigma = Sigma(Pnum);
    %curr_Alp = meanalph(Pnum);

    % Apply forward variance stabilizing transformation
    fz = GenAnscombe_forward(z,curr_Sigma, curr_Alp, curr_Mu); % Generalized Anscombe VST (J.L. Starck, F. Murtagh, and A. Bijaoui, Image  Processing  and  Data Analysis, Cambridge University Press, Cambridge, 1998)

    
    % DENOISING

    % Scale the image (BM3D processes inputs in [0,1] range)    
    scale_shift = (1-scale_range)/2;

    maxzans = max(fz(:));
    minzans = min(fz(:));
    fz = (fz-minzans)/(maxzans-minzans);
    curr_sigma_den = 0.96/(maxzans-minzans);
    fz = fz*scale_range+scale_shift;     

    curr_sigma_den = curr_sigma_den*scale_range;
%     if sum(Pnum==DiffInd_hard)>0
%         curr_sigma_den =curr_sigma_den *ScaleFac_hard;
%     elseif sum(Pnum==DiffInd_med)>0
%         curr_sigma_den =curr_sigma_den *ScaleFac_med;
%     end

    % run BM3D
    [~, D] = BM3D_optpara_20160404YY(fz,fz,curr_sigma_den*255,'np',0); 
    % Scale back to the initial VST range
    D = (D-scale_shift)/scale_range;
    D = D*(maxzans-minzans)+minzans;

    % Apply the inverse transformation
    yhat0 = GenAnscombe_inverse_exact_unbiased(D,curr_Sigma,curr_Alp,curr_Mu); 
    yhat0 = (yhat0 - curr_Mu)/curr_Alp;    
    
    SF = sum(z(:).*yhat0(:))/sum(yhat0(:).^2);
    Dset(:,:,Pnum) = fliplr(yhat0*SF);

    Pnum
    
end
%%
save DENOISED_dataset.mat Dset

