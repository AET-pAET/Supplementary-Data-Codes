% My_1Dcurvefit_comp_initguess_fixBG
% Author: Yongsoo Yang, UCLA Physics and Astronomy
%         yongsoo.ysyang@gmail.com

function [p, fminres, fitresult, eflag] = My_1Dcurvefit_comp_initguess_fixBG(refdata, testdata,start_ind,end_ind,init_guess)

    %opts = optimset('lsqcurvefit');
    opts = optimset('Display','off');
    % fitting for x axis projection
    %init_guess = [0 1 0];
    fun = @(p,xdata) real(My_FourierShift_1D_sidecut(xdata,p(1),start_ind,end_ind))*p(2);
    [p,fminres,~,eflag] = lsqcurvefit(fun,init_guess,testdata,refdata(start_ind:end_ind),[],[],opts);
    fitresult = real(My_FourierShift_1D_sidecut(testdata,p(1),start_ind,end_ind))*p(2);  
    
end