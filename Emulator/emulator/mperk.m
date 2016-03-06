function output = mperk(varargin);
% SYNTAX    output = mperk(varargin)
%
% Purposes  
%     - Estimate correlation parameters as MLE / REML  
%             for the Cubic, Gaussian, PowerExponential correlation functions
%     - Estimate the mean and variance of a Gaussian 
%                       Stochastic process
%     - Aquire the Best Linear unbiased predictions
%                      and prediction errors with the estimations
%     - Perform cross validation estimation on the training data. 
% 
% Authors         Gang Han and Tom Santner  
% Reference      : Santner, T., Williams, B., Notz, I.(2003)
%                   "The Design and Analysis of Computer Experiments"
%                   Page 64 -- Page 82 
% Organization   : The Ohio State Univ., Dept. of Statistics
% Updated        : 24 Jan 2008 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%  Input values, their possible values, and defaults 
%  
%   'X' an n by p matrix of inputs for training input values
%           -  n == number of training data cases 
%           -  p == Number of the computer code input for 
%                  each experiment
%           -  Default: none--IS A REQUIRED INPUT
%
%    'Y' an n by 1 vector of real-valued outputs of computer code 
%             at the training data sites.
%            -Default: none--IS A REQUIRED INPUT
%
%
%  'CorrelationEstimationMethod' one of the strings {'MLE', 'REML'} 
%        where
%        - 'MLE' states that the program to estimate the roughness 
%              parmeters by ML, and
%        - 'REML' states that the program to estimate the roughness 
%              parmeters by REML 
%        - Default == 'REML'
%
%  'CorrelationFamily' one of the strings {'PowerExponential', 
%                     'Gaussian', 'Cubic'} 
%        where    
%        - 'PowerExponential' means that the product power
%            exponential family is to be used with powers estimated,
%        - 'Gaussian' means that the product power
%            exponential family is to be used with powers set equal to 2
%        - 'Cubic' means that the product cubic correlation family is
%            to be used 
%        - Default == 'PowerExponential'
%  
%  'XPred'  npred by 1 matrix of new inputs at which predictions are
%         to be computed 
%         - npred is the number of the desired predictions.
%         - Default == no predictions requested
%
%  'RegressionModel' an  n by q  matrix of regression 
%         functions, F, describing the mean of the fitted GaSP.
%                   (E(GASP) = F * Beta)
%         - q is the number of regression parameters ('Beta'). 
%         - F contain the regression coefficients for the mean of
%            the GASP model   
%         - Default == ones(n,1)
%
%  'PredRegressionModel', an npred by q matrix(Fpred) of 
%         regression prediction inputs 'Xpred'.  
%         - Default == ones(npred,1)
% 
%  'CrossValidation' one of the strings {'Yes' , 'No'}
%         - 'Yes', mperk performs  cross validation and save the 
%                   results in output.cv. 
%         - 'No', mperk does NOT perform cross validation. 
%         - Defaut = 'No'
%
%  'SensitivityAnalysis', one of the strings {'Yes' , 'No'}
%         - 'Yes', indicates that the program produces sensitivity
%                   coefficients in output.sens and
%                   creates main effect plots (mainplot.eps)
%         WARNING sensitivity will only be computed for the Gaussian and
%                   Cubic Correlationfamily and the constant mean model.  
%                   If sensitivity analysis is requested for PowerExponential
%                   correlatilon families, the request is ignored. 
%         - 'No', mperk does NOT compute sensitivity coefficients 
%         -  Default = 'No'
%
%  'InitNumber' is a positive integer number \in {1,2,...,200}.
%         This value specifies the number of the starting points in the 
%         correlation parameter estimation search as follows. We first 
%         generate    (200 * p) starting points by creating a 
%         Latin Hypercube design with (200 * p) points where p
%         is the number(or dimension) of the correlation parameters. 
%         We compute the likelihood for each of these 200 * p parameter 
%         values points and then we pick the (InitNumber*p) 
%         parameter sets with largest likelihood values to run 
%         the optimization routine  fmincon.  The highest final
%         likelihood from these optimization runs yields the estimated 
%         correlation parameters used by the program.
%         (Default ==  5)
%
%  'UseCorModel' is a structure of the estimated model parameters. 
%          Using these parameters, one can predict new outputs and cross 
%          validate the training data. One could first estimate parameters using 
%	   'mperk.m' to get an output structure 'output1', then use 'mperk.m' 
%          again with 'output1' as an input to 'UseCorModel' for
%          prediction.
%
%  WARNING: mperk REQUIREs 'X' and 'Y' but all other inputs are optional.
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Output: The output of the Matlab mperk program is a MATLAB structure
%   with five parts, each part itself a MATLAB structure.  If 'output'
%   is the output, these parts are
%       a): output.job  
%       b): output.corparms
%       c): output.est
%       d): output.preds
%       e): output.cv  
%       f): output.sens
%   The description of each part is given as follows.
%
% 'output.job' is a MATLAB structure that describes the data and other job 
%        information such as the fitted model and method of 
%        fitting.  The components are :
%
%   (1) 'output.job.CorrelationFamily' describes the type of 
%          correlation function. (one of {'PowerExponential',
%          'Gaussian', 'Cubic'}) 
%                 
%   (2) 'output.job.CorrelationEstimationMethod' is the estimation
%          method. Its value is either 'MLE' or 'REML.'
%    
%   (3) 'output.job.X' is the n by p matrix of the training data
%          inputs.
%                              
%   (4) 'output.job.Y' is an n by 1 vector of the training data
%          responses.
%
%   (5) 'output.job.XPred' is the npred by p matrix of the
%          prediction inputs.
% 
%   (6) 'output.job.RegressionModel' is the matrix(n by q)
%          describing the mean the training data responses.
%
%   (7) 'output.job.PredRegressionModel' is the matrix(npred by q)
%          describing the mean of the prediction data outputs.
%
% 'output.corparms' is a MATLAB structure that gives the estimated 
%         correlation parameters.
%           
%   (1)  'output.corparms.theta' is a p by 1 vector.
%         If the correlation matrix is 'Gaussian' or
%             'PowerExponential', 'output.corparms.theta' is the
%              scale paramater.
%         If the correlation matrix is 'Cubic', 'output.corparms.theta'
%              contains all the correlationparameters
%
%   (2) 'output.corparms.power' is a p by 1 vector. 
%               If the correlation matrix is 'Gaussian' or
%                 'PowerExponential', 'output.corparms.power' is 
%                 the power.
%               If the correlation matrix is 'Cubic',
%                 'output.corparms.power' has no meaning.
%
% 'output.est'  is a MATLAB structure that gives the following 
%            estimated parameters and quantities.
% 
%   (1) 'output.est.loglik' is the maximized likelihood value with
%          the estimated correlation parameters.
%
%   (2) 'output.est.beta' is the q by 1 estimated regression
%          parameters. 
%       
%   (3) 'output.est.invR' is an n by n inverse correlation matrix
%          corresponding to the training data inputs output.job.X and
%          correlation parameters.
%
%   (4) 'output.est.sigma2' is the estimated process
%          variance(\sigma^2).
%
% 'output.preds' is a MATLAB structure that gives the predictions 
%          and the prediction uncertainties of those predictions
%
%   (1) 'output.preds.ypreds' is a npred by 1 matrix containing the
%          predictions.
%
%   (2) 'output.preds.se' is a npred by 1 matrix containing the 
%          standard error for each prediction in.
%
% 'output.cv' is a MATLAB structure that gives results related 
%           for with the leave-one-out cross validation.
%
%   (1) 'output.cv.ypreds' is an n by 1 matrix with the
%          leave one out prediction.
%
%   (2) 'output.cv.se' is an n by 1 matrix with the
%          standard errors of the predictions.
%
%   (3) 'output.cv.resids' is a n by 1 residule matrix with
%          each element equals the true response value minus the
%          leave one out predictions.
%
% 'output.sens' is a MATLAB structure that gives the estimated sensitivity
%           coefficients. 
%
%   (1) 'output.sens.sme' is a 1 by p vector having estimated main effect 
%           sensitivity indices.
%
%   (2) 'output.sens.ste' is a 1 by p vector having estimated total effect 
%           sensitivity indices.
%
%   (3) 'output.totalvar' is an estimated total variance.
%
%   (4) 'output.ngrid' is the number of grids of input region 
%           for evaluating main effect values.
%
%   (5) 'output.maineffect' is an ngrid by p matrix to contain 
%           the main effect values evaluated at the grids of the inputs 
%           which are used to create main effect plots.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% MATLAB programs used by mperk.m
%  Programs used for estimating the correlation parameters
%  
%  (1) 'estgaussiancor.m' : For Gaussian correlation function, it estimates
%  the scale parameters with power parameters fixed at 2.
%
%  (2) estexponentcor.m   : estimates the correlaton parameters, both the
%  scale parameter and the power parameters.
%
%  (3) 'estcubiccor.m'    : estimate the cubic correlation parameters.
%
%  (4) 'calicatMLE.m'     : calculates the likelihood.
%     
%  (5) 'calicatREML.m'    : calculate the restricted likelihood.
%     
%  (6) 'llik1d.m'         : calculate the likelihood.
%          
%  (7) 'cormatexp.m'      : returns the correlation function for Gaussian
%  or Power Exponential correlation function.
%  
%  (8) 'cubiccorrfn.m'    : returns the correlation for a pair of inputs
%  for cubic correlation function.
%
%  The files call each other as follows for estimating the correlation
%  parameters:
%
%     |-------->(1)('Gaussian')---|
%     |                           |---->(4)('MLE')----|--->(7)('Gaussian'
%  mperk.m--- ->(2)('PowerExp')---|                   |       /'PowerExp')
%     |                           |---->(5)('REML')---|--->(8)('Cubic')
%     |-------->(3)( 'Cubic'  )---|
%
%  
% Programs used for Prediction and Cross Validation
%    
%  (1) 'cormatexp.m'            
%  
%  (2) 'cubiccorrfn.m'          
%
%  (3) 'correxp.m' : returns the correlation matrix between the
%  predictions and the training data.
%  
% Programs used for Computation of the Sensitivity Analysis
%    
%  (1) 'sensitivity.m'   : returns the sensitivity coefficients and main effect plots.
%
%  (2) 'sgint.m'         : calculates an integration of correlation function.
%
%  (3) 'dbint.m'         : calculates a double integration of correlation function.
%
%  (4) 'mxint.m'         : calculates an integration of the product of
%  correlation functions.
%
%  (5) 'Cfunc.m'         : calculates a matrix C by calling the functions
%  sgint.m and mxint.m.
%
%Disclaimer: This software is provided as is without express or implied 
%            warranty that the output is errorfree or reliable.
%            The author accepts no responsibility  or liability for any
%            loss or damage occasioned by its use. 
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (nargin < 2)
  error('Too few input arguments');
end


% Set default inputs (to be modified as needed)
% 
output.job.CorrelationFamily            = 'PowerExponential';
output.job.CorrelationEstimationMethod  = 'REML';
input.xpred                             = Inf;
input.corparms.scale                    = Inf; 

input.corparms.smoothness               = Inf; 
%  default == estimate roughness parameters

input.corparms.ctype                    = 1; 
% default == estimate both smoothness and roughness parameters

input.corparms.fittype                  = 1; 
%  default == roughness parameters are to be estimated by REML

input.F                                 = Inf;   
% default == vector of ones

input.Fpred                             = Inf; 
% default == vector of ones

use_default_F                           =   1;
use_default_Fpred                       =   1;
use_default_Xpred                       =   1;
use_default_cross_valid                 =   0;
use_default_initnum                     =   5;
use_default_output                      =   0;
use_default_sens_coef                   =   0;

%%%% reset input parameters according to user input
paramPairs  = varargin(1:end);
for k       = 1:2:length(paramPairs)
  param     = lower(paramPairs{k});
  
%%%  disp(param)
  if (~ischar(param))
    error('Optional parameter names must be strings');
  end
  value     = paramPairs{k+1};

%%%  disp(value) 
  switch (param)
  case 'x'
    input.x      = value;
    output.job.X = value; 
    
  case 'y'
     input.y        = value;
     output.job.Y   = value;  
     
  case 'xpred'
     input.xpred        = value;
     output.job.XPred   = value;
     use_default_Xpred  = 0;
     
  case 'correlationestimationmethod'
      if (~strcmp(value,{'MLE','REML'}))
      error(['CorrelationEstimationMethod must be ''MLE'' or ''REML''.']);
      end
      if (strcmp(value,'MLE'))    
            input.corparms.fittype  = 0;
      else 
            input.corparms.fittype  = 1;
      end
      output.job.CorrelationEstimationMethod    = value;
      
  case 'correlationfamily'
      if (~strcmp(value,{'Gaussian','PowerExponential','Cubic'}))
          error...
          (['CorrelationFamily must be ''Gaussian'',...''PowerExponential'', or ''Cubic''.']);
      end
      if (strcmp(value,'Gaussian'))    
         input.corparms.ctype   = 0;
      elseif (strcmp(value,'PowerExponential'))   
         input.corparms.ctype   = 1;
      elseif (strcmp(value,'Cubic')) 
         input.corparms.ctype   = 2;
      end
      output.job.CorrelationFamily  = value;
      
  case 'regressionmodel'
      input.F                       = value;
      output.job.RegressionModel    = value;
      use_default_F                 = 0;
      
  case 'predregressionmodel'
      input.Fpred                    = value;
      output.job.PredRegressionModel = value;
      use_default_Fpred              = 0;
      
  case 'crossvalidation'
      if (strcmp(value,'Yes'))    
         use_default_cross_valid    = 1;
      elseif (strcmp(value,'No'))   
         use_default_cross_valid    = 0;
      else 
         error(['CrossValidation must be ''Yes'' or ''No''.']);
      end
      
  case 'sensitivityanalysis'
      if (strcmp(value,'Yes'))
        use_default_sens_coef     = 1;
      elseif (strcmp(value,'No'))
        use_default_sens_coef    = 0;
      else
        error(['SensitivityAnalysis must be ''Yes'' or ''No''.']);
      end
      
  case 'initnumber'
      use_default_initnum           = value;
      
  case 'usecormodel'
      output                        = value;
      use_default_output            = 1;
      
  otherwise
      error(['Unrecognized option ' param '.']);
  end % end switch
end % end k loop on data setup


%%%%%  Construct data structure for other functions
x       =   input.x;
xpred   =   input.xpred;
[n,p]   =   size(input.x);

% Set minimum and maximum values of input.x matrix
input.ranges            = zeros(size(input.x,2),2);
for i   =   1:p,        
    input.ranges(i,1)   = min(input.x(:,i)); % 0; 
    input.ranges(i,2)   = max(input.x(:,i)); % 1; 
end

% Transform x and xpred to the scale [0,1]
% and conduct the search in that range
for i = 1:p,
    x(:,i) = (x(:,i)-input.ranges(i,1))/(input.ranges(i,2)-input.ranges(i,1));
end
output.job.X    = input.x;
input.x = x;


if (use_default_F       == 1)
  input.F                       = Inf;
  output.job.RegressionModel    =   ones(size(input.x,1),1); 
end 

if (use_default_Fpred   == 1) 
  input.Fpred                   = Inf;   
  output.job.PredRegressionModel= ones(size(input.xpred,1),1);  
end

if (use_default_Xpred   == 1) 
  input.xpred                   = Inf;
  output.job.XPred              = 'None';
end    
n       =   size(input.x,1);
[p,k]   =   size(input.ranges);

% Make the design matrix F  for the training data 
if input.F  ==  Inf,
  input.F =   ones(n,1);
end

%%%%% Main program

if use_default_output == 0,
% Check to see if we want to estimate the correlation parameters
    %% We need to estimate the ...correlation parameters.

    if input.corparms.ctype ==  0,         
        %% Gaussian correlation: We do not need to estimate the smoothness parameter
        [corparms, loglik]  = estgaussiancor...
            (input.x,input.y,input.F,input.corparms,use_default_initnum); 
        %% estimate the Gaussian correlation function
        output.corparms.theta  = ones(p,1);
        output.corparms.power  = ones(p,1)*2; 
        for i=1:p,
            output.corparms.theta(i)  = corparms.scale(i);
        end  % if      
    elseif input.corparms.ctype ==1, 
        %% Power Exponential correlation: estimate the smoothness parameters
        [corparms1, loglik1]    = estgaussiancor...
            (input.x,input.y,input.F,input.corparms,use_default_initnum); 
        %% also search the boundary
            my_fittype = input.corparms.fittype;
        [corparms2, loglik2]    = estexponentcor...
            (input.x,input.y,input.F,my_fittype,use_default_initnum); 
        %% estimate the power exp correlation function
        %%% Note: The program favors the infinitely differentiable 
        %%% sample path. we compare loglik1-1 and loglik2 and choose 
        %%% the smaller one.
        if (loglik1-1)  > loglik2,  
            corparms    =   corparms2;
            loglik      =   loglik2;
        else
            corparms    =   corparms1;
            loglik      =   loglik1;
        end
        output.corparms.theta  = 1 * ones(p,1);
        output.corparms.power  = 2 * ones(p,1);  
        for i=1:p,
            output.corparms.theta(i)  = corparms.scale(i);
            output.corparms.power(i)  = corparms.smoothness(i);        
        end        
    elseif input.corparms.ctype ==  2,
        %%% cubic function 
        [corparms, loglik]      = estcubiccor(input.x,input.y,...
            input.F,input.corparms,use_default_initnum);
        output.corparms.theta   = 1 * ones(p,1);
        for i=1:p,
            output.corparms.theta(i) = corparms.scale(i);
        end    
        output.corparms.power  = Inf;
    end
else %%% For this case, the input has the correlation parameters.
    if input.corparms.ctype ==  2,
        cor     = output.corparms.theta;
        loglik  = llik1d(cor, output.job.Y, output.job.X,...
            input.corparms.ctype, input.corparms.fittype,...
            output.job.RegressionModel, 2*length(output.corparms.theta));
    else 
        cor     = [output.corparms.theta;output.corparms.power];
        loglik  = llik1d(cor, output.job.Y, output.job.X,...
            input.corparms.ctype, input.corparms.fittype,...
            output.job.RegressionModel, 2*length(output.corparms.theta));
    end
    %%% scale 'output.corparms.theta'
    if input.corparms.ctype ~=  2, 
    		for i=1:p,
        		output.corparms.theta(i)  = output.corparms.theta(i)*(...
            	(input.ranges(i,2)-input.ranges(i,1))^output.corparms.power(i));
    		end
	else 
    		for i = 1:p, 
        		output.corparms.theta(i)  = output.corparms.theta(i) / ...
            	(input.ranges(i,2)-input.ranges(i,1));
    		end
    end
end
output.est.loglik=-loglik/2;

%%%%%%  Compute covariance matrix %%%%%
if input.corparms.ctype == 2,
 %%% i.e., use cubic correlation function %%%
    [n,p]   =   size(input.x);
    RR      =   zeros(n,n);
    for i   =   1:n
      for j     =   1:n
         RR(i,j)   =   cubiccorrfn(input.x(i,:)-...
                input.x(j,:),output.corparms.theta);
      end
    end 
% Use the commented part if one wants to fit a model with measurement error.
%   RR    = RR; % * 0.99 + 0.01*eye(n,n); 
else
    R       = cormatexp(input.x,output.corparms.theta,...
                output.corparms.power);
    RR      = (R + R')/2;
% Use the commented part if one wants to fit a model with measurement error.
%   RR      = RR;%*0.99 + 0.01*eye(n,n);
end


%%% Calculate Betahat %%%
RinvF           =   RR\input.F;
tFRinvF         =   input.F'*RinvF;
RinvY           =   RR\input.y;
tFRinvFinvFt    =   tFRinvF\input.F';
Betahat         =   tFRinvFinvFt*RinvY;
output.est.beta =   Betahat;
output.est.inv_R=   inv(RR);


%%% Calculate sigma_squre_hat %%%
muhat               =  input.F*Betahat;
%%% MLE estimate of sigma2 the denominator is n %%%
sigma2hat           =   (1/n)*((input.y-muhat)'*output.est.inv_R*(input.y-muhat));
output.est.sigma2   =  sigma2hat;


%%% calculate t - statistics %%%
resid = input.y - input.F*Betahat;
Rinvresid = RR\resid;
if (input.corparms.fittype == 0)
    %%% MLE %%%
    T1hat2 = (resid'*Rinvresid)/n;
elseif (input.corparms.fittype == 1)
    %%% REML %%%
    T1hat2 = (resid'*Rinvresid)/(n - size(input.F,2));
end


%%% Predict outputs from computer experiment %%%
%%% Check to see if we want to do prediction %%%
if input.xpred  ~=Inf, %% In case we need to do prediction
    [npred,p]   =size(input.xpred);
    output.job.XPred= input.xpred;
    for i=1:p,
        xpred(:,i) = (xpred(:,i)-input.ranges(i,1))/(input.ranges(i,2)-input.ranges(i,1));
    end
    input.xpred = xpred;
    
    %%% check the constant mean condition%%%
    if input.Fpred  == Inf, 
        input.Fpred = ones(npred,1);
        npar        = 1;
    else
        npar        = size(input.F,2);
    end

    %%% Make the correlation matrices for Xpred
    if input.corparms.ctype ==  2,
        %%% calculate correlation function %%
        %%% for prediction inputs
        %%% !!!! We use the standardized data     !!!! %%%
        [npred,p]   =   size(input.xpred);
        V11         =   zeros(npred,npred);
        for i   =   1:npred
            for j        =   1:npred
                V11(i,j)   =   cubiccorrfn(input.xpred(i,:)-...
                    input.xpred(j,:),output.corparms.theta);
            end
        end
        V11     = V11; %*0.99 + 0.01*eye(npred,npred);
        %%% One can use the commented part if one wants to fit a model with
        %%% measurement error.
        
        %%% for training data and predictions
        [n,p]   =   size(output.job.X);
        V12     =   zeros(npred,n);
        for i   =   1:npred
          for j     =   1:n
            V12(i,j)  =  cubiccorrfn(input.xpred(i,:)-...
                         input.x(j,:),output.corparms.theta);
          end
        end            
    else     
        V11     = cormatexp(input.xpred, output.corparms.theta, ...
                    output.corparms.power);
        V11     = V11; %*0.99 + 0.01*eye(npred,npred);
        %%% One can use the commented part if one wants to fit a model with
        %%% measurement error.
        V12     = correxp(input.xpred, input.x, output.corparms.theta, ...
                    output.corparms.power);
    end
    % Compute preds and their SE's
    preds       = input.Fpred*Betahat + V12*Rinvresid;
    RinvV12t    = RR\V12';
    mse         = V11-V12*RinvV12t+(input.Fpred-V12*RinvF)*inv(tFRinvF)...
                  *(input.Fpred-V12*RinvF)';
   
    for i = 1:npred,
        if (mse(i,i) < eps)
            se(i,1) =0.0;
        elseif input.corparms.fittype == 0,
            se(i,1) = sqrt(T1hat2*mse(i,i));            
        elseif input.corparms.fittype ==1,
            se(i,1) = sqrt(((n-npar)/(n-npar-2))*T1hat2*mse(i,i));
        end % end if
    end
% replaced on 22May 09 by TJS    
%     for i = 1:npred,
%         if input.corparms.fittype == 0,
%             se(i,1) = sqrt(T1hat2*mse(i,i));            
%         elseif input.corparms.fittype ==1,
%             se(i,1) = sqrt(((n-npar)/(n-npar-2))*T1hat2*mse(i,i));
%         end    
%     end
    output.preds.ypreds  =  preds;
    output.preds.se      =  se;
end


%%% Cross - Validation for computer experiment %%%
%%% Predict for y(x_i) with the correlation parameters, save output %%%  
if use_default_cross_valid == 1,
    %%% check the constant mean condition%%%
    if input.Fpred==Inf, 
        npar=1;
    else
        npar=size(input.F,2);
    end   
    q               =   size(input.F,2);
    valid.x         =   zeros(n-1,p);
    valid.y         =   zeros(n-1,1);
    valid.F         =   zeros(n-1,q);
    valid.Fpred     =   zeros(1,q);
    valid.preds     =   zeros(n,1);
    valid.predsse   =   zeros(n,1);
    valid.resid     =   zeros(n,1);
    
    %%% Creat the input and output for cross validation %%%
    for i_vali   =    1 : n,
    
        %%% For x,y, and F
        if i_vali   ==  1,
            for j   =   1:(n-1),
                valid.x(j,:)    =    input.x(j+1,:);
                valid.y(j)      =    input.y(j+1);
                valid.F(j,:)    =    input.F(j+1,:);
            end
        else
            for j   =   1:(i_vali-1),
                valid.x(j,:)    =    input.x(j,:);
                valid.y(j)      =    input.y(j);
                valid.F(j,:)    =    input.F(j,:);
            end
    
            if  j   <   n-1,
                for j   =   (i_vali+1):n,
                    valid.x(j-1,:)  =   input.x(j,:);
                    valid.y(j-1)    =   input.y(j);
                    valid.F(j-1,:)  =   input.F(j-1,:);
                end
            end
        end
    
        %%% For the prediction inputs
        valid.Fpred  =    input.F(i_vali,:);
        valid.predin =    input.x(i_vali,:);
        valid.true   =    input.y(i_vali);    
    
        %%% Make the correlation matrices for Xpred %%%
        if input.corparms.ctype==2
            %%% correlation matrix for prediction inputs
            npred     =   1;
            V11       =   zeros(npred,npred);
            for i   =   1:npred,
                for j        =   1:npred,
                    V11(i,j)   =   cubiccorrfn(valid.predin(i,:)-...
                        valid.predin(j,:),output.corparms.theta);
                end
            end
            V11 = V11;%*0.99 + 0.01*eye(npred,npred);
            %%% One can use the commented part if one wants to fit a model
            %%% with measurement error.
            
            %%% for training data and predictions
            V12     =   zeros(npred,n-1);
            for i   =   1:npred,
                for j     =   1:(n-1)
                    V12(i,j)   =   cubiccorrfn(valid.predin(i,:)-...
                    valid.x(j,:),output.corparms.theta);
                end
            end  
        
            %%% for training data and predictions        
            RR     =   zeros(n-1,n-1);
            for i   =   1:(n-1),
                for j     =   1:(n-1)
                     RR(i,j)   =   cubiccorrfn(valid.x(i,:)-...
                        valid.x(j,:),output.corparms.theta);
                end
            end          
   
        else        
            %%% correlation matrix for prediction inputs
            V11 = cormatexp(valid.predin, output.corparms.theta, ...
                output.corparms.power);
            V11 = V11;%*0.99 + 0.01*eye(npred,npred);
            %%% for training data and predictions
            V12 = correxp(valid.predin, valid.x, output.corparms.theta,...
               output.corparms.power);
            %%% for training data and predictions        
            RR  =   zeros(n-1,n-1);
            output.corparms.theta;
            output.corparms.power;
            valid.x;
            RR  =   cormatexp(valid.x, output.corparms.theta,...
                output.corparms.power);   
        end
        RR;
        %%% Compute some quantities for the cross validation
        %%% Calculate Betahat %%%
        valid.RinvF         =   RR\valid.F;
        valid.tFRinvF       =   valid.F'*valid.RinvF;
        valid.RinvY         =   RR\valid.y;
        valid.tFRinvFinvFt  =   valid.tFRinvF\valid.F';
        valid.Betahat       =   valid.tFRinvFinvFt*valid.RinvY;
        %%% Calculate sigma_squre_hat %%%
        valid.muhat         =   valid.F*Betahat;
               %%% MLE estimate of sigma2 the denominator is n %%%
        valid.sigma2hat     =   (1/(n-1)) * ((valid.y-valid.muhat)'...
                        * inv(RR) * (valid.y-valid.muhat));
        %%% calculate t - statistics %%%
        valid.resid         =   valid.y - valid.F * valid.Betahat;
        valid.Rinvresid     =   RR\valid.resid;
        if (input.corparms.fittype == 0)
                  %%% MLE %%%
            valid.T1hat2    =   (valid.resid' * valid.Rinvresid)/(n-1);
    	elseif (input.corparms.fittype == 1)
                  %%% REML %%%
            valid.T1hat2    =   (valid.resid' * valid.Rinvresid)/((n-1)...
                                - size(valid.F,2));
        end
        
        %%% Compute preds and their SE's
        valid.predout       =   valid.Fpred * valid.Betahat + ...
                                V12 * valid.Rinvresid;
        valid.RinvV12t      =   RR\V12';
        valid.mse           =   V11 - V12 * valid.RinvV12t + ...
                                (valid.Fpred - V12 * valid.RinvF)...
                * inv(valid.tFRinvF) * (valid.Fpred - V12 * valid.RinvF)';
        if input.corparms.fittype == 0,
            valid.se = sqrt(valid.T1hat2 * valid.mse);            
        elseif input.corparms.fittype==1,
            valid.se = sqrt((((n-1)-npar)/((n-1)-npar-2))...
                * valid.T1hat2 * valid.mse);
        end    
        valid.preds(i_vali)   =   valid.predout;
        valid.predsse(i_vali) =   valid.se;    
    end
    
%%% Save the result of the validation
output.cv.ypreds        =   valid.preds;
output.cv.se            =   valid.predsse;
    %%% Residule %%%
output.cv.resids        =   input.y-valid.preds;
end


%%% Do sensitivity analysis %%%%
%%% Check to see if we want to do sensitivity analysis %%%%
if use_default_sens_coef == 1,
    
    %%% check the constant mean model condition
    if size(input.F,2) == 1 & sum(input.F==1)==n
        constant = 1;
    else
        constant = 0;
    end
    
    %%% do sensitivity analysis only for the Gaussian correlation family and constant mean model
    if (input.corparms.ctype ~=1 & constant == 1),  
        output.sens = sensitivity(input.x, input.y, output.est.beta, output.est.sigma2,...
                                 (output.est.inv_R/output.est.sigma2), output.corparms.theta, input.corparms.ctype);    
    else
        output.sens = 'MPERK only does sensitivity analysis for Gaussian or Cubic correlation family and constant mean.';
    end
end


%%% Transform the estimations of the correlation parameters to the scales 
%%% according to the support of the data. (The estimations were for the 
%%% support [0,1].)
if input.corparms.ctype ~=  2, 
	for i=1:p,
        output.corparms.theta(i)  = output.corparms.theta(i)/(...
        (input.ranges(i,2)-input.ranges(i,1))^output.corparms.power(i));
    end
else 
    for i = 1:p, 
        output.corparms.theta(i)  = output.corparms.theta(i) * ...
        (input.ranges(i,2)-input.ranges(i,1));
    end
end


%%% Display outputs
fprintf...
    ('Model settings (saved as output.job)\n');
disp(output.job);
fprintf...
    ('\nEstimated correlation parameters (saved as output.corparms)\n');
disp(output.corparms);
fprintf...
    ('\nEstimation of mean and variance (saved as output.est)\n');
disp(output.est);
if use_default_Xpred == 0,
	fprintf...
        ('\nPrediction and estimated standard errors (saved as output.preds)\n');
	disp(output.preds);
end
if use_default_cross_valid == 1,
   fprintf...
     ('\nCross Validation predictions, standard errors, and residules\n');
   fprintf...
       ('(saved as output.cv)\n');
   disp(output.cv);
end
if (use_default_sens_coef == 1 & input.corparms.ctype ~= 1 & constant == 1),
   fprintf...
     ('\nSensitivity indices, total variance, and main effect values over grids\n');
   fprintf...
       ('(saved as output.sens)\n');
   disp(output.sens);
end

function value = LocalToNum(value);

if ischar(value),
  value = str2num(value);
end


