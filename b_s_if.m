%% b_s
% Collection of methods for binocular summation analysis
%
% Used for analysis in Meier, Tarczy-Hornoch, Boynton, & Fine (under
% revision), first version up at bioRxiv preprint doi:
% https://doi.org/10.1101/2022.10.10.511635)

% Description of the four stages for the  'softmax' model:

% Stage 1: Joystick vs. contrast calibration
%
% J-hat(t) = a + b • J(t-d)
% J is the joystick lever position (0-1), d reflects a time delay between
% stimulus presentation and observer response, and a and b represent a
% linear scaling between joystick location and calibrated position J-hat.

% Stage 2: Monocular attenuation
%
% C-hat(t) = k_right•C_right(t) + k_left•C_left(t)
%            divide by max (k_r, k_l) to normalize
%            max k (i.e. k=1) is "fellow" eye FE, k < 1 is "amblyopic" AE
% C-hat(t) is the model prediction of perceived contrast, C_right and
% C_left are the contrasts presented to left and right eyes at time t
%
% uses parameters p.k(1) and p.k(2)
%
% out(:,i) = p.k(i)*in(:,i)
%
% set p.k = [1,1] by default for no effect

% Stage 3: Binocular interactions (normalization)
%
% equation for each eye:
% C-hat_AE(t) = k_AE•C_AE(t) / mu_AE•C_FE(t)+sigma
% C-hat_FE(t) = C_FE(t) / mu_FE•k_AE•C_AE(t)+sigma
%
% uses parameters p.sigma, p.U(1), p.U(2), p.U(3) and p.U(4).  Note, U used
% to be a 2x2 Mat but 'fitcon' had trouble with this.  So U is now a 1x4
% vector.
%
%  out(:,1) = in(:,1)./(p.U(1)*in(:,1) + p.U(2)*in(:,2)+p.sigma)   % LE
%  out(:,2) = in(:,2)./(p.U(3)*in(:,1) + p.U(4)*in(:,2)+p.sigma)   % RE
%
% Set U's to 0 and gamma to 1 by default for no effect of normalization.

% stage 4, softmax
%
% uses p.smax
% out = sum(out.^p.smax,2)'.^(1/p.smax)
% In paper, we assume the final perceived contrast is simply the sum of the
% outputs of the left and right eyes (no free parameters):
% C-hat = C-hat_AE + C-hat_FE


%% Inputs
% Expected input structure, functions below use the following fields
%
% data
%      .binocular
%               .t  (vector) timestamp for binocular trial phase
%      .experiment
%                .binoResponse  time x trials matrix of response (0-1)
%                               during binocular trial phases
%                .binoS         time x trials matrix of contrast (0-1)
%                               presented during binocular trial phases
%                .response      time x trials matrix of response (0-1)
%                               presented during dichoptic/monoptic phases
%                .LEcontrast    time x trials matrix, left eye contrast (0-1)
%                               presented during dichoptic/monoptic phases
%                .REcontrast    time x trials matrix, right eye contrast (0-1)
%                               presented during dichoptic/monoptic phases
%      .sID   string, subject ID
%      .t     vector timestamp for dichoptic/monotopic phase corresponding
%             with time matrix above (e.g. 0:1/samplingrate:totalduration)
%
% p
%  paramters used in calibration (step 1 above):
%  .joystickfunction
%  .clean_range
%  .delay       estimated delay (in sec) btwn presented contrast & response
%  .intercept   estimate of intercept (a) for calibration
%  .slope       estimate of slope (b) for calibration
%
%  parameters used in steps 2 & 3:
%  .costflag
%  .k
%  .U
%  .sigma
%
%  other stuff that gets used:
%  .dt              ∆t (eg diff(t(1:2))
%  .startT
%  .model
% There are additional flags here for fitting alternate models - see
% relevant sections below


%%
classdef b_s
    methods(Static)

        %% Joystick calibration using binocular trial data
        function [err, err_noCost, data] = getErrBinoMean(p,data)
            % Calculates MSE between binocularly-presented data (two eyes
            % get the same input) and model prediction, based on the mean
            % joystick response across multiple trials
            % relevant inputs:
            % data.binoResponse is a time x trials matrix of response
            % data.binoS is presented stimulus
            % relevant outputs:
            % err is error including any cost penalizations (eg max delay)
            % err_nocost is simple error between data and prediction
            data.experiment.binoMeanResp = nanmean(data.experiment.binoResponse, 1); % both the stimulus and the original response need to be scaled from 0-1 to -0.5-0.5.
            b_calib = b_s.joyFunction(p,data.experiment.binoMeanResp); % apply calibration function
            ind = ~isnan(b_calib) .* ~isnan(data.experiment.binoS);
            ind(1:round(p.startT/p.dt)) = 0; % time pts that don't get included in the err calc
            data.experiment.binogvals = find(ind);

            % penalize negative delays or delays longer than penalizeDelay
            cost = max(p.delay - p.penalizeDelay, 0)^4 + abs(min([p.slope, 0]))^10;
            err_noCost = (sum((b_calib(find(ind))-data.experiment.binoS(find(ind))).^2))/sum(ind);
            err = cost + sum((b_calib(find(ind))-data.experiment.binoS(find(ind))).^2); % add SSE to err
            err = err + b_s.belowZeroCost(p.delay); % penalize negative delays
            data.experiment.binoMeanRespCalib = b_calib;
        end

        function [err,err_noCost, data] = getErrBinoInd(p,data)
            % Calculates MSE between binocularly-presented data (two eyes
            % get the same input) and model prediction, similar to
            % gerErrBinoMean but uses individual trials

            err = 0;  err_noCost = 0; n = 0;
            for i = 1:size(data.experiment.binoResponse,1)
                b_calib = b_s.joyFunction(p,data.experiment.binoResponse(i,:));
                ind = ~isnan(b_calib) .* ~isnan(data.experiment.binoS) ;
                cost = max(p.delay - p.penalizeDelay, 0)^4 + abs(min([p.slope, 0]))^10;
                err = cost + err + sum((b_calib(find(ind))-data.experiment.binoS(find(ind))).^2); % add SSE to err
                err_noCost =  err_noCost + sum((b_calib(find(ind))-data.experiment.binoS(find(ind))).^2); % add SSE to err
                n = n+sum(ind);
                data.experiment.binoIndRespCalib{i} = b_calib;
                data.experiment.binoNgood(i)  = 1;
                if n ~= 0
                    err_noCost = err_noCost/n;
                end
            end
        end
        function [bR_Calib, bR_Delay] = joyFunction(p,bR)
            % binocular response calibrated is a function of the binocular
            % response
            % takes in:
            %   p: containing what sort of joystick calibration
            %   p.joystickfunction can be 'delay' and/or scale {'delay', 'scale'}
            %   bR: binocular response
            % returns:
            %       bR_Calib,  the full calibration
            %       bR_Delay: delay only
            bR_Delay = bR;
            if contains(p.joystickfunction,'delay')
                tmp = round(p.delay/p.dt); tmp=max(tmp, 1);
                tc = NaN(size(bR));
                tc(1:end-tmp+1) = bR(tmp:end);
                bR_Delay= tc; % shifting the response backwards in time to match the stim
            end
            bR_Calib = bR_Delay;     % now do a linear calibration
            if contains(p.joystickfunction,'scale')
                bR_Calib = p.intercept+(p.slope*bR_Delay);
            end
        end

        function [data_out, n_good]  = cleanData(data_in, p)
            % p.clean_range = (0-1); only calibrate or fit data where
            % there's this much range in the data (on a given trial)
            % data_out has 'clean' data only
            % n_good is number of trials kept
            data_out = data_in;
            for r  = 1:size(data_in, 1)
                if (max(data_in(r,:))-min(data_in(r, :)))<=p.clean_range  || length(find(isnan(data_in(r, :))))==size(data_in, 2)
                    data_in(r, :) = NaN;
                    n_range(r) = 0;  % number of runs where calibration used an adequate range
                else
                    n_range(r) = 1;  % number of runs where calibration used an adequate range
                end
            end
            for t = 1:size(data_in, 2)
                TF = isoutlier(data_in(:, t), 'mean');
                data_out(TF, t) = NaN;
            end
            n_good = sum(n_range);
        end

        function plotJoystickCalibration(data,outputdir)
            % plots the joystick calibration (no additional model fits)
            sIDtxt = strrep(data.sID,'_',' '); % for titles
            % Fit the joystick data to obtain calibration data

            % Mean response and prediction over time stimulus time-course.
            h(1) = plot(data.binocular.t,data.experiment.binoS, 'k--', 'LineWidth',1);   hold on% plot presented contrast
            h(2) = plot(data.binocular.t,data.experiment.binoMeanResp,'r--', 'LineWidth', 1);  % plot  joystick response with delay
            h(3) = plot(data.binocular.t,data.experiment.binoMeanRespCalib,'r-', 'LineWidth', 1);  % plot calibrated joystick response

            xlabel('Time (s)'); ylabel('Contrast');
            legend(h, {'presented', 'response', 'calibration'});  grid

            if exist('outputdir','var')
                saveas(gcf, [outputdir filesep data.sID filesep ...
                    data.sID '-Joystick-Calibration-respFunc.fig']);
            end
        end

        %% Model fits

        function [err,predModel,respJoy,respJoyCalib,t,stim,n] = getErr(p, data, varargin)
            if ~isfield(p, 'mid_range_flag')
                p.mid_range_flag = 0; % this is a weird condition where we constrain the fits to places where the signal isn't too different in the two eyes
            end
            if ~isfield(p, 'abs'), p.abs = 0; end
            if nargin == 2 % all data
                subset = 1:size(data.experiment.response,1);
            else % subset
                subset = varargin{1};
            end
            n = 0;
            err = 0;  t = {}; stim = {};

            for i = 1:length(subset)% for each run
                clear tmp_stim;
                t{i} = data.t;
                stim{i}(:, 1)=data.experiment.LEcontrast(subset(i), :); % left eye
                stim{i}(:, 2)=data.experiment.REcontrast(subset(i), :); % right eye
                if p.mid_range_flag == 1
                    ind = zeros(size(stim{i}(:, 1)));
                    ind((stim{i}(:, 1)-stim{i}(:, 2))<0.5) = 1;
                    ii = find(ind);
                else
                    ii = 1:length(stim{i});
                end
                tmp_stim(:, 1)= stim{i}(ii, 1);
                tmp_stim(:, 2)= stim{i}(ii, 2);
                stim{i} = tmp_stim;

                % Get model prediction for this data set
                eval(['[predModel{i},S{i}] = ', p.model, '(p, stim{i});']);
                respJoy{i} = data.experiment.response(subset(i),ii); % data.experiment.response(i,:) has the response for that particular run
                respJoyCalib{i}= b_s.joyFunction(p,respJoy{i}); % pass through calibration function
                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                
                tmp = sum((predModel{i}(id)-respJoyCalib{i}(id)).^2);
                err = err + tmp;
                n = n+sum(id); % count the number of data points, don't divide by this if costflag so the two types of error don't get confused
            end

            % cost for U<0, sigma<0 and k<0
            if strcmp(p.model, 'b_s.softmax') && p.costflag == 1

                if p.abs == 1
                    err = err +  b_s.aboveLimCost(p.U(2)) + ...
                        b_s.aboveLimCost(p.U(3));

                elseif p.abs == 0 % we now use abs(), keeping for historical reasons
                    err = err + b_s.belowZeroCost(p.U(2)) + ...
                        b_s.belowZeroCost(p.U(3)) + ...
                        b_s.belowZeroCost(p.sigma) + ...
                        b_s.belowZeroCost(p.k(1)) + ...
                        b_s.belowZeroCost(p.k(2)) + ...
                        b_s.aboveLimCost(p.U(2)) + ...
                        b_s.aboveLimCost(p.U(3));
                end
            end
            if p.costflag == 0
                err = err/n;
            end % don't divide by n if there's a cost flag. This means the fitting units are different between the two conditions and makes it harder to confuse them
        end

        function out = belowZeroCost(in)
            out = 10000*min(in,0)^2;
        end
        function out = aboveLimCost(in)
            out = 10000*max(in,1000)^2;
        end

        function plotContrastVsResponse(contrastVector, responseVector)
            plot([0,1],[0,1],':', 'LineWidth', 2, 'Color', [.6 .6 .6]); % Plot equity line
            hold on
            plot(contrastVector, responseVector,'b-', 'LineWidth', 2); % Plot response function
            xlabel('presented contrast');
            ylabel('predicted response');

            axis equal
            axis tight
            set(gca,'XLim',[0,1]);
            set(gca,'YLim',[0,1]);
            set(gca,'XTick',0:.2:1);
            set(gca,'YTick',0:.2:1);
            grid
        end

        %% helper functions

        function [out, p] = softmax(p, S)
            % Run through the four stages of the 'softmax' model
            S = S + 0.5; % move to 0-1 units
            [out, p] = b_s.fixed_nonlinearity(p, S);
            [out, p] = b_s.linear_attenuation(p, out);
            if isnan(p.tau)
                [out, p] = b_s.normalization(p,out);
            else
                [out, p] = b_s.tau_normalization(p,out);
            end

            if ~isfield(p, 'smax'); p.smax = 1; end
            % raise to smax then sum (then raise to 1/smax)
            % linear summation if smax=1 and acts like max rule if smax >> 5 or so.
            % If smax=2, it's similar to Shroedinger's model (See Wikipedia page
            % on 'Binocular Summation')

            out = sum(out.^p.smax,2).^(1/p.smax);
            out = out(:)';  % turn output into a row vector
            out = out -0.5; % go back into -0.5 - 0.5 units
        end
        function y=gamma(n,theta,t)
            % GAMMA
            %	y=Gamma(n,theta,t)
            %	returns a gamma function from [0:t];
            %	y=(t/theta).^(n-1).*exp(-t/theta)/(theta*factorial(n-1));
            %	which is the result of an n stage leaky integrator.
            %
            %	6/27/95 gmb

            flag=0;
            if t(1)==0
                t=t(2:length(t));
                flag=1;
            end
            id=find(t<=0);
            t(id)=ones(size(id));
            y = (  (theta'*(1./t)).^(1-n).*exp(-(1./(theta'*(1./t)))))./(theta'*ones(size(t))*factorial(n-1));
            y(id)=zeros(size(id));
            if flag==1
                y=[0;y']';
            end
        end
        function val = mse(x,y)
            val = nanmean((x-y).^2)
        end

        %% fixed nonlinearity, Meier et al model
        % set p.m = [1 1] for no effect [left right]
        % not implemented in current modeling
        function [out, p] = fixed_nonlinearity(p, in)
            if ~isfield(p, 'm'); p.m = [1 1]; end
            for i = 1:length(p.m)
                out(:, i) = in(:,i).^p.m(i);
            end
        end

        %% Linear attenuation, Meier et al model
        % k(1) = left, k(2) = right
        % set p.k = [1 1] for no effect
        function [out, p] = linear_attenuation(p, in)
        if ~isfield(p, 'offset'), p.offset = 0; end
            if ~isfield(p, 'k'), p.k = [1 1]; end
            if ~isfield(p, 'abs'), p.abs = 0; end
            for i = 1:length(p.k)
                if p.abs == 1
                    out(:, i) = abs(p.k(i))*in(:,i);
                elseif p.abs == 0
                    out(:, i) = p.k(i)*in(:,i);
                end
            end
            out = out+p.offset;
        end
        %% Normalization, Meier model
        % U(1) = weight on left input in demoninator for left output
        % U(2) = weight on right input in denom for left
        % U(3) = weight on left in denom for right
        % U(4) = weight on right in denom for right
        % current implementation fits U(2) and U(3)
        % set all to 0 for no effect
        function [out, p] = normalization(p, in)
            if ~isfield(p, 'U'); p.U = [0 0 ; 0 0]; end
            if ~isfield(p, 'sigma'); p.sigma = 1; end
            if ~isfield(p, 'abs'), p.abs = 0; end
            if p.abs == 1
                  out(:,1) = in(:,1)./(abs(p.U(1))*in(:,1) + abs(p.U(2))*in(:,2)+abs(p.sigma));
                  out(:,2) = in(:,2)./(abs(p.U(3))*in(:,1) + abs(p.U(4))*in(:,2)+abs(p.sigma));

            elseif p.abs == 0
                out(:,1) = in(:,1)./(p.U(1)*in(:,1) + p.U(2)*in(:,2)+p.sigma);
                out(:,2) = in(:,2)./(p.U(3)*in(:,1) + p.U(4)*in(:,2)+p.sigma);
            end
        end

        % for alternative models, incorporate time
        function [out, p] = tau_normalization(p, in)
            if ~isfield(p, 'U'); p.U = [0 0 ; 0 0]; end
            if ~isfield(p, 'sigma'); p.sigma = 1; end
            if ~isfield(p, 'abs'), p.abs = 0; end
            y = b_s.gamma(5, p.tau/1000,0:p.dt:p.dt*length(in));
            y = reshape(y, length(y), 1);

            for i = 1:2
                in_conv(:, i)  =p.dt.*conv(  in(:, i), y, 'same');
            end
            if p.abs == 1
                out(:,1) = in(:,1)./(abs(p.U(1)) * in_conv(:,1) + abs(p.U(2)) * in_conv(:,2) + abs(p.sigma));
                out(:,2) = in(:,2)./(abs(p.U(3)) * in_conv(:,1) + abs(p.U(4)) * in_conv(:,2)+abs(p.sigma));
            elseif p.abs == 0
                out(:,1) = in(:,1)./(p.U(1)*in_conv(:,1) + p.U(2)*in_conv(:,2)+p.sigma);
                out(:,2) = in(:,2)./(p.U(3)*in_conv(:,1) + p.U(4)*in_conv(:,2)+p.sigma);
            end
        end

        %% alternative models functions
        function kfoldErr = cross_calibrate(p,data, freeList)
            if ~isfield(p, 'Kfold')
                p.Kfold = 10;
            end
            p.costflag = 0;
            tt = randi( p.Kfold, size(data.experiment.response,1), 1);
            for f = 1:p.Kfold
                n = 0;
                err = 0;  t = {}; S = {}; stim = {};
                % do the training
                train_ind =find(tt~=f);
                test_ind =find(tt==f);

                % do the training
              %  p.costflag = 1; p = fitcon('b_s.getErr',p, freeList, data, train_ind);
                p.costflag = 1; p = fit('b_s.getErr',p, freeList, data, train_ind);
                % grab error for this model fit, don't include the limit
                % costs
                p.costflag = 0; [err,~,~,~,~,~,~] = b_s.getErr(p, data, train_ind);
                % measure the error on the test
                kfoldErr(f) = err;
                % do the test
            end % end of cross-validatation
        end
        function [err, predModel] = simpleAverage(p, data)
            % calculates the expected output and err for a simple average
            % across the two eyes, no free parameters so fitting and
            % cross-validation unnecessary
            err = 0; n = 0;
            for i = 1:size(data.experiment.response,1) % for each run
                t{i} = data.t;
                % begin by getting the calibrated joystick response
                respJoy{i}  = data.experiment.response(i,:); % data.experiment.response(i,:) has the response for that particular run
                respJoyCalib{i}= b_s.joyFunction(p, respJoy{i} );
                predModel{i} = 0.5*data.experiment.LEcontrast(i, :) + 0.5*data.experiment.REcontrast(i, :);
                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                tmp = sum((predModel{i}(id)-respJoyCalib{i}(id)).^2);
                err = err + tmp;
                n = n+sum(id); % count the number of data points
            end
            err = err/n;
        end
        function [err, predModel] = simpleMax(p, data)
            % calculates the expected output and err for a simple max
            % across the two eyes, no free parameters so fitting and
            % cross-validation unnecessary
            err = 0; n = 0;
            for i = 1:size(data.experiment.response,1) % for each run
                t{i} = data.t;
                % begin by getting the calibrated joystick response
                respJoy{i}  = data.experiment.response(i,:); % data.experiment.response(i,:) has the response for that particular run
                respJoyCalib{i}= b_s.joyFunction(p, respJoy{i} );
                predModel{i} = max(data.experiment.LEcontrast(i, :), data.experiment.REcontrast(i, :));
                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                tmp = sum((predModel{i}(id)-respJoyCalib{i}(id)).^2);
                err = err + tmp;
                n = n+sum(id); % count the number of data points
            end
            err = err/n;
        end
        function [err, predModel] = rivalry(p, data)
            % calculates the expected output and err for simple rivalry, no
            % free parameters so fitting and cross-validation unnecessary.
            % Works differently from the other models because no free
            % parameters to fit

            err = 0; n = 0;
            for i = 1:size(data.experiment.response,1) % for each run

                % begin by getting the calibrated joystick response
                respJoy{i} = data.experiment.response(i,:); % data.experiment.response(i,:) has the response for that particular run
                respJoyCalib{i}= b_s.joyFunction(p, respJoy{i} ) ;

                diffval = abs([data.experiment.LEcontrast(i, :)-respJoyCalib{i}; data.experiment.REcontrast(i, :)-respJoyCalib{i}]) ; % diff eye-contrast and joystick
                for  j = 1:length(diffval)
                    tmp = find(min(diffval(:, j)) ==diffval(:, j));
                    if ~isempty(tmp)
                        whicheye(j) = tmp(1);
                    else
                        whicheye(j) =NaN;
                    end
                end
                predModel{i} = NaN(size(data.experiment.LEcontrast(i, :)));
                predModel{i}(whicheye ==1) = data.experiment.LEcontrast(i, whicheye ==1);
                predModel{i}(whicheye ==2) = data.experiment.REcontrast(i, whicheye ==2);

                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                tmp = sum((predModel{i}(id)-respJoyCalib{i}(id)).^2);
                err = err + tmp;
                pred_i1 = mean(data.experiment.LEcontrast(i,id) + data.experiment.REcontrast(i,id))/2;
                n = n+sum(id); % count the number of data points
            end
            err = err/n;
        end
        function [err, predModel] = meanJoystick(p, data)
            % calculates the expected output and err for a simple mean. The second allows for two
            % intercepts, placed to minimize the mse.
            err = 0; n = 0; allJoy = [];
            for i = 1:size(data.experiment.response,1) % for each run
                respJoy{i} = data.experiment.response(i,:); % data.experiment.response(i,:) has the response for that particular run
                respJoyCalib{i}= b_s.joyFunction(p, respJoy{i} );
                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                tmp= respJoyCalib{i};
                allJoy = [allJoy tmp(id)];
            end
            for i = 1:size(data.experiment.response,1)
                predModel{i} = mean(allJoy); % calculate mse for simple mean
                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                err = err + sum((predModel{i} -respJoyCalib{i}(id)).^2);
                n = n+sum(id); % count the number of data points
            end
            err = err/n;
        end
        function [err, predModel] = DualMeanJoystick(p, data)
            % calculates the expected output and err for two intercept
            % models. The first is a simple mean. The second allows for two
            % intercepts, placed to minimize the mse.
            err = 0; n = 0; allJoy = [];
            for i = 1:size(data.experiment.response,1) % for each run
                respJoy{i} = data.experiment.response(i,:); % data.experiment.response(i,:) has the response for that particular run
                respJoyCalib{i}= b_s.joyFunction(p,respJoy{i});
                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                tmp= respJoyCalib{i};
                allJoy = [allJoy tmp(id)];
            end
            allJoyMean = nanmean(allJoy);
            allJoyUpper = nanmean(allJoy(allJoy>0));
            allJoyLower = nanmean(allJoy(allJoy<0));
            for i = 1:size(data.experiment.response,1)
                diffval = abs([allJoyUpper-respJoyCalib{i}; allJoyLower-respJoyCalib{i}]) ; % difference between each eye's  contrast and joystick position (2 x time)
                for  j = 1:length(diffval) % for each length in time
                    tmp = find(min(diffval(:, j)) ==diffval(:, j)); % which eye's contrast is closest to the joystick position
                    if ~isempty(tmp)
                        which_i(j) = tmp(1); % which intercept is closest
                    else
                        which_i(j) =NaN; end
                end
                id = ~isnan(respJoy{i}) & ~isnan(respJoyCalib{i});
                predModel{i} = NaN(size(data.experiment.LEcontrast(i, :)));
                predModel{i}(which_i ==1) =allJoyUpper;
                predModel{i}(which_i ==2) = allJoyLower;
                err = err + sum((predModel{i}(id)-respJoyCalib{i}(id)).^2);
                n = n+sum(id); % count the number of data points
            end
            err = err/n;
        end
        function [out, p] = weightedAverage(p, S)
            % weighted average across the two eyes, with a free parameter
            % determining the relative weights
            if ~isfield(p, 'wa'); p.wa = .5; end
            out = p.wa*S(:,1) + (1-p.wa)*S(:, 2);
            out = reshape(out, 1, length(out)); % turn into a row vector
        end
        function [out, p] = ds2006(p, S)
            Sconv = S+0.5; % move to  0 -1 units.
            y = b_s.gamma(5, p.tau/1000,0:p.dt:p.dt*length(S));
            y = reshape(y, length(y), 1);
            Sconv(:, 1)  =p.e(1)*p.dt.* conv( Sconv(:, 1), y, 'same');
            Sconv(:, 2)=  p.e(2)*p.dt.* conv( Sconv(:, 2), y, 'same');
            % A gain-control theory of binocular combinaton, Ding &
            % Sperling 2006, implenting eq 6
            numLE= 1+Sconv(:, 1); numRE= 1+Sconv(:, 2);
            denom= 1+Sconv(:, 1) + Sconv(:, 2) ;
            out = (numLE./denom).*p.e(3).*S(:, 1) + (numRE./denom).*p.e(4).*S(:, 2);
            out = reshape(out, 1, length(out)); % turn into a row vector
            out  = out -0.5; % move back into -0.5-0.5 units
        end
        function [out, p] = dskl(p,S)
            warning('calling b_s.dskl,  many params for our data'); return
            %             g(1, :) = (S(:, 1)./p.ge(1)).^p.gamma_ge(1);
            %             g(2, :) = (S(:, 2)./p.ge(2)).^p.gamma_ge(2);
            %             ge(1, :) = g(1).*p.mu(1).*S(:, 1);
            %             ge(2, :) = g(2).*p.mu(2).*S(:, 2);
            %
            %             gc(1, :) = p.gc(1)*p.mu(1)*S(:, 1);
            %             gc(2, :) = p.gc(2)*p.mu(2)*S(:, 2);
            %             b(1) = 1+p.beta(1).^p.gamma(1);
            %             b(2) = 1+p.beta(2).^p.gamma(2);
            %             a(1) = 1+p.alpha(1).^p.gamma(1);
            %             a(2) = 1+p.alpha(2).^p.gamma(2);
            %
            %             num(1, :) = 1+(ge(2, :)/(1+b(1)*gc(1, :))); % dominant eye in equation A7 2013 paper
            %             num(2, :) = 1+(ge(1, :)/(1+b(2)*gc(2, :))); % nondominant eye
            %             denom(1, :) = 1+(gc(2, :)/(1+a(1)*gc(1, :)));
            %             denom(2, :) = 1+(gc(1, :)/(1+a(2)*gc(2, :)));
            %
            %             out = ( S(:, 1).*p.mu(1).*(num(1, :)/denom(1, :)) ) + ( S(:, 2).*p.mu(2).*num(2, :)/denom(2, :) ); % tada!
            %             out = reshape(out, 1, length(out)); % turn into a row vector
        end
        function [out, p] = bmg2007(p,S )
            % Baker, Meese & Georgeson (2007) Spatial Vision, 20: 397‐413doi:10.1163/156856807781503622
            % also referenced in Hess, 2008
            % the amblyopia model is Baker et al. Contrast masking in
            % strabismic amblyopia: Attenuation, noise, interocular
            % suppression and binocular summation, 2008
            % these models work in percent contrast so need to do some
            % conversions
            S_tmp = (S+.5);
            out = (S_tmp(:, 1).^p.m)./(abs(p.S)+S_tmp(:, 1)+abs(p.w(2)).*S_tmp(:, 2)) +  ...
                (S_tmp(:, 2).^p.m)./(abs(p.S)+abs(p.w(1)).*S_tmp(:, 1)+S_tmp(:, 2));
            out = reshape(out, 1, length(out));
            out = (out.^p.pq(1))./(p.Z + out.^p.pq(2));
            out = out-.5;
            p.S = abs(p.S); p.w = abs(p.w); % used in their abs
        end
        function [err, p, out] = minkowski(p, S, calibrated_data)
            out = ((S(:, 1).^abs(p.n)+S(:, 2).^abs(p.n))/2).^(1/abs(p.n));
            err = sum((out-calibrated_data).^2);
        end

        function plot_alt_models(p, plotStr, varargin)
            if nargin ==2;      normalize_flag = 0;
            else    normalize_flag = varargin{1};       end
            for s = 1:length(plotStr)
                switch plotStr{s}
                    case {'softmax', 'softmax_tau','ds2006', 'ds2006_tau','bmg2007', 'bmh2008', 'weightedAverage'}
                        y(s) = eval(['p.', [plotStr{s}, '.kfoldErr']]);
                        yerr(s) = eval(['p.', [plotStr{s}, '.kfoldStd']]);
                    case {'rivalry', 'meanJoystick', 'DualMeanJoystick', 'simpleAverage', 'simpleMax'}
                        y(s) = eval(['p.', [plotStr{s}, '.Err']]);
                        yerr(s) = 0;
                    otherwise
                        warning(['var p has no Err for this model ... ' , plotStr{s}]);
                end
            end
            if normalize_flag;      yplot = y-y(1);
            else    yplot = y;     end
            if contains(p.sID(1:4) ,'AM');          clr = 'r'; offset = -.2;
            elseif contains(p.sID(1:4) ,'NS');      clr = 'b'; offset = 0.2;
            elseif contains(p.sID(1:4) ,'BD');  clr = 'g';      offset  = 0 ; end
            eh=errorbar([1:length(plotStr)]+offset+[.02*randn(size(yplot))], yplot, yerr, 'CapSize', 0); hold on
            set(eh, 'Marker','o', 'MarkerFaceColor', clr,'MarkerEdgeColor', clr, 'LineStyle', 'none', 'Color', clr, 'MarkerSize', 6);
            set(gca, 'XTick', 1:length(plotStr));
            set(gca, 'XTickLabel',plotStr);
            set(gca, 'XLim', [0 length(plotStr)+1])
            set(gca, 'YLim', [0 0.6])
        end
        function plot_parameters(p, plotStr, varargin)
            if nargin ==2
                logflag  = 0;
            else
                logflag = varargin{1};
            end
            for s = 1:length(plotStr)
                y = eval(['p.', plotStr{s}]);
                if contains(p.sID(1:4) ,'AM');          clr = 'r'; offset = -.2;
                elseif contains(p.sID(1:4) ,'NS');      clr = 'b'; offset = 0.2;
                elseif contains(p.sID(1:4) ,'BD');  clr = 'g';      offset  = 0 ; end
                if logflag  == 0
                    eh=plot(s+offset+[.02*randn(size(y))], y); hold on
                else
                    tol = 0.001;
                    y(y<=0) = tol;
                    eh=plot(s+offset+[.02*randn(size(y))], log(y)); hold on
                end
                set(eh, 'Marker','o', 'MarkerFaceColor', clr,'MarkerEdgeColor', clr, 'LineStyle', 'none', 'Color', clr, 'MarkerSize', 6);
            end
            set(gca, 'XTick', 1:length(plotStr));
            set(gca, 'XTickLabel',plotStr);
            set(gca, 'XLim', [0 length(plotStr)+1]);
        end
        function plot_modelresponses(eyeL, eyeR, resp, pred, clr)
            % plots timecourses of model responses
            h(1) = plot(eyeL, 'k--'); hold on
            h(2)  = plot(eyeR, 'k-.');
            h(3)  = plot(resp , 'r-', 'LineWidth', 1); hold on
            h(4) = plot(pred, clr); hold on
        end
        % legend(h, {'calib response', 'LE stim', 'RE stim', p.model});

        %% gridsearch function
        function [pBest,errBest] = gridsearch(funName,p,gridParams,gridList, varargin)
            % [pBest,errBest] = gridsearch(funName,p,gridParams,gridList, var1, var2, ..)
            %
            % Grid search to find best initial parameters for optimization.
            %
            % Inputs:
            %   funName            name of error function. Must have form compatible
            %                      with 'fit' and 'fitcon' :
            %                       [err] = <funName>(params, var1, var2, ...)
            %   p                  Structure with initial parameters
            %   gridParams         List of parameters (fields of p) to be gridded
            %   gridList           cell array of 1-d grid values, in order
            %                      corresponding to gridParams
            %   var1, etc.         additional variables to be passed in to funName
            %
            % Outputs:
            %   pBest              parameter structure with lowest err in the grid
            %   errBest            corresponding error value
            %
            % Note, parameters in 'p' that are not in 'gridParams' are kept at these
            % initial values for all function evaluations.
            %
            % Note also: gridParams can contain elements of an array, e.g. 'a(1)'.

            nParams = length(gridParams);

            % Generate a string to be evaluated that has the form (depending on
            % nParams): '[X{1},X{2},X{3}] = ndgrid(gridList{1},gridList{2},gridList{3});'

            leftStr = '[';
            rightStr = ' = ndgrid(';
            for i=1:nParams
                leftStr = strcat(leftStr,sprintf('X{%d},',i));
                rightStr = strcat(rightStr,sprintf('gridList{%d},',i));
            end
            str = sprintf('%s] %s);',leftStr(1:end-1),rightStr(1:end-1));
            eval(str);

            errBest = 1e10;
            pBest = p;

            % Loop through elements of the grid, saving the current best fit
            for i=1:length(X{1}(:))
                for j=1:nParams
                    % generate a string to be evaluated of the form (for example):
                    % 'p.a(1) = X{1}(1);'
                    str = sprintf('p.%s = X{%d}(%d);',gridParams{j},j,i);
                    eval(str)
                end

                % evaluate the function with these parameters
                err = feval(funName, p,varargin{:});

                % save the fit if it's currently the best
                if err<errBest
                    pBest = p;
                    errBest = err;
                end
            end
        end

        %% second stage analyses
        function gatherTable(datatype, varargin)
            if nargin == 1
                condition = NaN;
            else
                condition = varargin{1};
            end
            resultsDataDir = [cd filesep ['fitdata_', datatype] filesep 'model_fits'];

            if ~isnan(condition)
                csvSaveName = [datatype, '_', condition, '_fits'];
                files = dir([resultsDataDir filesep '*' condition '*.mat']);
            else
                csvSaveName = [datatype, '_fits'];
                files = dir([resultsDataDir filesep '*regular.mat']);
            end

            for i = 1:length(files)
                fileLoad = [files(i).folder filesep files(i).name];
                pdata = load(fileLoad);
                p = pdata.p; clear pdata;

                if strcmpi(condition, 'altmodels')
                    % collect error values for the alternative models
                    datatable.sID{i,1} = p.sID;
                    datatable.group{i,1} = p.sID(1:2);
                    datatable.softmax{i,1} = p.softmax.kfoldErr;
                    datatable.meanJoystick{i,1} = p.meanJoystick.Err;
                    datatable.simpleAverage{i,1} =p.simpleAverage.Err;
                    datatable.weightedAverage{i,1} = p.weightedAverage.kfoldErr;
                    datatable.softmax_tau{i,1} = p.softmax_tau.kfoldErr;
                    datatable.ds2006{i,1} = p.ds2006.kfoldErr;
                    datatable.ds2006_tau{i,1} = p.ds2006_tau.kfoldErr;
                    datatable.bmg2007{i,1} = p.bmg2007.kfoldErr;
                    datatable.simpleMax{i,1} = p.simpleMax.Err;
                    datatable.rivalry{i,1} = p.rivalry.Err;
                    datatable.dualMeanNull{i,1} = p.DualMeanJoystick.Err;

                else % end catch for alternative models

                    % in p:
                    % k(1) is left eye, k(2) is right eye
                    % U(2) = weight on right eye in denom for left eye
                    % U(3) = weight on left eye in denom for right eye

                    % whichever k == 1 is the fellow eye
                    if all(p.k == 1) % catch for ppl with both
                        disp([p.sID ': Both k == 1'])
                        FE = 1; AE = 2;
                        sAE = 'right'; sFE = 'left';

                        if strcmpi(p.sID(1:2), 'am') || strcmpi(p.sID(1:2), 'bd')

                            % amblyopia or binocular disorder
                            if strcmpi(p.sID(4), 'L')
                                % right eye is fellow
                                FE = 2; AE = 1;
                            elseif strcmpi(p.sID(4), 'R')
                                % left eye is fellow
                                FE = 1; AE = 2;
                            else
                                disp('  (sanity check, should not get to case)')
                            end
                        else
                            disp('.. control subject with k1=k2')
                        end

                    elseif p.k(1) == 1 % left eye is fellow
                        FE = 1; AE = 2;
                        sAE = 'right'; sFE = 'left';
                    elseif p.k(2) == 1 % right eye is fellow
                        FE = 2; AE = 1;
                        sAE = 'left'; sFE = 'right';
                    else
                        disp('sanity check: this should not happen')
                    end

                    tempU = p.U(2:3);
                    % tempU(1) = weight on RE in equation for LE
                    % tempU(2) = weight on LE in equation for RE

                    datatable.sID{i,1} = p.sID;
                    datatable.group{i,1} =  p.sID(1:2);
                    if ~isnan(condition)
                        datatable.condition{i,:} = condition;
                    end

                    % recording calibration parameters
                    datatable.delay{i,1} = p.delay;
                    datatable.intercept{i,1} = p.intercept;
                    datatable.slope{i,1} = p.slope;
                    datatable.clean_range{i,1} = p.clean_range;
                    datatable.n_good{i,1} = p.n_good;


                    if ~isnan(condition)
                        switch lower(condition)
                            case 'reduced'
                                datatable.calibErr{i,1} = p.calibErr;
                            otherwise
                                datatable.joyCalErrMean{i,1} = p.errMean;
                                datatable.joyCalErrInd{i,1} = p.errInd;
                        end
                    else
                        datatable.calibErr{i,1} = p.calibErr;
                    end

                    datatable.kLE{i,1} = p.k(1);
                    datatable.kRE{i,1} = p.k(2);
                    datatable.beforeRefit_kLE{i,1} = p.kBeforeRefit(1);
                    datatable.beforeRefit_kRE{i,1} = p.kBeforeRefit(2);

                    datatable.U2{i,1} = p.U(2);
                    datatable.U3{i,1} = p.U(3);

                    datatable.sigma{i,1} = p.sigma;

                    datatable.AE{i,1} = sAE; % should be less than 1
                    datatable.FE{i,1} = sFE; % should be equal to 1

                    datatable.kAE{i,1} = p.k(AE);
                    datatable.kFE{i,1} = p.k(FE);

                    datatable.uAE{i,1} = tempU(AE);
                    datatable.uFE{i,1} = tempU(FE);

                    datatable.step1attenuationErr{i,1} = p.step1attenuationErr;
                    datatable.step2normalizationErr{i,1} = p.step2normalizationErr;
                    datatable.softmaxErr{i,1} = p.softmaxErr;
                end % end if altmodels else everything else

                clear pdata p
            end

            T = struct2table(datatable);

            writetable(T, ['fitdata_', datatype filesep csvSaveName '.csv']);
            save(['fitdata_', datatype filesep csvSaveName '.mat'], 'T')
            disp('done')
        end % end gatherTable



    end % end methods
end % end class def