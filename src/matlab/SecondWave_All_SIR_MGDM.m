%clear all
close all
clf
% Just the normal stuff of clearing useless data and closing alll graphs.

addpath(genpath('mEpiTools'))
addpath('data')
addpath('subroutines')
addpath('align_Ylabels')
% Adding the folder paths which will be needed for this code to run.

ifplotPR = 1;
dists = 1;
debug = 0;
ifprint  = 1;
span = 7;  % 7 = weekly 14 = biweekly
% Defining some constants

global h1
global h2
global sb_number
global ax2
format long
global pos1
global mean_posterior
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
% Setting some global constands and some parameters. 

%%

%debug = 0;

if(dists) % This dists appears to be a boolean which says if we want to analyze districts or states.
% This part of the if should execute now as the dists is 1/True as per line 13.
    fullName = 'https://api.covid19india.org/csv/latest/districts.csv';
    name = getDistrictNames(); % This gets the names of districts to fetch from the API csv file.
    % This function is made on the line 865 which returns the list of the districts.
	tt1 = datestr(datenum(2020,4,26)); % Converts the numbers 26/4/2020 to datatype.
	sb_number = 4; % Setting value of a global constant
    % Use of sb_number is not clear
    
else
    fullName = 'https://api.covid19india.org/csv/latest/states.csv';
    name = getStateNames(); 
    % The stuff done here is similar, we analyze data of states here.
    % Declaration at SecondWave_All_SIR_Univ_MGDM line 743
	tt1 = datestr(datenum(2020,3,14)); % Getting date 14/3/2020
	sb_number = 5; % Still not clear what this value does.
    % Maybe this value just finds if we want to do analysis of states of districts or it can also have other cases.
    % If it was limited to districts or states, the work would have been done by only dists bool value.
end

%fullName = 'https://api.covid19india.org/csv/latest/districts.csv';
%fullName = 'https://api.covid19india.org/csv/latest/states.csv';

urlwrite(fullName,'dummy.csv');
% urlwrite documentation: https://in.mathworks.com/help/matlab/ref/urlwrite.html
% webwrite is preferred option now in R2021
% This file has date and district wise info of covid data and maybe write the data in our dummy.csv file.
tableAll =readtable('dummy.csv');
% We are reading the dummy.csv

span1 = 14;  % 7 = weekly 14 = biweekly

mean_si = 4.7;
min_mean_si = 3.7; 
max_mean_si = 6.0;
    
std_si = 2.9; 
min_std_si = 1.9; 
max_std_si = 4.9; 
% Now we have defined some new constants to scan for 14 days.

%%
std_std_si = (max_std_si-std_si)/2; 
std_mean_si = (max_mean_si-mean_si)/2; 
% this seems like comparing current mean and standard deviation with max one.

for n = 1:length(name)
    
    
    
    Location = name{n};
    
    if(dists)
        % This parts is executed if we are analyzing districts. 
        indR = find(contains(tableAll(:,3).District,Location)==1)
        % We are filtering out the rows with our district of interest.
        % Because they mean proper data for that time isn't avialable.
        % not clear: why are we not just directly extracting the districts
        % from this table as it can be easier to do so.
        % Why are we making a table with one table column.
        % 3rd column is for the districts. 
        % We are trying the get the array of some rows from the table.
        % We are extracting columns in vector from the tableAll data of
        % dummy.csv.
        R = table2array(tableAll(indR,5)); % Recovered column
        D = table2array(tableAll(indR,6)); % Deceased column
        C = table2array(tableAll(indR,4)); % Confirmed people.
        O = table2array(tableAll(indR,7)); % Other
        T = table2array(tableAll(indR,8)); % Number of tests 
        % Not clear here. 
        
        date = table2array(tableAll(indR,1)); % date colum. 
		states = table2array(tableAll(indR,2)); % the state column
		myFolder1 = states{1};		
		
		dat = datenum(2020,2,30); % converting 30/2/2020 to datetime.
        I = find((datenum(date) >= dat)); % Getting dates where date greater then dat in above line.

        R = R(I);
        D = D(I);
        C = C(I);
        O = O(I);
        T = T(I);		
        date = date(I);
        % Getting all the column output with date filtering.ate = date(I);
        
		
    else
        indR = find(contains(tableAll(:,2).State,Location)==1);
       
        
        C = table2array(tableAll(indR,3));
        R = table2array(tableAll(indR,4));
        D = table2array(tableAll(indR,5));
        O = table2array(tableAll(indR,6));
        T = table2array(tableAll(indR,7));
        date = table2array(tableAll(indR,1));
		
		dat = datenum(2020,2,30);
        I = find((datenum(date) >= dat));
        
        R = R(I);
        D = D(I);
        C = C(I);
        O = O(I);
        T = T(I);		
        date = date(I);	
        
    end
    %{
         C(end) = [];
         R(end) = [];
         D(end) = [];
         O(end) = [];
         T(end) = [];
         date(end) = [];
    %}
    
    A = C - (D+R+O); %active cases
    %span = 7;  % 7 = weekly 14 = biweekly
    
    dailyC = diff(C);
    dailyCmean = movmean(dailyC,span);
    
    if(debug) % This output is shown when debugging the code.        
        figure
        subplot(2,1,1)
        stem(date,C,'-o','linewidth',2)
        datetick('x','dd/mm','keeplimits','keepticks')
        title(Location)
        xlabel('Date')
        ylabel('Confirmed cases')
        xtickangle(45)
        subplot(2,1,2)
        stem(date(2:end),diff(C),'-o','linewidth',2)
        datetick('x','ddmmm','keeplimits','keepticks')
        xlabel('Date')
        ylabel('Infection rate')
        xtickangle(45)
    end
    
    
    
    PR = diff(movmean(C,span)) ./ diff(movmean(T,span)) *100;
    CFR = diff(movmean(D,span)) ./ diff(movmean(C,span)) *100;
    Testing = diff(T);
    
    %%
    
    date0 = date;
    %date1 = datenum(2020,4,26);
	if(dists)
		tt1 = datestr(datenum(date0(1)));		
		sb_number = 4;
		
	else
		tt1 = datestr(datenum(date0(1)));
		sb_number = 5  ;
	end
    
    date = datenum(date0);
    
    
    
    % date1= datenum(2021,3,5); %in paper
    %date1= datenum(2021,3,6); %india
    %  date1= datenum(2021,3,5); %in paper
    %   date1= datenum(2021,3,15); %in paper
    date1= datenum(2021,4,1); %in paper
    
    datestr(date1)
    
    %dateend= datenum(2021,5,1); %in manuscript
    dateend= date(end); %in manuscript
    % end is a matlab defined variable which means end of the array.
    
    I = find((date >= date1) & (date <= dateend));
    % So here we find the dates in the range starting date to ending date.
    
    init.country = Location;
    init.date = date(I);
    init.C = C(I)-C(I(1)-1);
    init.D = D(I) -D(I(1)-1) ;
    init.R = R(I)-R(I(1)-1);
    
    % init.C = C(I);
    % init.D = D(I);
    % init.R = R(I);
    
    init.Npop = 80e6; % Does not matter much beyond certain minimum value
    
    
    %% doubling time
    
    diffC = diff(init.C);
    % q1 = init.C(1);
    % q2 = init.C(end);
    q1 = diffC(1);
    q2 = diffC(end);
    t = length(diffC);
    
    
    td = t*log(2)/log(q2/q1);
    
    %% Using exponential model
    Cexp = init.C(1:end);
    texp = (1:length(Cexp))';
    [f,gof,out]=fit(texp,Cexp,'exp1');
    
    td_2 = log(2)/f.b; % doubling time using exponential model
    
    tnew = 1:texp(end)+120;
    date0_model = init.date(1)-1;
    
    EXP.date = tnew+date0_model;
    EXP.C = f(tnew);
    
    %% Using Logistic model
    fun = @(A1,A2,A3,x)(A1./(1+A2*exp(-A3*x)));
    a0 = iniGuess(init.C); % get better initial guess using code by Milan Batista
    
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' ); % Non-linear regression of square order polynomial.
    opts.Display = 'off';
    opts.Algorithm = 'Trust-Region';
    opts.Robust = 'Off';
    
    opts.Lower = [0 0 0];
    opts.StartPoint = a0;
    
    [flog, gof] = fit( texp, Cexp, fun, opts ); % Fitting the curve with given options.
    
    % For co-effcients
    fun = @(a,x)(a(1)./(1+a(2)*exp(-a(3)*x)));
    logmdl = NonLinearModel.fit(texp,Cexp,fun,a0);
    
    LOG.date = EXP.date;
    LOG.C = flog(tnew);
    
    
    %% SIR model - Original Code by Milan Batista
    [SIR] = SIR_RR(init,'prn','off', 'plt','off'); -2 /2
    
    
    SIR.R0
    SIR_Ca = SIR.Ca + C(I(1)-1);
    
    
    % first fit Gaussian for Decay data
    
    dt = SIR.t(2)-SIR.t(1);
    dI = diff(SIR.Ca)/dt; %diff(res.I)/dt;
    dI = dI.';
    [val, idx] = max(dI);
    
    peakdate = datestr(SIR.date(idx+2),'dd-mmm-yy');
    fprintf(['Date of Peak is : ',peakdate,'\n'])
    
    dI_data =  dI(idx+1:end);
    Sn = (0:length(dI_data)-1);
    [fI,gof,out]=fit(Sn',dI_data','gauss1');
    
    coeff = coeffvalues(fI);
    a = coeff(1);
    mu = coeff(2);
    sig = coeff(3);
    
    fIfit = @(a,mu,sig,x)(a.*exp(-((x-mu)/sig).^2));
    %a = dI_data(1);
    
    % correction for asymmetric gaussian
    al = 1.1; %1.3; % Correction factor
    sig = al*sig;
    x0= fIfit(a,mu,sig,Sn);
    scalex = dI_data(1)/ x0(1) % scale it so that the peak is same
    xnew = x0*scalex;
    
    %
    if(debug)
        figure
        plot(fI,Sn, dI_data )
        hold on
        plot(Sn,xnew,'.k')
        legend('fit','Actual','Corrected')
        title('Fitted vs Modified Gaussian')
    end
    
    % SEIC.date = SEIR.date;
    SIR_Cnew = SIR.Ca;
    dCnew = dI;
    dCnew(idx+1:end) = xnew;
    %
    dCnew(dCnew < 0) = 0;
    
    for i=idx+1:length(SIR.Ca)
        SIR_Cnew(i) = SIR_Cnew(i-1)+dt*dCnew(i-1);
    end
    
    SIR_Cnew(idx) = 0.5*(SIR_Cnew(idx-1) + SIR_Cnew(idx+1));
    
    
    %%
    
    
    %     x=x.';
    %     y=y.';
    x = SIR.date(2:end);
    %  y = diff(SIR.Ca)/(SIR.t(2)-SIR.t(1));
    y = diff(SIR_Cnew)/(dt);
    
    x = x.';
    y = y.';
    slope = diff(y)/(SIR.t(2)-SIR.t(1));
    err = [1.1*y; 0.9*y];
    
    
    %tstart = datenum(datetime(2021,02,16));
    tstart = datenum(datetime(2021,03,1));
    
    tend = datenum(datetime(2021,07,1));
    
    
    figure
    %set(gcf,'Position',[50 50 3.5*624/3 4*832/3])
	sz = 2/3;
    set(gcf,'Position',[50 50 676*sz 1040*sz*sb_number/5])
    h1 = subplot(sb_number,1,1);
	ax1 = gca;

    hold on

    h111 = bar(date(2:end),dailyC);
	set( get( get( h111, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    plot(date(2:end),dailyCmean,'k','LineWidth',2*sz)
    %    plot(SIR.date(2:end),diff(SIR.Ca)/(SIR.t(2)-SIR.t(1)),'R','linewidth',3)
    plot(SIR.date(2:end),diff(SIR_Cnew)/dt,'r','linewidth',3*sz)
    
    patch([x, fliplr(x)], [err(1,:), fliplr(err(2,:))], 'm', 'FaceAlpha',0.1,'LineStyle','--','LineWidth',1)   % Shaded Confidence Intervals
    %  plot(SIR.date(2:end),diff(SIR.Ca)/(SIR.t(2)-SIR.t(1)),'R','linewidth',3)
    plot(SIR.date(2:end),diff(SIR_Cnew)/dt,'r','linewidth',3*sz)

    xlim([tstart tend])
    
    ylim([0,1.2*max(diff(SIR_Cnew)/dt)]);
    
    
    xticks([ datenum(datetime(2021,03,01)) datenum(datetime(2021,04,01)) datenum(datetime(2021,05,01)) datenum(datetime(2021,06,01)) datenum(datetime(2021,07,01)) datenum(datetime(2021,08,01))])
    % xticklabels({'' '' '' '' '' ''})
    datetick('x','ddmmmyy','keeplimits','keepticks')
    xax = get(gca,'XAxis');
    set(xax,'TickDirection','out')
    %tname = sprintf('%s, %s, TPR = %s, CFR = %s',Location,datestr(dateend),num2str(PR(end),'%.0f'),num2str(CFR(end),'%.2f'));
    tname = sprintf('%s',Location);
    
    %title(Location)
    title(tname)
    
    %xlabel('Date')
    yl = ylabel('Daily Cases');
	set(yl,'Units','normalized')
	set(yl,'Pos',[  -0.1  0.50 0])
    %xtickangle(45)
    grid on
    box on	
	


	
	
	
	
    % plot and see full data
    
    date = datenum(date0);
    h3 = subplot(sb_number,1,4);
	ax3 = gca;

    bar(date(2:end),diff(D))
    hold on
    plot(date(2:end),movmean(diff(D),7),'k','LineWidth',2*sz)
    %datetick('x','dd/mm','keeplimits','keepticks')
    %title(Location)
    %xlabel('Date')
    yl = ylabel('Daily Fatalities');
	set(yl,'Units','normalized')
	set(yl,'Pos',[  -0.1  0.50 0])
    xticks([ datenum(datetime(2021,03,01)) datenum(datetime(2021,04,01)) datenum(datetime(2021,05,01)) datenum(datetime(2021,06,01)) datenum(datetime(2021,07,01)) datenum(datetime(2021,08,01))])
    % xticklabels({'' '' '' '' '' ''})
    %ylim([0, 1.33*max(movmean(diff(D(end-100:end)),14))])
    datetick('x','ddmmmyy','keeplimits','keepticks')
	
    ylim([0, 1.2*max(movmean(diff(D(end-150:end)),7))])
    if(strcmpi(Location,'India' ))
        h3.YAxis.Exponent = 3;
    else
        h3.YAxis.Exponent = 2;
    end
    xax = get(gca,'XAxis');
    set(xax,'TickDirection','out')
    % 	tstart = datenum(datetime(2021,02,16));
    %     tend = datenum(datetime(2021,05,04));
    xlim([tstart tend])
    %xlabel('Date')
    %xtickangle(45)
	
    grid on
    hold off
	
	yyaxis(ax3,'right')
    plot(date(2:end),CFR,'color',[1 0.58 0],'linewidth',3*sz)	    
    yl = ylabel('Fatality Rate')	;
	set(yl,'Units','normalized')
	set(yl,'Pos',[  1.065  0.50 0])	
    set(ax3,{'ycolor'},{[1 0.58 0]}) 

			
    h4 = subplot(sb_number,1,3);
	ax4 = gca;	

    %stem(date,C-R,'-o','linewidth',2)
    %    bar(date,C-R-D)
    bar(date,A)
    
    hold on
    %    plot(date,movmean(C-R-D,7),'k','LineWidth',2)
    plot(date,movmean(A,7),'k','LineWidth',2*sz)
    
    %datetick('x','ddmmm','keeplimits','keepticks')
    h4label = ylabel('Active Cases');
	set(h4label,'Units','normalized')
	set(h4label,'Pos',[  -0.1  0.50 0])

    %ylim([0, 1.33*max(movmean(C(end-100:end)-R(end-100:end)-D(end-100:end),14))])
    ylim([0, 1.2*max(movmean(A(end-150:end),7))])
    xticks([ datenum(datetime(2021,03,01)) datenum(datetime(2021,04,01)) datenum(datetime(2021,05,01)) datenum(datetime(2021,06,01)) datenum(datetime(2021,07,01)) datenum(datetime(2021,08,01))])
    xax = get(gca,'XAxis');
    set(xax,'TickDirection','out')
	% xticklabels({'01/03' '01/04' '01/05' '01/06' '01/07' '01/08'})
    datetick('x','ddmmmyy','keeplimits','keepticks')
    %  h4label.Position(2) = -1;
    
    % 	tstart = datenum(datetime(2021,02,16));
    %     tend = datenum(datetime(2021,05,04));
    xlim([tstart tend])
    %xtickangle(45)
    grid on
    hold off
    
    
	

	
	
	
	
	
	
	
	

	
	oldFolder = cd('data');
    pwd;
    path = pwd;
    
    % check country name
    nname = Location;
    if strcmp("",nname)
        continue
    end
	
	% get data
    try
    
    I = size(C,1);
    I(1) = C(1);
	
    for i = 2:length(C)
        I(i) = C(i) - C(i-1);            
        if I(i) < 0
            I(i) = 0;
        end
    end
	
    catch
        disp(nname)
        continue
    end
    
    % check data
    if isempty(C)
        disp(nname)
        continue
    end


    % setup file name
    nname = strrep(nname,'/','_');
    nname = strrep(nname,' ','_');
    nname = strrep(nname,'-','_');
    nname = strrep(nname,'''','_');
    nname = strrep(nname,'(','_');
    nname = strrep(nname,')','_');
    fname = sprintf('getIncid%s.m',nname);
    
    % open new file
    fid = fopen(fullfile(path,fname),'w');
    if fid < 0
        fprintf('***Fail to open %s\n',fname);
        continue
    end

    fprintf(fid,'function [incid,district] = getIncid%s()\n',nname);
    fprintf(fid,'%%GETINCID%s Coronavirus data for %s\n',upper(nname),nname);
    fprintf(fid,'district = ''%s'';\n',strrep(nname,'_',' '));
    fprintf(fid,'incid.I = [...\n');
    for m = 1:length(I)
    %         if T + m - 1 < T0
    %             continue
    %         end
        fprintf(fid,'%d \n',I(m));
    end
    fprintf(fid,'%%<-------------- add new data here\n');
    fprintf(fid,']'';\n');

        fprintf(fid,'incid.start_date = datetime(''%s'');\n',tt1);
    fprintf(fid,'end\n');
    
    % close file
    fclose(fid);
    cd(oldFolder);
	
    start_date = [];
	incid = importData(Location,start_date);
	
	if isempty(incid)
        error('No data for %s.',country)
    end
    T = length(incid.I)
    t_start = 2:T-span1+1; % starting at 2 as conditional on the past observations
    t_end = t_start + span1-1;
    res = estimate_R(incid, ...
        'uncertain_si',...
        'config', make_config(...
        't_start',t_start,'t_end',t_end,...
        'mean_si', mean_si,...
        'std_mean_si', std_mean_si,...
        'min_mean_si', min_mean_si,...
        'max_mean_si', max_mean_si,...
        'std_si', std_si,...
        'std_std_si', std_std_si,...
        'min_std_si', min_std_si,...
        'max_std_si', max_std_si));

    plot1(res,'what',['R'],'title',sprintf('%s',incid.name),...
        'opt_I',options(...
        'plot_type','bar'))	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
    
    
    %%
    tend2 = date(end);
    if(ifplotPR)
        tname = sprintf('%s, %s',Location,datestr(date(end)));
        
        if(dists == 0)
            h5 = subplot(sb_number,1,5);
			ax5 = gca;
            %   subplot(2,1,1)
			hold on
			bar(date(2:end),Testing)
            plot(date(2:end),movmean(Testing,7),'k','LineWidth',2*sz)	 
            xticks([ datenum(datetime(2021,03,01)) datenum(datetime(2021,04,01)) datenum(datetime(2021,05,01)) datenum(datetime(2021,06,01)) datenum(datetime(2021,07,01)) datenum(datetime(2021,08,01))])
            % xticklabels({'01/03' '01/04' '01/05' '01/06' '01/07' '01/08'})
            datetick('x','ddmmmyy','keeplimits','keepticks')
            %  h4label.Position(2) = -1;
            
            % 	tstart = datenum(datetime(2021,02,16));
            %     tend = datenum(datetime(2021,05,04));
            %xlim([tstart tend2])
            datetick('x','ddmmmyy','keeplimits','keepticks')
			ylim([0, 1.2*max(movmean(Testing(end-150:end),7))])
            %title(tname)
            xlabel('Date')
	        yl = ylabel('Daily Tests');			
			set(yl,'Units','normalized')
			set(yl,'Pos',[  -0.1  0.50 0])

            %xtickangle(45)
            grid on 
			box on
		    hold off
			
			yyaxis(ax5,'right')
            %plot(date(2:end),PR,'-o','linewidth',2*sz,'MarkerFaceColor', 'r', 'MarkerSize', 5*sz )			
            plot(date(2:end),PR,'color',[0.2 0.7 0],'linewidth',3*sz)			
            yl = ylabel('Positivity Rate')	;
			set(yl,'Units','normalized')
			%pos = get(yl,'position')
			set(yl,'Pos',[  1.065  0.50 0])

			set(ax5,{'ycolor'},{[0.2 0.7 0]}) 

		end

        
		if (dists==0)
		    linkaxes([h1,h2,h3,h4,h5],'x')
		else
		    linkaxes([h1,h2,h3,h4],'x')
	    end
		h1.XLim = [datenum(datetime(2021,03,01))  datenum(datetime(2021,07,01))];	
		set(gcf, 'visible', 'on');
		set(findall(gcf,'-property','FontSize'),'FontSize',10)
		%yax = gca.YAxis; %
		xax = get(gca,'XAxis');
		set(xax,'TickDirection','out')
		%    set(gca,'TickDir','out');
		
		
        ax1.XAxis.FontSize = 9;
        ax2.XAxis.FontSize = 9;
        ax3.XAxis.FontSize = 9;
        ax4.XAxis.FontSize = 9;
		
		if (dists==0)
            ax5.XAxis.FontSize = 9;
 	    end

		
		myFolder = datestr(dateend);
		myFolder_data = sprintf('%s_data',datestr(dateend))

		
		if not(isfolder(myFolder))
			mkdir(myFolder)
		end		

		if not(isfolder(myFolder_data))
			mkdir(myFolder_data)
		end	
		
		if (dists==1)
		oldFolder = cd(myFolder);
		if not(isfolder(myFolder1))
			mkdir(myFolder1)
		end
	    cd(oldFolder);	
		oldFolder = cd(myFolder_data);
		if not(isfolder(myFolder1))
			mkdir(myFolder1)
		end
	    cd(oldFolder);			
		end
		
		if(ifprint)
			fname = sprintf('%s.png',Location)
			if (dists==1)
			    picturename = fullfile(myFolder, myFolder1, fname)
			else
    		    picturename = fullfile(myFolder, fname)
			end
			saveas (gcf, picturename)
			print(gcf,picturename,'-dpng','-r300');

		end
		dailyD = diff(D);
		fname = sprintf('%s.txt',Location);
		if (dists==1)
		fid = fopen(fullfile(myFolder_data, myFolder1, fname),'w');
		else
		fid = fopen(fullfile(myFolder_data, fname),'w');
		end
		%{
		fprintf(fid,'On %s, %s reported %d new cases, %d deaths and %d active cases. ',datestr(dateend), Location, dailyC(end), dailyD(end), A(end));
		
		if(dists== 0 && strcmp(Location,'India')==0)
		    fprintf(fid,"On this date, %d tests were conducted. This region has positivity rate of %0.1f %% and fatality rate of %0.1f %%. ", Testing(end), PR(end), CFR(end)) ;
			if(PR(end) > 10)
		        fprintf(fid,"Such high positivity rate remains a concern. ") ;
			end	
		elseif(dists== 0 && strcmp(Location,'India')==1)
		    fprintf(fid,"On this date, %d tests were conducted. This region has positivity rate of %0.1f %% and fatality rate of %0.1f %%. ", Testing(end-1), PR(end), CFR(end)) ;	
			if(PR(end) > 10)
		        fprintf(fid,"Such high positivity rate remains a concern. ") ;
			end			
		else 
         	fprintf(fid,"On this date, fatality rate was %0.1f %%. ",CFR(end)) ;	
		end
		fprintf(fid,"The effective reproduction number based on latest data is %0.1f. ", mean_posterior(end));
        today = now
		if(datenum(peakdate) < today-1)
		    fprintf(fid,"According to the model, the daily cases peaked on %s. ", peakdate) ;
		elseif (today-1  < datenum(peakdate)) & (datenum(peakdate)< datenum(datetime(2021,06,07))) 
		    fprintf(fid,"According to the model, the daily cases are expected to peak on %s. ", peakdate) ;
		else
            disp('error')	
		
		end
        %}
		dailyDmean = movmean(dailyD,7);
		casepeak = max(dailyCmean);
		deathpeak = max(dailyDmean);
		I11 = find((dailyCmean == casepeak));
		I12 = find((dailyDmean == deathpeak));
		fprintf(fid,'On %s, %s reported %d new cases and %d new deaths. The number of active cases is %d. Based on a 7-day average, this region reported a maximum of %d cases on %s, while a maximum of %d deaths occured on %s. ',datestr(dateend), Location, dailyC(end), dailyD(end), A(end), floor(dailyCmean(I11(1))), datestr(date(I11(1))), floor(dailyDmean(I12(1))), datestr(date(I12(1))));
		
		today = floor(now);
		Amean = movmean(A,7);
		activepeak = max(Amean);		
		I13 = find((Amean == activepeak));		
		if(date(I11) < today-1)
		    fprintf(fid,"The daily count of infections has reduced by %0.1f %% compared to its value at its peak. ", ((casepeak-dailyCmean(end))*100)/casepeak ) ;
			if(date(I13) < today-1)
				fprintf(fid,"The active case has also gone down by %0.1f %% from its peak value of %d on date %s. The current case fatality rate based on average daily data is %0.1f %%. ", ((activepeak-Amean(end))*100)/activepeak, floor(activepeak(1)), datestr(date(I13(1))), CFR(end) ) ;
			end
	    
		elseif(date(I13) < today-1)
		    fprintf(fid,"The current case fatality rate based on average daily data is %0.1f %%. ", CFR(end) ) ;
        end
        
		I14 = find((SIR.date >= today-1-0.05) & (SIR.date <= today-1+0.05))	;
		modelcases = diff(SIR_Cnew)/dt;
		realC = dailyCmean(end)
		modelC = floor(modelcases(I14(1)));
		modelcases1 = modelcases(I14(1):1/dt:end);
		%mc2 = floor(modelcases1(1:14))
		%datess = SIR.date(I14(1):1/dt:end);
        %datess1 = datess(1:14) 		
        size12 = size(modelcases1)		
		if(datenum(peakdate) < today-1)
		    fprintf(fid,"\n According to the model, the daily cases are currently in decline phase. The effective reproduction number based on the latest data is %0.1f. ", mean_posterior(end)) ;
			if size12(1) >= 13
			    fprintf(fid,"With the current decay rate, %s is expected to have about %d daily cases in two weeks. ", Location, floor(modelcases1(14-1))) ;    
		    end
		else
		    fprintf(fid,"\n According to the model, the daily cases are still in the growth phase and are expected to peak on %s. The effective reproduction number based on the latest data is %0.1f. ", peakdate,  mean_posterior(end)) ;
		end		
		
		if(dists== 0 && strcmp(Location,'India')==0)
		    fprintf(fid,"\n %s conducted %d of tests on %s. The test positivity rate (TPR) is %0.1f %%. World Health Organization (WHO) recommends TPR to be lower than 5 %% for at least two weeks before the outbreak can be considered as under control. ", Location, Testing(end), datestr(dateend), PR(end)) ;
			if(PR(end) > 5)
		        fprintf(fid,"Relatively high value indicates that it may take some time for the interventions to be eased. ") ;
			else
			    fprintf(fid," Low TPR value indicates that the interventions are expected to be eased soon. ") ;
			end	
		elseif(dists== 0 && strcmp(Location,'India')==1)
		    fprintf(fid,"\n %s conducted %d tests on %s. The test positivity rate (TPR) is %0.1f %%. World Health Organization (WHO) recommends TPR to be lower than 5 %% for at least two weeks before the outbreak can be considered as under control. ", Location, Testing(end-1), datestr(dateend), PR(end)) ;
			if(PR(end) > 5)
		        fprintf(fid,"Relatively high value indicates that it may take some time for the interventions to be eased. ") ;
			else
			    fprintf(fid," Low TPR value indicates that the interventions are expected to be eased soon. ") ;
			end			
		end		

		%SIR.date(I14:1/dt:I14+14)
		%datestr(today-1)
		%SIR.date(end-5:end)
		%date(end)
		I15 = find((date == datenum(datetime(2021,02,13))))	
		cumc = cumsum(dailyC(I15-1:end));
		cumd = cumsum(dailyD(I15-1:end));
		fprintf(fid,"\n In the second wave, a total of %d has been reported in %s, along with %d deaths. However, the current trend indicates that the situation is improving. People are advised to strictly follow COVID protocols and guidelines even after the interventions are eased in order to avoid any fresh wave of infections (third wave). ", cumc(end), Location, cumd(end)) ;

        fclose(fid);

        
    end
    

	
	
    prend = num2str(PR(end),'%.0f');
    cfrend = num2str(CFR(end),'%.2f');
    diary Karnataka_May10.txt
    disp([Location, ' ', peakdate, ' ', prend, ' ', cfrend])
    diary off
    
    
end
%delete(findall(0))


function names = getStateNames()
names = {

    %'Andhra Pradesh' %20/4 1.2
	%'Manipur'
	%'Jammu and Kashmir'	
	%'Karnataka'
    %'Kerala'
	%'Telangana'  %1/4 1.2
    %%'Assam'   %20/4 1.4
    %%'Mizoram'
    %%'Odisha'
    %%'Sikkim'
    %%'West Bengal'	
    %%'Arunachal Pradesh'	
	'Nagaland'	%20/4 1.4 ylimit on rt
    %{
    'India'	%1/4 1.1
		
    %'Andhra Pradesh'
    %'Arunachal Pradesh'
    %%'Assam'
    'Bihar'
    'Chhattisgarh'
    'Goa'
    'Gujarat'
    'Haryana'
    'Himachal Pradesh'
    %'Jammu Kashmir'
	
    'Jharkhand'
    %%'Karnataka'
    %%'Kerala'
    'Madhya Pradesh'
    'Maharashtra'
    %%'Manipur'
    'Meghalaya'
	
    %%'Mizoram'
    %%'Odisha
    'Rajasthan'
    %%'Sikkim'
    'Tamil Nadu'
    %%'Telangana'
    'Uttar Pradesh'
    'Uttarakhand'
    %%'West Bengal'
    'Delhi'
    'Puducherry'
	
	%'Ladakh'	
	%}
	%'Chandigarh' %20/4 1.1
    %'Punjab'	
    %'West Bengal'	
    };
end


function names = getDistrictNames()
names = {


'Bengaluru Urban'
% 'Pune'
% 'Mumbai'
% 'Thane'
% 'Nagpur'
% 'Ernakulam'
% 'Ahmedabad'
% 'Lucknow'
% 'Nashik'
% 'Kozhikode'
% 'Jaipur'
% 'Chennai'
% 'Gurugram'
% 'Kannur'
% 'Raipur'
% 'Ahmednagar'
% 'Chandrapur'
% 'Dehradun'
% 'Tumakuru'
% 'Jodhpur'
% 'Solapur'
% 'Surat'
% 'Patna'
% 'Chittoor'
% 'Sangli'
% 'Durg'
% 'Visakhapatnam'
% 'Latur'
% 'Ranchi'
% 'Ballari'
% 'Palghar'
% 'Kasaragod'
% 'Guntur'%
% 'Beed'%
% 'Indore'%
% 'Kanpur Nagar'%
% %'S.P.S. Nellore' 'Prakasam'
% 'Varanasi'%
% 'Prayagraj'%
% %%%%%%%%%'Mysuru'
% 'Pathanamthitta'
% 'Bhopal'
% 'Wayanad'
% 'Kalaburagi'
% 'Nanded'
% 'Bhandara'
% 'Raigad'
% 'Buldhana'
% 'Anantapur'
% 'Bengaluru Rural'
% 'Faridabad'
% 'Parbhani'
% 'Jalgaon'
% 'Meerut'
% 'Chengalpattu'
% 'Haridwar'
% %'S.A.S. Nagar'
% 'Ludhiana'
% 'West Godavari'
% 'Kangra'
% 'Yavatmal'
% 'Rajnandgaon'
% %'Udaipur'
% 'Raigarh'
% 'Kollam'
% 'Gondia'
% 'Srinagar'
% 'Alwar'
% 'Amravati'
% 'Vadodara'
% 'Sundargarh'
% 'Bhilwara'
% 'Janjgir Champa'
% 
% 'Gwalior'
% %'Y.S.R. Kadapa'
% 'Bareilly'
%  
% 'Gorakhpur'
% 
% 'Baloda Bazar'
% 
% 'Moradabad'
% 'Jalna'
% %%%%%%%%%%%'Kota'
% 'Hisar'
% 'Korba'
% 'Osmanabad'
% 'Vizianagaram'
% 'Gautam Buddha Nagar'
% 'Mandya'
% 'Gaya'
% 'Sikar'
% 'Thiruvallur'
% %'Dhule'
% 'Jhansi'
% 'Jammu'
% 'Saharanpur'
% 
% 'Wardha'
% 'Bathinda'
% 
% 'Nainital'
% 
% 'Shivamogga'
% 
% 'Jalandhar'
% 'Nandurbar'
% 'Bagalkote'
% 'Hooghly'
% 'Churu'
% 'Ghaziabad'
% 
% 'Yadgir'
% 
% 'Udham Singh Nagar'
% 'Sonipat'
% %%%%%%%%%%%%%%'Akola'
% 'Bikaner'
% 
% 'East Singhbhum'
% 
% 'Panipat'
% 
% 'Tiruchirappalli'
% 'Ganganagar'
% 
% 'Amritsar'
% 'Ajmer'
% 'Karnal'
% 
% 'Mungeli'
% 'Jamnagar'
% 'Jabalpur'
% 'Muzaffarnagar'
% 'Tirunelveli'
% 'Kodagu'
% 'Lakhimpur Kheri'
% 'Pauri Garhwal'
% 'Ghazipur'
% 'Sirsa'
% 'Rajkot'
% 
% 'Nalanda'
% 'Purnia'
% 'Saran'
% %'Bilaspur'
% 'Mandi'
% 'Shimla'
% 'Ambala'
% 'Yamunanagar'
% 
% 'Bokaro'
% 'Hazaribagh'
% 'Chamarajanagara'
% 'Ujjain'
% 
% %{
% %%%%%%%%%%% 20/5 1.3 'Thrissur'
% 
% 'Thiruvananthapuram' 'Palakkad' 'East Godavari'
% 
% 'Alappuzha'
% 
% %%%%%%%%%%%%%%%'Satara' %%%%%%%%%%%%%%'Kolhapur' 'Srikakulam' 'Idukki'
% 'Kottayam'
% 
% 'Hassan' 'Khordha' 'Kurnool' 'Dakshina Kannada'
% 
% 'Cuttack' 'Udupi'
% 
% 'Chikkaballapura'
% 
% 'Uttara Kannada' 'Koppal'
% 
% 'Purba Medinipur' %%%%%%%%%%%%%%%
% %}
% %{
% 'Malappuram'  % 20/5 1.1 'Kolkata' 'North 24 Parganas' 'South 24 Parganas'
% 'Howrah' 'Madurai' 'Raichur' 'Nadia' 'Paschim Bardhaman' 'Dharwad'
% %}
% 
% 'West Champaran'
% 'Mehsana'
% 'Koriya'
% 'Muzaffarpur'  
% 'Sawai Madhopur'
% 'Barmer'
% 'Kanyakumari'
% 'Mahasamund'
% 'Bhavnagar'
% 'Kabeerdham'
% 'Chikkamagaluru'
% 'Jaunpur'
% 'Begusarai'
% 'Birbhum'
% 'Erode'
% 'Pali'
% 'Jashpur'
% 'Uttar Bastar Kanker'
% %%%%%%%%%%%%%%%%%%'Sindhudurg' 'Tiruppur'
% 'Rajsamand'
% 'Surajpur'
% 'Mahendragarh'
% %%%%%%%%%%%%%%%%%%%%'Kolar'
% 'Thoothukkudi'
% 'Agra'
% %%%%%%%%%%%%%%%%%%%%'Bharatpur' 'Thanjavur'
% 'Kancheepuram'
% 'Dhamtari'
% 'Tehri Garhwal'
% 'Chittorgarh'
% 'Purba Bardhaman'
% 'Surguja'
% 'Patiala'
% 'Bametara'
% 'Amroha'
% 'Azamgarh'
% 'Krishnagiri'
% 'Angul'
% 'Ballia'
% 'Sambalpur'
% 'Fazilka'
% 'Washim'
% 'Vijayapura'
% 'Shahjahanpur'
% 'Bidar'
% 'Sultanpur'
% 'Ramanagara'
% %%%%%%%%%%%%%%%'Salem'
% 'Vaishali'
% 'Paschim Medinipur'
% 'Gadchiroli'
% 'Malda'
% 'Anantnag'
% 'Rae Bareli'
% 'Budgam'
% 'Solan'
% 'Baramulla'
% 'Bargarh'
% 'Baran'
% 'Darjeeling'
% 'Vellore'
% 'Gariaband'
% 'Jhunjhunu'
% 'Bhagalpur'
% 'Barabanki'
% 'Sirohi'
% %%%%%%%%%%%%%'Virudhunagar'
% 'Sonbhadra'
% 'Junagadh'
% %%%%%%%%%%%%%'Cuddalore'
% 'Samastipur'
% 'Tiruvannamalai'
% %%%%%%%%%%%%'Ratlam'
% 'Mathura'
% 'Kalahandi'
% 'Chandauli'
% 'Puri'
% 'Bhiwani'
% 'Murshidabad'
% 'Supaul'
% 'Jind'
% 'Budaun'
% 'Theni'
% 'Mansa'
% 'Unnao'
% 'Aligarh'
% 'Bijnor'
% 'Kulgam'
% 'Bulandshahr'
% 'Sri Muktsar Sahib'
% 'Jhalawar'
% 'Jajpur'
% 'Gaurela Pendra Marwahi'
% 'Lalitpur'
% 'Madhubani'
% 'Sirmaur'
% 'Katihar'
% 'Deoria'
% 'Davanagere'
% 'Ramgarh'
% 'Fatehabad'
% 'Rajouri'
% 'Mirzapur'
% 'Gandhinagar'
% %%%%%%%%%%%%%'Viluppuram'
% 'Latehar'
% 'Kutch'
% 'Jharsuguda'
% 'Hoshiarpur'
% 'Panchkula'
% %%%%%%%%%%%%%%'East Champaran'
% 'Nabarangapur'
% 'Hardoi'
% 'Bankura'
% 'Balasore'
% 'Saharsa'
% 'Ranipet'
% 'Balangir'
% 'Dholpur'
% 'Mainpuri'
% 'Pathankot'
% 'Ayodhya'
% 'Etah'
% 'Shamli'
% %%%%%%%%%%%%%%'Munger'
% 'Kurukshetra'
% 'Una'
% 'Pulwama'
% 'Nuapada'
% 'Sambhal'
% 'Hapur'
% 'Kushinagar'
% 'Gopalganj'
% 'Namakkal'
% 'Jalpaiguri'
% 'Chamoli'
% 'Rohtas'
% 'Siwan'
% 'Madhepura'
% 'Etawah'
% 'Pilibhit'
% 'Koderma'
% 'Dimapur'
% 'Shivpuri'
% 'Sitapur'


    };
end
