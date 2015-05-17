function chain3(filename)
% SOLUTION OF 3 COMPONENT DECAY
% Computer Project 1 of NPRE 247
% by Hao, Xiong
%   OCTOBER 1, 2013
%   UNIVERSITY OF ILLINOIS AT URBANA-CHAMPAIGN
%
%   filename is the input data file, containing values of
%     half life of A B
%     initial number of A B C
%     Delta t
%     final time
%     unit (optional, default 's')

% Read input file with the input data
[TAhalf,TBhalf,NA0,NB0,NC0,dt,tend,unit] = read(filename);
% Numerical solutions of the equations for 0 < t < t_final
[t_c,NA_c,NB_c,NC_c] = numerical_solution(TAhalf,TBhalf,NA0,NB0,NC0,dt,tend);
[dt_s,t_s,NA_s,NB_s,NC_s] = numerical_solution_stable(TAhalf,TBhalf,NA0,NB0,NC0,dt,tend);
dt_m = (dt+dt_s)/2;
[t_m,NA_m,NB_m,NC_m] = numerical_solution(TAhalf,TBhalf,NA0,NB0,NC0,dt_m,tend);
% Analytical solutions of the equations for 0 < t < t_final
t_a = t_s;
[NA_a,NB_a,NC_a] = analytical_solution(t_a,TAhalf,TBhalf,NA0,NB0,NC0);

% t_max vs. 1/dt
N     = 200;
dtmax = linspace(dt,dt_s,N);
tmax  = 1:N;
for j = 1:length(dtmax)
  [t,NA,NB,NC] = numerical_solution(TAhalf,TBhalf,NA0,NB0,NC0,dtmax(j),tend);
  [C,I] = max(NB);
  tmax(j) = t(I);
end

% t_max calculated from analytical solution
lambdaA1 = log(2)/TAhalf;
lambdaB1 = log(2)/TBhalf;
tmax_a = 1/(lambdaB1-lambdaA1)*log((lambdaA1*lambdaB1*NA0+lambdaA1*lambdaB1*NB0-lambdaB1^2*NB0)/(lambdaA1^2*NA0));

% Write data to output file
filename(end-3:end) = [];
write([filename,'_output_coarse.txt'],TAhalf,TBhalf,NA0,NB0,NC0,[dt dt],tend,unit,t_c,NA_c,NB_c,NC_c);
write([filename,'_output_medium.txt'],TAhalf,TBhalf,NA0,NB0,NC0,[dt dt_m],tend,unit,t_m,NA_m,NB_m,NC_m);
write([filename,'_output_reliable.txt'],TAhalf,TBhalf,NA0,NB0,NC0,[dt dt_s],tend,unit,t_s,NA_s,NB_s,NC_s);
write([filename,'_output_analytical.txt'],TAhalf,TBhalf,NA0,NB0,NC0,[dt dt_s],tend,unit,t_a,NA_a,NB_a,NC_a);
write_tmax([filename,'_output_tmax.txt'],TAhalf,TBhalf,NA0,NB0,NC0,dt,tend,unit,dtmax,tmax,tmax_a);

% Plot --------------------------------------------------------------------
figure(1);clf;set(gca,'FontSize',14);set(gcf,'Color','w');hold on;box on;
plot(t_c,NB_c,'r--',t_m,NB_m,'b-.',t_s,NB_s,'k-',t_a,NB_a,'g-');
xlabel(['t ( ',unit,' )']);
ylabel('N (number)');
legend(['\Delta t = ',num2str(dt),unit,' (coarse N_B)'],...
  ['\Delta t = ',sprintf('%5.5f',dt_m),unit,' (medium N_B)'],...
  ['\Delta t = ',sprintf('%5.5f',dt_s),unit,' (stable N_B)'],'analytical N_B');

figure(2);clf;set(gca,'FontSize',14);set(gcf,'Color','w');hold on;box on;
plot(t_s,NA_s,'r-',t_s,NB_s,'g--',t_s,NC_s,'b-.',t_s,NA_s+NB_s+NC_s,'k--');
xlabel(['t ( ',unit,' )']);
ylabel('N (number)');
legend('N_A','N_B','N_C','N_A+N_B+N_C');
title(['stable N_B, \Delta t = ',num2str(dt),unit]);

figure(3);clf;set(gca,'FontSize',14);set(gcf,'Color','w');hold on;box on;
plot(1./dtmax,tmax,'-');
xlabel(['1/\Delta t ( ',unit,'^{-1} )']);
ylabel(['t_{max} ( ',unit,' )']);
title(['t_{max} = ',num2str(tmax_a),unit]);
limx = [0 200];
xlim(limx);
plot([limx(1) limx(2)],[tmax_a,tmax_a],'r-');
legend('numerical t_{max}','analytical value');

% Read input file with the input data -------------------------------------
  function [TAhalf,TBhalf,NA0,NB0,NC0,dt,tend,unit] = read(filename)
    fni = fopen(filename);
    while ~feof(fni)
      str = fgetl(fni);
      if ~isempty(str) && str(1) ~= '%'
        str = strsplit(str);
        TAhalf  = str2double(str{1});   % half life of A
        TBhalf  = str2double(str{2});   % half life of B
        NA0     = str2double(str{3});   % initial number of A
        NB0     = str2double(str{4});   % initial number of B
        NC0     = str2double(str{5});   % initial number of C
        dt      = str2double(str{6});   % initial Delta t
        tend    = str2double(str{7});   % final time
        % The 8th value is the unit used, default is 's'.
        if length(str) > 7
          unit = str{8};
        else
          unit = 's';
        end
        break;
      end
    end
    fclose(fni);
  end

% Calculates the number according to the given Delta t --------------------
  function [t,NA,NB,NC] = numerical_solution(TAhalf,TBhalf,NA0,NB0,NC0,dt,tend)
    lambdaA = log(2)/TAhalf;        % decay constant of A
    lambdaB = log(2)/TBhalf;        % decay constant of B
    t = 0:dt:tend;
    NA = t;NA(1) = NA0;
    NB = t;NB(1) = NB0;
    NC = t;NC(1) = NC0;
    
    for i = 1:length(t)-1
      NA(i+1) = (1-lambdaA*dt)*NA(i);
      NB(i+1) = lambdaA*dt*NA(i)+(1-lambdaB*dt)*NB(i);
      NC(i+1) = lambdaB*dt*NB(i)+NC(i);
    end
  end % cal

% Reliable numerical solution ---------------------------------------------
  function [dt_s,t,NA,NB,NC] = numerical_solution_stable(TAhalf,TBhalf,NA0,NB0,NC0,dt,tend)
    [t,NAt,NBt,NCt] = numerical_solution(TAhalf,TBhalf,NA0,NB0,NC0,dt,tend);
    dt_s = dt;
    error = 0.01; % value to check if solution is stable
    while 1
      dt_s = dt_s/2;
      [t,NA,NB,NC] = numerical_solution(TAhalf,TBhalf,NA0,NB0,NC0,dt_s,tend);
      buf = 1:length(NAt); % compare the corresponding results of dt and dt/2
      buf = 2*buf-1;
      % Check if solution is stable, if not, make dt half and calculate again
      if max([abs(max(NAt-NA(buf))),abs(max(NBt-NB(buf))),abs(max(NCt-NC(buf)))]) <= error
        break;
      else
        NAt = NA;
        NBt = NB;
        NCt = NC;
      end
    end
  end % numerical_solution

% analytical  solutions of the equations for 0 < t < t_final --------------
  function [NA,NB,NC] = analytical_solution(t,TAhalf,TBhalf,NA0,NB0,NC0)
    lambdaA = log(2)/TAhalf;        % decay constant of A
    lambdaB = log(2)/TBhalf;        % decay constant of B
    NA = NA0*exp(-lambdaA*t);
    NB = NB0*exp(-lambdaB*t)+lambdaA*NA0/(lambdaB-lambdaA)*(exp(-lambdaA*t)-exp(-lambdaB*t));
    NC = NC0+NB0*(1-exp(-lambdaB*t))+NA0/(lambdaB-lambdaA)*(lambdaB*(1-exp(-lambdaA*t))-lambdaA*(1-exp(-lambdaB*t)));
  end

% Output ------------------------------------------------------------------
  function write(filename,TAhalf,TBhalf,NA0,NB0,NC0,dt,tend,unit,t,NA,NB,NC)
    fno = fopen(filename,'w');
    fprintf(fno,'%% First line:\r\n');
    fprintf(fno,'%%     half life of A B\r\n');
    fprintf(fno,'%%     initial number of A B C\r\n');
    fprintf(fno,'%%     Delta t\r\n');
    fprintf(fno,'%%     final time\r\n');
    fprintf(fno,'%%     unit (optional, default ''s''\r\n');
    fprintf(fno,'%% Second Line: ''Delta t used'', ''data length''\r\n');
    fprintf(fno,'%% Third Line:  data\r\n');
    fprintf(fno,'%f %f %f %f %f %f %f ',TAhalf,TBhalf,NA0,NB0,NC0,dt(1),tend);
    fprintf(fno,'%s\r\n',unit);
    fprintf(fno,'%f %f\r\n',dt(2),length(NA));
    fprintf(fno,'%f\t%f\t%f\t%f\r\n',[t;NA;NB;NC]);
    fclose(fno);
  end

  function write_tmax(filename,TAhalf,TBhalf,NA0,NB0,NC0,dt,tend,unit,dtmax,tmax,tmax_a)
    fno = fopen(filename,'w');
    fprintf(fno,'%f %f %f %f %f %f %f ',TAhalf,TBhalf,NA0,NB0,NC0,dt,tend);
    fprintf(fno,'%s\r\n',unit);
    fprintf(fno,'%f %f\r\n',tmax_a,length(dtmax));
    fprintf(fno,'%f\t%f\r\n',[dtmax;tmax]);
    fclose(fno);
  end
end