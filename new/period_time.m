function T=period_time(s,interval)

    l=length(s);
    extremum=[];
    for i=2:l-1
       if s(i)>=s(i-1) && s(i)>s(i+1)
          extremum=[extremum,i];
       end
    end
    T=median(extremum(2:end)-extremum(1:end-1))*interval;

    