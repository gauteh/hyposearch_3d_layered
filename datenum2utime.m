function u = datenum2utime (d, print)
  if (~exist('print', 'var'))
    print = false;
  end
  
  u = 86400 * (d - datenum('01-Jan-1970'));
  
  if (print)
    fprintf ('%f\n', u);
  end
end
