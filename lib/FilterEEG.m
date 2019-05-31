function [output] = FilterEEG(EEG, cutoff, sample_rate, high_or_low, order)
%simple filter
  ny=sample_rate/2;

  [b,a] = butter(order, cutoff/ny, high_or_low);


  EEG_f = filter(b,a,EEG);

  EEG_f_reverse = EEG_f(end:-1:1, :);

  EEG_f_reverse_f = filter(b,a,EEG_f_reverse);

  output = EEG_f_reverse_f(end:-1:1, :);

end