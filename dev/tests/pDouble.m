classdef pDouble < cachedOwner
%pDouble  sonda: DOS handlers incrementales en la misma CP
  properties (SetAccess = private)
    N = 0
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'up' , 'down' }
  end
  methods
    function obj = pDouble()
      obj = obj.Define( 'val' , @(o) o.N , ...
                        'up'   , @(v,o,d) v + d , ...   % incremental #1
                        'down' , @(v,o,d) v - d );      % incremental #2 (¿se usa?)
    end
    function obj = Up( obj , d ),   obj.N = obj.N + d; obj = obj.fire({'up'},d);   end
    function obj = Down( obj , d ), obj.N = obj.N - d; obj = obj.fire({'down'},d); end
  end
end
