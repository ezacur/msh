classdef pTypo < cachedOwner
%pTypo  sonda: un fire() con el nombre de evento MAL ESCRITO
  properties (SetAccess = private)
    N = 0
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'bump' , 'general' }
  end
  methods
    function obj = pTypo()
      obj = obj.Define( 'dbl' , @(o) 2*o.N , 'bump' , @(v,o,d) v + 2*d );
    end
    function obj = BumpOK( obj , d )
      obj.N = obj.N + d;
      obj = obj.fire( { 'bump' } , d );          % correcto
    end
    function obj = BumpTypo( obj , d )
      obj.N = obj.N + d;
      obj = obj.fire( { 'bmup' } , d );          % ERRATA: 'bmup'
    end
  end
end
