classdef pInner < cachedOwner
%pInner  sonda: una CP cuyo VALOR es otro objeto con subsref sobrecargado
  properties (Constant, Hidden)
    CP_EVENTS = { 'chg' }
  end
  methods
    function obj = pInner()
      obj = obj.Define( 'inner' , @(o) cart().Demo() );   % el valor es un cart
    end
  end
end
