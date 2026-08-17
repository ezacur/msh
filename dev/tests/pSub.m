classdef pSub < cachedOwner
%pSub  sonda: una subclase que SOBREESCRIBE subsref y llama a los internos,
%      que es el camino de extension que el tutorial recomienda para msh
  properties (SetAccess = private)
    N = 7
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'chg' }
  end
  methods
    function obj = pSub()
      %el 'chg',[] no es decorativo: dbl LEE datos (o.N), y una CP que lee datos
      %sin declarar ningun evento sobrevive a todas las ediciones -- Define avisa
      %(cachedOwner:noEvents). Esta sonda no edita nunca, pero se declara bien
      obj = obj.Define( 'dbl' , @(o) 2*o.N , 'chg' , [] );
    end
    function varargout = subsref( obj , s )
      %imita lo que haria msh: resolver la CP llamando al motor de la base
      if strcmp( s(1).type , '.' ) && isscalar( s ) && ischar( s(1).subs ) && ...
         strcmp( s(1).subs , 'dbl' )
        varargout = { obj.access( 'dbl' ) };     % <-- ¿puede la subclase?
        return;
      end
      [ varargout{ 1:max(nargout,1) } ] = subsref@cachedOwner( obj , s );
    end
  end
end
