classdef pDiamond < cachedOwner
%pDiamond  sonda: DIAMANTE de dependencias.
%
%       base            top sale de left Y de right, y las dos salen de base.
%      /    \           Es el caso que el cierre transitivo dice que funciona
%   left    right       y que ningun ejemplo recorria: hasta ahora las deps se
%      \    /           probaban en CADENA (base -> mid -> top), no confluyendo.
%       top
%
%   Con N=1: base=10, left=11, right=20, top=31.
%   Al disparar 'chg' cae base, y con ella left, right Y top.
%
%   Sirve tambien para el grupo de save/load: al cargar, la cache llega vacia
%   (es Transient) y las dependencias se REDESCUBREN en el primer recalculo --
%   que es lo que hay que comprobar, porque las deps NO viajan en el .mat.
%
% See also test_cache, pComp, pIncDep, cachedOwner.
  properties (SetAccess = private)
    N = 1
  end
  properties (Constant, Hidden)
    CP_EVENTS = { 'chg' }
  end
  methods
    function obj = pDiamond()
      obj.DEBUG = false;
      obj = obj.Define( 'base'  , @(o) 10 * o.N            , 'chg' , [] );
      obj = obj.Define( 'left'  , @(o) o.access('base') + 1 );
      obj = obj.Define( 'right' , @(o) o.access('base') * 2 );
      obj = obj.Define( 'top'   , @(o) o.access('left') + o.access('right') );
    end
    function obj = Bump( obj )
      obj.N = obj.N + 1;
      obj = obj.fire( 'chg' );
    end
  end
end
