classdef pRefresh < cachedOwner
%pRefresh  sonda: el gancho refreshDefs -- al cargar un .mat, el CODIGO gana.
%
%   Sin el gancho, un objeto cargado usa los compute y handlers DEL FICHERO
%   para siempre. pRefresh re-registra 'val' en su refreshDefs con el compute
%   "de ahora" (@() 2); el test (grupo 25 de test_cache) guarda un objeto cuya
%   'val' venia del "codigo viejo" (@() 1) y comprueba que tras load sirve la
%   del codigo -- y que una CP NO re-registrada ('extra') conserva la suya.
%
% See also test_cache, cachedOwner.
  properties (Constant, Hidden)
    CP_EVENTS = { 'chg' }
  end
  methods (Access = protected)
    function obj = refreshDefs( obj )
      %la definicion DEL CODIGO: lo que un .mat viejo traiga para 'val' se pisa
      obj = obj.Define( 'val' , @(~) 2 , 'chg' , [] );
    end
  end
end
