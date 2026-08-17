classdef cacheProxy
%cacheProxy  Proxy del almacen cacheado de un objeto: su plano de CONTROL.
%
%   No tiene estado propio: envuelve a un objeto DUENO y reenvia toda la cadena
%   de subscripts a su despachador `cachedAccess`. El dueno devuelve una
%   cacheProxy en alguna de sus propiedades; llamando C a ese proxy y `foo` a una
%   clave cacheada, el llamante escribe:
%       C                     == imprime la tabla de definiciones y estados
%       C.foo                 == el valor (computa o hace replay si hace falta)
%       C.foo.campo           == indexar DENTRO del valor
%       C.foo.info            == informe de esa CP (estado, eventos, de quien
%                                sale, quien la usa, historial y coste)
%       C.foo.delete          == borrar el VALOR (la definicion queda)
%       o = C.foo.removeCP    == borrar definicion Y valor
%       o = C.foo.set( x )    == sembrar un valor a mano (aislado, COW)
%       C.foo.<evento>        == el handler de ese evento (o [] si invalida)
%   Las tres formas que devuelven `o` producen un dueno NUEVO (semantica de
%   valor): hay que reasignarlo.
%
%   CONTRATO DEL DUENO (duck typing: sirve cualquier clase que lo cumpla):
%     out = cachedAccess( owner , s )   recibe la cadena de subscripts SIN su
%                                       primer eslabon y devuelve un CELL de
%                                       salidas -- {} para las operaciones que no
%                                       devuelven nada, como .delete
%     displayCachedView( owner )        imprime la tabla
%     n = cachedNumArgs( owner , s , ctx )   OPCIONAL: cuantas salidas da esa
%                                       cadena (para no truncar las listas
%                                       separadas por comas). Si el dueno no lo
%                                       implementa se asume 1 -- que es lo que
%                                       habia --, asi que no rompe a nadie.
%   Todos pueden ser metodos Hidden: el proxy los invoca como funciones.
%
%   Esto es el plano de CONTROL. La lectura corriente no necesita el proxy: el
%   dueno puede (y conviene que lo haga) exponer un acceso directo mas barato.
%
% See also cacheHandle, cachedOwner.

  properties (Access = private)
    OWNER = []   % el objeto dueno (por valor; debe implementar cachedAccess y
                 % displayCachedView)
  end

  methods
    function obj = cacheProxy( M )
      if nargin, obj.OWNER = M; end
    end

    function varargout = subsref( obj , s )
      %la misma disciplina escalar que el dueno: sin esto, un array de proxies
      %([P,P] se puede formar, el proxy no bloquea concatenacion) moria en un
      %"Too many input arguments" criptico al leer obj.OWNER como cs-list
      if ~isscalar( obj )
        error('cacheProxy:notScalar', ...
              'cacheProxy es escalar: no forme arrays de proxies, obtenga uno por objeto (M.CP).');
      end
      if isempty( obj.OWNER )
        error('cacheProxy:orphan', ...
              'cacheProxy sin dueno: obtengalo del objeto (p.ej. M.CP), no lo construya a mano.');
      end
      out = cachedAccess( obj.OWNER , s );
      if isempty( out )
        varargout = {};
      else
        varargout = out( 1:max( min( nargout , numel(out) ) , 1 ) );
      end
    end

    function n = numArgumentsFromSubscript( obj , s , ctx )
      %igual que en cachedOwner: un statement puede devolver CERO salidas
      %(P.x.delete) y el default de MATLAB pediria 1.
      %Y si la cadena indexa DENTRO del valor, el numero lo decide ese valor
      %(P.foo.campo sobre un struct array es una lista separada por comas): se
      %le pregunta al DUENO, que es quien sabe distinguir una operacion del
      %plano de control de una indexacion. Es opcional en el contrato: un dueno
      %que no lo implemente (el motor propio de msh) se queda como estaba.
      if ctx == matlab.mixin.util.IndexingContext.Statement, n = 0; else, n = 1; end
      %se PRUEBA a llamarlo en vez de preguntar con ismethod: ismethod NO VE los
      %metodos Hidden -- comprobado: dice 0 hasta para cachedAccess, que es el
      %nucleo de este contrato --, la misma trampa que methods() en la guarda de
      %sombras de Define. Si el dueno no lo tiene, el catch deja el default.
      if isscalar( obj ) && ~isempty( obj.OWNER )
        try, n = cachedNumArgs( obj.OWNER , s , ctx ); catch, end
      end
    end

    function disp( obj )
      if ~isscalar( obj )
        fprintf('  (array %s de cacheProxy: use un proxy por objeto)\n' , mat2str(size(obj)) );
        return;
      end
      if isempty( obj.OWNER ), fprintf('  (cacheProxy sin dueno)\n'); return; end
      displayCachedView( obj.OWNER );
    end
  end
end
