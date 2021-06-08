classdef IVB < handle
  %IVB SLC Current vs B data and lookup methods
  properties
    Types={'0.813q17.7' '1.625q27.3' '2.13q38.31' '1.625s29.2' '1.625s9.06' '2.375SX11.22'  '2.13q38.31' '1.625q27.3'} ;
    Poly=[-0.50684	4.7485	-0.00859	0.000138	3.69E-05	-2.41E-06	4.86E-08 ; % 0.813q17.7
      5.5156	4.1166	0.6208	-0.0534	0.002322	-4.97E-05	4.21E-07 ; % 1.625q27.3
      -2.4238	1.23E+01	-1.83E-01	1.36E-02	-4.47E-04	4.83E-06	2.24E-08 ; % 2.13q38.31
      0	0.049286	4.47E-05	-7.93E-08	6.94E-11	-2.91E-14	4.73E-18 ; % 1.625s29.2
      -1.6675	0.22135	-0.0004	2.60E-06	-7.37E-09	5.39E-12	1.04E-14 ; % 1.625s9.06
      -0.33711	0.4335	1.22E-03	-1.51E-05	8.20E-08	-2.03E-10	1.88E-13 ; % 2.375SX11.22
      -2.4238	1.23E+01	-1.83E-01	1.36E-02	-4.47E-04	4.83E-06	2.24E-08 ; % 2.13q38.31
      5.5156	4.1166	0.6208	-0.0534	0.002322	-4.97E-05	4.21E-07 ; % 1.625q27.3
      ] ;
  end
  
  methods
    function I = GetI(obj,Type,BDES)
      %GetI Get current for magnet Type with given BDES in T
      idat=find(ismember(obj.Types,Type),1);
      if isempty(idat)
        error('No data stored for type: %s',Type)
      end
      I = polyval(flip(obj.Poly(idat,:)),abs(BDES));
    end
    function BDES = GetB(obj,Type,I)
      %GETB Get integrated B field from current lookup
      idat=find(ismember(obj.Types,Type),1);
      if isempty(idat)
        error('No data stored for type: %s',Type)
      end
      P=flip(obj.Poly(idat,:));
      yi=linspace(0,500,10000);
      xi=polyval(P,yi);
      BDES=interp1(xi,yi,I);
    end
  end
end