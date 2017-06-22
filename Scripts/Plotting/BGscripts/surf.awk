function asin(x) { return atan2(x, sqrt(1-x*x)) }
function acos(x) { return atan2(sqrt(1-x*x), x) }
function atan(x) { return atan2(x,1) }
BEGIN{
  xref = 0.0 ;
  pi   = acos(-1.0) ;
  d120 = 2.0*pi/3.0 ;
  d60  =     pi/3.0 ;
  xf   = sqrt(5.0/(16.0*pi)) ;
  gf   = pi/180.0 ;
}
{
  #mass = nn + zz ;
  xe   = 3.0*1.2*1.2*(mass^(5.0/3.0)) ;
   
  if ( $1 != "!" && NF > 0 ) { 
    q1    = $1 ; 
    q2    = $2 ;
    q0    = sqrt(q1*q1 + q2*q2 + q1*q2) ;
    if ( q1 < 0.1 && q2 < 0.1 ) {
      gamma = 0.0 ;
    }
    else {
      gamma = atan((q2*sqrt(3.0))/(2.0*q1 + q2)) ;
    }
    ee    = $3 ;
    iqq   = q0 ;
    iqb   = xf*((4.0*pi*iqq)/(xe));
    iqg   = gamma ;   
    
                                           printf("%12.6f %12.6f %12.3f\n",iqb*cos(       iqg), iqb*sin(       iqg), ee-xref) >> "ESPDATA1.bg";
                                           printf("%12.6f %12.6f %12.3f\n",iqb*cos(       iqg), iqb*sin(       iqg), ee-xref) >> "ESP.xflr6"   ;
    if ( iqq > 0.01 && iqg < 0.999*d60 ) { printf("%12.6f %12.6f %12.3f\n",iqb*cos( d120 -iqg), iqb*sin( d120 -iqg), ee-xref) >> "ESPDATA2.bg"; }
    if ( iqq > 0.01 && iqg > 0.001     ) { printf("%12.6f %12.6f %12.3f\n",iqb*cos( d120 +iqg), iqb*sin( d120 +iqg), ee-xref) >> "ESPDATA3.bg"; }
    if ( iqq > 0.01 && iqg < 0.999*d60 ) { printf("%12.6f %12.6f %12.3f\n",iqb*cos(-d120 -iqg), iqb*sin(-d120 -iqg), ee-xref) >> "ESPDATA4.bg"; }
    if ( iqq > 0.01 && iqg > 0.001     ) { printf("%12.6f %12.6f %12.3f\n",iqb*cos(-d120 +iqg), iqb*sin(-d120 +iqg), ee-xref) >> "ESPDATA5.bg"; }
    if ( iqq > 0.01 && iqg > 0.001     && iqg < 0.999*d60 ) 
    { 
                                             printf("%12.6f %12.6f %12.3f\n",iqb*cos(      -iqg), iqb*sin(      -iqg), ee-xref) >> "ESPDATA6.bg"; 
    }
  }
}
