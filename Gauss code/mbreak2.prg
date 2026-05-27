@
This version: May 13, 2008.
This program is distributed freely for non-profit academic purposes only. for other uses, please contact Zhongjun Qu (zqu at uiuc.edu).
A lot of effort has been put to construct this program and we would appreciate that
you acknowledge using this code in your research and cite the relevant paper on which it is based: 
Qu, Z. and P. Perron (2007): "Estimating and Testing for Structural Changes in Multivariate Regressions", Econometrica.
Although a lot of efforts have been put in constructing the program, we cannot be
held responsible for any consequences that could result from remaining errors.
Copyright: Zhongjun Qu and Pierre Perron (2007)
@


new;
library pgraph;
format /ld 6,3;
cls;

@select your output file@

output file = mbr-out.out reset; 

@load the data@
         
load data[150,5]=data1.dat;                  

@Set your T by n matrix of dependent variables@                                  
y=data[.,1:2];  

@Set your T by q matrix of regressors@  
                        
z=data[.,3:5];   

@Specify the number of observations@

T=rows(y);    

@Specify the matrix _S (here the first and second regressors are used in the 
first equation and the first and third regressors are used in the second equation)@                      
                                                                      
_S={1 0 0 0,                              
    0 1 0 0,                              
    0 0 0 0,  
    0 0 1 0, 
    0 0 0 0, 
    0 0 0 1};

@set the maximum number of breaks allowed. If use the WDmax test set m=5, since the critical values reported
correspond to this case@

m=2;

@If you have restrictions, specify the matrix R (here only the constants (first
regressor) is allowed to change in each equation). 
if you do not have restrictions, set R=eye(cols(_S)*(m+1))@
                                                                              
R={1 0 0 0 0 0 0 0,
0 1 0 0 0 0 0 0,                         
0 0 1 0 0 0 0 0,                        
0 0 0 1 0 0 0 0, 
0 0 0 0 1 0 0 0, 
0 1 0 0 0 0 0 0, 
0 0 0 0 0 1 0 0, 
0 0 0 1 0 0 0 0,
0 0 0 0 0 0 1 0,
0 1 0 0 0 0 0 0,
0 0 0 0 0 0 0 1,
0 0 0 1 0 0 0 0
};

        
@set the trimming parameter that specifies the minimal length of a segment as a proportion of the sample size.@

              
trm=0.2;                  

@set brv=1 when allowing breaks in the covariance matrix of the errors; otherwise, set brv=0@
                            

brv=0;   

@set brbeta=1 when allowing breaks in regression coefficients; otherwise, set brbeta=0.@                  
                            
brbeta=1;    

@set vauto=1 when applying a correction for serial correlation in the errors
(in which case Andrews' (1991) method is used to construct the robust 
covariance matrix.@ 
         
vauto=0;   

@set to prewhit=1 if you want to apply an AR(1) pre-whitening when estimating the 
robust covariance matrix (see Andrews and Monahan (1992).@
                
prewhit=0;                 

@Set hetq=1 if you want to allow the distribution of the regressors to change across
regimes (this is used only when constructing confidence intervals for the estimates 
of the break dates. For the construction of the test hetq=1 always.@

hetq=0;                    

/*Call the main procedure*/

call mainp(m,cols(_S),z,y,cols(y),trm,T,brv,brbeta,vauto,hetq,prewhit);

#include mbreak.src;

end;

