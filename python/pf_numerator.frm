#--
*off statistics;
CF ncmd, kval, num, a;


Auto S invd, E, p;
Auto S k;
Auto S f,z,r,x,s,c;
*#include- pf_2l_pentabox.frm
#include- pf_4l_2x2_fishnet.frm
*#include- pf_2l_sunrise.frm

#$LOOPS=4;

* numerator function in LMB
*#$numF = k0^2*k1^2+c1*k0^3*k1+c2*k0^2*k1^3;
#$numF = 1;

* Multiply numerator functions
Multiply $numF;
*print;
.sort:start;

* Start unfoding the numerator instructions loop by loop
#do i = 0,{`$LOOPS'-1}
	id once num(ncmd(?z),?x) =ncmd(?z,0)*num(?x);
	B+ ncmd, k`i';
	.sort:collect-ncmd;
	keep brackets;
	
	id k`i'^r?*ncmd(?z,z1?,0) = sum_(s,0, r - nargs_(?z) + 0, a(r,s,?z)*z1^s);
	B+ a;
	.sort:collect-a;
	keep brackets;

	repeat id a(r?,s?,?z,z1?) = -sum_(k, s+1,r-nargs_(?z) + 0, a(r,k,?z)*z1^(k-s-1));
	id a(r?,s?) = delta_(r,s);
	.sort:energy-k`i';
#enddo
id num=1;

* check that the substitution is complete
if (count(ncmd, 1,num,1,a, 1));
    Print "Unsubstituted ncmd: %t";
    exit "Critical error";
endif;
.sort:rm-num;

*Format 255;
Format C;
Format O4,stats=on, saIter=1000;
#Optimize F
#write<pf_form.out> "%O"
#write<pf_form.out> "      F = %e",F
*print +s;
.end


