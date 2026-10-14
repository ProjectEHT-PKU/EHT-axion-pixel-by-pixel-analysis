You can set these parameters(l=1 m=1 a=0.5,0.8,0.9,0.9375,0.95,0.99   or   l=2 m=2 a=0.95,0.99) in the Beginning of main.c(for a, Mu_min, Mu_max) or in the beginning of Axion_Radial.c=h(for l,m,FALG_L_M,LAMBDA_L_M)

IF YOU WANT TO CALCULATE ONE SPECIFIC RADIAL_FUNCTION: 

	For l=1 m=1 case:	First, you need set Mu_min in Axion_Radial.c line 214 : Mu_min = 0.1

		a = 0.99   Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.446
		a = 0.95   Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.367
		a = 0.9375 Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.316
		a = 0.9    Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.316
		a = 0.8    Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.249
		a = 0.5    Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.132

	For l=2 m=2 case:	First, you need set Mu_min in Axion_Radial.c line 214 : Mu_min = 0.320

		a = 0.95   Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.752         
		a = 0.99   Set Mu_max in Axion_Radial.c line 215 : Mu_max = 0.928

IF YOU WANT TO CALCULATE ALL THE POSSIBLE RADIAL_FUNCTION:

	Unnote the print_all_cable function in main.c line 19. Run the program and change the l,m,flag_l_m,lambda_l_m in 
	Radial_funtion.h, run again.

Some notices :
I used Interpolation.nb to do the interpolation work, used the Plot_Radial_funtion.nb to see the differences between different radial functions.

The program can output *.txt file since C become complicated when you want to output *.csv. So you need to change the filename manually to *.csv since Mathematica cannot recieve *.txt file properly.

In the l=1,m=1 case, a=0.95's real part of \omega is an interpolation of a=0.99 and a=0.9
In the l=1,m=1 case, a=0.9375's real and image part of \omega is an interpolation of a=0.95 and a=0.9
In the l=1,m=1 case, a=0.8's real part of \omega is an interpolation of a=0.7 and a=0.9
In the l=2,m=2 case, we use all the real data from l=1,m=1 case, since there is no data from Dolan.