! Amir Shahmoradi, Wednesday 5:43 PM, December 25, 2013, IFS, UT Austin
function delayed_rate_Belz_Li(z)

	implicit none
	real*8, intent(in) :: z
	real*8 :: delayed_rate_Belz_Li
	
	if (z>2.5d0 .and. z<=6.501d0) then
	
		delayed_rate_Belz_Li = -2.09118024744342 + 5.15382361299299*z - 5.46442271664195*z**2 + 3.29445310883082*z**3 - &
		1.24547016168265*z**4 + 0.306288936905084*z**5 - 0.0490440324964182*z**6 + 0.00493757380504717*z**7 - &
		2.8406197192875d-4*z**8 + 7.1267413875775d-6*z**9
			
	elseif (z>1.0d0 .and. z<=2.5d0) then
	
		delayed_rate_Belz_Li = -0.860225762659041 + 4.22669545558817*z - 8.8608672853467*z**2 + &
		10.4863792284648*z**3 - 7.64722909221129*z**4 + 3.51616699500767*z**5 - 0.99555474471022*z**6 + &
		0.158768937543719*z**7 - 0.0109254199773642*z**8
	
	elseif (z<=1.d0 .and. z>0.09d0) then
	
		delayed_rate_Belz_Li = 1.9259529998937d-4 - 0.00345273599582578*z + 0.0315750061532092*z**2 - &
		0.0447054552119846*z**3 + 0.0681248152128166*z**4 - 0.0384603341625357*z**5
	
	else
	
		write(*,*) 'Inappropriate redshift range of integration for delayed_rate_Belz_Li: ',z
		stop
	
	end if

end function delayed_rate_Belz_Li