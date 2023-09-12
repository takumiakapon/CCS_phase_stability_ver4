module mod_injection_data
    use mod_condition
    use mod_cubic_eq_d
    use mod_autodiff
    use mod_input

    implicit none

    contains

    subroutine injection_data_d(pp,MD_injection_d)
        implicit none
        type(diffs),intent(in)::pp
        type(diffs),intent(out)::MD_injection_d

        integer::i,j
        real(8)::a_mix_V,b_mix_V
        real(8),dimension(com_2phase)::Tc,Pc,acentric,a,b,x,y
        real(8),dimension(com_2phase,com_2phase)::bic
        type(diffs)::A_V,B_V,z_factor,a_coef(3)

        Tc(1) = 647.30d0
        Tc(2) = 304.2d0
        Tc(3) = 190.6d0
        Tc(4) = 373.2d0
        Pc(1) = 217.6d0 * 101325.d0
        Pc(2) = 72.8d0 * 101325.d0
        Pc(3) = 45.4d0 * 101325.d0
        Pc(4) = 88.2d0 * 101325.d0
        acentric(1) = 0.344d0
        acentric(2) = 0.225d0
        acentric(3) = 0.008d0
        acentric(4) = 0.1d0
        bic(1,1)=0.0d0
        bic(1,2)=0.2d0
        bic(1,3)=0.4907d0
        bic(1,4)=0.12d0
        bic(2,1)=bic(1,2)
        bic(2,2)=0.0d0
        bic(2,3)=0.105d0
        bic(2,4)=0.135d0
        bic(3,1)=bic(1,3)
        bic(3,2)=bic(2,3)
        bic(3,3)=0.0d0
        bic(3,4)=0.07d0
        bic(4,1)=bic(1,4)
        bic(4,2)=bic(2,4)
        bic(4,3)=bic(3,4)
        bic(4,4)=0.0d0

        !とりあえず気相でのCO2のみの圧入
        x=0.0d0
        y=0.0d0
        y(2)=1.0d0

        do j=1,com_2phase
            a(j) = (0.457235528921d0*(R**2.0d0)*(tc(j)**2.0d0)/pc(j)) * ((1.0d0+(0.37464d0+1.54226d0*acentric(j)-0.26992*&
            (acentric(j)**2.0d0))*(1.0d0-sqrt(temp/tc(j))) ))**2.0d0
            b(j) = 0.0777960739039d0*R*tc(j)/pc(j)
            b_mix_V=b_mix_V+y(j)*b(j)
        end do
        !write(*,*) a(4)
    
    
        do j=1,com_2phase
            do i=1,com_2phase
                a_mix_V=a_mix_V+y(j)*y(i)*(1.0d0-bic(j,i))*sqrt(a(j)*a(i))
            end do
        end do

        A_V          = (pp/((R*temp)**2.0d0)) * a_mix_V!!
        B_V          = (pp/(R*temp))        * b_mix_V!!
        a_coef(1)       =(B_V-1.0d0)!!
        a_coef(2)       = A_V - 2.0d0*B_V - 3.0d0*B_V**2.0d0!!
        a_coef(3)       =(B_V**2.0d0 + B_V**3.0d0-A_V*B_V)!!

        call cubic_eq_d(a_coef(1),a_coef(2),a_coef(3),z_factor)

        MD_injection_d=pp/z_factor/R/temp

    end subroutine
end module
