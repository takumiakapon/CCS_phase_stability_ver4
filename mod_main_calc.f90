module mod_main_calc
    use mod_autodiff
    use mod_condition
    use mod_cubic_eq_d
    use mod_input
    use mod_fugacity
    use mod_injection_data

    implicit none

    integer,private::i,j,jj

    contains
    subroutine main_calc(V0,lnk0,Nc0,Nc0old,Nm0,Nm0old,Nmini,P0,P0old,Pb0,fai00,fai000,q_judge,phase_judge,phase,&
                        Swd,krgd,krwd,Sw00,wc0,chemi_mat,theta0,g)
        implicit none
        integer,intent(inout)::q_judge,phase_judge(n),phase(n)
        real(8),intent(inout)::Pb0,Nmini(com_mine),chemi_mat(chemi+mine,com_all),theta0(chemi+mine,n)
        real(8),intent(inout),dimension(n)::V0,P0,P0old,fai00,Sw00,fai000
        real(8),intent(inout),dimension(com_2phase,n)::lnk0
        real(8),intent(inout),dimension(com_2phase+com_ion,n)::Nc0,Nc0old,wc0
        real(8),intent(inout),dimension(com_mine,n)::Nm0,Nm0old
        real(8),intent(in),dimension(21)::Swd,krgd,krwd
        type(diffs),allocatable,intent(out)::g(:)
        type(diffs),allocatable,target::xd(:)
        type(diffs),pointer::P(:),V(:),lnk1(:),lnk2(:),lnk3(:),lnk4(:),Nc1(:),Nc2(:),Nc3(:),Nc4(:),Nc5(:),Nc6(:),Nc7(:)
        type(diffs),pointer::Nc8(:),Nc9(:),Nc10(:),Nc11(:),Nc12(:),Nc13(:),Nc14(:),Nm1(:),Nm2(:),Nm3(:),Nm4(:),Nm5(:),Pb
        real(8),allocatable,dimension(:)::x0,kakuninn
        real(8)::kaku

        !!
        real(8)::z_factor0(n),v_L2(n),v_L3(n),v_L4(n),Mw(com_2phase+com_ion)
        type(diffs)::lnk(com_2phase,n),Nc(com_2phase+com_ion,n),Nm(com_mine,n),Nt(n),z(com_2phase+com_ion,n),rach(n)
        type(diffs)::kakuninnnnnn(com_2phase+com_ion),x(com_2phase+com_ion,n),y(com_2phase,n),k(com_2phase,n)
        type(diffs)::lnfai_L(com_2phase,n),lnfai_V(com_2phase,n),z_factor(n),v_L1(n),MV_V(n),MD_V(n)
        type(diffs)::v_L(com_2phase+com_ion,n),MV_L(n),MD_L(n),Sw(n),Sg(n),L(n),Mw_ave_L(n),Mw_ave_V(n)
        type(diffs)::phase_d_L(n),phase_d_V(n)

        real(8)::beta_v,myu_H2O_20,sa,sisuu,myu_L_normal,Vc(com_2phase),Tc(com_2phase),Pc(com_2phase),myu_c(com_2phase)
        real(8)::av_0,av_1,av_2,av_3,av_4
        type(diffs)::myu_L(n),Pc_ave_V(n),Tc_ave_V(n),Vc_ave_V(n),zeta(n),ro_r_V(n),myu_V_SC(n),myu_V(n),krg(n),krw(n)
        real(8),allocatable,dimension(:)::Sw0,fai0
        type(diffs)::faimine(n),fai(n),kk(n),ramda
        real(8)::faimineold(n),faiold(n),NmMD(com_mine)
        
        real(8)::re,W(n)
        type(diffs)::q,MD_injection_d,T_L(com_2phase+com_ion,n),T_V(com_2phase,n)

        !化学反応
        real(8)::ks(chemi+mine),ke(chemi),km(mine),min0
        real(8),allocatable,dimension(:)::w10,w20,w30,w40,w50,w60,w70,w80,w90,w100,w110,w120,w130,w140
        !real(8),allocatable,dimension(:)::Q10,Q20,Q30,Q40,Q50,Q60,Q70,Q80,Q90,Q100
        real(8),allocatable,dimension(:)::theta10,theta20,theta30,theta40,theta50,theta60,theta70,theta80,theta90,theta100
        type(diffs)::wc(com_2phase+com_ion,n),Ac(com_mine,n),Q_chemi(chemi+mine,n),theta(chemi+mine,n),min
        type(diffs),dimension(n)::w1,w2,w3,w4,w5,w6,w7,w8,w9,wten,w11,w12,w13,w14
        !type(diffs),dimension(n)::Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Qten
        type(diffs),dimension(n)::theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8,theta9,thetaten
        type(diffs)::rs(chemi+mine,n),rs_sum(com_all,n),rs1(n)
        real(8)::AA(mine)

        allocate(x0(n*eq+q_judge))

        !!自動微分の下準備
        do i=1,n
            !?相平衡
            x0(i*eq-24)=lnk0(1,i)
            x0(i*eq-23)=lnk0(2,i)
            x0(i*eq-22)=lnk0(3,i)
            x0(i*eq-21)=lnk0(4,i)
            !?物質収支
            x0(i*eq-20)=Nc0(1,i)
            x0(i*eq-19)=Nc0(2,i)
            x0(i*eq-18)=Nc0(3,i)
            x0(i*eq-17)=Nc0(4,i)
            x0(i*eq-16)=Nc0(5,i)
            x0(i*eq-15)=Nc0(6,i)
            x0(i*eq-14)=Nc0(7,i)
            x0(i*eq-13)=Nc0(8,i)
            x0(i*eq-12)=Nc0(9,i)
            x0(i*eq-11)=Nc0(10,i)
            x0(i*eq-10)=Nc0(11,i)
            x0(i*eq-9)=Nc0(12,i)
            x0(i*eq-8)=Nc0(13,i)
            x0(i*eq-7)=Nc0(14,i)
            x0(i*eq-6)=Nm0(1,i)
            x0(i*eq-5)=Nm0(2,i)
            x0(i*eq-4)=Nm0(3,i)
            x0(i*eq-3)=Nm0(4,i)
            x0(i*eq-2)=Nm0(5,i)
            !?飽和率の制約式
            x0(i*eq-1)=P0(i)
            !?rachford-rice
            x0(i*eq)=V0(i)
        end do
        if (q_judge == 1) then
            x0(n*eq+q_judge)=Pb0
        end if

        call diffsset1(x0,xd)
        call sizeset(x0,g)

        lnk1 => xd(eq-24:n*eq-24:eq)
        lnk2 => xd(eq-23:n*eq-23:eq)
        lnk3 => xd(eq-22:n*eq-22:eq)
        lnk4 => xd(eq-21:n*eq-21:eq)
        Nc1 => xd(eq-20:n*eq-20:eq)
        Nc2 => xd(eq-19:n*eq-19:eq)
        Nc3 => xd(eq-18:n*eq-18:eq)
        Nc4 => xd(eq-17:n*eq-17:eq)
        Nc5 => xd(eq-16:n*eq-16:eq)
        Nc6 => xd(eq-15:n*eq-15:eq)
        Nc7 => xd(eq-14:n*eq-14:eq)
        Nc8 => xd(eq-13:n*eq-13:eq)
        Nc9 => xd(eq-12:n*eq-12:eq)
        Nc10 => xd(eq-11:n*eq-11:eq)
        Nc11 => xd(eq-10:n*eq-10:eq)
        Nc12 => xd(eq-9:n*eq-9:eq)
        Nc13 => xd(eq-8:n*eq-8:eq)
        Nc14 => xd(eq-7:n*eq-7:eq)
        Nm1 => xd(eq-6:n*eq-6:eq)
        Nm2 => xd(eq-5:n*eq-5:eq)
        Nm3 => xd(eq-4:n*eq-4:eq)
        Nm4 => xd(eq-3:n*eq-3:eq)
        Nm5 => xd(eq-2:n*eq-2:eq)
        P => xd(eq-1:n*eq-1:eq)
        V => xd(eq:n*eq:eq)
        if (q_judge == 1) then
            Pb => xd(n*eq+q_judge)
        end if
        !?-------------------------------

        do i=1,n
            lnk(1,i)=lnk1(i)
            lnk(2,i)=lnk2(i)
            lnk(3,i)=lnk3(i)
            lnk(4,i)=lnk4(i)
            Nc(1,i)=Nc1(i)
            Nc(2,i)=Nc2(i)
            Nc(3,i)=Nc3(i)
            Nc(4,i)=Nc4(i)
            Nc(5,i)=Nc5(i)
            Nc(6,i)=Nc6(i)
            Nc(7,i)=Nc7(i)
            Nc(8,i)=Nc8(i)
            Nc(9,i)=Nc9(i)
            Nc(10,i)=Nc10(i)
            Nc(11,i)=Nc11(i)
            Nc(12,i)=Nc12(i)
            Nc(13,i)=Nc13(i)
            Nc(14,i)=Nc14(i)
            Nm(1,i)=Nm1(i)
            Nm(2,i)=Nm2(i)
            Nm(3,i)=Nm3(i)
            Nm(4,i)=Nm4(i)
            Nm(5,i)=Nm5(i)
            do j=1,com_2phase
                k(j,i) = exp(lnk(j,i))
            end do
        end do

        !!rachford-rice
        !#TODO相の数で判断するか否か
        do i=1,n
            call residualvectorset3(n*eq+q_judge,Nt(i))
            do j=1,com_2phase+com_ion
                Nt(i) = Nt(i) + Nc(j,i)
            end do
            do j=1,com_2phase+com_ion
                z(j,i) = Nc(j,i) /Nt(i)
            end do
            call residualvectorset3(n*eq+q_judge,rach(i))
            do j=1,com_2phase
                rach(i) = rach(i) +(1.0d0-k(j,i))*z(j,i)/(1.0d0-V(i)+V(i)*k(j,i))
            end do
            do j=com_2phase+1,com_2phase+com_ion
                rach(i) = rach(i) +z(j,i)/(1.0d0-V(i))
            end do
        end do
        !do j=1,com_2phase+com_ion
        !    kakuninnnnnn(j) = z(j,1)
        !end do
        !call outxs(kakuninnnnnn,kakuninn)
        !write(*,*) kakuninn

        !!モル分率----------------------
        do i=1,n
            if (phase_judge(i) == 2) then
                do j=1,com_2phase
                    x(j,i) = z(j,i)/(1.0d0-V(i)+V(i)*k(j,i))
                    y(j,i) = x(j,i)*k(j,i)
                end do
                do j=com_2phase+1,com_2phase+com_ion
                    x(j,i) = z(j,i)/(1.0d0-V(i))
                end do
            elseif (phase(i) == 1) then
                do j=1,com_2phase
                    call residualvectorset3(n*eq+q_judge,y(j,i))
                end do
                do j=1,com_2phase+com_ion
                    x(j,i) = z(j,i)
                end do
            end if
        end do
        !do j=1,com_2phase+com_ion
        !    kakuninnnnnn(j) = y(j,1)
        !end do
        !call outxs(kakuninnnnnn,kakuninn)
        !write(*,*) kakuninn

        !!フガシティ
        call vapor_fugacity_main(y,lnfai_V,P,z_factor0,z_factor,q_judge,phase_judge)
        call liquid_fugacity_main(lnfai_L,P,v_L1,v_L2,v_L3,v_L4)

        !do j=1,com_2phase!+com_ion
        !    kakuninnnnnn(j) = lnfai_L(j,1)
        !end do
        !call outxs(kakuninnnnnn,kakuninn)
        !write(*,*) kakuninn

        
        !!モル密度、体積
        !?気相
        do i=1,n
            if (phase_judge(i) == 2) then
                MD_V(i) = P(i)/z_factor(i)/R/temp
                MV_V(i) = 1.0d0/MD_V(i)
            else
                call residualvectorset3(n*eq+q_judge,MD_V(i))
                call residualvectorset3(n*eq+q_judge,MV_V(i))
            end if
        end do
        !?液相
        do i=1,n
            v_L(1,i) = v_L1(i)
            call residualvectorset4(v_L2(i),n*eq+q_judge,v_L(2,i))
            call residualvectorset4(v_L3(i),n*eq+q_judge,v_L(3,i))
            call residualvectorset4(v_L4(i),n*eq+q_judge,v_L(4,i)) 
            do j=com_2phase+1,com_2phase+com_ion
                call residualvectorset3(n*eq+q_judge,v_L(j,i))!!要検討====================
            end do
            call residualvectorset3(n*eq+q_judge,MV_L(i))
            do j=1,com_2phase+com_ion
                MV_L(i) = MV_L(i) + x(j,i)*v_L(j,i)
            end do
            MD_L(i) = 1.0d0/MV_L(i)
        end do

        !call outxs(MD_L,kakuninn)
        !write(*,*) kakuninn


        !!飽和率の算出
        do i=1,n
            L(i)=1.0d0-V(i)
            Sw(i)=Nt(i)*L(i)*MV_L(i)
            if(phase_judge(i) == 2) then
                Sg(i)=Nt(i)*V(i)*MV_V(i)
            else
                call residualvectorset3(n*eq+q_judge,Sg(i))
            end if
        end do
        call outxs(Sw,Sw0)
        !write(*,*) Sw0
        
        
        !!相質量密度[kg/m^3]
        Mw(1)=18.01528d0 !モル質量[g/mol]
        Mw(2)=44.01d0
        Mw(3)=16.04d0
        Mw(4)=32.082d0
        Mw(5)=1.0079d0
        Mw(6)=61.0171d0
        Mw(7)=60.00092d0
        Mw(8)=17.0073d0
        Mw(9)=40.078d0
        Mw(10)=24.0d0
        Mw(11)=60.0d0
        Mw(12)=28.9799d0
        Mw(13)=33.0735d0
        Mw(14)=32.07d0

        do i=1,n
            call residualvectorset3(n*eq+q_judge,Mw_ave_L(i))
            call residualvectorset3(n*eq+q_judge,Mw_ave_V(i))
            do j=1,com_2phase+com_ion
                Mw_ave_L(i)=Mw_ave_L(i)+Mw(i)*x(j,i)
            end do
            do j=1,com_2phase
                Mw_ave_V(i)=Mw_ave_V(i)+Mw(i)+y(j,i)
            end do
            phase_d_L(i)=Mw_ave_L(i)*MD_L(i)*10.0d0**(-3.0d0) ![kg/m^3]
            phase_d_V(i)=Mw_ave_V(i)*MD_V(i)*10.0d0**(-3.0d0)
        end do
        call outxs(phase_d_V,kakuninn)
        !write(*,*) kakuninn

        !!相粘度
        !?液相
        beta_v=-1.297d0+(0.574*10.0d0**(-1.0d0))*(temp-273.15)-(0.697d0*10.0d0**(-3.0d0))*(temp-273.15)**2.0d0+&
            (0.447d0*10.0d0**(-5.0d0))*(temp-273.15)**3.0d0-(0.105d0*10.0d0**(-7.0d0))*(temp-273.15)**4.0d0 ![1/GPa]
        myu_H2O_20=1002.d0 !20[℃] ![μPa・s]
        sa=20.0d0+273.15d0-temp
        sisuu=(sa*(1.2378d0-(1.303d0*10.0d0**(-3.0d0))*sa+&
            (3.06d0*10.0d0**(-6.0d0))*sa**2.0d0+&
            (2.55d0*10.0d0**(-8.0d0))*sa**3.0d0)/(96.0d0+temp-273.15d0))
        myu_L_normal=myu_H2O_20*10.0d0**sisuu ![μPa・s]
        do i=1,n
            myu_L(i)=(myu_L_normal*10.0d0**(-6.0d0))*(1.0d0+beta_v*(P(i)*10.0d0**(-9.0d0))) ![Pa・s]
        end do

        !call outxs(myu_L,kakuninn)
        !write(*,*) kakuninn
        !?気相
        !vc臨界体積[m^3/mol],1:H2O
        Vc(1)=56.3d0*10.0d0**(-6.0d0)
        Vc(2)=94.0d0*10.0d0**(-6.0d0)
        Vc(3)=99.0d0*10.0d0**(-6.0d0)
        Vc(4)=98.5d0*10.0d0**(-6.0d0)
        Tc(1) = 647.30d0
        Tc(2) = 304.2d0
        Tc(3) = 190.6d0
        Tc(4) = 373.2d0
        Pc(1) = 217.6d0 * 101325.d0
        Pc(2) = 72.8d0 * 101325.d0
        Pc(3) = 45.4d0 * 101325.d0
        Pc(4) = 88.2d0 * 101325.d0
        
        do i=1,n
            call residualvectorset3(n*eq+q_judge,Pc_ave_V(i))
            call residualvectorset3(n*eq+q_judge,Tc_ave_V(i))
            call residualvectorset3(n*eq+q_judge,Vc_ave_V(i))
            if (phase_judge(i) == 2) then
                do j=1,com_2phase
                    Pc_ave_V(i)=Pc_ave_V(i)+Pc(j)/101325.0d0*y(j,i) ![Pa→atm]
                    Tc_ave_V(i)=Tc_ave_V(i)+Tc(j)*y(j,i) ![K]
                    Vc_ave_V(i)=Vc_ave_V(i)+Vc(j)*y(j,i) ![m^3/mol]
                end do
                zeta(i)=Tc_ave_V(i)**(1.0d0/6.0d0)/((Mw_ave_V(i)**(1.0d0/2.0d0))*(Pc_ave_V(i)**(2.0d0/3.0d0)))
                ro_r_V(i)=MD_V(i)*Vc_ave_V(i)
            elseif (phase(i) == 1) then
                call residualvectorset3(n*eq+q_judge,zeta(i))
                call residualvectorset3(n*eq+q_judge,ro_r_V(i))
            end if
        end do

        do i=1,com_2phase
            myu_c(i)=sqrt(Mw(i))*(Pc(i)/101325.d0)**(2.0d0/3.0d0)/(Tc(i)**(1.0d0/6.0d0))*(4.610d0*(temp/Tc(i))**0.618d0-&
                    2.04d0*exp(-0.449d0*(temp/Tc(i)))+1.94d0*exp(-4.058d0*(temp/Tc(i)))+0.1d0)*10.0d0**(-4.0d0) ![cP]
        end do
        !write(*,*) myu_c

        do i=1,n
            if (phase_judge(i) == 2) then
                !myu_V_SC(i)=(y(1,i)*myu_c(1)*sqrt(Mw(1))+y(2,i)*myu_c(2)*sqrt(Mw(2)))/&
                !            (y(1,i)*sqrt(Mw(1))+y(2,i)*sqrt(Mw(2))) ![cP]
                myu_V_SC(i)=(y(1,i)*myu_c(1)*sqrt(Mw(1))+y(2,i)*myu_c(2)*sqrt(Mw(2))&
                            +y(3,i)*myu_c(3)*sqrt(Mw(3))+y(4,i)*myu_c(4)*sqrt(Mw(4)))/&
                            (y(1,i)*sqrt(Mw(1))+y(2,i)*sqrt(Mw(2))+y(3,i)*sqrt(Mw(3))+y(4,i)*sqrt(Mw(4))) ![cP]
            elseif (phase(i) == 1) then
                call residualvectorset3(n*eq+q_judge,myu_V_SC(i))
            end if
        end do

        av_0=1.0230d0*10.0d0**(-1.0d0)
        av_1=2.3364d0*10.0d0**(-2.0d0)
        av_2=5.8533d0*10.0d0**(-2.0d0)
        av_3=-4.0758d0*10.0d0**(-2.0d0)
        av_4=9.3324d0*10.0d0**(-3.0d0)

        do i=1,n
            if (phase_judge(i) == 2) then
                myu_V(i)=(myu_V_SC(i)+(((av_0+av_1*ro_r_V(i)+av_2*ro_r_V(i)**2.0d0+av_3*ro_r_V(i)**3.0d0+&
                         av_4*ro_r_V(i)**4.0d0)**4.0d0-10.0d0**(-4.0d0))/zeta(i)))*10.0d0**(-3.0d0) !cP→Pa・s
            elseif (phase(i) == 1) then
                call residualvectorset3(n*eq+q_judge,myu_V(i))
            end if
        end do
        
        call outxs(myu_V,kakuninn)
        !write(*,*) kakuninn!!粘度よさそう

        !!相対浸透率について
        do i=1,n
            !do j=1,20
            !    if (Swd(j) < Sw0(i) .and. Sw0(i) <= Swd(j+1)) then
            !        krg(i)=(krgd(j+1)-krgd(j))/(Swd(j+1)-Swd(j))*(Sw(i)-Swd(j))+krgd(j)
            !        krw(i)=(krwd(j+1)-krwd(j))/(Swd(j+1)-Swd(j))*(Sw(i)-Swd(j))+krwd(j)
            !    end if
            !end do
            if (Sw0(i) < Swc) then
                call residualvectorset4(0.0d0,n*eq+q_judge,krw(i))
                call residualvectorset4(krg_swc,n*eq+q_judge,krg(i))
            else if (1.0d0-Sgr < Sw0(i)) then
                call residualvectorset4(krw_sgr,n*eq+q_judge,krw(i))
                call residualvectorset4(0.0d0,n*eq+q_judge,krg(i))
            else
                krw(i)=krw_sgr*((Sw(i)-Swc)/(1.0d0-Swc-Sgr))**nw
                krg(i)=krg_swc*((1.0d0-Sw(i)-Sgr)/(1.0d0-Swc-Sgr))**ng
            end if
        end do
        !call outxs(krg,kakuninn)
        !write(*,*) kakuninn

        !!孔隙率について
        NmMD(1)=Nm1_MD
        NmMD(2)=Nm2_MD
        NmMD(3)=Nm3_MD
        NmMD(4)=Nm4_MD
        NmMD(5)=Nm5_MD
        do i=1,n
            call residualvectorset4(fai00(i),n*eq+q_judge,faimine(i))
            faimineold(i)=fai00(i)
            do j=1,com_mine
                faimine(i)=faimine(i)-(Nm(j,i)-Nmini(j))/NmMD(j)
                faimineold(i)=faimineold(i)-(Nm0old(j,i)-Nmini(j))/NmMD(j)
            end do
            fai(i)=faimine(i)*(1.0d0+Cr*(P(i)-iniPressure))
            faiold(i)=faimineold(i)*(1.0d0+Cr*(P0old(i)-iniPressure))            
        end do
        call outxs(fai,fai0)
        
        fai000=fai0
        !write(*,*) fai000

        !!絶対浸透率
        do i=1,n
            !kk(i)=k_ini*(fai(i)/faiini)**3.0d0*((1.0d0-faiini)/(1.0d0-fai(i)))**2.0d0
            call residualvectorset4(k_ini,n*eq+q_judge,kk(i))
        end do


        !!injection
        if (phase_judge(1) == 2) then !?2相の時
            ramda = krg(1)/myu_V(1) + krw(1)/myu_L(1)
        elseif (phase(1) == 1) then !?液相のみ
            ramda = 1.0d0/myu_L(1)
        else !?気相のみ 
            ramda = 1.0d0/myu_V(1)
        end if
        re=0.14d0*sqrt(dx**2.0d0+dy**2.0d0)
        if (q_judge == 1) then
            q=2.0d0*pi*kk(1)*dz*ramda*(Pb-P(1))/(log(re/rw)+skin)/dx/dy/dz !?流量制御
        else
            !?q=2.0d0*pi*kk(1)*dz*ramda*(Pbh0-P(1))/(log(re/rw)+skin)/dx/dy/dz !?圧力制御
            !?call residualvectorset3(n*eq+q_judge,q) !?圧力制御 !!ここ覚えていない、上の方が正しそうだけど、なんだっけ？
        end if

        call injection_data_d(P(1),MD_injection_d)

        !!weighting factor------
        do i=1,n-1
            if (P0(i+1) >= P0(i)) then
                W(i)=0.0d0
            else 
                w(i)=1.0d0
            end if
        end do
        w(:) = 1.0d0
        !write(1,*)w
        write(10,*) W
        !?----------------------

        !!トランスミッシビリティー-----
        do i=1,n
            if (phase_judge(i) == 2) then !?2相の時
                do j=1,com_2phase+com_ion
                    T_L(j,i)=kk(i)*krw(i)*MD_L(i)*x(j,i)/myu_L(i)
                end do
                do j=1,com_2phase
                    T_V(j,i)=kk(i)*krg(i)*MD_V(i)*y(j,i)/myu_V(i)
                end do
            elseif (phase(i) == 1) then !?液相だけの時
                do j=1,com_2phase+com_ion
                    T_L(j,i)=kk(i)*MD_L(i)*x(j,i)/myu_L(i)
                end do
                do j=1,com_2phase
                    call residualvectorset3(n*eq+q_judge,T_V(j,i))
                end do
            else !?気相だけの時
                do j=1,com_2phase+com_ion
                    call residualvectorset3(n*eq+q_judge,T_L(j,i))
                end do
                do j=1,com_2phase
                    T_V(j,i)=kk(i)*MD_V(i)*y(j,i)/myu_V(i)
                end do
            end if
        end do


        !!化学反応======================================
        !!速度定数[mol/m^3/s]
        do i=1,chemi
            ks(i) =1.0d0*10.0d0**(-4.0d0) !化学反応
        end do
        ks(6)=10.0d0**(-9.12d0) !アノーサイト
        ks(7)=10.0d0**(-12.0d0)!10.0d0**(-12.7d0) !エンスタタイト
        ks(8)=10.0d0**(-10.64d0) !フォルステライト
        ks(9)=10.0d0**(-5.81d0) !カルサイト
        ks(10)=10.0d0**(-12.0d0)!10.0d0**(-9.34d0) !マグネサイト

        !!貯留層の温度における鉱物反応の速度定数--------
        !ks(6)=ks(6)*exp(-Ea1/R*(1.0d0/temp-1.0d0/(273.15d0+25.0d0)))
        !ks(7)=ks(7)*exp(-Ea2/R*(1.0d0/temp-1.0d0/(273.15d0+25.0d0)))
        !ks(8)=ks(8)*exp(-Ea3/R*(1.0d0/temp-1.0d0/(273.15d0+25.0d0)))
        !ks(9)=ks(9)*exp(-Ea4/R*(1.0d0/temp-1.0d0/(273.15d0+25.0d0)))
        !ks(10)=ks(10)*exp(-Ea5/R*(1.0d0/temp-1.0d0/(273.15d0+25.0d0)))
        
        
        !!平衡定数
        ke(1) = ke1
        ke(2) = ke2
        ke(3) = ke3
        ke(4) = ke4
        ke(5) = ke5
        km(1) = km1
        km(2) = km2
        km(3) = km3
        km(4) = km4
        km(5) = km5

        do i=1,n
            if (phase(i) == 1 .or. phase_judge(i) ==2) then
            do j=1,com_2phase+com_ion
                wc(j,i)=MD_L(i)*x(j,i)/(MD_L(i)*x(1,i)*18.015d0*10.0d0**(-3.0d0))
            end do
            call residualvectorset4(1.0d0,n*eq+q_judge,wc(1,i))
            w1(i)=wc(1,i)
            w2(i)=wc(2,i)
            w3(i)=wc(3,i)
            w4(i)=wc(4,i)
            w5(i)=wc(5,i)
            w6(i)=wc(6,i)
            w7(i)=wc(7,i)
            w8(i)=wc(8,i)
            w9(i)=wc(9,i)
            wten(i)=wc(10,i)
            w11(i)=wc(11,i)
            w12(i)=wc(12,i)
            w13(i)=wc(13,i)
            w14(i)=wc(14,i)
        end if
        end do
        call outxs(w1,w10)
        call outxs(w2,w20)
        call outxs(w3,w30)
        call outxs(w4,w40)
        call outxs(w5,w50)
        call outxs(w6,w60)
        call outxs(w7,w70)
        call outxs(w8,w80)
        call outxs(w9,w90)
        call outxs(wten,w100)
        call outxs(w11,w110)
        call outxs(w12,w120)
        call outxs(w13,w130)
        call outxs(w14,w140)
        do i=1,n
            wc0(1,i)=w10(i)
            wc0(2,i)=w20(i)
            wc0(3,i)=w30(i)
            wc0(4,i)=w40(i)
            wc0(5,i)=w50(i)
            wc0(6,i)=w60(i)
            wc0(7,i)=w70(i)
            wc0(8,i)=w80(i)
            wc0(9,i)=w90(i)
            wc0(10,i)=w100(i)
            wc0(11,i)=w110(i)
            wc0(12,i)=w120(i)
            wc0(13,i)=w130(i)
            wc0(14,i)=w140(i)
        end do

        !!比表面積
        AA(1)=AA1
        AA(2)=AA2
        AA(3)=AA3
        AA(4)=AA4
        AA(5)=AA5
        do i=1,n
            do j=1,com_mine
                Ac(j,i)=AA(j)*Nm(j,i)/Nmini(j)
            end do
        end do

        !!活量積
        do i=1,n
            do j=1,chemi
                call residualvectorset4(1.0d0,n*eq+q_judge,Q_chemi(j,i))
                do jj=1,com_2phase+com_ion
                    if (Nc0(jj,i) /= 0) then
                        Q_chemi(j,i)=Q_chemi(j,i)*wc(jj,i)**chemi_mat(j,jj)
                    end if
                    theta(j,i)=Q_chemi(j,i)/ke(j)
                enddo
            end do
            do j=1,mine
                call residualvectorset4(1.0d0,n*eq+q_judge,Q_chemi(j+chemi,i))
                do jj=1,com_2phase+com_ion
                    if (Nc0(jj,i) /= 0) then
                        Q_chemi(j+chemi,i)=Q_chemi(j+chemi,i)*wc(jj,i)**chemi_mat(j+chemi,jj)
                    end if
                    theta(j+chemi,i)=Q_chemi(j+chemi,i)/km(j)
                enddo
            end do
            theta1(i)=theta(1,i)
            theta2(i)=theta(2,i)
            theta3(i)=theta(3,i)
            theta4(i)=theta(4,i)
            theta5(i)=theta(5,i)
            theta6(i)=theta(6,i)
            theta7(i)=theta(7,i)
            theta8(i)=theta(8,i)
            theta9(i)=theta(9,i)
            thetaten(i)=theta(10,i)
        end do
        call outxs(theta1,theta10)
        call outxs(theta2,theta20)
        call outxs(theta3,theta30)
        call outxs(theta4,theta40)
        call outxs(theta5,theta50)
        call outxs(theta6,theta60)
        call outxs(theta7,theta70)
        call outxs(theta8,theta80)
        call outxs(theta9,theta90)
        call outxs(thetaten,theta100)
        !write(*,*) theta10
        !write(*,*) theta20
        !write(*,*) theta30
        !write(*,*) theta40
        !write(*,*) theta50
        !write(*,*) theta60
        !write(*,*) theta70
        !write(*,*) theta80
        !write(*,*) theta90
        !write(*,*) theta100
        do i=1,n
            theta0(1,i)=theta10(i)
            theta0(2,i)=theta20(i)
            theta0(3,i)=theta30(i)
            theta0(4,i)=theta40(i)
            theta0(5,i)=theta50(i)
            theta0(6,i)=theta60(i)
            theta0(7,i)=theta70(i)
            theta0(8,i)=theta80(i)
            theta0(9,i)=theta90(i)
            theta0(10,i)=theta100(i)
        end do

        

        do i=1,n
            !!化学反応の反応速度
            !H2O + CO2 = HCO3- + H+
            !HCO3- = H+ + CO32-
            !H2O = H+ + OH-
            
            do j=1,chemi
                min0=10000000000000000.0d0
                call residualvectorset4(min0,n*eq+q_judge,min)
                if (theta0(j,i) < 1.0d0) then !順反応
                    do jj=1,com_2phase+com_ion
                        if (chemi_mat(j,jj) < 0.0d0) then !減少成分
                            if (min0 > wc0(jj,i)) then !最小成分
                                min0=wc0(jj,i)
                                min=wc(jj,i)
                            end if
                        end if
                    end do
                else !逆反応
                    do jj=1,com_2phase+com_ion
                        if(chemi_mat(j,jj) > 0.0d0) then
                            if (min0 > wc0(jj,i)) then
                                min0=wc0(jj,i)
                                min=wc(jj,i)
                            end if
                        end if
                    end do
                end if
                rs(j,i)=ks(j)*(1.0d0-theta(j,i))*fai(i)*Sw(i)*min
            end do


            !!鉱物反応の反応速度
            !αCO3 + H+ = α2+ + HCO3-
            !αSiO3 + 2H+ = α2+ + H2O + SiO2
            do j=chemi+1,chemi+mine
                rs(j,i)=Ac(j-chemi,i)*ks(j)*(1.0d0-theta(j,i))*Sw(i)!*fai(i)
            end do
        end do
        
        do i=1,n
        !    rs1(i)=rs(3,i)
        end do
        !call outxs(rs1,kakuninn)
        !write(*,*) kakuninn

        !!反応速度の合計
        do i=1,n
            do j=1,com_all
                call residualvectorset3(n*eq+q_judge,rs_sum(j,i))
                do jj=1,chemi+mine
                    rs_sum(j,i)=rs_sum(j,i)+rs(jj,i)*chemi_mat(jj,j)
                end do
            end do
        end do

        do i=1,n
        !    rs1(i) = rs_sum(9,i)
        end do
        !call outxs(rs1,kakuninn)
        !write(*,*) kakuninn


        !!物質収支式
        !?grid1
        do j=1,com_2phase+com_ion
            g(1*eq-eq+com_2phase+j)=(-1.0d0/dt)*(Nc(j,1)*fai(1)-Nc0old(j,1)*faiold(1))
            if (phase(1) == 1 .or. phase(1) == 2) then !?液相あり
                g(1*eq-eq+com_2phase+j)=&
                g(1*eq-eq+com_2phase+j)+(1.0d0/dx**2.0d0)*((w(1)*T_L(j,1)+(1.0d0-W(1))*T_L(j,1))*(P(2)-P(1)))&    
                                                +rs_sum(j,1)
            end if
        end do
        if (phase(1) == 3 .or. phase(1) == 2) then !?気相あり
            do j=1,com_2phase
                g(1*eq-eq+com_2phase+j)=&
                g(1*eq-eq+com_2phase+j)+(1.0d0/dx**2.0d0)*((w(1)*T_V(j,1)+(1.0d0-W(1))*T_V(j,1))*(P(2)-P(1)))
            end do
        end if

        g(1*eq-eq+com_2phase+2)=g(1*eq-eq+com_2phase+2)+q*MD_injection_d




        !?境界以外
        do i=2,n-1
            do j=1,com_2phase+com_ion
                g(i*eq-eq+j+com_2phase)=(-1.0d0/dt)*(Nc(j,i)*fai(i)-Nc0old(j,i)*faiold(i))
                if (phase(i) == 1 .or. phase(i) == 2) then !?液相あり
                    g(i*eq-eq+j+com_2phase)=&
                    g(i*eq-eq+j+com_2phase)+(1.0d0/dx**2.0d0)*((W(i)*T_L(j,i)+(1.0d0-W(i))*T_L(j,i+1))*(P(i+1)-P(i))-&
                                                               (W(i-1)*T_L(j,i-1)+(1.0d0-W(i-1))*T_L(j,i))*(P(i)-P(i-1)))&
                                                            +rs_sum(j,i)
                end if
            end do
            if (phase(i) == 3 .or. phase(i) == 2) then !?気相あり
                do j=1,com_2phase
                    g(i*eq-eq+j+com_2phase)=&
                    g(i*eq-eq+j+com_2phase)+(1.0d0/dx**2.0d0)*((W(i)*T_V(j,i)+(1.0d0-W(i))*T_V(j,i+1))*(P(i+1)-P(i))-&
                                                               (W(i-1)*T_V(j,i-1)+(1.0d0-W(i-1))*T_V(j,i))*(P(i)-P(i-1)))
                end do
            end if            
        end do


        !?gridn
        do j=1,com_2phase+com_ion
            g(n*eq-eq+com_2phase+j)=(-1.0d0/dt)*(Nc(j,n)*fai(n)-Nc0old(j,n)*faiold(n))
            if (phase(n) == 1 .or. phase(n) == 2) then !?液相あり
                g(n*eq-eq+com_2phase+j)=&
                g(n*eq-eq+com_2phase+j)+(-1.0d0/dx**2.0d0)*((W(n-1)*T_L(j,n-1)+(1.0d0-W(n-1))*T_L(j,n))*(P(n)-P(n-1)))&
                                                +rs_sum(j,n)
            end if
        end do
        if (phase(n) == 3 .or. phase(n) == 2) then !?気相あり
            do j=1,com_2phase
                g(n*eq-eq+com_2phase+j)=&
                g(n*eq-eq+com_2phase+j)+(-1.0d0/dx**2.0d0)*((W(n-1)*T_V(j,n-1)+(1.0d0-W(n-1))*T_V(j,n))*(P(n)-P(n-1)))
            end do
        end if



        do i =1,n
            do j=1,com_mine
                g(i*eq-eq+j+com_2phase+com_2phase+com_ion)=&
                    (-1.0d0)*(Nm(j,i)*(1.0d0-fai(i))-Nm0old(j,i)*(1.0d0-faiold(i)))/dt +rs_sum(j+com_2phase+com_ion,i)
            end do
        end do







        
        do i=1,n
            do j=1,com_2phase
                g(i*eq-eq+j)=lnk(j,i)+lnfai_V(j,i)-lnfai_L(j,i) !!熱力学的条件式
            end do
            g(i*eq-1)=Sw(i)+Sg(i)-1.0d0 !!飽和率の制約式
            g(i*eq)=rach(i) !!rachford-rice
        end do

        !!流量制御
        if (q_judge == 1) then
            g(n*eq+q_judge) = q_input - q*MD_injection_d !!坑井式
        end if
        Sw00=Sw0

        !call outxs(rach,kakuninn)
        !write(*,*) kakuninn
        
    end subroutine

end module
