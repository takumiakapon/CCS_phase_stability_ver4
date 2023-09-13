!#TODOメモリリーク説があるので、初期計算は別でやって、無理やりインプットさせる！！
program  main
    
    use mod_autodiff
    use mod_condition
    use mod_cubic_eq_d
    use mod_cubic_eq
    use mod_gauss
    use mod_ini_N
    use mod_input
    use mod_phase_stability_analysis
    use mod_fugacity
    use mod_ini_flash
    use mod_ini_chemical
    use mod_main_calc
    implicit none


    
    !!初期計算用
    integer::i,j,jj,iteration,time,phase_judge0,phase0,ii,year,q_judge,hour,day
    real(8),allocatable,dimension(:,:)::amat,cmat,emat,gmat
    real(8),allocatable,dimension(:)::bmat,z0,k0,lnk0,x0,y0,w0,alpha0,dmat,fmat,hmat
    real(8)::V0,L0,P0,P0old,error,wt,lumda,z_factor0
    type(diffs),allocatable::fxs(:)
    
    !!
    real(8),dimension(com_2phase+com_ion)::Nc0,Nc0old
    real(8),dimension(com_mine)::Nm0,Nm0old
    real(8)::Nt,Sw0,lnk1,lnk2,lnk3,lnk4
    real(8),allocatable,dimension(:,:)::chemi_mat
    
    !main計算
    integer,dimension(n)::phase_judge,phase
    real(8),dimension(n)::V,L,P,Pold,fai000
    real(8),dimension(com_2phase,n)::lnk
    real(8),dimension(com_2phase+com_ion,n)::Nc,Ncold,z
    real(8),dimension(com_mine,n)::Nm,Nmold
    real(8)::Pb0,fai0,Nmini(com_mine),theta0(chemi+mine,n)
    real(8),dimension(21)::Swd,krgd,krwd
    real(8),allocatable,dimension(:)::Sw,fai
    real(8),allocatable,dimension(:,:)::wc
    
    
    
    call file_open
    
    allocate(z0(com_2phase+com_ion),k0(com_2phase),lnk0(com_2phase)&
        ,x0(com_2phase),y0(com_2phase+com_ion),w0(com_2phase),alpha0(com_2phase),chemi_mat(chemi+mine,com_all)&
        ,Sw(n),wc(com_2phase+com_ion,n),fai(n))
     
    !!相対浸透率のテーブルデータのインプット
    do i=1,21
        read(101,*) Swd(i),krgd(i),krwd(i)
    end do

    !!化学反応の係数のインプット
    do i=1,chemi+mine
        read(1,*) (chemi_mat(i,j),j=1,com_all)
    end do
    
    !初期設定
    
    
    P0 =iniPressure
    
    
    !平衡定数の初期予測
    !lnk0(1) = -5.0d0
    !lnk0(2) = 3.0d0
    !lnk0(3) = 5.0d0
    !lnk0(4) = 3.0d0
    !k0(1) = exp(lnk0(1))!1.0d0/Pc(1)*exp(5.37*(1.0d0+acentric(1)*(1.0d0-1.0d0/tc(1))))!exp(lnk0(1))!=1.0d0/pr(i)*exp(5.37*(1.0d0+omega(i)*(1.0d0-1.0d0/tr(i)))
    !k0(2) = exp(lnk0(2))!1.0d0/Pc(2)*exp(5.37*(1.0d0+acentric(2)*(1.0d0-1.0d0/tc(2))))!exp(lnk0(2))
    !k0(3) = exp(lnk0(3))
    !k0(4) = exp(lnk0(4))
    
    
    !初期の全モル分率
    !z0(1) = 0.95d0 !H2O
    !z0(2) = 0.01d0 !CO2
    !z0(3) = 0.0d0 !CH4
    !z0(4) = 0.01d0 !H2S
    !do i =com_2phase+1,com_2phase+com_ion
    !    z0(i) = (1.0d0 - (z0(1)+z0(2)+z0(3)+z0(4)))/(com_2phase+com_ion-com_2phase) !そのほかのイオン
    !end do
    
    !!!検証 これ
    !z0=0.0000d0
    !z0(1)=0.982!0.98295661537557144!0.98294928708382090!0.98294  !H2O
    !z0(2)=2.34d0*10.0d0**(-7.0d0)!3.4817505806941896E-007!1.8992194875568811E-005!0.00001 !CO2
    !z0(3)=1.70d0*10.0d0**(-2.0d0)!1.6990282021069454E-002!1.6990001958143271E-002!0.01699 !CH4
    !z0(5)=6.57d0*10.0d0**(-11.0d0)!6.1817158796360923E-011!6.3674605894641180E-007!1E-5 !H+
    !z0(6)=1.76d0*10.0d0**(-5.0d0)!2.8331115169197601E-005!1.1005716001412855E-006!1E-5 !HCO3-
    !z0(7)=2.75d0*10.0d0**(-7.0d0)!1.3212068937893931E-006!9.9072355493762180E-006!1E-5 !CO32-
    !z0(8)=2.60d0*10.0d0**(-7.0d0)!1.0885035523759696E-006!9.9072356593756316E-006!1E-5 !OH-
    !z0(10)=4.23d0*10.0d0**(-6.0d0)!1.1030901971414173E-005!1.0092765812147074E-005!1E-5 !Mg2+
    !z0(11)=1.0d0-z0(1)-z0(2)-z0(3)-z0(5)-z0(6)-z0(7)-z0(8)-z0(10) !SiO2
    
    !!ケースこれ
    !z0=0.0d0
    !z0(1)=0.990029871
    !z0(2)=1.68415E-05
    !z0(3)=0.008000107!0.0d0
    !z0(4)=0.0d0
    !z0(5)=9.43763E-10
    !z0(6)=0.000118031
    !z0(7)=1.3007E-07
    !z0(8)=1.86134E-08
    !z0(10)=9.13985E-06
    !z0(11)=1.0d0-z0(1)-z0(2)-z0(3)-z0(5)-z0(6)-z0(7)-z0(8)-z0(10) !SiO2
    
    
    !相安定解析の主要変数の初期値
    !if (z0(1) >= z0(2)+z0(3)+z0(4)) then !液相から気相が出てくるパターン
    !    do i =1,com_2phase
    !        w0(i) = z0(i)*k0(i)
    !    end do
    !else !気相から液相が出てくるパターン
    !    do i=1,com_2phase
    !        w0(i) = z0(i)/k0(i)
    !    end do
    !end if
    !V0 = 1.0d0 - z0(1)
    !L0 = 1.0d0 -V0
    !do i=1,com_2phase
    !    alpha0(i) = 2.0d0 *sqrt(w0(i))
    !end do
    
    !write(*,*) lnk0
    
    
    !do iteration =1,100
    !    if ((z0(1) >= z0(2)+z0(3)+z0(4))) then
    !        call phase_stability_liquid(alpha0,P0,z0,fxs)
    !    else
    !        call phase_stability_vapor(alpha0,P0,z0,fxs)
    !    end if
        
    !    call jacobian_mat(fxs,amat)
    !    call outxs(fxs,bmat)
    !    bmat = -bmat
    !    do i=1,com_2phase
    !        if (z0(i) == 0.0d0) then !存在しない成分は計算しないよ
    !            do j =1,com_2phase
    !                amat(i,j) = 0.0d0
    !                amat(j,i) = 0.0d0
    !            end do
    !            amat(i,i) = 1.0d0
    !            bmat(i) =0.0d0
    !        end if
    !    end do
            
    !    call pvgauss(com_2phase,amat,bmat)
        
    !    do i=1,com_2phase
    !        alpha0(i) =alpha0(i) +bmat(i)
    !        w0(i) = (alpha0(i) / 2.0d0) ** 2.0d0
    !    end do
        
    !    error = sqrt(dot_product(bmat,bmat))
        !write(*,*) error,iteration
    !    if (error < 10.0d0**(-5.0d0)) then
    !        exit
    !    end if
    !end do
    
    
    !wt =0.0d0
    !do i=1,com_2phase
    !    wt =wt+w0(i)
    !end do
    
    !write(*,*)(log(wt))
    !lumda =1.0d0 -log(wt)
    !if (lumda >= 1.0d0) then
    !    phase_judge0 = 1
    !    if (z0(1) > z0(2) +z0(3)+z0(4)) then
    !        V0 = 0.0d0
    !        L0 = 1.0d0 - V0
    !        phase0 = 1 !液相のみ
    !    else
    !        V0 = 1.0d0
    !        L0 = 1.0d0 - V0
    !        phase0 = 3 !気相のみ
            !write(*,*) V0
    !    end if
    !else
    !    phase_judge0 = 2
    !    phase0 = 2 !2相
    !    v0 = z0(2) + z0(3) + z0(4)
    !    L0 = 1.0d0 - v0
    !    if (z0(1) >= (z0(2)+z0(3)+z0(4))) then !液相から気相が出てくるパターン
    !        do i=1,com_2phase
    !            if (z0(i) == 0d0) then
    !                k0(i) = 1.0d0
    !                lnk0(i) = log(k0(I))
    !            else
    !                k0(i) = w0(i) / z0(i)
    !                lnk0(i) = log(k0(i))
    !            end if
    !        end do
    !    else !気相から液相が出てくるパターン
    !        do i=1,com_2phase
    !            if (z0(i)== 0.0d0) then
    !                k0(i) = 5.0
    !                lnk0(i) = log(k0(i))
    !            else
    !                k0(i) = z0(i) / w0(i)
    !                lnk0(i) = log(k0(i))
    !            end if
    !        end do
    !    end if
    !end if
    !deallocate(amat,bmat)
    !write(*,*) lnk0,lumda,v0
    !do i=1,com_2phase
    !    write(*,*) z0(i)
    !end do
    !write(*,*) phase_judge0,'phase'
    
    
    !allocate(cmat(com_2phase+1,com_2phase+1),dmat(com_2phase+1))
    !if (phase_judge0 == 2) then
    !    do iteration = 1,100
    !        call ini_fla(P0,z0,lnk0,v0,z_factor0,fxs)
    !        call jacobian_mat(fxs,cmat)
    !        call outxs(fxs,dmat)
    !        dmat = -dmat
    !        do i=1,com_2phase
    !            if (z0(i) == 0.0d0) then !存在しない成分は計算しないよ
    !                do j =1,com_2phase
    !                    cmat(i,j) = 0.0d0
    !                    cmat(j,i) = 0.0d0
    !                end do
    !                cmat(i,i) = 1.0d0
    !                dmat(i) =0.0d0
    !            end if
    !        end do
            !do i=1,com_2phase+1
                !write(*,*) (cmat(i,j),j=1,com_2phase+1)
            !end do
    !        call pvgauss(com_2phase+1,cmat,dmat)
    !        do i=1,com_2phase
    !            lnk0(i) = lnk0(i) + dmat(i)
    !            k0(i) = exp(lnk0(i))
    !        end do
    !        v0 = v0 +dmat(com_2phase+1)
    !        
    !        error = sqrt(dot_product(dmat,dmat))
    !        if (error < 10.0d0**(-5.0d0)) then
    !            L0 = 1.0d0 - V0
    !            exit
    !        end if
    !    end do
    !end if
    !deallocate(cmat,dmat)
    !lnk1=lnk0(1) !chemicalの前に謎にlnkの値が0になるバグが10回に1回くらいおこるので、一旦保存して、あとで　代入する。まじでなぞ
    !lnk2=lnk0(2)
    !lnk3=lnk0(3)
    !lnk4=lnk0(4)
    
    !if (phase_judge0 == 2) then !2相
    !    call ini_N_2phase(z_factor0,V0,z0,x0,y0,k0,Nc0,Nt,Sw0)
    !else if (phase0 == 1) then !液相のみ
    !    call ini_N_liquid(V0,z0,x0,y0,Nc0,Nt,Sw0)
    !else !気相のみ
    !    call ini_N_vapor(V0,z0,x0,y0,Nc0,Nt,Sw0)
    !end if
    
    !!岩石の初期モル数
    !Nm0(1)=Nm1_ini
    !Nm0(2)=Nm2_ini
    !Nm0(3)=Nm3_ini
    !Nm0(4)=Nm4_ini
    !Nm0(5)=Nm5_ini
    !do i=1,com_2phase+com_ion
    !    Nc0old(i) =Nc0(i)
    !end do
    !do i=1,com_mine
    !    Nm0old(i) =Nm0(i)
    !end do
    !P0old =p0
    
    
    !lnk0(1)=lnk1
    !lnk0(2)=lnk2
    !lnk0(3)=lnk3
    !lnk0(4)=lnk4
    
    !write(*,*) Nc0
    
    !液相がある場合、電離
    !if (phase_judge0 == 2 .or. phase0 == 1) then
    !    allocate(emat(eq,eq),fmat(eq))
    !    do time =1,100!00 !time iteration
            !write(11,*) 'day',time
    !        do iteration =1,100
                !write(11,*) iteration
    !            call ini_chemi(P0,P0old,Nc0,Nc0old,Nm0,Nm0old,lnk0,V0,Sw0,fai0,chemi_mat,phase_judge0,z0,theta0,z_factor0,fxs)
    !            call jacobian_mat(fxs,emat)
    !            call outxs(fxs,fmat)
    !            fmat = -fmat
                !do i=1,eq
                    !write(10,*) (emat(i,j),j=1,eq)
                !end do
    !            do i=1,com_2phase
    !                if (Nc0(i) == 0.0d0) then !存在しない成分は計算しないよ
            !            write(10,*) 'a'
    !                    do j =1,eq
    !                        emat(i,j) = 0.0d0
    !                        emat(j,i) = 0.0d0
    !                    end do
    !                    emat(i,i) = 1.0d0
    !                    fmat(i) =0.0d0
    !                end if
    !            end do
    !            do i=1,com_2phase+com_ion
    !                if (Nc0(i) == 0.0d0) then !存在しない成分は計算しないよ
            !            write(10,*) 'i'
    !                    do j =1,eq
    !                        emat(i+chemi,j) = 0.0d0
    !                        emat(j,i+chemi) = 0.0d0
    !                    end do
    !                    emat(i+chemi,i+chemi) = 1.0d0
    !                    fmat(i+chemi) =0.0d0
    !                end if
    !            end do
    !            if (phase0 == 1) then !水相だけの時はrachford解かないよ
            !        write(10,*) 'u'
    !                do i =1,eq
    !                    emat(i,eq) = 0.0d0
    !                    emat(eq,i) = 0.0d0
    !                end do
    !                emat(eq,eq) = 1.0d0
    !                fmat(eq) = 0.0d0
    !            end if
                !do i=1,eq
                    !write(10,*) (emat(i,j),j=1,eq)
                !end do
            
    !            call pvgauss(eq,emat,fmat)
                
                
    !            do i=1,com_2phase
    !                lnk0(i) = lnk0(i) + fmat(i)
    !                k0(i) = exp(lnk0(i))
    !            end do
    !            do i=1,com_2phase+com_ion
    !                Nc0(i) = Nc0(i) + fmat(i+com_2phase)
    !            end do
    !            do i=1,com_mine
    !                Nm0(i) = Nm0(i) + fmat(i+com_2phase+com_2phase+com_ion)
    !            end do
    !            P0 = P0 + fmat(eq-1)
    !            v0 = v0 + fmat(eq)
    !            error = sqrt(dot_product(fmat,fmat))
                !write(10,*) time,iteration,error,'---------------------------------------------------'
                 !write(*,*) 'V:',V0
                 !write(*,*) p0
                !do i=1,eq
                    !write(10,*) (emat(i,j),j=1,eq)
                !end do
    !            if (error < 10.0d0**(-5.0d0)) then
    !                goto 10
    !            end if
                
    !        end do
!10          &
    !        do j=1,com_2phase+com_ion
    !            Nc0old(j) = Nc0(j)
    !        end do
    !        do j=1,com_mine
    !            Nm0old(j) = Nm0(j)
    !        end do
    !        P0 = inipressure
    !        P0old = P0
            
            
    !        if (iteration == 100) then
    !            write(*,*) '収束せず'
    !        end if
            
    !    end do
    !    deallocate(emat,fmat)
    !end if
    
    
    !write(*,*) Nc0
    !write(14,*) phase_judge0,'phases'
    !write(14,*) P0,'pressure'
    !write(14,*) 'z_factor:',z_factor0
    !write(14,*) 'V0',V0,'Sw',Sw0
    !write(14,*) 'Nc','z'
    !do i=1,com_2phase+com_ion
    !    write(14,*) Nc0(i),z0(i)
    !end do
    !do i=1,com_mine
    !    write(14,*) Nm0(i)
    !end do
    !write(14,*) 'theta'
    !do i=1,chemi+mine
    !    write(14,*) theta0(i)
    !end do
    read(11,*) V0
    do j=1,com_2phase+com_ion
        read(11,*) Nc0(j),z0(j)
    end do
    do j=1,com_mine
        read(11,*) Nm0(j)
    end do
    do j=1,com_2phase
        read(11,*) lnk0(j)
    end do
    read(11,*) fai0
    read(11,*) P0
    read(11,*) phase_judge0
    read(11,*) phase0

    !!初期の結果の引継ぎ(grid1→5)
    do i=1,n
        V(i)=v0
        L(i)=1.0d0-V(i)
        do j=1,com_2phase+com_ion
            Nc(j,i)=Nc0(j)
            Ncold(j,i)=Nc(j,i)
            z(j,i)=z0(j)
        end do
        do j=1,com_2phase
            lnk(j,i)=lnk0(j)
        end do
        do j=1,com_mine
            Nm(j,i)=Nm0(j)
            Nmini(j)=Nm0(j)
            Nmold(j,i)=Nm(j,i)
        end do
        fai(i)=fai0
        P(i)=P0
        Pold(i)=P(i)
        phase_judge(i)=phase_judge0
        phase(i)=phase0
    end do
    Pb0 = P(1)
    
    
    !!ようやくメイン計算！！！
    allocate(amat(com_2phase,com_2phase),bmat(com_2phase),gmat(n*eq,n*eq),hmat(n*eq))
    do year=1,1!3!50!000
        do day =1,3!10!50!0!150!0
        do hour =1,1!24    
        !    !!相安定解析
            do ii=1,n !gridごとに相安定性解析するよ
                do i=1,com_2phase+com_ion
                    z0(i)=z(i,ii)
                end do
                do i=1,com_2phase
                    k0(i)=exp(lnk(i,ii))
                end do
                
            
                !相安定解析の主要変数の初期値
                if (z0(1) >= z0(2)+z0(3)+z0(4)) then !液相から気相が出てくるパターン
                    do i =1,com_2phase
                        w0(i) = z0(i)*k0(i)
                    end do
                else !気相から液相が出てくるパターン
                    do i=1,com_2phase
                        w0(i) = z0(i)/k0(i)
                    end do
                end if
                !write(*,*) w0
                V0 = 1.0d0 - z0(1)
                L0 = 1.0d0 -V0
                do i=1,com_2phase
                    alpha0(i) = 2.0d0 *sqrt(w0(i))
                end do
                !write(*,*) alpha0
            
                !write(*,*) P0
            
                do iteration =1,100
                    if ((z0(1) >= z0(2)+z0(3)+z0(4))) then
                        
                        !write(*,*) 'main',ii
                        call phase_stability_liquid2(alpha0,P0,z0,fxs) !?こいつがなんか悪そう
                        !write(*,*) 'a'
                    else
                        !write(*,*) 'main',ii
                        !write(*,*) 'vapor'
                        call phase_stability_vapor(alpha0,P0,z0,fxs)
                    end if
            !write(*,*) 'main',ii
                    
                    call jacobian_mat(fxs,amat)
                    call outxs(fxs,bmat)
                    bmat = -bmat
                              
                    
                    do i=1,com_2phase
                        if (z0(i) == 0.0d0) then !存在しない成分は計算しないよ
                            do j =1,com_2phase
                                amat(i,j) = 0.0d0
                                amat(j,i) = 0.0d0
                            end do
                            amat(i,i) = 1.0d0
                            bmat(i) =0.0d0
                        end if
                    end do
                    !do i=1,com_2phase
                    !    write(*,*) (amat(i,j),j=1,com_2phase)
                    !end do
                    call pvgauss(com_2phase,amat,bmat)
                    
                    do i=1,com_2phase
                        alpha0(i) =alpha0(i) +bmat(i)
                        w0(i) = (alpha0(i) / 2.0d0) ** 2.0d0
                    end do
                
                    error = sqrt(dot_product(bmat,bmat))
                   !write(*,*) error,iteration
                    if (error < 10.0d0**(-3.0d0)) then
                        exit
                    end if
                end do
                
                wt =0.0d0
                do i=1,com_2phase
                    wt =wt+w0(i)
                end do
                
                !!write(*,*)(log(wt))
                lumda =1.0d0 -log(wt)
                if (lumda >= 1.0d0) then
                    phase_judge(ii) = 1
                    if (z0(1) > z0(2) +z0(3)+z0(4)) then
                        V0 = 0.0d0
                        L0 = 1.0d0 - V0
                        phase(ii) = 1 !液相のみ
                    else
                        V0 = 1.0d0
                        L0 = 1.0d0 - V0
                        phase(ii) = 3 !気相のみ
                !        !write(*,*) V0
                    end if
                else
                    phase_judge(ii) = 2
                    phase(ii) = 2 !2相
                    v0 = z0(2) + z0(3) + z0(4)
                    L0 = 1.0d0 - v0
                    if (z0(1) >= (z0(2)+z0(3)+z0(4))) then !液相から気相が出てくるパターン
                       do i=1,com_2phase
                            if (z0(i) == 0d0) then
                                k0(i) = 1.0d0
                                lnk0(i) = log(k0(i))
                            else
                                k0(i) = w0(i) / z0(i)
                                lnk0(i) = log(k0(i))
                            end if
                        end do
                    else !気相から液相が出てくるパターン
                        do i=1,com_2phase
                            if (z0(i)== 0.0d0) then
                                k0(i) = 5.0
                                lnk0(i) = log(k0(i))
                            else
                                k0(i) = z0(i) / w0(i)
                                lnk0(i) = log(k0(i))
                            end if
                        end do
                    end if
                end if
                do i=1,com_2phase
                    lnk(i,ii) = lnk0(i)
                end do
                V(ii) = V0
                L(ii) = 1.0d0-V(ii)
                !write(*,*) '------------'
                !write(*,*) day,ii
                !write(*,*) phase
                !write(*,*) 'VVVVV',V
    
                !write(*,*) lumda
                    !write(*,*) ii,'grid',phase_judge(ii),'phase'
            end do
            

        !!mainの流動計算
        write(30,*) phase    
        if (year < 4) then
            q_judge = 1 !流量制御
        else
            q_judge = 0 !孔底圧力制御
        end if
        q_judge = 1 !!とりあえず流量制御で固定

        do iteration=1,100
            call main_calc(V,lnk,Nc,Ncold,Nm,Nmold,Nmini,P,Pold,Pb0,fai,fai000,q_judge,phase_judge,phase,&
                            Swd,krgd,krwd,Sw,wc,chemi_mat,theta0,fxs)

            deallocate(gmat,hmat)
            allocate(gmat(n*eq+q_judge,n*eq+q_judge),hmat(n*eq+q_judge))
            call jacobian_mat(fxs,gmat)
            call outxs(fxs,hmat)
            hmat=-hmat


    !        do i=1,n*eq+q_judge
            !    write(13,*) (gmat(i,j),j=1,n*eq+q_judge),hmat(i)
    !        end do
            do i=1,n

                do j=1,com_2phase
                    if (Nc(j,i) == 0.0d0) then
                        do jj=1,eq
                            gmat(i*eq-eq+j,jj) = 0.0d0
                            gmat(jj,i*eq-eq+j) = 0.0d0
                        end do
                        gmat(i*eq-eq+j,i*eq-eq+j) = 1.0d0 !?存在しない成分の熱力学平衡条件式は計算しない？
                        hmat(i*eq-eq+j) = 0.0d0
                    end if
                end do

                do j=1,com_2phase+com_ion
                    !write(*,*) Nc(j,i),j,i
                    if (Nc(j,i) == 0.0d0) then
                        !write(*,*) j,i
                        do jj=1,eq
                            gmat(i*eq-eq+j+com_2phase,jj) = 0.0d0
                            gmat(jj,i*eq-eq+j+com_2phase) = 0.0d0
                        end do
                        gmat(i*eq-eq+j+com_2phase,i*eq-eq+j+com_2phase) = 1.0d0 !?存在しない成分の物質収支は計算しない？
                        hmat(i*eq-eq+j+com_2phase) = 0.0d0
                    end if
                end do
                !write(*,*) day,i,phase_judge(i),hmat(i*eq)

                if (phase_judge(i) /= 2) then
                    do jj=1,eq
                        gmat(i*eq,jj) = 0.0d0
                        gmat(jj,i*eq) = 0.0d0
                    end do
                    gmat(i*eq,i*eq) = 1.0d0 !?1相の時は、rachford解かない
                    hmat(i*eq) = 0.0d0


                    do j=1,com_2phase!+com_ion
                        do jj=1,eq
                            gmat(i*eq-eq+j,jj) = 0.0d0
                            gmat(jj,i*eq-eq+j) = 0.0d0
                        end do
                        gmat(i*eq-eq+j,i*eq-eq+j) = 1.0d0 !?1相の時は、lnk解かない
                        hmat(i*eq-eq+j) = 0.0d0
                    end do
                end if 
                !write(*,*) day,i,phase_judge(i),hmat(i*eq)
            end do

            do i=1,n*eq+q_judge
                !write(13,*) !(gmat(i,j),j=1,n*eq+q_judge),hmat(i)
            end do
            

            call pvgauss(n*eq+q_judge,gmat,hmat)
            
    !        do i=1,n*eq+q_judge
            !    write(13,*) (gmat(i,j),j=1,n*eq+q_judge),hmat(i)
    !        end do
            

            !!更新
            do i=1,n
                do j=1,com_2phase
                    lnk(j,i) = lnk(j,i) + hmat(i*eq-eq+j)
                end do
                do j=1,com_2phase+com_ion
                    Nc(j,i) = Nc(j,i) + hmat(i*eq-eq+com_2phase+j)
                end do
                do j=1,com_mine
                    Nm(j,i) = Nm(j,i) + hmat(i*eq-eq+com_2phase+com_2phase+com_ion+j)
                end do
                P(i) = P(i) + hmat(i*eq-1)
                V(i) = V(i) + hmat(i*eq)
                write(*,*) day,'day',i,phase(i),V(i)
            end do
            if (q_judge == 1) then
                Pb0 = Pb0 + hmat(n*eq+q_judge)
            end if 
            error = sqrt(dot_product(hmat,hmat))

            !write(*,*) iteration,error
            if (error < 10.0d0**(-3.0d0)) then
                exit
            end if

            

        end do !iteration loop
        
        

        Pold = P
        Ncold = Nc
        Nmold = Nm
        do i =1,n
            Nt = 0.0d0
            do j =1,com_2phase+com_ion
                Nt = Nt + Nc(j,i)
            end do
            do j=1,com_2phase+com_ion
                z(j,i) = Nc(j,i) / Nt 
            end do
        end do

        !?タイムループ回った！次回からは検証に向けて条件整える
    end do !hour loop
    !do i=1,n
    !        write(*,*) day,iteration!,P(i)
    !end do
    
    !write(*,*) day,error,iteration,phase!,Nc(1,1),Nc(2,1),P(1)
    
    
    do i=1,n
        write(30+i,*) P(i)
        write(35+i,*) Sw(i)
        write(40+i,*) wc(1,i)
        write(45+i,*) wc(2,i)
        write(50+i,*) wc(3,i)
        write(55+i,*) wc(4,i)
        write(60+i,*) wc(5,i)
        write(65+i,*) wc(6,i)
        write(70+i,*) wc(7,i)
        write(75+i,*) wc(8,i)
        write(80+i,*) wc(10,i)
        write(85+i,*) wc(11,i)
        write(90+i,*) Nm(2,i)*dx*dy*dz*(1.0d0-fai000(i))
        write(95+i,*) Nm(4,i)*dx*dy*dz*(1.0d0-fai000(i))
        !write(*,*) fai000(i)!#TODOSiO2おかしいかも 
        write(200+i,*) theta0(7,i)

    end do


    



        
        !do i=1,n
        !        write(*,*) day,iteration!,P(i)
        !end do
        
        
        
        
            
        end do !day loop
        
        end do !year loop
        deallocate(amat,bmat,gmat,hmat)
    
    
    write(*,*) 'finish!!!'
end program  main
