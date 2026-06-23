!========================================!
!      CHANNEL STAGGERED BI_PARALLEL     !
!========================================!
! 程序目的：使用合成湍流方法（STG）生成入口边界速度场，支持MPI并行。
! 说明：本程序基于von Kármán-Pao能谱生成各向同性湍流，并通过插值得到入口面的速度时间序列。

    PROGRAM MAIN
    
        !----------------------------------------
            
            INCLUDE 'head.fi'          ! 包含头文件，定义全局变量和数组尺寸
            INCLUDE 'parameter.h'      ! 包含参数文件，定义NPR, NPC, NPP等并行参数
        
            
            REAL*8 T0,T1               ! 用于CPU计时
        !----------------------------------------
        !  execute statement
        !----------------------------------------
        ! 并行初始化并划分笛卡尔拓扑：
            CALL MPI_Init( iMPI_ErrorInfo )                     ! 初始化MPI环境
            CALL MPI_Comm_Rank(MPI_COMM_WORLD, iMPI_MyID , iMPI_ErrorInfo)   ! 获取当前进程ID
            CALL MPI_Comm_Size(MPI_COMM_WORLD, iMPI_NumProcs , iMPI_ErrorInfo) ! 获取总进程数
            DIMS(0)= NPR                                        ! 笛卡尔网格的行数（x方向）
            DIMS(1)= NPC                                        ! 笛卡尔网格的列数（y方向）
            PERIODIC(0) = .FALSE.                               ! x方向非周期性（入口出口边界）
            PERIODIC(1) = .FALSE.                               ! y方向非周期性
            CALL MPI_CART_CREATE(MPI_COMM_WORLD, 2, DIMS, PERIODIC,.FALSE.,MPI_COMM_CART,iMPI_ErrorInfo) ! 创建2D笛卡尔通信子
            CALL MPI_CART_COORDS(MPI_COMM_CART, iMPI_MyID, 2, COORD,iMPI_ErrorInfo) ! 获取当前进程在笛卡尔网格中的坐标
        
            ! 提取列通信子（沿x方向，即每列共享y,z平面）
            CALL MPI_CART_SUB(MPI_COMM_CART, (/.TRUE.,.FALSE./), MPI_COMM_COL, iMPI_ErrorInfo)
            ! 提取行通信子（沿y方向，即每行共享x,z平面）
            CALL MPI_CART_SUB(MPI_COMM_CART, (/.FALSE.,.TRUE./), MPI_COMM_ROW, iMPI_ErrorInfo)
            ! 检查进程数是否与预期一致（NPP = NPR * NPC）
            !!!NP=NPROW*NPCOL
            IF ( NPP .EQ. iMPI_NumProcs ) THEN
        
            !   进程数正确，继续执行
        
                CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )   ! 全局同步
        
            !   主程序调用
                CALL CPU_TIME(T0)                                   ! 开始计时
                CALL MAINPROGRAM
                CALL CPU_TIME(T1)                                   ! 结束计时
                WRITE(*,*) iMPI_MyID,T1-T0                          ! 输出每个进程的运行时间
        
            ENDIF
        
            !Print*, 'FINISHED'
        
            CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )       ! 同步后结束
            CALL MPI_Finalize( iMPI_ErrorInfo )                     ! 终止MPI环境
        
        !----------------------------------------
        
            END PROGRAM
        
        !===== END PROGRAM ======================
        
        !-------------------------------------------------------------------------
        !  MAIN PROGRAM
        !-------------------------------------------------------------------------
        ! 主程序子程序：生成湍流速度场，计算入口速度时间序列，并输出结果文件
        
            SUBROUTINE MAINPROGRAM
            
        !----------------------------------------
            
            INCLUDE 'head.fi'          ! 包含全局变量定义
            INCLUDE 'parameter.h'      ! 包含参数定义
        
            Real*8 M_PI                ! 圆周率
            Real*8 km0,lx,ly,lz,dx,hdx,dy,hdy,dz,hdz   ! 几何参数：计算域尺寸，网格步长及半长
            Real*8 phi(0:nModes),tnu(0:nModes),theta(0:nModes),psi(0:nModes) ! 随机角度和相位
            integer i,j,k,m,nall,nx_st,nx_ed,ip,NNT    ! 循环索引和局部变量
            Real*8 Ran,kmmax,dk,km(0:nModes)           ! 波数相关变量
            Real*8 kx(0:nModes),ky(0:nModes),kz(0:nModes)        ! 波数向量分量
            Real*8 ktx(0:nModes),kty(0:nModes),ktz(0:nModes)     ! 修正后的波数（交错网格）
            Real*8 sxm(0:nModes),sym(0:nModes),szm(0:nModes)     ! 速度幅值向量（含能谱缩放）
            Real*8 zetax,zetay,zetaz,smag,Ekm,espec    ! 临时变量：单位方向，能谱值，幅值缩放
            Real*8 xc(nx),yc(ny),zc(nz)                ! 网格中心坐标（交错网格）
            Real*8 u(nxp,ny,nz),v(nxp,ny,nz),w(nxp,ny,nz)  ! 速度场（每个进程存储nxp个x平面）
            Real*8 arg                            ! 相位参数


            REAL*8 XREF0,YREF0,ZREF0,XREF1              ! 参考坐标，用于定位入口平面
            REAL*8 INLET_X(NXALL),INLET_Y(NXALL),INLET_Z(NXALL)   ! 入口网格点的全局坐标（读取）
            REAL*8 INLET_XP(NXALL),INLET_YP(NXALL),INLET_ZP(NXALL) ! 移动后的入口点坐标
            REAL*8 TIME_ST
            REAL*8 UIN(3,NXALL),R
            REAL*8 TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6     ! 插值临时变量
            REAL*8 UTIMEH(3,NTALL),UTIME(3,NTALL)          ! 时间序列（局部和全局归约）
            REAL*8 USUM0(3,NXALL),USUM(3,NXALL),UIN_NEW(3,NXALL) ! 平均速度和修正后的入口速度
            INTEGER IA,JA,KA                               ! 网格索引
            CHARACTER*20 CHAR1                             ! 输出文件名序号
            CHARACTER*10 CHAR3
          
          
        
            M_PI=4.d0*atan(1.d0)       ! 定义圆周率
            nall=nx*ny*nz              ! 总网格点数（单个进程）
            
            ! 根据N（无量纲化参数）和 L(2L*2L包围INLEI1) 指定计算域尺寸 
            ! 湍流块 总体积 lx*ly*lz
            lx = dble(nx/N)*L ! 25*L
            ly = dble(ny/N)*L ! 2*L
            lz = dble(nz/N)*L ! 2*L
           
            ! 最小波数(最大波长L)，一般设为2π/L（取三个方向的最大值）
            km0 = max(2.0*M_PI/lx,2.0*M_PI/ly,2.0*M_PI/lz)
        
            ! 计算网格步长及半长 (L/N)
            dx  = lx/nx
            hdx = dx/2.0
            dy  = ly/ny
            hdy = dy/2.0
            dz  = lz/nz
            hdz = dz/2.0
           
            ! 生成网格中心坐标
            Do i=1,nx
            xc(i) = hdx + (i-1)*dx
            End do
        
            Do j=1,ny
            yc(j) = hdy + (j-1)*dy
            End do
        
            Do k=1,nz
            zc(k) = hdz + (k-1)*dz
            End do
          
            ! 读取入口网格点坐标（文件：ReadData/INLET1.plt）
            OPEN(1,FILE='ReadData/INLET1.plt',STATUS='OLD')   
            DO I=1,NXALL 
                READ(1,*) INLET_X(I),INLET_Y(I),INLET_Z(I)
            ENDDO  
            CLOSE(1)
            
            ! 定义参考点：出口附近某点（nx-5），y和z方向中间点
            XREF0=xc(nx-5)
            YREF0=0.5*(yc(ny/2)+yc(ny/2+1))
            ZREF0=0.5*(zc(nz/2)+zc(nz/2+1))
            


            
            ! 计算每个进程的起始x索引：根据时间推进，入口平面随时间移动（泰勒冻结假设）
            ! INLET_X(1)是max(INLET_X(:))
            XREF1=XREF0+INLET_X(1)-iMPI_MyID*DT*NTP*UREF
          
            ! 检查XREF1是否小于0（越界错误）
            IF (XREF1 < 0.0d0) THEN
                WRITE(*,*) 'ERROR: Process', iMPI_MyID, 'XREF1 =', XREF1, '< 0'
                CALL MPI_ABORT(MPI_COMM_WORLD, 1, iMPI_ErrorInfo)
            END IF

            ! 计算每个进程负责的x范围（nx_st到nx_ed），确保覆盖入口平面移动后的区域
            nx_ed= FLOOR((XREF1+hdx)/dx)+3
            nx_st=nx_ed-nxp+1
        
        
          
          
          
            
            ! 生成随机数种子，仅0号进程生成随机相位phi, theta, psi，然后广播
            IF (iMPI_MyID.EQ.0) THEN
            call random_seed() 
            Do i=0,nModes
            call random_number(ran) 
            phi(i)=2.0*M_PI*ran
            call random_number(ran)
            tnu(i)=ran
            theta(i)=acos(2.0*tnu(i) -1.0)
            call random_number(ran)
            psi(i)   = M_PI*  ran - M_PI/2.0
            End do    
            ENDIF
            call mpi_bcast( phi, nModes+1, mpi_real8, 0, mpi_comm_world, impi_errorinfo )
            call mpi_bcast( theta, nModes+1, mpi_real8, 0, mpi_comm_world, impi_errorinfo )
            call mpi_bcast( psi, nModes+1, mpi_real8, 0, mpi_comm_world, impi_errorinfo )
            CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
            
            ! 最大波数由网格分辨率决定（Nyquist频率）
            kmmax = M_PI/dx;
            IF (iMPI_MyID.EQ.0) THEN
            write(*,*) "I will generate data up to wave number: ", kmmax
            ENDIF
            dk = (kmmax-km0)/nModes          ! 波数间隔
        
            ! 创建波数数组km(i)
            Do i=0,nModes
            km(i) = km0 + i*dk
            End do
        
            ! 根据球面角生成波数向量分量
            Do i=0,nModes
            kx(i) = sin(theta(i))*cos(phi(i))*km(i)
            ky(i) = sin(theta(i))*sin(phi(i))*km(i)
            kz(i) = cos(theta(i))*km(i)
            End do
        
            ! 修正波数（交错网格效应，sin(k*dx/2)/dx）
            Do i=0,nModes
            ktx(i) = sin(kx(i)*hdx)/dx
            kty(i) = sin(ky(i)*hdy)/dy
            ktz(i) = sin(kz(i)*hdz)/dz
            End do
             
            ! 重新生成随机角度phi和theta（用于生成涡量方向）
            IF (iMPI_MyID.EQ.0) THEN
            Do i=0,nModes
            call random_number(ran)
            phi(i) = 2.0*M_PI* ran
            call random_number(ran)
            tnu(i)    = ran
            theta(i) = acos(2.0*tnu(i) -1.0);
            End do    
            ENDIF
            call mpi_bcast( phi, nModes+1, mpi_real8, 0, mpi_comm_world, impi_errorinfo )
            call mpi_bcast( theta, nModes+1, mpi_real8, 0, mpi_comm_world, impi_errorinfo )
            CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
            
            ! 计算每个模式的涡量方向单位向量，并取与修正波数的叉积得到速度方向
            Do i=0,nModes
            zetax = sin(theta(i))*cos(phi(i));
            zetay = sin(theta(i))*sin(phi(i));
            zetaz = cos(theta(i));
            ! 计算叉积：zeta × k_tilde
            sxm(i) =  (zetay*ktz(i) - zetaz*kty(i))
            sym(i) = -(zetax*ktz(i) - zetaz*ktx(i))
            szm(i) =  (zetax*kty(i) - zetay*ktx(i))
            ! 单位化
            smag = 1.0/sqrt(sxm(i)*sxm(i) + sym(i)*sym(i) + szm(i)*szm(i))
            sxm(i) = sxm(i)*smag
            sym(i) = sym(i)*smag
            szm(i) = szm(i)*smag
            End do
        
        
            ! 输出von Kármán-Pao能谱曲线到Tecplot文件
            IF (iMPI_MyID.EQ.0) OPEN(70,FILE='Tecplot_InputFiles/von Karman-Pao.plt',STATUS='UNKNOWN')
            WRITE(70, '(A)') 'TITLE = "von Karman-Pao"'
            WRITE(70, '(A)') 'VARIABLES = "k", "E"'
            WRITE(70, '(A, I0, A)') 'ZONE T="von Karman-Pao", I=', nModes+1, ', F=POINT'
            ! 更新速度幅值向量：根据能谱值缩放
            Do i=0,nModes
            call karman(km(i),Ekm)     ! 计算能谱值E(k)
        
            IF (iMPI_MyID.EQ.0) WRITE(70,101) km(i),Ekm
        
            espec = 2.0*sqrt(Ekm * dk)  ! 幅值 = 2*sqrt(E(k)*Δk)
            sxm(i) =sxm(i)* espec
            sym(i) =sym(i)* espec
            szm(i) =szm(i)* espec
            End do
            IF (iMPI_MyID.EQ.0) close(70)
        
          
        
            u=0.0
            v=0.0
            w=0.0
        
            ! 合成湍流速度场：对所有模式求和，每个模式贡献一个傅里叶模态
            Do k=1,nz
            Do j=1,ny
            Do i=nx_st,nx_ed
            ip=i-nx_st+1
            Do m=0,nModes
               arg=kx(m)*xc(i) + ky(m)*yc(j) + kz(m)*zc(k) - psi(m) ! 相位
               u(ip,j,k) =u(ip,j,k)+ cos(arg)*sxm(m)
               v(ip,j,k) =v(ip,j,k)+ cos(arg)*sym(m)
               w(ip,j,k) =w(ip,j,k)+ cos(arg)*szm(m)  
            end do
            end do
            end do
            end do
           
        
        
            ! 计算第一段速度时间序列并输出平均场
            USUM0=0.0
            DO NNT=1,NTP
            DO K=1,N1
            DO J=1,N1
            I=N1*(K-1)+J
            ! 根据泰勒冻结假设移动入口平面
            INLET_XP(I)=INLET_X(I)+XREF0-(iMPI_MyID*DT*NTP+NNT*DT)*UREF
            INLET_YP(I)=INLET_Y(I)+YREF0
            INLET_ZP(I)=INLET_Z(I)+ZREF0
            ! 找到移动后的点所在的网格单元（三线性插值）
            IA=FLOOR((INLET_XP(I)+hdx)/dx)
            ip=IA-nx_st+1
            IF (ip < 1 .OR. ip+1 > nxp) THEN
                WRITE(*,*) 'ERROR ip out of range:', iMPI_MyID, ip, IA, nx_st, nx_ed
                CALL MPI_ABORT(MPI_COMM_WORLD, 1, iMPI_ErrorInfo)
            ENDIF
            JA=FLOOR((INLET_YP(I)+hdy)/dy)
            KA=FLOOR((INLET_ZP(I)+hdz)/dz)
            ! 插值u分量
            TEMP1=( (INLET_XP(I)-xc(IA))*u(ip+1,JA  ,KA)+(xc(IA+1)-INLET_XP(I))*u(ip,JA  ,KA) )/( xc(IA+1)-xc(IA) )
            TEMP2=( (INLET_XP(I)-xc(IA))*u(ip+1,JA+1,KA)+(xc(IA+1)-INLET_XP(I))*u(ip,JA+1,KA) )/( xc(IA+1)-xc(IA) )
            TEMP3=( (INLET_YP(I)-yc(JA))*TEMP2+(yc(JA+1)-INLET_YP(I))*TEMP1 )/( yc(JA+1)-yc(JA) )
            TEMP4=( (INLET_XP(I)-xc(IA))*u(ip+1,JA  ,KA+1)+(xc(IA+1)-INLET_XP(I))*u(ip,JA  ,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP5=( (INLET_XP(I)-xc(IA))*u(ip+1,JA+1,KA+1)+(xc(IA+1)-INLET_XP(I))*u(ip,JA+1,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP6=( (INLET_YP(I)-yc(JA))*TEMP5+(yc(JA+1)-INLET_YP(I))*TEMP4 )/( yc(JA+1)-yc(JA) )	
            UIN(1,I)=( (INLET_ZP(I)-zc(KA))*TEMP6+(zc(KA+1)-INLET_ZP(I))*TEMP3 )/( zc(KA+1)-zc(KA) )
            ! 插值v分量
            TEMP1=( (INLET_XP(I)-xc(IA))*v(ip+1,JA  ,KA)+(xc(IA+1)-INLET_XP(I))*v(ip,JA  ,KA) )/( xc(IA+1)-xc(IA) )
            TEMP2=( (INLET_XP(I)-xc(IA))*v(ip+1,JA+1,KA)+(xc(IA+1)-INLET_XP(I))*v(ip,JA+1,KA) )/( xc(IA+1)-xc(IA) )
            TEMP3=( (INLET_YP(I)-yc(JA))*TEMP2+(yc(JA+1)-INLET_YP(I))*TEMP1 )/( yc(JA+1)-yc(JA) )
            TEMP4=( (INLET_XP(I)-xc(IA))*v(ip+1,JA  ,KA+1)+(xc(IA+1)-INLET_XP(I))*v(ip,JA  ,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP5=( (INLET_XP(I)-xc(IA))*v(ip+1,JA+1,KA+1)+(xc(IA+1)-INLET_XP(I))*v(ip,JA+1,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP6=( (INLET_YP(I)-yc(JA))*TEMP5+(yc(JA+1)-INLET_YP(I))*TEMP4 )/( yc(JA+1)-yc(JA) )	
            UIN(2,I)=( (INLET_ZP(I)-zc(KA))*TEMP6+(zc(KA+1)-INLET_ZP(I))*TEMP3 )/( zc(KA+1)-zc(KA) )
            ! 插值w分量
            TEMP1=( (INLET_XP(I)-xc(IA))*w(ip+1,JA  ,KA)+(xc(IA+1)-INLET_XP(I))*w(ip,JA  ,KA) )/( xc(IA+1)-xc(IA) )
            TEMP2=( (INLET_XP(I)-xc(IA))*w(ip+1,JA+1,KA)+(xc(IA+1)-INLET_XP(I))*w(ip,JA+1,KA) )/( xc(IA+1)-xc(IA) )
            TEMP3=( (INLET_YP(I)-yc(JA))*TEMP2+(yc(JA+1)-INLET_YP(I))*TEMP1 )/( yc(JA+1)-yc(JA) )
            TEMP4=( (INLET_XP(I)-xc(IA))*w(ip+1,JA  ,KA+1)+(xc(IA+1)-INLET_XP(I))*w(ip,JA  ,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP5=( (INLET_XP(I)-xc(IA))*w(ip+1,JA+1,KA+1)+(xc(IA+1)-INLET_XP(I))*w(ip,JA+1,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP6=( (INLET_YP(I)-yc(JA))*TEMP5+(yc(JA+1)-INLET_YP(I))*TEMP4 )/( yc(JA+1)-yc(JA) )	
            UIN(3,I)=( (INLET_ZP(I)-zc(KA))*TEMP6+(zc(KA+1)-INLET_ZP(I))*TEMP3 )/( zc(KA+1)-zc(KA) )
            ENDDO
            ENDDO
            
            ! 存储每个时刻中心点速度（用于时间序列分析）
            UTIMEH(:,iMPI_MyID*NTP+NNT)=UIN(:,NXALL/2)
            
            USUM0=USUM0+UIN   ! 累加时间平均
            
            ENDDO
            USUM0=USUM0/DBLE(NTP)   ! 计算时间平均（湍流部分）
          
            ! 归约所有进程的UTIMEH和USUM0到0号进程
            CALL MPI_REDUCE(UTIMEH,UTIME,3*NTALL,mpi_real8,MPI_SUM,0,MPI_COMM_WORLD,iMPI_ErrorInfo)
            CALL MPI_REDUCE(USUM0,USUM,3*NXALL,mpi_real8,MPI_SUM,0,MPI_COMM_WORLD,iMPI_ErrorInfo) 
            CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )
            call mpi_bcast( USUM,3*NXALL, mpi_real8, 0, MPI_COMM_WORLD, iMPI_ErrorInfo ) ! 广播平均场给所有进程

            USUM = USUM/DBLE(NPP)
            
            ! ! 0号进程输出平均场（湍流部分）
            ! IF (iMPI_MyID.EQ.0) THEN
            ! OPEN(1,FILE='Tecplot_InputFiles/UMEAN0_INLET.plt',STATUS='UNKNOWN')  
            ! WRITE(1, '(A)') 'TITLE = "UMEAN0_INLET"'
            ! WRITE(1, '(A)') 'VARIABLES = "X", "Y", "Z","U","V","W"'
            ! WRITE(1, '(A, I0, A, I0, A)') 'ZONE T="UMEAN0_INLET", I=', N1, ', J=', N1, ', F=POINT'
            ! DO I=1,NXALL
            !      WRITE(1,101) INLET_X(I),INLET_Y(I),INLET_Z(I),USUM(1,I),USUM(2,I),USUM(3,I)
            ! ENDDO
            ! CLOSE(1)
            ! ENDIF
            
            ! ! 0号进程输出时间序列（湍流部分，中心点）
            ! IF (iMPI_MyID.EQ.0) THEN  
            ! OPEN(71,FILE='Tecplot_InputFiles/UTIME0_INLET.plt',STATUS='UNKNOWN')
            ! WRITE(71, '(A)') 'TITLE = "UTIME0_INLET"'
            ! WRITE(71, '(A)') 'VARIABLES = "T","u","v","w"'
            ! WRITE(71, '(A, I0, A)') 'ZONE T="UTIME0_INLET", I=', NTALL, ', F=POINT'
            ! DO i=1,NTALL
            ! WRITE(71,101) DT*DBLE(i),UTIME(1,i),UTIME(2,i),UTIME(3,i)
            ! end do
            ! close(71)
            ! ENDIF
            
            
   
            

            ! *****************************************************************
            USUM0=0.0
            ! 第二段：生成修正后的入口速度（加上平均流U_RECT，并减去第一段的平均场）
            DO NNT=1,NTP
            DO K=1,N1
            DO J=1,N1
            I=N1*(K-1)+J
            INLET_XP(I)=INLET_X(I)+XREF0-(iMPI_MyID*DT*NTP+NNT*DT)*UREF
            INLET_YP(I)=INLET_Y(I)+YREF0
            INLET_ZP(I)=INLET_Z(I)+ZREF0
            ! 同样的三线性插值，但此时速度场为湍流部分u,v,w
            IA=FLOOR((INLET_XP(I)+hdx)/dx)
            ip=IA-nx_st+1
            JA=FLOOR((INLET_YP(I)+hdy)/dy)
            KA=FLOOR((INLET_ZP(I)+hdz)/dz)
            TEMP1=( (INLET_XP(I)-xc(IA))*u(ip+1,JA  ,KA)+(xc(IA+1)-INLET_XP(I))*u(ip,JA  ,KA) )/( xc(IA+1)-xc(IA) )
            TEMP2=( (INLET_XP(I)-xc(IA))*u(ip+1,JA+1,KA)+(xc(IA+1)-INLET_XP(I))*u(ip,JA+1,KA) )/( xc(IA+1)-xc(IA) )
            TEMP3=( (INLET_YP(I)-yc(JA))*TEMP2+(yc(JA+1)-INLET_YP(I))*TEMP1 )/( yc(JA+1)-yc(JA) )
            TEMP4=( (INLET_XP(I)-xc(IA))*u(ip+1,JA  ,KA+1)+(xc(IA+1)-INLET_XP(I))*u(ip,JA  ,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP5=( (INLET_XP(I)-xc(IA))*u(ip+1,JA+1,KA+1)+(xc(IA+1)-INLET_XP(I))*u(ip,JA+1,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP6=( (INLET_YP(I)-yc(JA))*TEMP5+(yc(JA+1)-INLET_YP(I))*TEMP4 )/( yc(JA+1)-yc(JA) )	
            UIN(1,I)=( (INLET_ZP(I)-zc(KA))*TEMP6+(zc(KA+1)-INLET_ZP(I))*TEMP3 )/( zc(KA+1)-zc(KA) )  + UREF
            TEMP1=( (INLET_XP(I)-xc(IA))*v(ip+1,JA  ,KA)+(xc(IA+1)-INLET_XP(I))*v(ip,JA  ,KA) )/( xc(IA+1)-xc(IA) )
            TEMP2=( (INLET_XP(I)-xc(IA))*v(ip+1,JA+1,KA)+(xc(IA+1)-INLET_XP(I))*v(ip,JA+1,KA) )/( xc(IA+1)-xc(IA) )
            TEMP3=( (INLET_YP(I)-yc(JA))*TEMP2+(yc(JA+1)-INLET_YP(I))*TEMP1 )/( yc(JA+1)-yc(JA) )
            TEMP4=( (INLET_XP(I)-xc(IA))*v(ip+1,JA  ,KA+1)+(xc(IA+1)-INLET_XP(I))*v(ip,JA  ,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP5=( (INLET_XP(I)-xc(IA))*v(ip+1,JA+1,KA+1)+(xc(IA+1)-INLET_XP(I))*v(ip,JA+1,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP6=( (INLET_YP(I)-yc(JA))*TEMP5+(yc(JA+1)-INLET_YP(I))*TEMP4 )/( yc(JA+1)-yc(JA) )	
            UIN(2,I)=( (INLET_ZP(I)-zc(KA))*TEMP6+(zc(KA+1)-INLET_ZP(I))*TEMP3 )/( zc(KA+1)-zc(KA) )
            TEMP1=( (INLET_XP(I)-xc(IA))*w(ip+1,JA  ,KA)+(xc(IA+1)-INLET_XP(I))*w(ip,JA  ,KA) )/( xc(IA+1)-xc(IA) )
            TEMP2=( (INLET_XP(I)-xc(IA))*w(ip+1,JA+1,KA)+(xc(IA+1)-INLET_XP(I))*w(ip,JA+1,KA) )/( xc(IA+1)-xc(IA) )
            TEMP3=( (INLET_YP(I)-yc(JA))*TEMP2+(yc(JA+1)-INLET_YP(I))*TEMP1 )/( yc(JA+1)-yc(JA) )
            TEMP4=( (INLET_XP(I)-xc(IA))*w(ip+1,JA  ,KA+1)+(xc(IA+1)-INLET_XP(I))*w(ip,JA  ,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP5=( (INLET_XP(I)-xc(IA))*w(ip+1,JA+1,KA+1)+(xc(IA+1)-INLET_XP(I))*w(ip,JA+1,KA+1) )/( xc(IA+1)-xc(IA) )
            TEMP6=( (INLET_YP(I)-yc(JA))*TEMP5+(yc(JA+1)-INLET_YP(I))*TEMP4 )/( yc(JA+1)-yc(JA) )	
            UIN(3,I)=( (INLET_ZP(I)-zc(KA))*TEMP6+(zc(KA+1)-INLET_ZP(I))*TEMP3 )/( zc(KA+1)-zc(KA) )
            ENDDO
            ENDDO
            
            ! 减去之前计算的平均湍流场，得到纯脉动速度（保证脉动时间均值为0）
            UIN_NEW=UIN-USUM
            
            ! 存储中心点时间序列（修正后）
            UTIMEH(:,iMPI_MyID*NTP+NNT)=UIN_NEW(:,NXALL/2)
            
            ! 输出每个时刻的入口速度场（所有点）
            WRITE(CHAR1,'(I6.6)') iMPI_MyID*NTP+NNT
            OPEN(1,FILE='vel/u'//TRIM(ADJUSTL(CHAR1))//'.plt',STATUS='UNKNOWN')  
            DO I=1,NXALL
                WRITE(1,230) UIN_NEW(1,I),UIN_NEW(2,I),UIN_NEW(3,I)
            ENDDO
            CLOSE(1)
        230 FORMAT(3(E22.15,2X))
        
            ! 如果是0号进程且第一个时刻，输出修正后的入口速度场（包含坐标）
            IF (iMPI_MyID.EQ.0.AND.NNT.EQ.1) THEN
                OPEN(1,FILE='Tecplot_InputFiles/UModify_INLET.plt',STATUS='UNKNOWN')  
                WRITE(1, '(A)') 'TITLE = "UModify_INLET"'
                WRITE(1, '(A)') 'VARIABLES = "X", "Y", "Z","u","v","w"'
                WRITE(1, '(A, I0, A, I0, A)') 'ZONE T="UModify_INLET", I=', N1, ', J=', N1, ', F=POINT'
                DO I=1,NXALL
                    WRITE(1,101) INLET_X(I),INLET_Y(I),INLET_Z(I),UIN_NEW(1,I),UIN_NEW(2,I),UIN_NEW(3,I)
                ENDDO
                CLOSE(1)
            ENDIF
            
        
            USUM0=USUM0+UIN_NEW   ! 累加修正后的速度（用于时间平均）
            ENDDO
            USUM0=USUM0/DBLE(NTP)   ! 修正后的时间平均

            
            ! 归约所有进程的时间序列和平均场
            CALL MPI_REDUCE(UTIMEH,UTIME,3*NTALL,mpi_real8,MPI_SUM,0,MPI_COMM_WORLD,iMPI_ErrorInfo)
            CALL MPI_REDUCE(USUM0,USUM,3*NXALL,mpi_real8,MPI_SUM,0,MPI_COMM_WORLD,iMPI_ErrorInfo) 
            CALL MPI_Barrier( MPI_COMM_WORLD, iMPI_ErrorInfo )

            
            USUM=USUM/DBLE(NPP)
            ! ! 0号进程输出修正后的平均场
            ! IF (iMPI_MyID.EQ.0) THEN
            !     OPEN(1,FILE='Tecplot_InputFiles/UMEAN_Modify_INLET.plt',STATUS='UNKNOWN')  
            !     WRITE(1, '(A)') 'TITLE = "UMEAN_Modify_INLET"'
            !     WRITE(1, '(A)') 'VARIABLES = "X", "Y", "Z","U","V","W"'
            !     WRITE(1, '(A, I0, A, I0, A)') 'ZONE T="UMEAN_Modify_INLET", I=', N1, ', J=', N1, ', F=POINT'
            !     DO I=1,NXALL
            !          WRITE(1,101) INLET_X(I),INLET_Y(I),INLET_Z(I),USUM(1,I),USUM(2,I),USUM(3,I)
            !     ENDDO
            !     CLOSE(1)
            ! ENDIF
            
        
            ! 0号进程输出修正后的时间序列（中心点）
            IF (iMPI_MyID.EQ.0) THEN  
                OPEN(71,FILE='Tecplot_InputFiles/UTIME_CenterLine_INLET.plt',STATUS='UNKNOWN')
                WRITE(71, '(A)') 'TITLE = "UTIME_CenterLine_INLET"'
                WRITE(71, '(A)') 'VARIABLES = "T","u","v","w"'
                WRITE(71, '(A, I0, A)') 'ZONE T="UTIME_CenterLine_INLET", I=', NTALL, ', F=POINT'
                DO i=1,NTALL
                WRITE(71,101) DT*DBLE(i),UTIME(1,i),UTIME(2,i),UTIME(3,i)
                end do
                close(71)
            ENDIF
            
            ! 输出每个进程的湍流速度场（三维，用于Tecplot可视化）
            IF (0.LE.iMPI_MyID.AND.iMPI_MyID.LE.5) THEN   ! 仅前6个进程输出（调试用）
            WRITE(CHAR1,'(I6.6)') iMPI_MyID
             OPEN(71,FILE='Tecplot_InputFiles/vel_3d/vel'//TRIM(ADJUSTL(CHAR1))//'.plt',STATUS='UNKNOWN')    
             WRITE(71, '(A)') 'TITLE = "Turb_block'//TRIM(ADJUSTL(CHAR1))//'"'
             WRITE(71, '(A)') 'VARIABLES = "X", "Y", "Z","u","v","w"'
             WRITE(71, '(A, I0, A, I0, A, I0, A)') 'ZONE T="Turb_block'//TRIM(ADJUSTL(CHAR1))//'", I=', nxp, ', J=', ny, ', K=',nz, ', F=POINT'
             DO k=1,nz
             DO j=1,ny
             DO ip=1,nxp
                i=ip+nx_st-1
             WRITE(71,101) xc(i),yc(j),zc(k),u(ip,j,k),v(ip,j,k),w(ip,j,k)
             end do
             end do
             end do
             close(71)
            ENDIF
            
        
        
        101 FORMAT(6(E22.15,2X))
        102 FORMAT('Zone T=vel i=',I0,' j=',I0,' k=',I0,' f=point')   
        103 FORMAT(4096(E22.15,2X)) 
        104 FORMAT('Zone T=inlet_t0 i=',I0,' j=',I0,' f=point') 
        105 FORMAT(4(E22.15,2X))  
        !-------------------------------------------------
        
            END SUBROUTINE
        
            
            !********output E(km)*****
            ! 子程序：计算von Kármán-Pao能谱值
            Subroutine karman(k,tke)
            INCLUDE 'parameter.h'
            Real*8 k,ke_,nu,urms
            Real*8 ke,Lin,alpha,epsilone,keta,r1,r2,tke
        
            !input--------------------------------------------------------
            ! 以下参数根据湍流强度、积分尺度和粘性系数计算能谱
            !integral length
            Lin=0.005 ! R/10
            nu=1.31e-6
            urms=UREF*0.6
           !--------------------------------------------------------------
           ! 根据输入参数计算能谱公式中的相关量
           ! 能谱峰值对应的波数
           ke= 0.746834/Lin
           ! 标度常数
           alpha = 1.452762113
           ! 湍流耗散率
           epsilone = urms*urms*urms/Lin
           ! Kolmogorov波数
           keta = epsilone**(0.25)*nu**(-3.0/4.0)
           
           r1 = k/ke
           r2 = k/keta
           ! von Kármán-Pao能谱公式
           tke = alpha*(urms*urms/ke)*(r1*r1*r1*r1 / ((1.0 + r1*r1)** (17.0/6.0)))*exp(-2.0*r2*r2)
           return 
           End