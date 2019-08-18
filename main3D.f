                        program main

      implicit none
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c                                                        
c  The driver program for the implicit quasi-static      
c  3D nonlinear finite element code.                     
c  Based on total lagrangian formulation                 
c  Book
c  1. Nonlinear finite elements for continua and structures 
c     by T. Belytschko, W.M. Liu and B. Moran
c  2. Nonlinear continuum mechanics for finite element analysis
c     by Bonet and Wood
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      






c=========================================================
c  Different material constitutive model
c  Generalized Money-Rivilin (compressible)    
c  Psi=c(J-1)^2-dlog(J)+c1(Ic-3)+c2(IIc-3)
c  Material Parameters::   c,c1,c2, d=2(c1+2c2)
c  MatType 11: plain strain   
    
c  Nearly incompressible Neo-Hookean Material 
c  Psi=mu/2(trace(C)-3)+1/2kappa(J-1)^2
c  Material Parameters :: kappa,mu
c  MatType 21: plain strain


c  Compressible Neo-Hookean Material
c  Psi=mu/2(Ic-3)-mulog(J)+lambda/2(log J)^2
c  Material Parameters:: lambda,mu
c  MatType 31: plain strain
c  MatType 32: plain stress


c  Linearly elasticity 
c  Material Parameters:: lambda,mu
c  MatType:40
c  MatType 41: plain strain
c  MatType 42: plain stress





c=========================================================


      include 'header2D.h'
c
c*********************************************************
c Open file static.inp, read the problem headers, and open
c other files also
c*********************************************************
c
c     write(*,*)'Opening static2D.inp'
      open(15, file = 'static2D.inp', status = 'old')
      read(15,20)
      read(15,40) infile

20    format(A100)      
40    format(A8)      

      inlen = 8
      do i = 1,8
          if (infile(i:i).ne.' ') then
              inlen = i
          endif
      enddo

c------------open basic input/output files--------------------------------

      open(12,file=infile(1:inlen)//'/'//infile(1:inlen)//'.mat',
     $                             status='old') !material data
c     write(*,*)'...Opened material data file'

      open(13,file=infile(1:inlen)//'/'//infile(1:inlen)//'.dat',
     $                             status='old') !mesh data
c        write(*,*)'...Opened mesh data file'

      open(14,file=infile(1:inlen)//'/'//infile(1:inlen)//'.load',
     $                             status='old') !load data
c     write(*,*)'...Opened load data file'

c      open(70,file=infile(1:inlen)//'/'//infile(1:inlen)//'.charge',
c     $                             status='old') !load data
c      write(*,*)'...Opened load data file'



      open(22,file='Output/'//infile(1:inlen)//'/Constv.'
     $                     //infile(1:inlen),status='unknown') !reaction

      open(23,file='Output/'//infile(1:inlen)//'/Cohe_flux.'
     $                     //infile(1:inlen),status='unknown') !reaction  


      open(123,file='Concentration.dat',status='unknown')

      open(220,file='Cohesive.dat',status='unknown')
      open(320,file='FracOutput.dat',status='unknown')
 
c*********************************************************
c Read the mesh data file and create required arrays
c*********************************************************


       read(13,20)
       read(13,*)numnp,numelv,numelc,numperiod,numimp,numatv,numatc,
     $           num_press,npmark





c*********************************************************
c    Control options
c*********************************************************
      read(15,20)
      read(15,*)outputStep         !output file steps
      read(15,20)
      read(15,*)CohMatchSwitch     !option for nonmatching cohesive
      read(15,20)
      read(15,20)
      read(15,*)LineSwitch         !option for linesearch
      read(15,20)
      read(15,*)nsear_max          !maximum linear search iteration
      read(15,20)
      read(15,*)isolve_switch      !Solver option ( symmetric isolve_switch=0)
                                   !              ( unsymmetric isolve_switch=1)
                                   !              (GMRES =2)

      read(15,20)                  !Solver Type (Implicit=1,Explicit=0) 
      read(15,*)solver_type
      read(15,20)                  !Tolerance value for L2 norm of disp
      read(15,*)epsilon_disp
      read(15,20)                  !Maximum iteration for convergenece
      read(15,*)icount_max
      read(15,20) 
      read(15,20) 
      read(15,*)IRVESwitch         !switch for RVE Analysis with Deformatiom
                                   ! Gradient as input
      read(15,20) 
      read(15,*)multiphysicsswitch   !switch for multiphysics        

      read(15,20) 
      read(15,*)cycloading        !Switch for multiphysics 
      read(15,20)
      read(15,*)Plasticity_switch !Switch if plasticity 


c*********************************************************
c Initialize all the arrays
c*********************************************************


c-----global variable allocation------------------

       allocate(x(1:3,1:numnp))
       allocate(lmv(1:np_elev_max+1,1:numelv))
       allocate(elem_type(1:numelv)) 
       allocate(lmc(1:np_elec_max+1,1:numelc))
       allocate(coh_elem_type(1:numelc)) 
       allocate(lmp(1:np_elep_max+1,1:num_press))
       allocate(pr_elem_type(1:num_press)) 
       allocate(node_imposed(1:2,1:numimp))
       allocate(MatType(1:numatv))  
       allocate(matvprop(1:numatv,1:6))
       allocate(matcprop(1:numatc,1:6))
       allocate(idof(1:2,1:3*numnp))
       allocate(xdof(1:3*numnp))
       allocate(dof_period(1:2,1:2*numperiod))
       allocate(bounmark(1:numnp))
       allocate(P(1:npmark)) 
       allocate(dP(1:npmark))





c global external load vector
       allocate(Rext(1:3*numnp))
       allocate(dRext(1:3*numnp))
       Rext = 0
       dRext = 0

      allocate(F_Macro_0(1:3,1:3))
      allocate(Fij_M(1:3,1:3))
      allocate(dFij_M(1:3,1:3))



c*************************************************************
       call readInput(numnp,numelv,numelc,numimp,numatv,numatc,
     $      num_press,lmp,np_elev_max,np_elec_max,np_elep_max,npmark,P,
     $     x,lmv,lmc,elem_type,coh_elem_type,pr_elem_type,node_imposed
     $     ,MatType,matvprop,matcprop,idof
     $     ,numperiod,dof_period,loading_disp_max_X,loading_disp_max_Y,
     $     loading_disp_max_Z,numloadSteps,
     $     bounmark,Rext,F_Macro_0)  

c************************************************************

      numeq = maxval(idof(1,:))
      numpres = 3*numnp-numeq

c----------Variable for displacement-------------------------------------
c global nodal displacement vector

       allocate(Uarr(1:3*numnp))
       allocate(Uarr_n(1:3*numnp))
       allocate(Uarr_k(1:3*numnp))
       allocate(dUarr_k(1:3*numnp))


      Uarr=0.0d0
      Uarr_n=0.0d0
      Uarr_k=0.0d0
      dUarr_k=0.0d0

c-------------------------------------------------------
      allocate(Fp_arr(1:numelv,1:max_gp,1:3,1:3))

      Fp_arr(1:numelv,1:max_gp,1,1)=1.0d0
      Fp_arr(1:numelv,1:max_gp,2,2)=1.0d0
      Fp_arr(1:numelv,1:max_gp,3,3)=1.0d0

      allocate(PlasticPara_arr(1:numelv,1:max_gp))
      PlasticPara_arr=0.0d0

c-----------------------------------------------------
c   Damage parameters
c-----------------------------------------------------

      allocate(Svar(1:numelv,1:max_gp,1:num_Svar_max))
      Svar=0.0d0


c------------------------------------------------------

      allocate(dUf(1:numeq))
      allocate(dUp(1:3*numpres))
      dUf=0.0d0


      allocate(Ru(1:3*numnp))


      allocate(Ru_f(1:numeq))
      allocate(Rf(1:numeq))

c----------------------------------------------------------------
c      cohesive parameters
c----------------------------------------------------------------

       allocate(CohesiveHistPara(1:numelc,1:max_gp,1:2))
       CohesiveHistPara=0.0d0
       Area_Fractured=0.0d0
       allocate(CohElemNodeIndex(1:numelc,1:2,1:np_elec_max+1))
       CohElemNodeIndex=0
       call CohesiveElementFaceConnectivity(numelv,numelc,lmv,lmc
     $      ,np_elec_max,np_elev_max,coh_elem_type,CohElemNodeIndex)
c============================================================
c       Data file for concentration           
c================================================================


        open(60,file=infile(1:inlen)//'/'//infile(1:inlen)//'.conc',
     $                             status='old') !read internal variable data

        read(60,20)
        read(60,*)numnp_q,numel_q,numimp_q,numat_q,numelvf,nfmark   ! number of imposed conc/temperature
        numeq_q=numnp_q-numimp_q                     ! number of equation for internal variable 

        allocate(node_q(1:numnp_q))                  ! node list for internal variable
        allocate(node_q_imp(1:numimp_q))             ! node list with imposed int. variable
        allocate(idof_q(1:numnp))                    ! dof for concentration
        allocate(Dq_coff(1:numat_q))                 ! diffusivity for each material
        allocate(flux(1:nfmark,1:3))
        allocate(flux_max(1:nfmark,1:3))        
        allocate(alpha(1:numat_q))
        allocate(Jq_vec(1:numel_q))                  ! flux vector 
        allocate(Qp(1:numimp_q))                     ! prescribed field var
        allocate(dQp(1:numimp_q))                    ! increment in prescribed field
        allocate(Qf(1:numeq_q))                      ! unprescribed field var
        allocate(dQf(1:numeq_q))                     ! increment in unprescribed field
        allocate(lmq(1:np_elev_max+1,1:numel_q))     ! connectivity for internal variable
        allocate(lmq_f(1:np_elec_max+1,1:numelvf))   ! Connectivity for flux boundary condition
        allocate(Rq(1:numnp))                        ! Global load vector for field only
        allocate(Rflux(1:numnp))
        allocate(elem_type_q(1:numel_q))
        allocate(elem_type_f(1:numelvf))  

        allocate(Qarr(1:numnp))                      ! Global array of field var  
        allocate(Qarr_n(1:numnp))                    ! Global array of field var previous step
        allocate(Qarr_k(1:numnp))                    ! Global array of field var previous iteration
        allocate(dQarr_k(1:numnp))                   ! Increment at k th iteration step
 


        allocate(Rq_f(1:numeq_q))
        allocate(Rq_fl(1:numeq_q)) 
        allocate(Rq_p(1:numimp_q))
    

        allocate(Grad_P(numel_q,1:max_gp,1:3))
        allocate(dp_dc(1:numnp_q,1))





      
        Qp=0.0d0
        dQp=0.0d0
        Qf=0.0d0
        dQf=0.0d0
        Dq_coff=0.d0
        Jq_vec=0.0d0
        Qarr=0.0d0
        idof_q=0




       call readInternalVar(numnp,numnp_q,numel_q,numelvf,
     $      nfmark,numimp_q,numat_q,node_q,node_q_imp,lmq,np_elev_max,
     $      elem_type_q,idof_q,Dq_coff,delta_t,alpha,flux,
     $      elem_type_f,lmq_f,np_elec_max,Jq_vec,Qp)

c      C_max=3.6547e-13 
c      write(*,*)flux,alpha,C_max
c      write(*,*)delta_t,Dq_coff

c       call IntercalationParameters(numat_q,nfmark,numloadSteps
c     $    ,ncycle,delta_t,C_max,Dq_coff,alpha,flux)

      write(*,*)flux,alpha,C_max
      write(*,*)delta_t,Dq_coff


c=====================================================================
c    SuperLU index for concentration
c=====================================================================

          allocate(Qval_index_ff(1:numeq_q,1:numeq_q))
          allocate(Qval_index_fp(1:numeq_q,1:numimp_q))

          allocate(Qrowind_temp(1:numeq_q*numeq_q))
          allocate(Qcolptr(1:numeq_q+1))

          Qrowind_temp=0
          Qcolptr=0
          Qval_index_ff=0
          Qval_index_fp=0

            call SuperLUMatrixIndex_conc(max_gp
     $      ,np_elev_max,numel_q,lmq
     $      ,elem_type_q,numnp_q,numeq_q,numimp_q,idof_q
     $      ,Qrowind_temp,Qcolptr,Qval_index_ff,Qval_index_fp,Qnnz
     $      ,Qnnzimp)



          allocate(Qrowind(1:Qnnz))
          allocate(Qvalues_ff(1:Qnnz))
          allocate(Qvalues_fp(1:Qnnzimp))
         Qrowind(1:Qnnz)=Qrowind_temp(1:Qnnz)
          deallocate(Qrowind_temp)

c==========================================================================
c         SuperLU index for displacement
c==========================================================================
          allocate(val_index_ff(1:numeq,1:numeq))
          allocate(val_index_fp(1:numeq,1:numpres))

          allocate(rowind_temp(1:numeq*numeq))
          allocate(colptr(1:numeq+1))

          rowind_temp=0
          colptr=0
          val_index_ff=0
          val_index_fp=0

          call SuperLUMatrixIndex(np_elev_max,numelv,lmv
     $      ,elem_type,numnp,numeq,numpres,idof
     $      ,rowind_temp,colptr,val_index_ff,val_index_fp,nnz,nnzimp)


          allocate(rowind(1:nnz))
          allocate(values_ff(1:nnz))
          allocate(values_fp(1:nnzimp))
          rowind(1:nnz)=rowind_temp(1:nnz)
          deallocate(rowind_temp)
c------------------------------------------------
  

c============================================================
c    Data file for charge discharge cycle
c============================================================


              

c====================================================================
c         Store prescribed concentration
c====================================================================
          do i=1,numnp_q
          if(idof_q(node_q(i)).lt.0)then
           Qarr(node_q(i))=Qp(-idof_q(node_q(i)))
          else
           Qarr(node_q(i))=0.0d0
          endif
          enddo



   

c=====================================================================
c  
       if(IRVESwitch.eq.1) then
       dFij_M=F_Macro_0/float(numloadSteps)
        Fij_M=0.0d0
        Fij_M(1,1)=1.0d0
        Fij_M(2,2)=1.0d0
        Fij_M(3,3)=1.0d0
       else
        dFij_M=0.0d0
       endif



      UimposedX = 0.0d0
      UimposedY = 0.0d0
      UimposedZ = 0.0d0
      dUimposedX = 0.0d0
      dUimposedY = 0.0d0
      dUimposedZ = 0.0d0
c=======================================================================

c============================================================================
c     Time loop
c============================================================================
  
      numstep = 0
      DO WHILE (numstep.le.numloadSteps-1)

 
       numstep = numstep + 1   ! time step counter


c      if (mod(numstep,numloadSteps).eq.0)then
c      flux(1,1)=-flux(1,1)
c      flux(1,1)=-flux(1,2) 
c      flux(1,3)=-flux(1,3)
c      icycle=-icycle                ! mark for lithiation/delithiation 
c      endif

c=================================================================
c   apply displacement, load, pressure and time for field variable 
c=================================================================

       time=time+delta_t       ! new field variable 


c=================================================================
c     Apply displacement Boundary condition
c=================================================================  

        xdof=0.0d0

       call ImposedDisplacement(numnp,numloadSteps,idof
     $    ,loading_disp_max_X,loading_disp_max_Y,loading_disp_max_Z
     $    ,dUimposedX,dUimposedY,dUimposedZ,xdof)


c          write(*,*)'flux=',flux
c          stop
c=================================================================
c Cyclic loading for flux  
c================================================================= 
c       if(cycloading.eq.1)then
c       do i=1,nfmark
c          do i1=1,3
c            if(flux_max(i,i1).ne.0.0)then
c       call cyc_loading(flux_max(i,i1),loadratio,frequency,numcyc,
c     $     steppercyc,numloadSteps,numstep,flux(i,i1),time,Delta_t) 
c            endif
c           enddo
c       enddo
c
c       elseif(mod(numstep-1,100/2).eq.0)then
c              flux=-flux           
c       endif !cycle loading part ends     
                
c=================================================================  



         Rq=0.0d0   
         Rflux=0.0d0

         norm_Q=1e6 
         dQf=0.d0

         dQarr_k=0.d0
         Qarr_k=0.d0

         sumQ=0.d0
         sumdQ=0.d0 





         dUarr_k=0.0d0

c=================================================================
c       iteration loop starts here
c=================================================================

        Qarr_n=Qarr           
        Uarr_n=Uarr

        Qarr_k=Qarr_n         ! Initial guess for concentration
        Uarr_k=Uarr_n         ! Initial guess for displacement

        icount=0              ! iteration counter
        epsilon_Q=0.005       ! tolerance for convergence


        DO WHILE (norm_Q.ge.epsilon_Q.and.icount.lt.icount_max)  ! convergence checking

        icount=icount+1         ! increment in iteration

c---------------initialization--------------------------------------------


c------------- calculate flux --------------------------------------------
         call flux_loading(numelvf,elem_type_f,numnp,lmq_f,x
     $                ,np_elec_max
     $                ,flux,nfmark,Rflux)  

c------------- calculate pressure-----------------------------------------
        Grad_P=0.0d0
        dp_dc=0.0d0
 
c        call PressureParameters(max_gp,numnp,numelv,numnp_q
c     $  ,numel_q,numimp_q,numat_q,elem_type_q,lmq,x
c     $  ,numatv,matvprop,MatType
c     $  ,np_elev_max,Qarr_k,Uarr_k,alpha,Grad_P,dp_dc)


c----------- calculate diffusion matrices---------------------------------


          Qvalues_ff=0.0d0 
          Qvalues_fp=0.0d0 
          Rq_f=0.0d0

         call CoupledFiledMatrices(numnp,numelv,numnp_q,numel_q
     $  ,numimp_q,numat_q,elem_type_q,alpha,max_gp,delta_t,Dq_coff,lmq
     $  ,x,np_elev_max,Qarr_k,Qarr_n,Grad_P,dp_dc,idof_q
     $  ,Qnnz,Qnnzimp,numeq_q,Qval_index_ff,Qval_index_fp
     $  ,Qvalues_ff,Qvalues_fp,Rq_f)



c          call InterfaceDiffusion(max_gp,numnp,numatc
c     $ ,numelc,numel_q,np_elec_max,matcprop,x,coh_elem_type,lmc
c     $ ,np_elev_max,CohElemNodeIndex,numat_q
c     $ ,Dq_coff,elem_type_q,lmq,Qarr_k,idof_q,numnp_q
c     $ ,Qnnz,Qnnzimp,numeq_q,numimp_q,Qval_index_ff,Qval_index_fp
c     $ ,Qvalues_ff,Qvalues_fp,Rq_f)

         


c==================================================================
c      impose field boundary condition
c      Separate free and prescribed dof
c      for both displacement and field
c===================================================================
       dQp=0.0d0
       call ImposedFieldBndry(numnp,numelv,numpres,numnp_q,numel_q
     $ ,numimp_q,numeq,numeq_q,np_elev_max
     $ ,lmq,idof,idof_q,node_q,dQp,Rflux 
     $ ,Qnnzimp,Qval_index_fp,Qvalues_fp,Rq_f) 


      

c=======================================================================
c        Solve for concentration
c=======================================================================
           
       Rq_fl=Rq_f
      call SuperLUSolver2(numeq_q,Qnnz,Qcolptr
     $         ,Qrowind,Qvalues_ff,Rq_fl)

c======================================================================
c         update concentration
c======================================================================
          dQf=Rq_fl(1:numeq_q) 
 
          do i=1,numnp_q
           if(idof_q(node_q(i)).gt.0)then
                dQarr_k(node_q(i))= dQf(idof_q(node_q(i))) 
                Qarr_k(node_q(i))= Qarr_k(node_q(i))+dQarr_k(node_q(i))
           endif
          enddo






c=======================================================================
c PART II: Solve equilibrium condition
c Generate Bulk Stiffness 
c Generate Cohesive Stiffness at the interface
c=======================================================================

c       Kuu=0.0d0
c       Ru=0.0d0

c      call Stiff_brick_gen(np_elev_max,numelv,numatv,numnp,lmv,x
c     $       ,Qarr_k,numat_q,alpha
c     $       ,elem_type,max_gp,matvprop,MatType,Uarr_k,Kuu,Ru)


       values_ff=0.0d0
       values_fp=0.0d0
              Rf=0.0d0
             dUp=0.0d0

       call Stiff_brick_gen(np_elev_max,numelv,numatv,numnp,lmv,x
     $       ,Qarr_k,numat_q,alpha,elem_type,max_gp,matvprop,MatType
     $       ,Fp_arr,PlasticPara_arr,Uarr_k
     $       ,nnz,nnzimp,numeq,numpres
     $       ,idof,val_index_ff,val_index_fp
     $       ,values_ff,values_fp,Rf)
 
c       call Stiff_brick_Damage(np_elev_max,numelv,numatv,numnp,lmv,x
c     $       ,Qarr_k,numat_q,alpha,elem_type,max_gp,matvprop,MatType
c     $       ,Fp_arr,PlasticPara_arr,Uarr_k
c     $       ,nnz,nnzimp,numeq,numpres
c     $       ,idof,val_index_ff,val_index_fp
c     $       ,values_ff,values_fp,Rf,num_Svar_max,Svar)


       call CohesiveLargeStrainOPModel3D(np_elec_max,max_gp
     $           ,numelc,numatc,numnp,numeq,numpres
     $           ,Uarr_k,coh_elem_type,matcprop,lmc,x
     $           ,CohesiveHistPara
     $           ,nnz,nnzimp,idof,val_index_ff
     $           ,val_index_fp,values_ff,values_fp,Rf)

 
 
       call ImposedDisplaceFieldMatrices(numnp,numimp,numeq
     $  ,numpres,idof,xdof,dUp,dRext
     $  ,nnzimp,val_index_fp,values_fp,Rf)




 
c     Kuu_ff=0.0d0
c      Kuu_fp=0.0d0
c      Ru_f=0.0d0
c      Ru_p=0.0d0
c      dUp=0.0d0

c      call ImposedDisplaceFieldMatrices(numnp,numimp,numeq
c     $  ,numpres,idof,xdof,dUp,dRext
c     $  ,Kuu,Ru,Kuu_ff,Kuu_fp,Ru_f,Ru_p)


c=======================================================================
c        Solve for displacement
c=======================================================================

          Ru_f=Rf
          call SuperLUSolver2(numeq,nnz,colptr,rowind,values_ff,Ru_f)



c========================================================================
c        Update nodal displacement (unprescribed only)
c========================================================================
         dUf=Ru_f


        do i = 1,3*numnp
            if(idof(1,i).gt.0) then
               dUarr_k(i)=dUf(idof(1,i))
               Uarr_k(i)=Uarr_k(i)+dUarr_k(i)          
            endif
         enddo


c======================================================================
c    Check for convergence criteria
c======================================================================

c------------- Concentration norm calculation-------------------------
         sumQ=0.d0
         sumdQ=0.d0  
         do i = 1,numnp_q
            if(idof_q(node_q(i)).gt.0) then
               sumQ=sumQ+Qarr_k(node_q(i))*Qarr_k(node_q(i))          ! L2 norm of concentration previous step
               sumdQ=sumdQ+dQarr_k(node_q(i))*dQarr_k(node_q(i))      ! L2 norm of increment current step
            endif
         enddo
   
         norm_Q=dsqrt(sumdQ/sumQ)   ! Concentration norm ratio    


 

       ENDDO         ! iteration loop ends                 

        if (icount.eq.icount_max)then
        write(*,*)'============================================'
        write(*,*)'Solution reached maximum number of iteration'
        write(*,*)'No convergence possible'
        write(*,*)'============================================'
        stop
        else
        write(*,*)'Timestep, Interation and Norm',numstep,icount,norm_Q 
        endif



c======================================================================
c update converged displacement 
c (converged and prescribed)
c======================================================================

           do i = 1,3*numnp
             if(idof(1,i).gt.0) then
               Uarr(i)=Uarr_k(i)                    ! update unprescribed displacement
             else
               Uarr(i)=Uarr_k(i)+dUp(-idof(1,i))      ! update prescribed displacement 
              endif
           enddo
       




c=========================================================================
c       Update converged concentration (main concentration array)
c=========================================================================
           do i=1,numnp_q
            if(idof_q(node_q(i)).lt.0)then
            Qarr(node_q(i))=Qp(-idof_q(node_q(i))) ! update prescribed concentration
            else
            Qarr(node_q(i))=Qarr_k(node_q(i))      ! update unprescribed concentration
            endif
          enddo


c          write(123,*)time,Qarr(1)



           UimposedX = UimposedX+dUimposedX    ! Ux=Ux+dUx
           UimposedY = UimposedY+dUimposedY    ! Uy=Uy+dUy
           UimposedZ = UimposedZ+dUimposedZ    ! Uz=Uz+dUx
 
c============================================================================
c        Tecplot and paraview output of Stress and Concentration
c============================================================================

       call Plot_tecplot_gen(np_elev_max,numnp,numelv,numatv
     $       ,lmv,x,Uarr,elem_type,max_gp,MatType,Fp_arr,PlasticPara_arr
     $       ,matvprop,numstep 
     $       ,infile,Qarr,C_max,numat_q,alpha 
     $       ,inlen,outputStep,UimposedX,UimposedY,UimposedZ)

c       call Plot_tecplot_Damage(np_elev_max,numnp,numelv,numatv
c     $       ,lmv,x,Uarr,elem_type,max_gp,MatType,Fp_arr,PlasticPara_arr
c     $       ,matvprop,numstep,num_Svar_max,Svar 
c     $       ,infile,Qarr,C_max,numat_q,alpha 
c     $       ,inlen,outputStep,UimposedX,UimposedY,UimposedZ)


c=============================================================================
c     Fractured surface calculation
c=============================================================================

       call FracturedSurfaceCalc(np_elec_max,max_gp
     $           ,numstep,outputStep,infile,inlen
     $           ,numelc,numatc,numnp,numeq,numpres
     $           ,Uarr,coh_elem_type,matcprop,lmc,x
     $           ,CohesiveHistPara,Area_Fractured)


c=============================================================================
c        Average surface concentration
c============================================================================= 
       call AverageSurfaceConc(max_gp,numelvf,elem_type_f,numnp
     $                ,lmq_f,x,np_elec_max
     $                ,nfmark,Qarr,Cbar)  


      write(320,*)numstep,Cbar,Area_Fractured
c============================================================================
c       Calculate voltage profile
c============================================================================


c============================================================================



          ENDDO ! loop for load step


          call cpu_time(c_time)
          write(*,*)'simulation time=',c_time

      end program main

