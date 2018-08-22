!> @file crs_hypre.f
!! @ingroup crs_hypre
!! @brief Coarse grid (CRS) setup and solver using the Hypre library
!! @author Nicolas Offermans
!! @date Jun 25, 2018
c-----------------------------------------------------------------------
!> @brief Wrapper for CRS setup given the field number
!! @param[in] ifield       field number
!! @param[in] comm         MPI communicator
!! @param[in] a            coarse grid matrix  
!! @param[in] nxc          pol. order of CRS problem = nxc-1
!! @param[in] np           number of MPI processes
!! @note The data for the CRS solver are stored in an include file
!!       called CRS_HYPRE.
      subroutine hypre_crs_setup(ifield,comm,a,nxc,np)
      include 'SIZE'
      include 'DOMAIN'
      include 'CRS_HYPRE'

c     Subroutine arguments
      integer ifield,comm,nxc,np
      real a(lcr,lcr,lelv)      

      call hypre_crs_setup_(hypre_solver(ifield),hypre_A(ifield),
     $     hypre_x(ifield),hypre_b(ifield),hypre_gsh(ifield),
     $     hypre_ilower(ifield),hypre_nrows(ifield),hypre_map(1,ifield),
     $     comm,a,nxc,np)
      
      return
      end
c-----------------------------------------------------------------------
!> @brief Wrapper for CRS setup in stress mode
!! @param[in] comm         MPI communicator
!! @param[in] a            coarse grid matrix  
!! @param[in] nxc          pol. order of CRS problem = nxc-1
!! @param[in] np           number of MPI processes
!! @note The data for the CRS solver are stored in an include file
!!       called CRS_HYPRE.      
      subroutine hypre_crs_setup_strs(comm,a,nxc,np)
      include 'SIZE'
      include 'DOMAIN'
      include 'CRS_HYPRE'
      
c     Subroutine arguments
      integer comm,nxc,np
      real a(lcr,lcr,lelv)

      call hypre_crs_setup_(hypre_solver_strs,hypre_A_strs,
     $     hypre_x_strs,hypre_b_strs,hypre_gsh_strs,hypre_ilower_strs,
     $     hypre_nrows_strs,hypre_map_strs,comm,a,nxc,np)
      
      return
      end
c-----------------------------------------------------------------------
!> @brief Wrapper for CRS solver given the field number
!! @param[in] ifield       field number
!! @param[out] x           solution  
!! @param[in] b            right hand side
      subroutine hypre_crs_solve(ifield,x,b)
      implicit none
      include 'SIZE'
      include 'DOMAIN'          ! lcr
      include 'CRS_HYPRE'

c     Arguments list
      integer ifield
      real x(lcr,lcr,lelv),b(lcr,lcr,lelv)

      call hypre_crs_solve_(hypre_solver(ifield),hypre_A(ifield),
     $     hypre_x(ifield),hypre_b(ifield),hypre_gsh(ifield),
     $     hypre_ilower(ifield),hypre_nrows(ifield),hypre_map(1,ifield),
     $     x,b)
      
      return
      end
c-----------------------------------------------------------------------
!> @brief Wrapper for CRS solver in stress mode
!! @param[out] x           solution  
!! @param[in] b            right hand side      
      subroutine hypre_crs_solve_strs(x,b)
      implicit none
      include 'SIZE'
      include 'DOMAIN'          ! lcr
      include 'CRS_HYPRE'

c     Arguments list
      real x(lcr,lcr,lelv),b(lcr,lcr,lelv)

      call hypre_crs_solve_(hypre_solver_strs,hypre_A_strs,
     $     hypre_x_strs,hypre_b_strs,hypre_gsh_strs,
     $     hypre_ilower_strs,hypre_nrows_strs,hypre_map_strs,x,b)
      
      return
      end      
c-----------------------------------------------------------------------
!> @brief Wrapper for free CRS data given the field number
!! @param[in] ifield       field number
      subroutine hypre_crs_free (ifield)
      implicit none
      include 'SIZE'
      include 'CRS_HYPRE'

c     Arguments list      
      integer ifield

      call hypre_crs_free_ (hypre_solver(ifield),hypre_A(ifield),
     $     hypre_x(ifield),hypre_b(ifield),hypre_gsh(ifield))
      
      return
      end
c-----------------------------------------------------------------------
!> @brief Wrapper for freeing CRS data in stress mode 
      subroutine hypre_crs_free_strs
      implicit none
      include 'SIZE'
      include 'CRS_HYPRE'

      call hypre_crs_free_ (hypre_solver_strs,hypre_A_strs,
     $     hypre_x_strs,hypre_b_strs,hypre_gsh_strs)
      
      return
      end
c-----------------------------------------------------------------------
#ifdef HYPRE
c-----------------------------------------------------------------------
!> @brief Actual routine performing the CRS setup
!! @param[out] solver      handle to the Hypre AMG solver
!! @param[out] ij_A        Hypre matrix in IJ format
!! @param[out] ij_x        Hypre solution vector in IJ format
!! @param[out] ij_b        Hypre right hand side vector in IJ format
!! @param[out] gsh         gather-scatter handle for the CRS problem
!! @param[out] ilower      lowest row number owned by the process
!! @param[out] un          number of CRS nodes owned by process
!! @param[out] crs2hypre   mapping from CRS to Hypre data layout
!! @param[in] comm         MPI communicator
!! @param[in] a            coarse grid matrix  
!! @param[in] nxc          pol. order of CRS problem = nxc-1
!! @param[in] np           number of MPI processes
!! @note The vectors ij_x and ij_b are only initialized.      
      subroutine hypre_crs_setup_(solver,ij_A,ij_x,ij_b,gsh,ilower,
     $     un,crs2hypre,comm,a,nxc,np)
      implicit none
      include 'SIZE'
      include 'DOMAIN'          ! lcr,se_to_gcrs
      include 'HYPREf.h'

c     Arguments list
      integer*8 solver,ij_A,ij_x,ij_b,un,ilower
      integer gsh,crs2hypre(lelt),comm,nxc,np
      real a(lcr,lcr,lelv)

c     Hypre data
      integer glhid(lcr,lelv) ! global Hypre id
      integer*8 par_A,par_x,par_b
      integer*8 iupper,jlower,jupper,nrows,ncols,rows,cols
      integer*8 ierr
      real values
      
c     GS data 
      integer*8 uid(lcr,lelv) ! unique id array
      
c     Crystal router data
      integer crh,nmax
      parameter (nmax=1.2*lcr*(lelv+20)) ! somewhat arbitrary value
      integer vi(lcr+2,nmax)
      integer*8 vl
      real vr(lcr,nmax)

c     Other
      integer ir,ic,iel,icr,idh,ncr,ncrv,n
      integer igl_running_sum
      integer pid_owner(lcr,lelv)
      integer*8 i8

      common /scrxxti/ glhid,pid_owner

      ncr = nxc**ldim
      ncrv = ncr*nelv

c     Unique handle to determine node ownership
      call i8copy(uid,se_to_gcrs,ncrv)
      call fgslib_gs_unique(uid,ncrv,comm,np)
      call fgslib_gs_setup(gsh,uid,ncrv,comm,np)

c     Identify process id of owner
      un = 0
      do iel=1,nelv
         do icr=1,ncr
            if (uid(icr,iel).gt.0) then
               pid_owner(icr,iel)=nid
               un=un+1
            else
               pid_owner(icr,iel)=0
            endif 
         enddo
      enddo
      call fgslib_gs_op(gsh,pid_owner,2,1,0) ! scatter

c     Renumber nodes according to Hypre's requirement
      iupper=igl_running_sum(un)
      ilower=iupper-un+1
      jlower=ilower
      jupper=iupper

      idh=ilower
      do iel=1,nelv
         do icr=1,ncr
            if (uid(icr,iel).gt.0) then
               glhid(icr,iel)=idh
               crs2hypre(idh-ilower+1)=icr+(iel-1)*ncr
               idh=idh+1
            else
               glhid(icr,iel)=0
            endif
         enddo
      enddo

      if ((idh-1).ne.iupper) then
         call exitti('Error when building CRS matrix: idh~=iupper.$',1)
      endif
      
      call fgslib_gs_op(gsh,glhid,2,1,0) ! scatter

c     Transfer local matrix to process owner
      call fgslib_crystal_setup(crh,comm,np)
      do iel=1,nelv
         do icr=1,ncr
            ic=(iel-1)*ncr+icr
            do ir=1,ncr
               vi(ir,ic)=glhid(ir,iel) ! column number
               vr(ir,ic)=a(icr,ir,iel) ! matrix entry
            enddo
            vi(ncr+1,ic)=glhid(icr,iel) ! row number
            vi(ncr+2,ic)=pid_owner(icr,iel) ! id of proc. owner
         enddo
      enddo

      n=ncrv
      call fgslib_crystal_tuple_transfer(crh,n,nmax,vi,ncr+2,
     $     vl,0,vr,ncr,ncr+2)

      call fgslib_crystal_free(crh)      

      if (n.gt.nmax) then
         call exitti('Error when building CRS matrix: n>nmax.
     $ Increase nmax.$',nmax)
      endif

c     Initialize Hypre matrix
      call HYPRE_IJMatrixCreate(comm,ilower,iupper,jlower,jupper,ij_A,
     $     ierr)
      call HYPRE_IJMatrixSetObjectType(ij_A,int(HYPRE_PARCSR,kind(i8)),
     $     ierr)
      call HYPRE_IJMatrixInitialize(ij_A,ierr)
      
c     Fill in Hypre matrix
      nrows=1
      ncols=1
      do ic=1,n
         rows=vi(ncr+1,ic)
         do ir=1,ncr
            cols=vi(ir,ic)
            values=vr(ir,ic)
            if (rows.ne.0.and.cols.ne.0.and.values.ne.0.) then
               call HYPRE_IJMatrixAddToValues(ij_A,nrows,ncols,
     $              rows,cols,values,ierr)
            endif
         enddo
      enddo

c     Assemble matrix
      call HYPRE_IJMatrixAssemble(ij_A,ierr)

c     Create AMG solver
      call HYPRE_BoomerAMGCreate(solver,ierr)

c     Set AMG parameters
      call HYPRE_BoomerAMGSetPrintLevel(solver,int(1,kind(i8)),ierr) ! Verbose level: nothing(->0)/setup info(->1)/solver info(->2)/both(->3)
      call HYPRE_BoomerAMGSetCoarsenType(solver,int(8,kind(i8)),ierr) ! PMIS
      call HYPRE_BoomerAMGSetInterpType(solver,int(14,kind(i8)),ierr) ! extended interpolation
      call HYPRE_BoomerAMGSetRelaxType(solver,int(8,kind(i8)),ierr) ! l1-scaled hybrid symmetric Gauss-Seidel
      call HYPRE_BoomerAMGSetMaxCoarseSize(solver,int(5,kind(i8)),ierr) 
      call HYPRE_BoomerAMGSetStrongThrshld(solver,0.5,ierr) ! Increase for better convergence. Decrease for faster solver time.
      call HYPRE_BoomerAMGSetMeasureType(solver,int(1,kind(i8)),ierr) ! Local(->0)/Global(->1) measure
      call HYPRE_BoomerAMGSetTol(solver,0.,ierr) ! Decrease for better convergence. Increase for faster solver time.
      call HYPRE_BoomerAMGSetMaxIter(solver,int(1,kind(i8)),ierr) ! Increase for better convergence. Decrease for faster solver time.

c     Create and initialize rhs and solution vectors
      call HYPRE_IJVectorCreate(comm,jlower,jupper,ij_b,ierr)
      call HYPRE_IJVectorSetObjectType(ij_b,int(HYPRE_PARCSR,kind(i8)),
     $     ierr)
      call HYPRE_IJVectorInitialize(ij_b,ierr)
      call HYPRE_IJVectorAssemble(ij_b,ierr)
      
      call HYPRE_IJVectorCreate(comm,jlower,jupper,ij_x,ierr)
      call HYPRE_IJVectorSetObjectType(ij_x,int(HYPRE_PARCSR,kind(i8)),
     $     ierr)
      call HYPRE_IJVectorInitialize(ij_x,ierr)
      call HYPRE_IJVectorAssemble(ij_x,ierr)

c     Perform AMG setup
      call HYPRE_IJMatrixGetObject(ij_A,par_A,ierr)
      call HYPRE_IJVectorGetObject(ij_b,par_b,ierr)
      call HYPRE_IJVectorGetObject(ij_x,par_x,ierr)      

      call HYPRE_BoomerAMGSetup(solver,par_A,par_b,par_x,ierr)

      return
      end
c-----------------------------------------------------------------------
!> @brief Actual routine performing the CRS solver
!! @param[in] solver      handle to the Hypre AMG solver
!! @param[in] ij_A         Hypre matrix in IJ format
!! @param[in] ij_x         Hypre solution vector in IJ format
!! @param[in] ij_b         Hypre right hand side vector in IJ format
!! @param[in] gsh          gather-scatter handle for the CRS problem
!! @param[in] hil          lowest row number owned by the process
!! @param[in] un           number of CRS nodes owned by process
!! @param[in] crs2hypre    mapping from CRS to Hypre data layout
!! @param[out] x           solution
!! @param[in] b            right hand side
!! @note Vectors ij_x and ij_b only need to be initialized at input,
!!       their entries are set in the routine.
      subroutine hypre_crs_solve_(solver,ij_A,ij_x,ij_b,gsh,hil,un,
     $     crs2hypre,x,b)
      implicit none
      include 'SIZE'
      include 'DOMAIN'

c     Arguments list      
      integer*8 solver,ij_A,ij_x,ij_b
      integer*8 hil,un
      integer crs2hypre(lelt)
      integer gsh
      real x(1),b(1)

c     Local variables      
      integer*8 par_A,par_x,par_b
      integer*8 ierr,i,ir,i8
      
      call fgslib_gs_op(gsh,b,1,1,1) ! gather

      do i=1,un
         ir=hil+i-1
         call HYPRE_IJVectorSetValues(ij_b,int(1,kind(i8)),ir,
     $        b(crs2hypre(i)),ierr)
         call HYPRE_IJVectorSetValues(ij_x,int(1,kind(i8)),ir,0.,ierr)
      enddo

      call HYPRE_IJVectorAssemble(ij_b,ierr)
      call HYPRE_IJVectorGetObject(ij_b,par_b,ierr)

      call HYPRE_IJVectorAssemble(ij_x,ierr)
      call HYPRE_IJVectorGetObject(ij_x,par_x,ierr)

      call HYPRE_IJMatrixGetObject(ij_A,par_A,ierr)
      call HYPRE_BoomerAMGSolve(solver,par_A,par_b,par_x,ierr)

      do i=1,un
         ir=hil+i-1
         call HYPRE_IJVectorGetValues(ij_x,int(1,kind(i8)),ir,
     $        x(crs2hypre(i)),ierr)
      enddo

      call fgslib_gs_op(gsh,x,1,1,0) ! scatter
      
      return
      end
c-----------------------------------------------------------------------
!> @brief Actual routine for freeing CRS data
!! @param[in] solver       handle to the Hypre AMG solver
!! @param[in] ij_A         Hypre matrix in IJ format
!! @param[in] ij_x         Hypre solution vector in IJ format
!! @param[in] ij_b         Hypre right hand side vector in IJ format
!! @param[in] gsh          gather-scatter handle for the CRS problem      
      subroutine hypre_crs_free_ (solver,ij_A,ij_x,ij_b,gsh)
      implicit none
      include 'SIZE'
      include 'CRS_HYPRE'
      
c     Arguments list      
      integer*8 solver,ij_A,ij_x,ij_b
      integer gsh

c     Local variables     
      integer*8 ierr

      call HYPRE_BoomerAMGDestroy(solver,ierr)
      call HYPRE_IJMatrixDestroy(ij_A,ierr)
      call HYPRE_IJVectorDestroy(ij_x,ierr)
      call HYPRE_IJVectorDestroy(ij_b,ierr)
      call fgslib_gs_free(gsh)
      
      return
      end
c-----------------------------------------------------------------------
#else
c-----------------------------------------------------------------------
      subroutine hypre_crs_solve_(solver,ij_A,ij_x,ij_b,gsh,hil,un,
     $     crs2hypre,x,b)
      call exitti("Please recompile with HYPRE support.$",1)
      return
      end
c-----------------------------------------------------------------------
      subroutine hypre_crs_setup_(solver,ij_A,ij_x,ij_b,gsh,ilower,
     $     un,crs2hypre,comm,a,nxc,np)
      call exitti("Please recompile with HYPRE support.$",1)
      return
      end
c-----------------------------------------------------------------------
      subroutine hypre_crs_free_ (solver,ij_A,ij_x,ij_b,gsh)
      call exitti("Please recompile with HYPRE support.$",1)
      return
      end
#endif
